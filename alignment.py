import json
import os
import logging
from reference import ReferenceDatabase
from typing import Dict, List, Set
from collections import defaultdict


logging.basicConfig(filename='alignment.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class PseudoAligner:
    """Class for performing pseudo-alignment."""

    def __init__(self, reference_file: str, reads_file: str, output_file: str,
                min_read_quality: int = None, min_kmer_quality: int = None, max_genomes: int = None,
                reverse_complement: bool = False, coverage: bool = False, genomes: List[str] = None,
                min_coverage: int = 1, full_coverage: bool = False, kmer_size: int = 31,
                m: int = 1, p: int = 1):
        self.reference_db = ReferenceDatabase.load(reference_file) if reference_file else None
        self.reads_file = reads_file
        self.output_file = output_file
        self.min_read_quality = min_read_quality
        self.min_kmer_quality = min_kmer_quality
        self.max_genomes = max_genomes
        self.reverse_complement = reverse_complement
        self.coverage = coverage
        self.genomes = genomes if genomes else []
        self.min_coverage = min_coverage
        self.full_coverage = full_coverage
        self.kmer_size = kmer_size
        self.coverage_data = {}
        self.match_counts = {}
        self.unique_counts = defaultdict(int)
        self.ambiguous_counts = defaultdict(int)
        self.unmapped_reads = 0
        self.filtered_quality_kmers = 0
        self.filtered_hr_kmers = 0
        self.m = m
        self.p = p

        logging.info("PseudoAligner initialized with given parameters.")


    def run(self) -> None:
        """Performs pseudo-alignment of reads to the reference database with quality filtering and coverage tracking."""
        filtered_quality_reads = 0
        filtered_quality_kmers = 0
        filtered_hr_kmers = 0

        logging.info("Starting pseudo-alignment.")

        reads = PseudoAligner._read_fastq(self.reads_file)
        results = {}

        if self.coverage:
            self._initialize_coverage()

        for read, quality in reads:
            if self.min_read_quality and self._calculate_quality_mean(quality) < self.min_read_quality:
                filtered_quality_reads += 1
                continue

            if self.reverse_complement:
                forward_matches = self._find_matches(read, quality)
                reverse_read = self.reverse_complement_seq(read)
                reverse_quality = quality[::-1]
                reverse_matches = self._find_matches(reverse_read, reverse_quality)

                forward_unique = len(forward_matches) if len(forward_matches) == 1 else 0
                reverse_unique = len(reverse_matches) if len(reverse_matches) == 1 else 0

                if reverse_unique > forward_unique:
                    selected_matches = reverse_matches
                    selected_read = reverse_read
                else:
                    selected_matches = forward_matches
                    selected_read = read
            else:
                selected_matches = self._find_matches(read, quality)
                selected_read = read

            if self.coverage:
                self._update_coverage(selected_matches, selected_read)

            results[read] = list(selected_matches) if selected_matches else ["Unmapped"]

        results["filtering_stats"] = {
            "filtered_quality_reads": filtered_quality_reads,
            "filtered_quality_kmers": filtered_quality_kmers,
            "filtered_hr_kmers": filtered_hr_kmers
        }

        with open(self.output_file, 'w') as file:
            json.dump(results, file, indent=4)

        logging.info(f"Alignment results saved to {self.output_file}")

        if not results:
            logging.warning("No alignment results were saved! Check filtering parameters.")

        if self.coverage:
            self._save_coverage_data()

        self.save_alignment_summary(self.output_file)


    def _find_matches(self, read: str, quality: str) -> Set[str]:
        """
        Performs pseudo-alignment using specific and total k-mer matching.
        Implements rules based on parameters m and p.
        Supports reverse complement alignment if --reverse-complement is enabled.
        """

        k = self.reference_db.k
        m = self.m
        p = self.p
        kmer_to_genomes = self.reference_db.reference

        def get_counts(sequence, qual):
            specific = defaultdict(int)
            total = defaultdict(int)

            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]

                if "N" in kmer:
                    continue

                if self.min_kmer_quality and self._calculate_quality_mean(qual[i:i + k]) < self.min_kmer_quality:
                    self.filtered_quality_kmers += 1
                    continue

                if kmer in kmer_to_genomes:
                    genome_ids = list(kmer_to_genomes[kmer].keys())

                    if self.max_genomes and len(genome_ids) > self.max_genomes:
                        self.filtered_hr_kmers += 1
                        continue

                    if len(genome_ids) == 1:
                        specific[genome_ids[0]] += 1

                    for genome in genome_ids:
                        total[genome] += 1

            return specific, total

        fwd_specific, fwd_total = get_counts(read, quality)

        if self.reverse_complement:
            rev_read = self.reverse_complement_seq(read)
            rev_quality = quality[::-1]
            rev_specific, rev_total = get_counts(rev_read, rev_quality)

            fwd_best = max(fwd_specific.values()) if fwd_specific else 0
            rev_best = max(rev_specific.values()) if rev_specific else 0

            if rev_best > fwd_best:
                specific_counts = rev_specific
                total_counts = rev_total
            else:
                specific_counts = fwd_specific
                total_counts = fwd_total
        else:
            specific_counts = fwd_specific
            total_counts = fwd_total

        if not specific_counts:
            self.unmapped_reads += 1
            return set()

        sorted_specific = sorted(specific_counts.items(), key=lambda x: x[1], reverse=True)
        top_genome, top_count = sorted_specific[0]
        second_count = sorted_specific[1][1] if len(sorted_specific) > 1 else 0
        diff = top_count - second_count

        if diff >= m:
            top_total = total_counts[top_genome]
            max_total = max(total_counts.values())
            if (max_total - top_total) > p:
                ambiguous = [genome for genome, spec in specific_counts.items()
                            if total_counts[genome] >= top_total and spec > 0]
                for genome in ambiguous:
                    self.ambiguous_counts[genome] += 1
                return set(ambiguous)
            else:
                self.unique_counts[top_genome] += 1
                return {top_genome}
        else:
            ambiguous = [genome for genome, count in specific_counts.items() if count >= second_count]
            for genome in ambiguous:
                self.ambiguous_counts[genome] += 1
            return set(ambiguous)

    def save_alignment_summary(self, alignment_filename: str):
        output_file = f"{os.path.splitext(alignment_filename)[0]}_summary.json"

        summary = {
            "Statistics": {
                "unique_mapped_reads": sum(self.unique_counts.values()),
                "ambiguous_mapped_reads": sum(self.ambiguous_counts.values()),
                "unmapped_reads": self.unmapped_reads
            },
            "Summary": {}
        }

        all_genomes = set(self.unique_counts.keys()) | set(self.ambiguous_counts.keys())
        for genome in all_genomes:
            summary["Summary"][genome] = {
                "unique_reads": self.unique_counts.get(genome, 0),
                "ambiguous_reads": self.ambiguous_counts.get(genome, 0)
            }

        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=4)
            

    
    def _initialize_coverage(self):
        """Initializes coverage tracking."""
        logging.info("Initializing coverage tracking.")
        genome_lengths = {}

        for kmer, genome_mapping in self.reference_db.reference.items():
            for genome_id, positions in genome_mapping.items():
                max_pos = (max(positions) + self.reference_db.k) if positions else self.reference_db.k
                genome_lengths[genome_id] = max(genome_lengths.get(genome_id, 0), max_pos)

        for genome_id, genome_length in genome_lengths.items():
            self.coverage_data[genome_id] = {
                "unique_cov": [0] * genome_length,
                "ambiguous_cov": [0] * genome_length
            }
            self.match_counts[genome_id] = 0

        logging.info("Coverage initialized.")


    def _update_coverage(self, matched_genomes: Set[str], read: str):
        """Updates coverage maps for uniquely and ambiguously mapped reads."""
        if not matched_genomes:
            return 

        k = self.reference_db.k
        read_kmers = set()

        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            read_kmers.add(kmer)

        for kmer in read_kmers:
            if kmer not in self.reference_db.reference:
                continue
            for genome_id, positions in self.reference_db.reference[kmer].items():
                if genome_id not in matched_genomes:
                    continue
                for pos in positions:
                    if len(matched_genomes) == 1:
                        if pos < len(self.coverage_data[genome_id]["unique_cov"]):
                            self.coverage_data[genome_id]["unique_cov"][pos] += 1
                    else:
                        if pos < len(self.coverage_data[genome_id]["ambiguous_cov"]):
                            self.coverage_data[genome_id]["ambiguous_cov"][pos] += 1



    def _save_coverage_data(self):
        """Save genome coverage data to the output file."""
        if not self.coverage:
            return

        coverage_summary = {}
        coverage_details = {}

        for genome_id, data in self.coverage_data.items():
            unique = data["unique_cov"]
            ambiguous = data["ambiguous_cov"]
            total_length = len(unique)

            covered_unique = sum(1 for x in unique if x >= self.min_coverage)
            covered_ambiguous = sum(1 for x in ambiguous if x >= self.min_coverage)
            mean_unique = round(sum(unique) / total_length, 1)
            mean_ambiguous = round(sum(ambiguous) / total_length, 1)

            coverage_summary[genome_id] = {
                "covered_bases_unique": covered_unique,
                "covered_bases_ambiguous": covered_ambiguous,
                "mean_coverage_unique": mean_unique,
                "mean_coverage_ambiguous": mean_ambiguous
            }

            if self.full_coverage:
                coverage_details[genome_id] = {
                    "unique_cov": unique,
                    "ambiguous_cov": ambiguous
                }

        result = {
            "Coverage": coverage_summary
        }

        if self.full_coverage:
            result["Details"] = coverage_details

        with open(self.output_file, 'w') as outfile:
            json.dump(result, outfile, indent=4)



    @staticmethod
    def _calculate_quality_mean(quality_string: str) -> float:
        """Calculates the mean quality score from a FASTQ quality string using Phred+33."""
        if not quality_string:
            return 0.0  

        return sum(ord(q) - 33 for q in quality_string) / len(quality_string)

    @staticmethod
    def reverse_complement_seq(seq):
        complement_map = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
            'N': 'N'
        }
        try:
            complement_seq = "".join(complement_map.get(base.upper(), 'N') for base in seq)
            return complement_seq[::-1]
        except Exception as e:
            logging.error(f"Error generating reverse complement: {e}")
            return seq[::-1]


    @staticmethod
    def _read_fastq(filename: str):
        """Reads a FASTQ file and returns a list of (read_sequence, quality_score)."""
        with open(filename, 'r') as file:
            lines = [line.strip() for line in file if line.strip()]  
        
        reads = []
        for i in range(0, len(lines), 4):
            if i + 1 < len(lines) and i + 3 < len(lines): 
                read_seq = lines[i + 1].strip()
                quality = lines[i + 3].strip() if (i + 3 < len(lines)) else ""
                if read_seq.startswith("@") or read_seq == "": 
                    continue  
                reads.append((read_seq, quality))
        
        return reads
    
    @staticmethod
    def dump_json(filename: str, output_path: str):
        """Dumps the alignment results in the desired summary JSON format."""
        try:
            with open(filename, 'r') as file:
                alignment_data = json.load(file)

            summary = {
                "Statistics": {
                    "unique_mapped_reads": 0,
                    "ambiguous_mapped_reads": 0,
                    "unmapped_reads": 0
                },
                "Summary": {}
            }

            for read, genomes in alignment_data.items():
                if read == "filtering_stats":
                    continue
                if genomes == ["Unmapped"]:
                    summary["Statistics"]["unmapped_reads"] += 1
                    continue

                if len(genomes) == 1:
                    summary["Statistics"]["unique_mapped_reads"] += 1
                    category = "unique_reads"
                else:
                    summary["Statistics"]["ambiguous_mapped_reads"] += 1
                    category = "ambiguous_reads"

                for genome in genomes:
                    if genome not in summary["Summary"]:
                        summary["Summary"][genome] = {
                            "unique_reads": 0,
                            "ambiguous_reads": 0
                        }
                    summary["Summary"][genome][category] += 1

            if isinstance(output_path, str):
                with open(output_path, 'w') as out_file:
                    json.dump(summary, out_file, indent=4)
            else:
                json.dump(summary, output_path, indent=4)


        except FileNotFoundError:
            print("Error: The file was not found. Please make sure the file exists at the specified path.")





