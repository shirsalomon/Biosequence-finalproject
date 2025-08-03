import gzip
import pickle
import json
import logging
from typing import Dict, List
from collections import defaultdict


logging.basicConfig(filename='reference_database.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class ReferenceDatabase:
    """Class for managing the k-mer reference database."""
    
    def __init__(self, k: int, filter_similar: bool = False, similarity_threshold: float = 0.95):
        """Initializes the reference database with the given k-mer size."""
        self.k: int = k
        self.reference: Dict[str, Dict[str, List[int]]] = {}
        self.genome_lengths: Dict[str, int] = {}  
        self.filter_similar = filter_similar  
        self.similarity_threshold = similarity_threshold  

    def build_from_fasta(self, filename_or_genomes) -> None:
        """Builds the k-mer reference database from a FASTA file or a dictionary of genomes."""
        logging.info(f"Building reference database")

        if isinstance(filename_or_genomes, dict): 
            genomes = filename_or_genomes
        else: 
            genomes: Dict[str, str] = self._read_fasta(filename_or_genomes)

        if not genomes:
            logging.error("No valid genomes found")
            raise ValueError("Error: No valid genomes found")
            
        if getattr(self, "filter_similar", False):
            genomes = self.filter_similar_genomes(genomes)

        for genome_id, sequence in genomes.items():
            self.genome_lengths[genome_id] = len(sequence)

        for genome_id, sequence in genomes.items():
            kmers = self._generate_kmers(sequence)
            if kmers:
                for kmer, positions in kmers.items():
                    if kmer not in self.reference:
                        self.reference[kmer] = {}
                    self.reference[kmer][genome_id] = positions

        logging.info(f"Successfully built reference database with {len(self.reference)} k-mers")





    def save(self, filename: str) -> None:
        """Saves the reference database to a compressed .kdb file."""
        with gzip.open(filename, 'wb') as file:
            pickle.dump(self, file, protocol=pickle.HIGHEST_PROTOCOL)
        logging.info(f"Reference database saved to {filename}")


    @classmethod
    def load(cls, filename: str):
        """Loads the reference database from a compressed .kdb file."""
        try:
            with gzip.open(filename, 'rb') as file:
                obj = pickle.load(file)
                if isinstance(obj, cls):
                    logging.info(f"Reference database loaded from {filename}")
                    return obj
                else:
                    logging.error(f"Loaded object is not of type ReferenceDatabase")
                    return None
        except Exception as e:
            logging.error(f"Error loading reference database: {e}")
            return None


    def dump_json(self, filename: str) -> None:
        """Dumps a detailed summary of the reference database in JSON format as required."""

        kmers_output = {}
        genome_kmer_count = {genome: {"unique": 0, "multi": 0} for genome in self.genome_lengths}
        kmer_occurrence = {}

        for kmer, genome_map in self.reference.items():
            kmers_output[kmer] = {}
            for genome_id, positions in genome_map.items():
                kmers_output[kmer][genome_id] = sorted(positions)

            kmer_occurrence[kmer] = len(genome_map)

        for kmer, genome_map in self.reference.items():
            if len(genome_map) == 1:
                genome_id = next(iter(genome_map))
                genome_kmer_count[genome_id]["unique"] += 1
            else:
                for genome_id in genome_map:
                    genome_kmer_count[genome_id]["multi"] += 1

        summary_output = {}
        for genome_id, counts in genome_kmer_count.items():
            summary_output[genome_id] = {
                "total_bases": self.genome_lengths.get(genome_id, 0),
                "unique_kmers": counts["unique"],
                "multi_mapping_kmers": counts["multi"]
            }

        output_data = {
            "Kmers": kmers_output,
            "Summary": summary_output
        }

        with open(filename, 'w') as f:
            json.dump(output_data, f, indent=4)
        logging.info(f"Reference summary saved to {filename}")
    
    def dump_json_stdout(self):
        kmers = {}
        genome_kmer_count = defaultdict(lambda: {"unique_kmers": 0, "multi_mapping_kmers": 0})
        genome_lengths = self.genome_lengths.copy()

        for kmer, genome_dict in self.reference.items():
            kmers[kmer] = {gid: sorted(pos) for gid, pos in genome_dict.items()}

        for kmer, genome_dict in self.reference.items():
            if len(genome_dict) == 1:
                genome_id = next(iter(genome_dict))
                genome_kmer_count[genome_id]["unique_kmers"] += 1
            else:
                for genome_id in genome_dict:
                    genome_kmer_count[genome_id]["multi_mapping_kmers"] += 1

        summary = {
            "Kmers": kmers,
            "Summary": {
                genome_id: {
                    "total_bases": genome_lengths.get(genome_id, 0),
                    "unique_kmers": genome_kmer_count[genome_id]["unique_kmers"],
                    "multi_mapping_kmers": genome_kmer_count[genome_id]["multi_mapping_kmers"]
                } for genome_id in genome_lengths
            }
        }
            
        if hasattr(self, 'similarity_report'):
            summary["Similarity"] = self.similarity_report

        print(json.dumps(summary, indent=4))


        

    def dump_json_to_file(self, output_path: str):
        kmers = {}
        genome_kmer_count = defaultdict(lambda: {"unique_kmers": 0, "multi_mapping_kmers": 0})
        genome_lengths = self.genome_lengths.copy()

        for kmer, genome_dict in self.reference.items():
            kmers[kmer] = {gid: sorted(pos) for gid, pos in genome_dict.items()}

        for kmer, genome_dict in self.reference.items():
            if len(genome_dict) == 1:
                genome_id = next(iter(genome_dict))
                genome_kmer_count[genome_id]["unique_kmers"] += 1
            else:
                for genome_id in genome_dict:
                    genome_kmer_count[genome_id]["multi_mapping_kmers"] += 1

        output = {
            "Kmers": kmers,
            "Summary": {
                genome_id: {
                    "total_bases": genome_lengths.get(genome_id, 0),
                    "unique_kmers": genome_kmer_count[genome_id]["unique_kmers"],
                    "multi_mapping_kmers": genome_kmer_count[genome_id]["multi_mapping_kmers"]
                } for genome_id in genome_lengths
            }
        }

        if hasattr(self, "similarity_report") and self.similarity_report:
            output["Similarity"] = self.similarity_report

        with open(output_path, 'w') as f:
            json.dump(output, f, indent=4)




    def _read_fasta(self, filename: str) -> Dict[str, str]:
        """Reads a FASTA file and returns a dictionary with genome IDs and their sequences."""
        genomes = {}
        genome_id = None
        sequence = []
        try:
            with open(filename, 'r') as file:
                for line in file:
                    if line.startswith(">"):
                        if genome_id:
                            genomes[genome_id] = ''.join(sequence)
                        genome_id = line[1:].strip()
                        sequence = []
                    else:
                        sequence.append(line.strip())
                if genome_id:
                    genomes[genome_id] = ''.join(sequence)
            return genomes
        except FileNotFoundError:
            logging.error(f"The file {filename} does not exist.")
            return {}
        

    def _generate_kmers(self, sequence: str) -> Dict[str, List[int]]:
        """Generates all k-mers of length k from a given DNA sequence."""
        if not sequence or len(sequence) < self.k:
            return {}

        kmers = {}
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            if 'N' not in kmer:
                if kmer not in kmers:
                    kmers[kmer] = []
                kmers[kmer].append(i)
        
        return kmers
    



    def filter_similar_genomes(self, genomes: Dict[str, str]) -> Dict[str, str]:
        logging.info("Filtering similar genomes")
        kept_genomes = {}
        removed_genomes = {}
        similarity_report = {}

        genome_kmer_sets = {
            gid: set(self._generate_kmers(seq).keys())
            for gid, seq in genomes.items()
        }

        genome_lengths = {gid: len(seq) for gid, seq in genomes.items()}

        genome_stats = []
        for gid, kmers in genome_kmer_sets.items():
            total_kmers = len(kmers)
            unique_kmers = sum(
                1 for k in kmers if sum(gid in self.reference.get(k, {}) for gid in genomes) == 1
            )
            genome_stats.append((gid, unique_kmers, total_kmers, genome_lengths[gid]))

        genome_stats.sort(key=lambda x: (x[1], x[2], x[3]))

        G = set(gid for gid, _, _, _ in genome_stats)

        for i, (gid1, uniq1, total1, len1) in enumerate(genome_stats):
            if gid1 not in G:
                continue
            K1 = genome_kmer_sets[gid1]
            similarity_report[gid1] = {
                "kept": "yes",
                "unique_kmers": uniq1,
                "total_kmers": total1,
                "genome_length": len1,
                "similar_to": "NA",
                "similarity_score": "NA"
            }
            for j in range(i + 1, len(genome_stats)):
                gid2, uniq2, total2, len2 = genome_stats[j]
                if gid2 not in G:
                    continue
                K2 = genome_kmer_sets[gid2]
                shared = len(K1 & K2)
                similarity = shared / min(len(K1), len(K2))
                if similarity > self.similarity_threshold:
                    G.remove(gid2)
                    similarity_report[gid2] = {
                        "kept": "no",
                        "unique_kmers": uniq2,
                        "total_kmers": total2,
                        "genome_length": len2,
                        "similar_to": gid1,
                        "similarity_score": round(similarity, 2)
                    }

        self.similarity_report = similarity_report
        return {g: genomes[g] for g in G}




    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement.get(base, base) for base in reversed(sequence))
