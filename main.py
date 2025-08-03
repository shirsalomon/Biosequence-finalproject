import argparse
import logging
import os
import sys
import json
from reference import ReferenceDatabase
from alignment import PseudoAligner


def check_file_exists(filepath, file_type="input"):
    """Check if a file exists and is not empty."""
    if not os.path.isfile(filepath) or os.stat(filepath).st_size == 0:
        print(f"Error: {file_type} file '{filepath}' is missing or empty.")
        sys.exit(1)


logging.basicConfig(filename='pseudo_aligner.log', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def main():
    parser = argparse.ArgumentParser(description="Shotgun Metagenomics Pseudo-aligner")
    parser.add_argument("-t", "--task", required=True,
                        choices=["reference", "dumpref", "align", "dumpalign"],
                        help="Task to perform")

    parser.add_argument("-g", "--genome", help="FASTA file with reference genomes")
    parser.add_argument("-r", "--reference", help="Reference database file (.kdb)")
    parser.add_argument("-k", "--kmer_size", type=int, default=31, help="K-mer size (default: 31)")
    parser.add_argument("-a", "--alignment", help="Output alignment file (.aln)")
    parser.add_argument("--reads", help="FASTQ file with sequencing reads")
    parser.add_argument("-m", type=int, default=1, help="Unique k-mers threshold (default: 1)")
    parser.add_argument("-p", type=int, default=1, help="Ambiguous k-mers threshold (default: 1)")


    # Parameters for filtering similar genomes
    parser.add_argument("--filter-similar", action="store_true", help="Remove highly similar genomes from reference.")
    parser.add_argument("--similarity-threshold", type=float, default=0.95, help="Similarity threshold for genome filtering.")


    # Quality filtering options
    parser.add_argument("--min-read-quality", type=int, default=None, help="Minimum mean quality for a read to be considered.")
    parser.add_argument("--min-kmer-quality", type=int, default=None, help="Minimum mean quality for an individual k-mer.")
    parser.add_argument("--max-genomes", type=int, default=None, help="Ignore k-mers that map to more than this number of genomes.")

    # Reverse complement option
    parser.add_argument("--reverse-complement", action="store_true",
                        help="Perform pseudo-alignment on both forward and reverse complement.")

    # Coverage tracking options
    parser.add_argument("--coverage", action="store_true", help="Track genome coverage during alignment.")
    parser.add_argument("--genomes", help="Comma-separated list of genomes to calculate coverage for (default: all).")
    parser.add_argument("--min-coverage", type=int, default=1, help="Minimum coverage depth to report (default: 1).")
    parser.add_argument("--full-coverage", action="store_true", help="Output detailed position coverage.")


    args = parser.parse_args()
    logging.info(f"Running task: {args.task}")

    if args.task == "reference":
        if not args.genome or not args.reference:
            raise ValueError("For reference building, both --genome and --reference are required.")

        logging.info(f"Building reference database from {args.genome}, k-mer size: {args.kmer_size}")
        check_file_exists(args.genome, "Reference genome")
        db = ReferenceDatabase(k=args.kmer_size)
        db.filter_similar = args.filter_similar
        db.similarity_threshold = args.similarity_threshold or 0.95

        genomes = db._read_fasta(args.genome)  
        if args.filter_similar:
            logging.info("Filtering similar genomes before building reference.")
            db.similarity_threshold = args.similarity_threshold or 0.95
            genomes = db.filter_similar_genomes(genomes)

        db.build_from_fasta(genomes)
        db.save(args.reference)
        print(f"Reference database saved to {args.reference}")
        logging.info(f"Reference database saved to {args.reference}")
        return
    
    if args.task == "dumpref":
        if args.reference:  
            logging.info(f"Dumping reference database from {args.reference}")
            check_file_exists(args.reference, "Reference database")

            db = ReferenceDatabase.load(args.reference)
            if db is None:
                raise ValueError("Error: Failed to load the reference database.")
            db.dump_json_stdout()

            return

        elif args.genome and args.kmer_size:
            logging.info(f"Building and dumping reference from {args.genome}")
            check_file_exists(args.genome, "Genome FASTA file")

            db = ReferenceDatabase(k=args.kmer_size)
            genomes = db._read_fasta(args.genome)

            if args.filter_similar:
                logging.info("Filtering similar genomes before building reference.")
                db.similarity_threshold = args.similarity_threshold or 0.95
                genomes = db.filter_similar_genomes(genomes)

            db.build_from_fasta(genomes)

            summary = db.dump_json_stdout()
            print(json.dumps(summary, indent=4))
            return


        else:
            raise ValueError("For dumping reference, provide either --reference or --genome.")
        




    if args.task == "align":
        if not all([args.reads, args.alignment]) or (not args.reference and not args.genome):
            raise ValueError("For alignment, --reads, --alignment, and one of --reference or --genome are required.")

        check_file_exists(args.reads, "Reads file")

        if args.reference:
            check_file_exists(args.reference, "Reference database")
            reference_db = ReferenceDatabase.load(args.reference)
            if reference_db is None:
                raise ValueError("Error loading reference database.")
        elif args.genome:
            check_file_exists(args.genome, "Reference genome")
            reference_db = ReferenceDatabase(k=args.kmer_size)
            genomes = reference_db._read_fasta(args.genome)
            reference_db.build_from_fasta(genomes)
        else:
            raise ValueError("Must provide either --reference or --genome.")

        aligner = PseudoAligner(
            reference_file=args.reference,
            reads_file=args.reads,
            output_file = args.alignment if args.alignment else "temp.json",
            min_read_quality=args.min_read_quality,
            min_kmer_quality=args.min_kmer_quality,
            max_genomes=args.max_genomes,
            reverse_complement=args.reverse_complement,
            coverage=args.coverage,
            genomes=args.genomes.split(",") if args.genomes else None,
            min_coverage=args.min_coverage,
            full_coverage=args.full_coverage,
            kmer_size=args.kmer_size,
            m=args.m,
            p=args.p
        )

        aligner.reference_db = reference_db

        aligner.run()
        print(f"Alignment results saved to {args.alignment}")
        return

    
    elif args.task == "dumpalign":
        if args.reference and args.reads:
            run_alignment_and_print_summary(
                reference_file=args.reference,
                reads_file=args.reads,
                m=args.m if args.m else 1,
                p=args.p if args.p else 1,
                kmer_size=args.kmer_size
            )
            return
        elif args.genome and args.kmer_size and args.reads:
            logging.info("Running full pipeline: build reference + align + dump summary")
            check_file_exists(args.genome, "Genome FASTA file")
            check_file_exists(args.reads, "Reads file")

            db = ReferenceDatabase(k=args.kmer_size)
            genomes = db._read_fasta(args.genome)
            db.build_from_fasta(genomes)

            aligner = PseudoAligner(
                reference_file=None,
                reads_file=args.reads,
                output_file="temp.json",  
                kmer_size=args.kmer_size,
                m=args.m if args.m else 1,
                p=args.p if args.p else 1
            )
            aligner.reference_db = db
            aligner.run()

            import sys
            PseudoAligner.dump_json("temp.json", sys.stdout)
            return
        if not args.alignment:
            raise ValueError("For dumping alignment results, --alignment is required.")

        logging.info(f"Dumping alignment results from {args.alignment}")
        check_file_exists(args.alignment, "Alignment results")

        output_path = f"{os.path.splitext(args.alignment)[0]}_summary.json"
        PseudoAligner.dump_json(args.alignment, output_path)
        print(f"Alignment summary saved to {output_path}")
        logging.info(f"Alignment summary saved to {output_path}")
        return



def run_alignment_and_print_summary(reference_file, reads_file, m=1, p=1, kmer_size=31):
    reference_db = ReferenceDatabase.load(reference_file)
    if not reference_db:
        raise ValueError("Failed to load reference database.")

    aligner = PseudoAligner(
        reference_file=reference_file,
        reads_file=reads_file,
        output_file=None,  
        kmer_size=kmer_size
    )
    
    aligner.m = m
    aligner.p = p

    reads = aligner._read_fastq(reads_file)
    for read, quality in reads:
        if aligner.min_read_quality and aligner._calculate_quality_mean(quality) < aligner.min_read_quality:
            aligner.unmapped_reads += 1
            continue
        aligner._find_matches(read, quality)

    summary = {
        "Statistics": {
            "unique_mapped_reads": sum(aligner.unique_counts.values()),
            "ambiguous_mapped_reads": sum(aligner.ambiguous_counts.values()),
            "unmapped_reads": aligner.unmapped_reads
        },
        "Summary": {}
    }

    all_genomes = set(aligner.unique_counts.keys()) | set(aligner.ambiguous_counts.keys())
    for genome in all_genomes:
        summary["Summary"][genome] = {
            "unique_reads": aligner.unique_counts.get(genome, 0),
            "ambiguous_reads": aligner.ambiguous_counts.get(genome, 0)
        }

    print(json.dumps(summary, indent=4))
    
if __name__ == "__main__":
        main()