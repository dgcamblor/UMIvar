#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# ██    ██ ███    ███ ██ ██    ██  █████  ██████  
# ██    ██ ████  ████ ██ ██    ██ ██   ██ ██   ██ 
# ██    ██ ██ ████ ██ ██ ██    ██ ███████ ██████  
# ██    ██ ██  ██  ██ ██  ██  ██  ██   ██ ██   ██ 
#  ██████  ██      ██ ██   ████   ██   ██ ██   ██ 
#-------------------------------------------------------------------------------
# This program introduces UMI-aware variants in a BAM file.
#-------------------------------------------------------------------------------

__author__ = "@dgcamblor"
__version__ = "0.9"


from subprocess import run
import argparse
import os
from timeit import default_timer as timer
import json
import sys

import pysam
import numpy as np

# Add src to the path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

from src.constants import *
from src.checkups import check_dependencies, check_dirs
from src.io_functions import create_ontarget, extract_variants, write_reads
from src.adding_variants import add_variants
from src.misc import ontarget_hash 

# Bash scripts
TO_FASTQ = "./src/to_fastq.sh"


def parse_args():
    parser = argparse.ArgumentParser(description="Introduce UMI-aware variants in a BAM file.")
    
    parser.add_argument("-i", "--input_bam", help="BAM file to introduce variants in.", required=True)
    parser.add_argument("-v", "--variants", help="CSV file with the variants to introduce.", required=True)

    parser.add_argument("-o", "--output_bam", help="Output BAM file.", required=False)
    parser.add_argument("-b", "--bed", help="BED file with the target regions to delimit the file.", required=False)
    parser.add_argument("-f", "--freq_mode", help="Mode to calculate the allelic frequency of the variants. Options: dedup | umi", default="dedup")
    parser.add_argument("-s", "--seed", help="Seed for the random number generators.", type=int, default=12)
    parser.add_argument("-c", "--sb_correction", help="Correct for strand bias.", action="store_true")

    return parser.parse_args()


def main():
    print(
    f"""
    {PURPLE}██    ██ ███    ███ ██ {END}██    ██  █████  ██████  
    {PURPLE}██    ██ ████  ████ ██ {END}██    ██ ██   ██ ██   ██ 
    {PURPLE}██    ██ ██ ████ ██ ██ {END}██    ██ ███████ ██████  
    {PURPLE}██    ██ ██  ██  ██ ██ {END} ██  ██  ██   ██ ██   ██ 
    {PURPLE} ██████  ██      ██ ██ {END}  ████   ██   ██ ██   ██ 

    Ver. {__version__}                                             
    By: {GREEN}{__author__}{END}                               
    """
    )

    np.random.seed(args.seed)

    stats = {
        "written_reads": 0,
        "variant_reads": 0,
        "added_variants": 0,
        "time": 0
    }

    t0 = timer()

    input_name = args.input_bam.split("/")[-1].split(".")[0]

    if not args.output_bam:
        output_path = "results"
        output_name = f"{input_name}_UV"
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        output_bam = f"{output_path}/{output_name}.bam"
    else:
        output_path = "/".join(args.output_bam.split("/")[:-1])
        output_name = args.output_bam.split("/")[-1].split(".")[0]
        output_bam = args.output_bam

    output_var = output_bam.replace(".bam", ".csv")

    print(f"{BLUE}Sample:{END} {input_name}")
    print(f"{BLUE}Output type:{END} {'on target' if args.bed else 'whole bam'}")
    print(f"{BLUE}Output path:{END} {output_bam}")
    print()

    check_dependencies()
    check_dirs()

    step = 1

    # Create on-target BAM file only if the user wants to
    if args.bed:
        print(f"{BLUE}[{step}] Creating on-target BAM file...{END}"); step += 1

        if os.path.exists("tmp/ontarget_hashes.json"):
            with open("tmp/ontarget_hashes.json", "r") as hash_file:
                try:
                    hashes = json.load(hash_file)
                except json.decoder.JSONDecodeError:
                    print(f"{YELLOW}WARNING: ontarget_hashes.json is corrupted, creating new file.{END}")
                    hashes = {}
        else:
            hashes = {}
        
        # If the ontarget file does not exist, create it
        if not os.path.exists(f"tmp/{input_name}_ontarget.sorted.bam"):
            create_ontarget(args.bed, args.input_bam)

            # Create a hash from the original BAM file and the BED to check if 
            # the ontarget file is up to date
            with open(f"tmp/ontarget_hashes.json", "w") as hash_file:
                hashes[input_name] = ontarget_hash(args.bed, args.input_bam)
                json.dump(hashes, hash_file)

        # If the ontarget file exists, check if it is up to date
        else:
            # The ontarget file is found but not the hash
            if input_name not in hashes:
                create_ontarget(args.bed, args.input_bam)

                with open(f"tmp/ontarget_hashes.json", "w") as hash_file:
                    hashes[input_name] = ontarget_hash(args.bed, args.input_bam)
                    json.dump(hashes, hash_file)

            # If the ontarget file is not up to date, create it again
            elif hashes[input_name] != ontarget_hash(args.bed, args.input_bam):
                create_ontarget(args.bed, args.input_bam)

                # Update the hash
                with open(f"tmp/ontarget_hashes.json", "w") as hash_file:
                    hashes[input_name] = ontarget_hash(args.bed, args.input_bam)
                    json.dump(hashes, hash_file)

            elif hashes[input_name] == ontarget_hash(args.bed, args.input_bam):
                print(f"{GREEN}On-target BAM file found and up to date.{END}")

        print()

    input_bam = f"tmp/{input_name}_ontarget.sorted.bam" if args.bed else args.input_bam

    with (
        pysam.AlignmentFile(input_bam, "rb") as in_bamfile,
        pysam.AlignmentFile(output_bam, "wb", template=in_bamfile) as out_bamfile,
        open(args.variants, "r") as in_varfile,
        open(output_var, "w") as out_varfile):

        out_varfile.write("chr,pos,ref,alt,af,umis,sb\n")

        # Obtain the variants from the CSV file 
        variants = extract_variants(in_varfile)

        # For each variant, introduce it in its covering reads from the BAM file
        print(f"{BLUE}[{step}] Introducing variants...{END}"); step += 1

        correct_sb = True if args.sb_correction else False
        record = add_variants(in_bamfile, variants, out_varfile, args.freq_mode, correct_sb, stats)

        print()

        # Write the reads to the output file   
        print(f"{BLUE}[{step}] Writing reads to output file...{END}"); step += 1

        write_reads(in_bamfile, out_bamfile, record, stats)

    print(); print()

    # Sort and index the output file
    print(f"{BLUE}[{step}] Sorting and indexing output file...{END}"); step += 1
    print()

    pysam.sort("-o", f"{output_path}/{output_name}.sorted.bam", output_bam)
    pysam.index(f"{output_path}/{output_name}.sorted.bam")

    # Produce FASTQ files for the reads
    print(f"{BLUE}[{step}] Producing FASTQ files...{END}"); step += 1
    print()

    run(f"{TO_FASTQ} -i {output_path}/{output_name}.sorted.bam", shell=True)

    t1 = timer()

    stats["time"] = t1 - t0

    # Print stats
    print()
    print(BLUE, end="")
    print("###############################")
    print("############ STATS ############")
    print("###############################")
    print(END, end="")
    print(f"{BLUE}- Total reads written:{END} {stats['written_reads']}")
    print(f"{BLUE}- Total reads with variants written:{END} {stats['variant_reads']}")
    all_vars = stats["added_variants"] == len(variants)
    print(f"{BLUE}- Number of variants written:{END} {YELLOW if not all_vars else END}{stats['added_variants']}{END} / {len(variants)}")
    print(f"{BLUE}- Time elapsed:{END} {stats['time']:.2f} s")
    print(f"{BLUE}Results stored in:{END} {output_path}")
    print(END, end="")
    print()


if __name__ == "__main__":
    args = parse_args()
    main()