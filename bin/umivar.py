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
__version__ = "1.0"

import sys, os
from subprocess import run
import argparse
from timeit import default_timer as timer

import pysam
import numpy as np

# Add src to the path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

from src.adding_variants import add_variants
from src.check_utils import check_dependencies, check_dirs
from src.constants import *
from src.io_functions import extract_variants, write_reads, handle_ontarget, print_stats, get_outputs
from src.misc_utils import handle_groups
from src.random_mode import random_variants


def parse_args():
    parser = argparse.ArgumentParser(description="Introduce UMI-aware variants in a BAM file.")

    # Required arguments    
    parser.add_argument("-i", "--input_bam", help="BAM file to introduce variants in.", required=True)
    parser.add_argument("-v", "--variants", help="CSV/VCF file with the variants to introduce. If a number is provided, that number of random variants will be simulated (requires -f).", required=True)

    # Optional arguments
    parser.add_argument("-o", "--output_bam", help="Output BAM file.", required=False)
    parser.add_argument("-b", "--bed", help="BED file with the target regions to delimit the file.", required=False)
    parser.add_argument("-f", "--ref_fasta", help="Reference FASTA file. Optional: only needed in the random generation of variants (numeric input in -v).", required=False)
    parser.add_argument("-e", "--edit_threshold", help="Edit distance threshold for the UMI-tools grouping algorithm.", type=int, required=False)
    parser.add_argument("-s", "--seed", help="Seed for the random number generators.", type=int, default=12)

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

    random_mode = True if args.variants.isdigit() else False  # Random variants mode

    if random_mode and not args.ref_fasta:
        print(f"{RED}ERROR:{END} Reference FASTA file is required in random mode. Use -f.")
        sys.exit(1)
    elif not random_mode and args.ref_fasta:
        print(f"{YELLOW}WARNING:{END} Reference FASTA file is not required in non-random mode. Ignoring -f."); print()

    # Check dependencies
    extra_depends = []
    if args.edit_threshold:
        extra_depends.append("umi_tools")

    check_dependencies(extra_depends)
    check_dirs()

    np.random.seed(args.seed)  # Fix the seed for the random number generators
    
    # Stats for the run
    stats = {
        "written_reads": 0,
        "variant_reads": 0,
        "added_variants": 0,
        "time": 0,
        "total_vars": 0,  # Total number of variants in the CSV/VCF file
        "all_vars": False  # Whether all the variants were added or not
    }

    t0 = timer()

    in_bam_name = args.input_bam.split("/")[-1].split(".")[0]

    # Handle output files
    out_path, out_bam_name, out_bam_path = get_outputs(in_bam_name)
    out_var_path = out_bam_path.replace(".bam", ".csv")  # Output CSV file with the introduced variants

    print(f"{BLUE}Sample:{END} {in_bam_name}")
    print(f"{BLUE}Variants:{END} {args.variants} {'(random)' if random_mode else ''}")
    print(f"{BLUE}Output type:{END} {'on target' if args.bed else 'whole BAM'}")
    print(f"{BLUE}Output path:{END} {out_bam_path}")
    print()

    step = 1

    # Create on-target BAM file only if the user wants to
    if args.bed:
        print(f"{BLUE}[{step}] Creating on-target BAM file...{END}"); step += 1
        handle_ontarget(in_bam_name, args.input_bam, args.bed)
        print()

    # The input BAM file will be the on-target BAM file if the user wants to
    in_bam_path = f"tmp/{in_bam_name}_ontarget.sorted.bam" if args.bed else args.input_bam

    # Obtain umi-tools groups
    if args.edit_threshold:
        print(f"{BLUE}[{step}] Obtaining UMI-tools groups...{END}"); step += 1

        umi_groups = handle_groups(in_bam_name, in_bam_path, args.edit_threshold, args.seed)

        print()
    else:
        umi_groups = None

    # Obtain the variants to introduce
    if random_mode:
        print(f"{BLUE}[{step}] Obtaining random variants...{END}"); step += 1
        variants = random_variants(in_bam_path, int(args.variants), args.ref_fasta)
        print()
    else:
        with open(args.variants, "r") as in_var_file:
            variants = extract_variants(in_var_file)

    # Variant simulation
    with (
        pysam.AlignmentFile(in_bam_path, "rb") as in_bam_file,
        pysam.AlignmentFile(out_bam_path, "wb", template=in_bam_file) as out_bam_file,
        open(out_var_path, "w") as out_var_file):

        out_var_file.write("chr,pos,ref,alt,af,cov,umis,sb\n")

        # For each variant, introduce it in its covering reads from the BAM file
        print(f"{BLUE}[{step}] Introducing variants...{END}"); step += 1

        record = add_variants(in_bam_file, variants, out_var_file, stats, umi_groups)

        print()

        # Write the reads to the output file   
        print(f"{BLUE}[{step}] Writing reads to output file...{END}"); step += 1

        write_reads(in_bam_file, out_bam_file, record, stats)

    print(); print()

    # Sort and index the output file
    print(f"{BLUE}[{step}] Sorting and indexing output file...{END}"); step += 1
    print()

    pysam.sort("-o", f"{out_path}/{out_bam_name}.sorted.bam", out_bam_path)
    pysam.index(f"{out_path}/{out_bam_name}.sorted.bam")

    # Produce FASTQ files for the reads
    print(f"{BLUE}[{step}] Producing FASTQ files...{END}"); step += 1
    print()

    run(f"{TO_FASTQ} -i {out_path}/{out_bam_name}.sorted.bam", shell=True)

    t1 = timer()

    stats["time"] = t1 - t0
    stats["total_vars"] = len(variants)
    stats["all_vars"] = stats["added_variants"] == stats["total_vars"]

    # Print stats
    print_stats(stats, out_path)


if __name__ == "__main__":
    args = parse_args()
    main()