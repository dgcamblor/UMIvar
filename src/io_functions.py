# io_functions.py | Functions for reading and writing in umivar

import sys, os
import json
from subprocess import run

from src.misc_utils import ontarget_hash
from src.variant import Variant
from src.adding_variants import get_read_id
from multiprocessing import cpu_count
from src.constants import *


def get_outputs(bam_name, output_bam=None):
    """Get the output path and name from the input BAM file.

    Args:
        bam_name (str): Name of the input BAM file.
        output_bam (str): Path to the output BAM file.

    Returns:
        tuple: Tuple with the output path, output name and output BAM file.
    """
    if not output_bam:  # Default output path
        output_path = "results"
        output_name = f"{bam_name}_UV"
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        output_bam = f"{output_path}/{output_name}.bam"
    else:  # Custom output path
        output_path = "/".join(output_bam.split("/")[:-1])
        output_name = output_bam.split("/")[-1].split(".")[0]
        output_bam = output_bam

    return output_path, output_name, output_bam


def create_ontarget(bed_file, bam_file):
    """Create a BAM file with only the reads that overlap the target regions.
    Shortens computation time and memory usage.
    
    Args:
        bed_file (file): path to BED file with the target regions.
        bam_file (file): path to BAM file with the reads.
    """
    in_bam_name = bam_file.split("/")[-1].split(".")[0]

    run(f"samtools view -b -L {bed_file} {bam_file} > tmp/{in_bam_name}_ontarget.bam", shell=True)
    run(f"samtools sort tmp/{in_bam_name}_ontarget.bam -o tmp/{in_bam_name}_ontarget.sorted.bam", shell=True)
    run(f"samtools index tmp/{in_bam_name}_ontarget.sorted.bam", shell=True)


def extract_variants(var_file):
    """Extract variants from a CSV/VCF file.
    
    Args:
        var_file (file): File with the variants. If CSV, each line must have the following format:
            chr,pos1,ref,alt,af
        If VCF, the af field must be in the INFO column, under the AF tag.
    Returns:
        list: List of Variant objects
    """
    variants = []

    if var_file.name.endswith(".csv"):
        for line in var_file:
            fields = line.rstrip().split(",")
            if not fields[1].isnumeric():  # Check if it is a header (the pos field is not a number)
                continue
            variants.append(Variant(chr = fields[0], 
                                    pos1 = int(fields[1]),
                                    ref = fields[2], 
                                    alt = fields[3], 
                                    af = float(fields[4])))
            
    elif var_file.name.endswith(".vcf"):
        for line in var_file:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            af = float(fields[7].split("AF=")[1].split(";")[0])
            variants.append(Variant(chr = fields[0], 
                                    pos1 = int(fields[1]),
                                    ref = fields[3], 
                                    alt = fields[4], 
                                    af = af))
            
    else:
        print(f"{RED}ERROR: Variants file must be a CSV or VCF file.{END}")
        sys.exit(1)

    return variants


def write_reads(in_bam_file, out_bam_file, record, stats):
    """Write the reads to the output file. If a read has a variant, write the mutated read instead.

    Args:
        in_bam_file (pysam.AlignmentFile): Input BAM file.
        out_bam_file (pysam.AlignmentFile): Output BAM file.
        record (dict): Dictionary with the reads that have a variant.
    Returns:
        dict: Dictionary with the number of reads and variants written.
    """
    cpus = cpu_count()  # Shortens computation time of counting reads
    total = int(run(f"samtools view -@ {cpus} -c {in_bam_file.filename.decode('utf-8')}", shell=True, capture_output=True).stdout.decode("utf-8").rstrip())
    written_reads, variant_reads = 0, 0
    
    for read in in_bam_file:
        if read.reference_name in record and read.reference_start in record[read.reference_name]:
            read_id = get_read_id(read)
            if read_id not in record[read.reference_name][read.reference_start]:
                out_bam_file.write(read)
                written_reads += 1
            else:
                mut_read = record[read.reference_name][read.reference_start][read_id]
                out_bam_file.write(mut_read)
                written_reads += 1
                variant_reads += 1
        else:
            out_bam_file.write(read)
            written_reads += 1

        print(f"Progress: {written_reads}/{total} ({written_reads/total*100:.2f}%)", end="\r")

    stats["written_reads"] = written_reads
    stats["variant_reads"] = variant_reads


def handle_ontarget(in_bam_name, in_bam_path, bed):
    """Handles the creation of the on-target BAM file. If the file does not exist, create it. 
    If it exists, check if it is up to date.

    Args:
        in_bam_name (str): Name of the input BAM file.
        in_bam_path (str): Path to the input BAM file.
        bed (str): Path to the BED file.
    """
    # Check if the hash list file exists, if so, load it
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
    if not os.path.exists(f"tmp/{in_bam_name}_ontarget.sorted.bam"):
        create_ontarget(bed, in_bam_path)

        # Create a hash from the original BAM file and the BED to check if 
        # the ontarget file is up to date
        with open(f"tmp/ontarget_hashes.json", "w") as hash_file:
            hashes[in_bam_name] = ontarget_hash(bed, in_bam_path)
            json.dump(hashes, hash_file)

    # If the ontarget file exists, check if it is up to date
    else:
        # The ontarget file is found but not the hash
        if in_bam_name not in hashes:
            create_ontarget(bed, in_bam_path)

            with open(f"tmp/ontarget_hashes.json", "w") as hash_file:
                hashes[in_bam_name] = ontarget_hash(bed, in_bam_path)
                json.dump(hashes, hash_file)

            # If the ontarget file is not up to date, create it again
        elif hashes[in_bam_name] != ontarget_hash(bed, in_bam_path):
            create_ontarget(bed, in_bam_path)

                # Update the hash
            with open(f"tmp/ontarget_hashes.json", "w") as hash_file:
                hashes[in_bam_name] = ontarget_hash(bed, in_bam_path)
                json.dump(hashes, hash_file)

        elif hashes[in_bam_name] == ontarget_hash(bed, in_bam_path):
            print(f"{GREEN}On-target BAM file found and up to date.{END}")


def print_stats(stats, out_path):
    """Print the statistics of the run."""
    print()
    print(BLUE, end="")
    print("###############################")
    print("############ STATS ############")
    print("###############################")
    print(END, end="")
    print(f"{BLUE}- Total reads written:{END} {stats['written_reads']}")
    print(f"{BLUE}- Total reads with variants written:{END} {stats['variant_reads']}")
    print(f"{BLUE}- Number of variants written:{END} {YELLOW if not stats['all_vars'] else END}{stats['added_variants']}{END} / {stats['total_vars']}")
    print(f"{BLUE}- Time elapsed:{END} {stats['time']:.2f} s")
    print(f"{BLUE}Results stored in:{END} {out_path}")
    print(END, end="")
    print()