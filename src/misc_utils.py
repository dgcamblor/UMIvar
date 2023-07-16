# misc_utils.py | Contains miscellaneous functions for the main program

import os
from subprocess import run

from src.constants import *


def ontarget_hash(bam, bed):
    """Create a simples hash from the original BAM file and the BED to 
    check if the on-target BAM file is up to date and avoid re-creating it.
    
    Args:
        bam (str): Path to the original BAM file
        bed (str): Path to the BED file
    
    Returns:
        str: Hash
    """
    bam_size = os.path.getsize(bam)
    bed_size = os.path.getsize(bed)

    bam_creation = os.path.getctime(bam)
    bed_creation = os.path.getctime(bed)

    bam_last_mod = os.path.getmtime(bam)
    bed_last_mod = os.path.getmtime(bed)

    return str(hash(bam_size + bed_size + bam_creation + bed_creation + bam_last_mod + bed_last_mod))


def handle_groups(in_bam_name, in_bam_path, edit_threshold, seed):
    """Handles the creation of the UMI-tools groups file. If the file does not exist, create it.

    Args:
        in_bam_name (str): Name of the input BAM file.
        in_bam_path (str): Path to the input BAM file.
        edit_threshold (int): Edit distance threshold for UMI-tools group.
        seed (int): Seed for the random number generator.

    Returns:
        dict: Dictionary with the UMI-tools groups.
    """
    groups_file = f"tmp/{in_bam_name}_groups.tsv"

    run(f"umi_tools group -I {in_bam_path} --group-out {groups_file} --edit-distance-threshold {edit_threshold} --random-seed {seed}", shell=True)

    with open(groups_file, "r") as gf:
        gf.readline()  # Skip header
        umi_groups = {line.split("\t")[4]: line.split("\t")[6] for line in gf}
        
    return umi_groups