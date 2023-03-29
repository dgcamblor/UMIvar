# misc.py | Contains miscellaneous functions for the main program

import os


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