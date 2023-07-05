# random_mode.py | Functions for the random variants mode

import pysam
from pyfaidx import Fasta

from src.constants import *
from src.variant import Variant


def get_positions(bam_path, n_pos, ref_fasta):
    """Obtain a set of positions from a BAM file where the variants will be introduced.
    These positions are selected among the positions with the highest coverage in the BAM file,
    to ensure the correct introduction of the variants.


    Args:
        bam_path (str): Path to BAM file
        n_pos (int): Number of positions to select
        ref_fasta (Fasta): Reference FASTA file opened as a pyfaidx Fasta
    Returns:
        List: list of 0-based positions where the variants are to be introduced
    """
    all_pos = []

    # Iterate over the pileup columns (positions)
    with pysam.AlignmentFile(bam_path, "r") as bam_file:
        for pileupcolumn in bam_file.pileup():
            nucleotides = [nucleotide.upper() for nucleotide in pileupcolumn.get_query_sequences(add_indels=True)]

            # Only consider positions with all the nucleotides being the reference nucleotide
            ref_nucleotide = ref_fasta[pileupcolumn.reference_name][pileupcolumn.reference_pos].seq.upper()
            if all(nucleotide == ref_nucleotide for nucleotide in nucleotides):
                all_pos.append((pileupcolumn.reference_name, pileupcolumn.reference_pos, pileupcolumn.n, ref_nucleotide))

    # Sort by the coverage of the position (third field)
    all_pos.sort(key=lambda x: x[2], reverse=True)

    # Select the positions
    selected_pos = []
    for pos in all_pos:
        if len(selected_pos) >= n_pos + EXTRA_POS:
            break
        else:
            # Check if the position is not too close to another chosen position
            if all(abs(pos[1] - x[1]) > POS_DIFF for x in selected_pos):
                selected_pos.append(pos)

    # Randomly choose n_pos positions
    random_indexes = np.random.choice(len(selected_pos), n_pos, replace=False)
    selected_pos = [selected_pos[i] for i in random_indexes]

    return selected_pos


def random_variants(bam_path, n_pos, ref_fasta_path):
    """Generate a set of random variants to be introduced in the BAM file.

    Args:
        bam_path (str): Path to BAM file
        n_pos (int): Number of positions to select
        ref_fasta_path (str): Path to the reference FASTA file
    Returns:
        list: List of Variant objects
    """
    variants = []

    ref_fasta = Fasta(ref_fasta_path)

    positions = get_positions(bam_path, n_pos, ref_fasta)  # 0-based positions

    for position in positions:
        chrom = position[0]
        pos1 = position[1] + 1  # 1-based position (input for Variant object)
        ref = position[3]

        # Randomly choose the type of variant
        var_type = np.random.choice(VAR_TYPES)

        if var_type == "SNV":
            alt = np.random.choice([n for n in NUCLEOTIDES if n != ref])
            
        elif var_type == "INS":
            indel_length = np.random.choice(INDEL_LENGTHS)
            alt = ref + "".join(np.random.choice(NUCLEOTIDES) for _ in range(indel_length))

        elif var_type == "DEL":
            indel_length = np.random.choice(INDEL_LENGTHS)

            # Access the reference FASTA file to obtain the deleted sequence
            ref = ref_fasta[chrom][pos1 - 1:pos1 + indel_length].seq.upper()
            alt = ref[0]

        af = np.random.choice(ALLELIC_FREQS)
        variants.append(Variant(chrom, pos1, ref, alt, af))

    variants.sort()

    """
    # Writing the random target variants to a file
    with open("tmp/target_variants.csv", "w") as tv:
        tv.write("chr,pos,ref,alt,af\n")
        for variant in variants:
            tv.write(f"{variant.chr},{variant.pos1},{variant.ref},{variant.alt},{variant.af}\n")
    """

    return variants