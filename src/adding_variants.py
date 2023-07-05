# adding_variants.py | Functions for inserting variants into a bam file

from collections import defaultdict
import itertools
import array
import re
import numpy as np

from src.constants import *


def rec_defaultdict():
    """Allows for a nested defaultdict using recursion."""
    return defaultdict(rec_defaultdict)


def unique_umis(umis):
    """Get the unique UMIs in order of appearance. Preserving the order allows
    the seed to remain useful for reproducibility."""
    seen = set()
    return [umi for umi in umis if not (umi in seen or seen.add(umi))]


def get_umi_from_qname(read_qname):
    """Get the UMI from the read query name.
    
    Args:
        read_qname (str): Read qname. The UMI must be the last element of the qname, 
        separated by a non-word character (commonly '_' or ':').

    Returns:
        str: UMI
    """
    return re.split("[^A-Z]", read_qname)[-1]


def get_read_id(read):
    """Get a unique identifier for a read, using all of its information"""
    return hash(f"{read}")


'''
def umi_freqs(umis):
    """Get the frequencies of the UMIs.

    Args:
        umis (list): List of UMIs covering a corresponding variant

    Returns:
        tuple: Tuple with the following elements:
            umi_absfreqs (dict): Dictionary with the absolute frequencies of the UMIs
            umi_relfreqs (dict): Dictionary with the relative frequencies of the UMIs
            freq_counts (dict): Dictionary with the count of each frequency appearing in the UMIs
    """
    # Get the count of reads with the same UMI (absolute frequencies)
    umi_absfreqs = defaultdict(int)
    for umi in umis:
        umi_absfreqs[umi] += 1

    # Get the relative frequencies of the UMIs
    umi_relfreqs = defaultdict(float)
    total_count = sum(umi_absfreqs.values())
    for umi in umi_absfreqs:
        umi_relfreqs[umi] = umi_absfreqs[umi] / total_count

    return umi_absfreqs, umi_relfreqs
'''


def umis_to_groups(umis, umi_groups):
    """Convert a list of UMIs to their corresponding groups.

    Args:
        umis (list): List of UMIs
        umi_groups (dict): Dictionary containing the UMIs (keys) and their corresponding groups (values).

    Returns:
        list: List of groups of UMIs
    """
    return unique_umis([umi_groups[umi] for umi in umis if umi in umi_groups])  #NOTE: Only use the UMIs that are in the dictionary


def groups_to_umis(groups, umi_groups):
    """For each UMI in groups, gets each UMI (key) that belongs to the group (value) 
    in the umi_groups dictionary and adds it to the list. 

    Args:
        groups (list): List of groups of UMIs
        umi_groups (dict): Dictionary containing the UMIs (keys) and their corresponding groups (values).
    
    Returns:
        list: List of UMIs
    """
    return unique_umis([umi for group in groups for umi in umi_groups if umi_groups[umi] == group])


def get_umi_combination(variant, umis, umis_fwd, umis_rev, umi_blacklist, umi_groups=None):
    """Get a combination of UMIs to cover a variant at the specified allelic frequency.
    
    Args:
        variant (Variant): Variant to cover
        umis (dict): Dictionary with the relative frequencies of the UMIs
        umis_fwd (list): List of UMIs covering the variant in the forward strand
        umis_rev (list): List of UMIs covering the variant in the reverse strand
        umi_blacklist (list): List of UMIs to ignore
        umi_groups (dict): Dictionary containing the groups of UMIs. If None, the UMIs will not be grouped.

    Returns:
        tuple: Tuple with the following elements:
            achieved_af (float): Best allelic frequency
            umi_combination (list): List of UMIs to use
    """
    if umi_groups:
        # Convert each UMI in umis, umis_fwd, umis_rev and umi_blacklist to its group
        umis = umis_to_groups(umis, umi_groups)
        umis_fwd = umis_to_groups(umis_fwd, umi_groups)
        umis_rev = umis_to_groups(umis_rev, umi_groups)
        umi_blacklist = umis_to_groups(umi_blacklist, umi_groups)

    if variant.af == 0:
        achieved_af = 0
        umi_combination = []
    
    elif variant.af == 1:
        umi_combination = [umi for umi in umis if umi not in umi_blacklist]
        achieved_af = 1 if not umi_blacklist else sum([umis[umi] for umi in umi_combination])
    
    else:
        ind_freq = 1/len(umis)
        achieved_af = 0
        available_umis = [umi for umi in umis if umi not in umi_blacklist]

        for i in range(1, len(available_umis)):
            achieved_af += ind_freq
            # If the current AF is greater than the target AF, evaluate it and the previous one
            # The closest will be selected
            if achieved_af >= variant.af:
                diff_curr = abs(achieved_af - variant.af)
                diff_prev = abs((achieved_af - ind_freq) - variant.af)
                if diff_curr < diff_prev: 
                    n_umis = i
                elif diff_curr > diff_prev:
                    n_umis = i - 1
            
                break
        # If the for loop finishes, use all the UMIs
        else:
            n_umis = len(available_umis)

        # Get a combination of UMIs that achieves the best target AF
        umi_combination = choose_umis(umis_fwd, umis_rev, n_umis, umi_blacklist)

    if umi_groups:  # Convert the UMIs back to their original form
        umi_combination = groups_to_umis(umi_combination, umi_groups)

    return umi_combination, achieved_af


def choose_umis(umis_fwd, umis_rev, n_umis, umi_blacklist):
    """Choose UMIs to use for the variant based on the UMIs alone. 

    Args:
        umis_fwd (list): List of UMIs that are on the forward strand
        umis_rev (list): List of UMIs that are on the reverse strand
        n_umis (int): Number of UMIs to use
        umi_blacklist (list): List of UMIs to ignore

    Returns:
        umi_combination (list): List of frequencies of the UMIs to use
    """
    possible_fwd = [umi for umi in umis_fwd if umi not in umi_blacklist]
    possible_rev = [umi for umi in umis_rev if umi not in umi_blacklist]
    umi_combination = []

    # Choose randomly between the forward and reverse strands until the number of UMIs is reached
    for _ in range(n_umis):
        if len(possible_fwd) > len(possible_rev):
            fwd_umi = np.random.choice(possible_fwd)
            umi_combination.append(fwd_umi)
            possible_fwd.remove(fwd_umi)
            if fwd_umi in possible_rev:
                possible_rev.remove(fwd_umi)

        elif len(possible_fwd) < len(possible_rev):
            rev_umi = np.random.choice(possible_rev)
            umi_combination.append(rev_umi)
            possible_rev.remove(rev_umi)
            if rev_umi in possible_fwd:
                possible_fwd.remove(rev_umi)

        else:
            umi = np.random.choice(unique_umis(possible_fwd + possible_rev))
            umi_combination.append(umi)
            if umi in possible_fwd:
                possible_fwd.remove(umi)
            if umi in possible_rev:
                possible_rev.remove(umi)
                
    return umi_combination


def add_variants(bam_file, variants, var_file, stats, umi_groups=None):
    """Adds variants in the reads of a BAM file. The variants are introduced in the reads
    covering the variant position. The reads are mutated in such a way that the allele frequency of the
    variant is as close as possible to the specified allele frequency.

    Variants are not introduced directly in the bam_file, but, instead, mutated reads are recorded
    in a dictionary. The dictionary is then used to create a new BAM file with the mutated reads.

    Args:
        bam_file (pysam.AlignmentFile): BAM file containing the reads to mutate
        variants (list): List of variants to introduce
        var_file (str): Path to the file where the mutated reads will be stored
        stats (dict): Dictionary containing the statistics of the program.
        umi_groups (dict): Dictionary containing the groups of UMIs. If None, the UMIs will not be grouped.

    Returns:
        dict: A dictionary containing, for each chromosome, the IDs of the reads that have been mutated and 
        the corresponding reads and their information.
    """
    var_record = rec_defaultdict()
    indel_record = rec_defaultdict()

    added_vars = 0

    for variant in variants:
        try:  # Failsafe in case the variant cannot be introduced
            # Get all the reads covering the variant
            covering_reads = [read for read in bam_file.fetch(variant.chr, variant.pos0, variant.pos0 + 1)]
            if not covering_reads:
                print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [no coverture]{END}")
                continue

            # Separate the reads by strand
            reads_fwd = [read for read in covering_reads if not read.is_reverse]
            reads_rev = [read for read in covering_reads if read.is_reverse]

            # Get the UMIs from the IDs of each read covering the variant
            umis = unique_umis([get_umi_from_qname(read.query_name) for read in covering_reads])
            umis_fwd = unique_umis([get_umi_from_qname(read.query_name) for read in reads_fwd])
            umis_rev = unique_umis([get_umi_from_qname(read.query_name) for read in reads_rev])

            var_coverage = len(umis)

            # Check the validity of each read
            read_blacklist, umi_blacklist = check_read_validity(var_record, variant, covering_reads)

            if len(covering_reads) == len(read_blacklist):
                print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [no valid reads]{END}")
                continue

            # Get the UMIs to be mutated
            umis_to_mutate, achieved_af = get_umi_combination(variant, umis, umis_fwd, umis_rev, umi_blacklist, umi_groups)

            if not umis_to_mutate and umis:
                print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [too low AF]{END}")
                continue
            elif not umis_to_mutate and not umis:
                print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [no reads]{END}")
                continue

            reads_to_mutate = [read for read in covering_reads if get_umi_from_qname(read.query_name) in umis_to_mutate]

            # Mutate the reads
            mut_strands = []  # Record the strands of the mutated reads 

            for read in reads_to_mutate:
                read_id = get_read_id(read)

                # If the read has already been previously mutated, get the mutated read to mutate it again
                if read_id in var_record[variant.chr][read.reference_start]:
                    read = var_record[variant.chr][read.reference_start][read_id]

                # Get the sequence of the read
                seq = read.query_sequence
                quals = read.query_qualities

                # Get the position of the variant in the read
                seq_pos = read.get_reference_positions().index(variant.pos0)

                if read_id in indel_record[variant.chr][read.reference_start] or check_indel_cigar(read.cigarstring):
                    cigar_read = indel_record[variant.chr][read.reference_start][read_id] if read_id in indel_record[variant.chr][read.reference_start] else read
                    seq_pos, cigar_pos = change_pos(seq_pos, cigar_read)
                else:
                    cigar_pos = seq_pos

                # Mutate the base
                if variant.type == "SNV" or variant.type == "INS":
                    seq = seq[:seq_pos] + variant.alt + seq[seq_pos + 1:]
                elif variant.type == "MNV":
                    seq = seq[:seq_pos] + variant.alt + seq[seq_pos + variant.length:]
                elif variant.type == "DEL":
                    seq = seq[:seq_pos+1] + seq[seq_pos + variant.length + 1:]

                # Change qualities for indels
                if variant.type == "INS":
                    new_quals = [np.random.choice(INS_QUALITIES) for _ in range(len(variant.alt) - 1)]
                    new_quals = array.array('B', new_quals)
                    quals = quals[:seq_pos + 1] + new_quals + quals[seq_pos + 1:]
                elif variant.type == "DEL":
                    quals = quals[:seq_pos + 1] + quals[seq_pos + variant.length + 1:]

                # Change the CIGAR string for indels (SNVs are still a match)
                if variant.type == "INS" or variant.type == "DEL":
                    read.cigarstring = change_cigar(read.cigarstring, variant, cigar_pos)  

                # NOTE: No need to modify the length of the read, as it is automatically
                # changed when the sequence is changed

                # Update the sequence and the qualities read
                read.query_sequence = seq
                read.query_qualities = quals

                # Check if the read is valid after the mutation
                if not valid_sim(read):
                    print(f"{YELLOW}WARNING: Read {read.query_name} has not been mutated correctly — skipped. {END}")

                    # Delete the read from the record if it has been mutated
                    if read_id in var_record[variant.chr][read.reference_start]:
                        del var_record[variant.chr][read.reference_start][read_id]

                    continue

                # Add the read to the record
                var_record[variant.chr][read.reference_start][read_id] = read

                # Add the read to the indel record if it contains an indel
                if check_indel_cigar(read.cigarstring):
                    indel_record[variant.chr][read.reference_start][read_id] = read

                # Get the strand of the read
                mut_strands += ["-"] if read.is_reverse else ["+"]

                bam_file.reset()

        except Exception as ex:
            print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [internal error: {ex.__class__.__name__}]{END}")
            continue

        added_vars += 1

        # Calculate strand bias
        b = mut_strands.count("+")
        d = mut_strands.count("-")
        a = len(reads_fwd) - b
        c = len(reads_rev) - d

        strand_bias = calculate_sb(a, b, c, d)

        bam_file.reset()  # Reset the BAM file to the beginning

        if added_vars <= 7:
            print(f"{GREEN}Variant {variant} introduced at frequency {achieved_af:3f} in {len(umis_to_mutate)} UMIs (SB: {strand_bias:3f}){END}")
        elif added_vars == 8:
            print(f"{GREEN}More variants introduced... (see {var_file.name}){END}")
        
        var_file.write(f"{variant.chr},{variant.pos1},{variant.ref},{variant.alt},{achieved_af},{var_coverage},{len(umis_to_mutate)},{strand_bias}\n")

    stats["added_variants"] = added_vars

    return var_record


def check_read_validity(var_record, variant, covering_reads):
    """Check if the reads are valid for mutation.

    Args:
        var_record (dict): Dictionary containing the previously mutated reads.
        variant (Variant): Variant to be introduced.
        covering_reads (list): List of reads covering the variant.

    Returns:
        tuple: Tuple containing the blacklisted reads and the blacklisted UMIs.
    """
    read_blacklist = []
    umi_blacklist = []

    for read in covering_reads:
        read_id = get_read_id(read)

        # If the read has already been mutated, use the mutated read instead
        if read_id in var_record[variant.chr][read.reference_start]:
            read = var_record[variant.chr][read.reference_start][read_id]

        umi = get_umi_from_qname(read.query_name)
        expanded_cigar = expand_cigar(read.cigarstring)

        # If the UMI has already been blacklisted, skip the read
        if umi in umi_blacklist:
            read_blacklist.append(read)
            continue

        # If the position of the variant is not in the read (e.g. deletion) blacklist the UMI
        if variant.pos0 not in read.get_reference_positions():
            read_blacklist.append(read)
            umi_blacklist.append(umi)
            continue

        # If the variant is a MNV or a DEL, check that the read covers the whole variant
        if variant.type == "MNV" or variant.type == "DEL":
            indel_pos = read.get_reference_positions().index(variant.pos0) + 1
            length = len(read.seq)

            if indel_pos + variant.length > length:
                read_blacklist.append(read)
                umi_blacklist.append(umi)
                continue
                
        seq_pos = read.get_reference_positions().index(variant.pos0)
        var_cigar = expanded_cigar[seq_pos+1:seq_pos+variant.length+1]

        # If the variant is inside a clipped region, blacklist the UMI
        if "S" in var_cigar or "H" in var_cigar:
            read_blacklist.append(read)
            umi_blacklist.append(umi)
            continue

    return read_blacklist, umi_blacklist


def valid_sim(read):
    """Check if the simulation has been done properly. This is done by checking
    if the length of the CIGAR string matches the length of the sequence."""
    expanded_cigar = expand_cigar(read.cigarstring)

    cigar_length = len(expanded_cigar) - expanded_cigar.count("D") - expanded_cigar.count("H")

    return cigar_length == len(read.seq)


def calculate_sb(a, b, c, d):
    """Calculate strand bias from the counts of the four possible combinations of strands and mutations.

    Args:
        a (int): Number of reads with the reference allele in the forward strand.
        b (int): Number of reads with the alternative allele in the forward strand.
        c (int): Number of reads with the reference allele in the reverse strand.
        d (int): Number of reads with the alternative allele in the reverse strand.

    Returns:
        float: Strand bias.
    """
    return abs((b / (a+b)) - (d / (c+d))) / ((b+d) / (a+b+c+d))


def change_pos(pos, read):
    """Change the position of a variant in a read if an indel has been introduced (or is
    already present in the read). Deletions decrease the position of the variant, 
    while insertions increase it.

    Position in sequence and position in CIGAR string are handled separately,
    to comply with the pysam API.
    
    Args:
        pos (int): Intended position of the variant in the read.
        read (pysam.AlignedSegment): Read containing the variant.
    
    Returns:
        tuple: Tuple with the following elements:
            int: New position of the variant in the read.
            int: New position of the variant in the CIGAR string.
    """
    new_pos, cigar_pos = pos, pos 

    expanded_cigar = expand_cigar(read.cigarstring)

    for elem in expanded_cigar:
        if elem == "D":
            new_pos += 0
            cigar_pos += 1
        elif elem == "I":
            new_pos += 1
            cigar_pos += 1

    return new_pos, cigar_pos


def check_indel_cigar(cigar):
    """Return True if the CIGAR string contains an indel"""
    return "I" in cigar or "D" in cigar


def change_cigar(cigar, variant, pos):
    """Change the CIGAR string of a read after an indel has been introduced.

    Args:
        cigar (str): CIGAR string of the read.
        variant (Variant): Variant to be introduced.
        pos (int): Position of the variant in the read.

    Returns:
        str: New CIGAR string of the read.
    """
    expanded_cigar = expand_cigar(cigar)

    if variant.type == "INS":
        new_cigar = expanded_cigar[:pos+1] + ("I" * variant.length) + expanded_cigar[pos+1:]

    elif variant.type == "DEL":
        new_cigar = expanded_cigar[:pos+1] + ("D" * variant.length) + expanded_cigar[pos + variant.length + 1:]

    new_cigar = "".join([str(len(list(group))) + key for key, group in itertools.groupby(new_cigar)])

    return new_cigar


def expand_cigar(cigar):
    """Expand a CIGAR string to a string of characters representing the identity of each base.

    Args:
        cigar (str): CIGAR string of the read to be expanded.

    Returns:
        str: Expanded CIGAR string.
    """
    cigar_groups = re.findall(CIGAR_PATTERN, cigar)
    expanded_cigar = "".join([elem[1] * int(elem[0]) for elem in cigar_groups])

    return expanded_cigar