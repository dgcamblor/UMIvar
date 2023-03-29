# adding_variants.py | Functions for inserting variants into a bam file

from collections import defaultdict
import itertools
import numpy as np
import array
import re

from src.constants import *


def rec_defaultdict():
    """Allows for a nested defaultdict using recursion."""
    return defaultdict(rec_defaultdict)


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
    """Get a unique identifier for a read, using allf of its information"""
    return hash(f"{read}")


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

    # Get the count of each frequency appearing in the UMIs
    freq_counts = defaultdict(int)
    for umi_relfreq in umi_relfreqs.values():
        freq_counts[umi_relfreq] += 1

    return umi_absfreqs, umi_relfreqs, freq_counts


def get_umi_combination(variant, umi_relfreqs, freq_counts, umis_fwd, umis_rev, umi_blacklist, sb_correct, mode):
    """Get the best combination of UMI to cover a variant at the specified allelic frequency.
    
    Args:
        variant (Variant): Variant to cover
        umi_relfreqs (dict): Dictionary with the relative frequencies of the UMIs
        freq_counts (dict): Dictionary with the count of each relative frequency appearing in the UMIs
        umis_fwd (list): List of UMIs covering the variant in the forward strand
        umis_rev (list): List of UMIs covering the variant in the reverse strand
        umi_blacklist (list): List of UMIs to ignore
        sb_correct (bool): Correct for strand bias or not
        mode (str): Mode to use. Can be "dedup" or "dup"

    Returns:
        tuple: Tuple with the following elements:
            best_af (float): Best allelic frequency
            best_combination (list): List of frequencies of the UMIs to use
    """
    if variant.af == 0:
        best_af = 0
        best_combination = []
        
        return best_af, best_combination
    
    elif variant.af == 1:
        best_combination = [umi for umi in umi_relfreqs if umi not in umi_blacklist]
        best_af = 1 if not umi_blacklist else sum([umi_relfreqs[umi] for umi in best_combination])

        return best_af, best_combination
    
    else:
        if mode == "dedup":
            ind_freq = 1/len(umi_relfreqs)
            best_af = 0

            available_umis = [umi for umi in umi_relfreqs if umi not in umi_blacklist]

            for i in range(1, len(available_umis)):
                best_af += ind_freq

                # If the current AF is greater than the target AF, evaluate it and the previous one
                if best_af >= variant.af:
                    diff_curr = abs(best_af - variant.af)
                    diff_prev = abs((best_af - ind_freq) - variant.af)

                    if diff_curr < diff_prev: 
                        n_umis = i
                    elif diff_curr > diff_prev:
                        n_umis = i - 1
                
                    break
            # If the for loop finishes, use all the UMIs
            else:
                n_umis = len(available_umis)

            # Get a best combination of UMIs
            best_combination = choose_umis_dedup(umi_relfreqs, umis_fwd, umis_rev, sb_correct, n_umis, umi_blacklist)

        elif mode == "dup":  
            # Invert frequencies for AF > 0.5 (speeds up the algorithm)
            if variant.af < 0.5:
                target_af = variant.af
            else:
                target_af = 1 - variant.af

            best_distances = []
            
            # Remove the UMIs in the blacklist from the counts of each frequency
            if umi_blacklist:
                for umi in umi_blacklist:
                    black_freq = umi_relfreqs[umi]
                    freq_counts[black_freq] -= 1

            for i in itertools.count(0):
                combinations = itertools.combinations_with_replacement(freq_counts, i)
                comb_sums = []
                for combination in combinations:
                    # Get the count of each frequency
                    comb_counts = defaultdict(int)
                    for freq in combination:
                        comb_counts[freq] += 1

                    # Check if the combination is valid (i.e. the number of UMIs with each 
                    # frequency is not greater than the number of UMIs with that frequency in the sample)
                    for freq in comb_counts:
                        if comb_counts[freq] > freq_counts[freq]:
                            break  # Not valid, skip to the next combination
                    else:
                        comb_sum = sum(combination)
                        diff = abs(comb_sum - target_af)
                        comb_sums.append((diff, combination, comb_sum))

                comb_sums.sort(key=lambda x: x[0])
                best_distances.append(comb_sums[0])

                if i > 0 and best_distances[i][0] > best_distances[i - 1][0]:
                    break

            best_af = best_distances[-2][2]
            freq_combination = best_distances[-2][1]

            # If the variant has AF > 0.5, the frequencies must be inverted
            if variant.af >= 0.5:
                best_af = 1 - best_af

                # Remove the determined best frequencies from the list of all frequencies
                all_freqs = []
                for frequency in freq_counts:
                    all_freqs += [frequency] * freq_counts[frequency]
                for frequency in freq_combination:
                    all_freqs.remove(frequency)

                freq_combination = all_freqs 

            best_combination = choose_umis_dup(umi_relfreqs, umis_fwd, umis_rev, sb_correct, freq_combination, umi_blacklist)

    return best_af, best_combination


def choose_umis_dedup(umi_relfreqs, umis_fwd, umis_rev, sb_correct, n_umis, umi_blacklist):
    """Choose UMIs to use for the variant based on the UMIs alone. Intended for
    the "dedup" mode, where frequencies are computed after removing duplicates.

    Args:
        umi_relfreqs (dict): Dictionary with the relative frequencies of the UMIs
        umis_fwd (list): List of UMIs that are on the forward strand
        umis_rev (list): List of UMIs that are on the reverse strand
        sb_correct (bool): Whether to perform a simple strand bias correction
        n_umis (int): Number of UMIs to use
        umi_blacklist (list): List of UMIs to ignore

    Returns:
        best_combination (list): List of frequencies of the UMIs to use
    """
    if not sb_correct:
        possible = [umi for umi in umi_relfreqs if umi not in umi_blacklist]
        best_combination = list(np.random.choice(possible, n_umis, replace=False))
    else:
        # Simple algorithm: try to get the same number of UMIs from each strand
        possible_fwd = [umi for umi in umis_fwd if umi not in umi_blacklist]
        possible_rev = [umi for umi in umis_rev if umi not in umi_blacklist]
        best_combination = []

        # For each UMI, choose a random UMI from the other strand
        for _ in range(n_umis):
            if len(possible_fwd) > len(possible_rev):
                fwd_umi = np.random.choice(possible_fwd)
                best_combination.append(fwd_umi)
                possible_fwd.remove(fwd_umi)
                if fwd_umi in possible_rev:
                    possible_rev.remove(fwd_umi)

            elif len(possible_fwd) < len(possible_rev):
                rev_umi = np.random.choice(possible_rev)
                best_combination.append(rev_umi)
                possible_rev.remove(rev_umi)
                if rev_umi in possible_fwd:
                    possible_fwd.remove(rev_umi)

            else:
                umi = np.random.choice(list(set(possible_fwd + possible_rev)))
                best_combination.append(umi)
                if umi in possible_fwd:
                    possible_fwd.remove(umi)
                if umi in possible_rev:
                    possible_rev.remove(umi)

    return best_combination


def choose_umis_dup(umi_relfreqs, umis_fwd, umis_rev, sb_correct, freq_combination, umi_blacklist):
    """Choose UMIs to use for the variant based on the frequencies of the UMIs.
    Intended for the "dup" mode, where frequencies are computed before removing
    duplicates.

    Args:
        umi_relfreqs (dict): Dictionary with the relative frequencies of the UMIs
        umis_fwd (list): List of UMIs that are on the forward strand
        umis_rev (list): List of UMIs that are on the reverse strand
        sb_correct (bool): Whether to perform a simple strand bias correction
        freq_combination (list): List of frequencies of the UMIs to use
        umi_blacklist (list): List of UMIs to ignore

    Returns:
        best_combination (list): List of frequencies of the UMIs to use
    """
    best_combination = []

    if not sb_correct:
        for freq in freq_combination:
            possible = [umi for umi in umi_relfreqs if umi_relfreqs[umi] == freq and umi not in best_combination and umi not in umi_blacklist]
            choice = np.random.choice(possible)
            best_combination.append(choice)

    else:
        # Simple algorithm: try to get the same number of UMIs from each strand
        strands = []
        for freq in freq_combination:
            possible_fwd = [umi for umi in umis_fwd if umi_relfreqs[umi] == freq and umi not in best_combination and umi not in umi_blacklist]
            possible_rev = [umi for umi in umis_rev if umi_relfreqs[umi] == freq and umi not in best_combination and umi not in umi_blacklist]

            if strands.count("+") > strands.count("-"):
                possible = possible_rev if possible_rev else possible_fwd
            elif strands.count("+") < strands.count("-"):
                possible = possible_fwd if possible_fwd else possible_rev
            else:
                possible = list(set(possible_fwd + possible_rev))

            choice = np.random.choice(possible)
            best_combination.append(choice)

            if choice in umis_fwd:
                strands.append("+")
            if choice in umis_rev:
                strands.append("-")

    return best_combination


def add_variants(bamfile, variants, out_varfile, af_mode, sb_correct, stats):
    """Adds variants in the reads of a BAM file. The variants are introduced in the reads
    covering the variant position. The reads are mutated in such a way that the allele frequency of the
    variant is as close as possible to the specified allele frequency.

    Variants are not introduced directly in the bamfile, but, instead, mutated reads are recorded
    in a dictionary. The dictionary is then used to create a new BAM file with the mutated reads.

    Args:
        bamfile (pysam.AlignmentFile): BAM file containing the reads to mutate
        variants (list): List of variants to introduce
        out_varfile (str): Path to the output file where the mutated reads will be stored
        af_mode (str): Mode to use for the variant introduction. Can be "dedup" or "dup".
        sb_correct (bool): Whether to correct for strand bias or not.
        stats (dict): Dictionary containing the statistics of the program.

    Returns:
        dict: A dictionary containing, for each chromosome, the IDs of the reads that have been mutated and 
        the corresponding reads and their information.
    """

    var_record = rec_defaultdict()
    indel_record = rec_defaultdict()

    added_vars = 0

    for variant in variants:
        try:
            # Get all the reads covering the variant
            covering_reads = [read for read in bamfile.fetch(variant.chr, variant.pos0, variant.pos0 + 1)]

            if not covering_reads:
                print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [no reads]{END}")
                continue

            # Check the validity of each read
            read_blacklist, umi_blacklist = check_read_validity(var_record, variant, covering_reads)

            if len(covering_reads) == len(read_blacklist):
                print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [no valid reads]{END}")
                continue

            forward_reads = [read for read in covering_reads if not read.is_reverse]
            reverse_reads = [read for read in covering_reads if read.is_reverse]

            umis = list(map(lambda x: get_umi_from_qname(x.query_name), covering_reads))  # Extract the UMIs from the read IDs
            umis_fwd = list(map(lambda x: get_umi_from_qname(x.query_name), forward_reads))
            umis_rev = list(map(lambda x: get_umi_from_qname(x.query_name), reverse_reads))

            # Get the frequencies of each UMI and their counts
            _, umi_relfreqs, freq_counts = umi_freqs(umis)

            # Get the UMIs to be mutated
            achieved_af, umis_to_mutate = get_umi_combination(variant, umi_relfreqs, freq_counts, umis_fwd, umis_rev, umi_blacklist, sb_correct, af_mode)
            print(umis_to_mutate)

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

                if read_id in indel_record[variant.chr][read.reference_start]:
                    seq_pos, cigar_pos = change_pos(seq_pos, indel_record[variant.chr][read.reference_start][read_id])
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
                    new_quals = [np.random.choice(QUALITIES) for _ in range(len(variant.alt) - 1)]
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

                # Add the read to the record
                var_record[variant.chr][read.reference_start][read_id] = read

                # Add the read to the indel record if it contains an indel
                if variant.type == "INS" or variant.type == "DEL":
                    indel_record[variant.chr][read.reference_start][read_id] = expand_cigar(read.cigarstring)

                # Get the strand of the read
                mut_strands += ["-"] if read.is_reverse else ["+"]

                bamfile.reset()

        except TypeError:
            print(f"{YELLOW}WARNING: Variant {variant} could not be introduced [internal error]{END}")
            continue

        added_vars += 1

        # Calculate strand bias
        b = mut_strands.count("+")
        d = mut_strands.count("-")
        a = len(forward_reads) - b
        c = len(reverse_reads) - d

        strand_bias = calculate_sb(a, b, c, d)

        bamfile.reset()

        if added_vars <= 7:
            print(f"{GREEN}Variant {variant} introduced at frequency {achieved_af:3f} in {len(umis_to_mutate)} UMIs (SB: {strand_bias:3f}){END}")
        elif added_vars == 8:
            print(f"{GREEN}More variants introduced... (see {out_varfile.name}){END}")
        out_varfile.write(f"{variant.chr},{variant.pos1},{variant.ref},{variant.alt},{achieved_af},{len(umis_to_mutate)},{strand_bias}\n")

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


def change_pos(pos, expanded_cigar):
    """Change the position of a variant in a read after an indel has been introduced.
    Deletions decrease the position of the variant, while insertions increase it.

    Position in sequence and position in CIGAR string are handled separately,
    to comply with the pysam API.
    
    Args:
        pos (int): Intended position of the variant in the read.
        expanded_cigar (str): Expanded CIGAR string of the read.
    
    Returns:
        tuple: Tuple with the following elements:
            int: New position of the variant in the read.
            int: New position of the variant in the CIGAR string.
    """
    new_pos, cigar_pos = pos, pos 

    for elem in expanded_cigar:
        if elem == "D":
            new_pos += 0
            cigar_pos += 1
        elif elem == "I":
            new_pos += 1
            cigar_pos += 1

    return new_pos, cigar_pos


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