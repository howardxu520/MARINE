import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import sparse
from copy import copy

complements = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N'
}

def reverse_complement(input_seq):
    """
    Get the reverse complement of a sequence
    """
    rev_input_seq = reversed(input_seq)
    new_seq = []
    for r in rev_input_seq:
        new_seq.append(complements.get(r))
        
    return ''.join(new_seq)


def incorporate_insertions_and_deletions(aligned_sequence, cigar_tuples, insertions=True, deletions=True, junctions=True):
    """
    Update an aligned sequence to reflect any insertions (take away those positions) such
    that it can be better compared base-to-base to a reference sequence.
    """
    new_seq = ''
    
    current_pos = 0
    
    positions_of_deletions = []
    
    for mod, num_bases in cigar_tuples:
        if mod in [0, 7, 8]:
            # 0 = alignment match, 7 = sequence match, 8 = sequence mismatch
            new_seq += aligned_sequence[current_pos:current_pos+num_bases]
            current_pos += num_bases
        if mod in [1]:
            # insertion -- only do this for aligned sequence, not reference
            if insertions:
                current_pos += num_bases
            
        if mod in [2]:
            # deletion -- only do this for aligned sequence, not reference
            if deletions:
                new_seq += ''.join(['*' for r in range(num_bases)])
        if mod in [3]:
            # N
            if junctions:
                new_seq += ''.join(['n' for r in range(num_bases)])
            
    return new_seq


def get_hamming_distance(str1, str2):
    assert(len(str1) == len(str2))
    distance = 0
    for i, v1 in enumerate(str1):
        v2 = str2[i]
        if v1 != v2:
            distance += 1
    return distance

    
    
def has_edits(read):
    # Are there any replacements? This will always return true if a read has any deletions,
    # as the deletions will also be followed by ACT or G...
    try:
        md_tag = read.get_tag('MD')
    except Exception as e:
        print("It seems like there is an MD tag missing", e)
        
    if ('G' in md_tag or 'A' in md_tag or 'T' in md_tag or 'C' in md_tag):
        # Edits present in this read, based on MD tag contents
        return True

def get_total_coverage_for_contig_at_position(r, coverage_dict):
    position = r.position
    contig = r.contig
    barcode = r.barcode
    return coverage_dict.get(contig).get(barcode)[position]


def print_read_info(read):
    md_tag = read.get_tag('MD')
    read_id = read.query_name
    cigar_string = read.cigarstring

    if read.has_tag('CB'):
        barcode = read.get_tag('CB')
        
    print('MD tag', md_tag)
    print("CIGAR tag", cigar_string)
    print("is_reverse", read.is_reverse)
    print("is_read1", read.is_read1)
    print("is_read2", read.is_read2)
    print("is_paired", read.is_paired)
    print("is_proper_pair", read.is_proper_pair)
    print("mate_is_reverse", read.mate_is_reverse)
    print("read id", read.query_name)

    print(str(read))
    
def get_read_information(read, contig, barcode_tag='CB', verbose=False, reverse_stranded=True, min_read_quality = 0):
    if barcode_tag is None:
        read_barcode = 'no_barcode'
    elif read.has_tag(barcode_tag):
        read_barcode = read.get_tag(barcode_tag)
    else:
        read_barcode = None
        
    if not read_barcode:
        return 'no_{}_tag'.format(barcode_tag), [], {}

    # For 10x data, exclude reads that are not counted towards cellranger UMI read counts
    # https://kb.10xgenomics.com/hc/en-us/articles/115003710383-Which-reads-are-considered-for-UMI-counting-by-Cell-Ranger
    if read.has_tag('xf'):
        if not read.get_tag('xf') == 25:
            return 'xf:{}'.format(read.get_tag('xf')), [], {}
    
    is_reverse = read.is_reverse
    reverse_or_forward = '+'

    if barcode_tag:
        # Single-cell (10x)
        
        if is_reverse:
            reverse_or_forward = '-'
            
    elif (read.is_read1 or read.is_read2):
        # Paired end
        if reverse_stranded:
            if (read.is_read1 and not is_reverse) or (read.is_read2 and is_reverse):
                reverse_or_forward = '-'
                
        if not reverse_stranded:
            if (read.is_read1 and is_reverse) or (read.is_read2 and not is_reverse):
                reverse_or_forward = '-'
    
    else:
        # Single end
        if is_reverse:
            if reverse_stranded:
                reverse_or_forward = '+'
            else:
                reverse_stranded = '-'
        else:
            if reverse_stranded:
                reverse_or_forward = '-'
            else:
                reverse_stranded = '+'
        
        
    reference_start = read.reference_start
    reference_end = read.reference_end
    read_id = read.query_name
    mapq = read.mapping_quality        
    cigarstring = read.cigarstring
    
    
    # ERROR CHECKS, WITH RETURN CODE SPECIFIED        
    if mapq < min_read_quality:
        return 'mapq_low', [], {}
        
    if not has_edits(read):
        return 'no_edits', [], {}

    # Defaults for coverage counting as well 
    # count_coverage function at: https://pysam.readthedocs.io/en/v0.16.0.1/api.html
    
    if read.is_secondary:
        return 'secondary', [], {}

    if read.is_unmapped:
        return 'is_unmapped', [], {}

    if read.is_qcfail:
        return 'is_qcfail', [], {}

    if read.is_duplicate:
        return 'is_duplicate', [], {}
    
    
    #if 'N' in cigarstring:
    #    return 'N', [], {}
    
    # PROCESS READ TO EXTRACT EDIT INFORMATION
    if verbose:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print_read_info(read)
        print('reverse_or_forward:', reverse_or_forward)
        #print("Read ID:", read_id)
        print("----------------------------")
        
    alt_bases, ref_bases, qualities, positions_replaced = get_edit_information_wrapper(read, verbose=verbose)
    if verbose:
        print("Successfully ran get_edit_information_wrapper\nalt bases: {}, ref bases: {}".format(alt_bases, ref_bases))
        
    if len(alt_bases) == 0:
        # These are reads that had deletions, and no edits.
        # They are categorized later because it is hard to tell from the MD tag if they have
        # edits at first when deletions are also indicated.
        return 'no_edits', [], {}
    
    num_edits_of_each_type = defaultdict(lambda:0)
    
    list_of_rows = []
    
    for alt, ref, qual, pos in zip(alt_bases, ref_bases, qualities, positions_replaced):
        if alt == "N" or ref == "N":
            continue

        if verbose:
            print("Getting info:", alt, ref, qual, pos)
            
        assert(alt != ref)
        updated_position = pos+reference_start
        
        distance_from_read_end = np.min([updated_position - reference_start, reference_end - updated_position])

        list_of_rows.append([
            read_barcode, str(contig), str(updated_position), ref, alt, read_id, reverse_or_forward, str(distance_from_read_end), str(qual), str(mapq)
        ])
        
        num_edits_of_each_type['{}>{}'.format(ref, alt)] += 1
    
    return None, list_of_rows, num_edits_of_each_type



def get_positions_from_md_tag(md_tag, verbose=False):
    """
    Figure out which positions are replaced, from the MD tag.
    """ 
    md_tag_parsed = []
    
    in_deletion = False
    
    for c in md_tag:
        if c == '^':
            in_deletion = True
            continue
        else:
            
            try:
                value = str(int(c))
                
                if in_deletion:
                    in_deletion = False
                    md_tag_parsed.append('-')
                    
                md_tag_parsed.append(value)
                
            except Exception as e:
                if not in_deletion:
                    md_tag_parsed.append('-')
                else:
                    md_tag_parsed.append('+1')

    positions = []

    try:
        position_splitters = [i for i in ''.join(md_tag_parsed).split('-')]
    except Exception as e:
        print("Failed splitting possition on {}, {}".format(md_tag_parsed, e))
        return None
    
    if verbose:
        print(position_splitters)
    
    for s in position_splitters:
        # account for plus signs
        if '+' in s:
            s_sum = np.sum([int(i) for i in s.split('+')])
            s = s_sum - 1
        else:
            s = int(s)
        if len(positions) == 0:
            positions.append(s)
        else:
            positions.append(positions[-1] + s + 1)
            
    if verbose:
        print(positions)
        
    return positions


def incorporate_replaced_pos_info(aligned_seq, positions_replaced, positions_deleted=[], qualities=False):
    """
    Return the aligned sequence string, with specified positions indicated as upper case
    and others as lower case. Also returns the bases at these positions themselves.
    """
    def upper(x): return x.upper()
    def lower(x): return x.lower()
    def nothing(x): return str(x)
    
    if not qualities:
        differences_function = upper
        others_function = lower
    else:
        differences_function = nothing
        others_function = nothing
        
    indicated_aligned_seq = []
    bases_at_pos = []
    for i, a in enumerate(aligned_seq):
        if a == '*':
            indicated_aligned_seq.append(a)
            continue
            
        if i in positions_replaced and i not in positions_deleted:
            indicated_aligned_seq.append(differences_function(a))
            if not qualities:
                bases_at_pos.append(a.upper())
            else:
                bases_at_pos.append(str(a))
        else:
            indicated_aligned_seq.append(others_function(a))
    return ''.join(indicated_aligned_seq), bases_at_pos

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


def remove_softclipped_bases(cigar_tuples, aligned_sequence):
    had_front_clipped = 0
    had_back_clipped = 0
    
    first_tuple = cigar_tuples[0]
    last_tuple = cigar_tuples[-1]
    
    to_clip_from_front = 0
    to_clip_from_back = 0
    
    if first_tuple[0] == 4:
        to_clip_from_front = first_tuple[1]
        had_front_clipped = 1
    if last_tuple[0] == 4:
        to_clip_from_back = last_tuple[1]
        had_back_clipped = 1
        
    cropped_sequence = aligned_sequence[to_clip_from_front:(len(aligned_sequence)-to_clip_from_back)] 
    cropped_tuples = cigar_tuples[had_front_clipped:len(cigar_tuples)-had_back_clipped]
    
    return cropped_sequence, cropped_tuples  


def get_edit_information(md_tag, cigar_tuples, aligned_seq, reference_seq, query_qualities, hamming_check=False, verbose=False): 
    if verbose:
            print('CIGAR tuples before clipping (if needed):\n', cigar_tuples)
            print('Aligned sequence before clipping (if needed):\n', aligned_seq)
            print("Qualities before clipping:\n", query_qualities)

    original_cigar_tuples = copy(cigar_tuples)
    aligned_seq, cigar_tuples = remove_softclipped_bases(original_cigar_tuples, aligned_seq)
    if verbose:
        print("Soft clipping quality scores ...")
    query_qualities, cigar_tuples = remove_softclipped_bases(original_cigar_tuples, query_qualities)

    if verbose:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('CIGAR tuples after clipping (if needed):\n', cigar_tuples)
        print('Aligned sequence after clipping (if needed):\n', aligned_seq)
        print("Qualities after clipping:\n", query_qualities)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


    positions_replaced = get_positions_from_md_tag(md_tag, verbose=verbose)
    # Incorporate insertions, deletions, and splice junctions (N)
    fixed_aligned_seq_for_deletion_locations = incorporate_insertions_and_deletions(aligned_seq, cigar_tuples, junctions=False)
    fixed_aligned_seq = incorporate_insertions_and_deletions(aligned_seq, cigar_tuples)

    # Account for deletions
    if '*' in fixed_aligned_seq_for_deletion_locations:
        # These are coordinates after already fixing the read with deletions, insertions and junctions
        positions_deleted = find(fixed_aligned_seq_for_deletion_locations, '*')
    else:
        positions_deleted = []

    indicated_reference_seq, ref_bases = incorporate_replaced_pos_info(reference_seq, positions_replaced, positions_deleted=positions_deleted)
    fixed_reference_seq = incorporate_insertions_and_deletions(indicated_reference_seq, cigar_tuples, insertions=False, deletions=False)
    indicated_qualities, qualities = incorporate_replaced_pos_info(query_qualities, positions_replaced, qualities=True)


    # Get global coordinates of positions replaced
    global_positions_replaced = []
    finalized_fixed_aligned_seq = ''
    alt_bases = []
    for i, character in enumerate(fixed_reference_seq):
        if character.isupper():
            global_positions_replaced.append(i)
            upper_char = fixed_aligned_seq[i].upper()
            finalized_fixed_aligned_seq += upper_char
            if upper_char != '*':
                alt_bases.append(upper_char)
        else:
            lower_char = fixed_aligned_seq[i].lower()
            finalized_fixed_aligned_seq += lower_char

    if verbose:
        if 'n' in fixed_aligned_seq:
            num_n = fixed_aligned_seq.count("n")
            n_to_replace = ''.join([i for i in fixed_aligned_seq if i == "n"])
        else:
            num_n = 1
            n_to_replace = 'n'
            
        print("Indicated reference seq:\n", indicated_reference_seq.replace(n_to_replace, "{}*n".format(num_n)))
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print("Fixed reference seq:\n", fixed_reference_seq.replace(n_to_replace, "{}*n".format(num_n)))
            
        print("Fixed aligned seq:\n", fixed_aligned_seq.replace(n_to_replace, "{}*n".format(num_n)))
        print("Finalized fixed aligned seq:\n", finalized_fixed_aligned_seq.replace(n_to_replace, "{}*n".format(num_n)))

        print("Indicated qualities:\n", indicated_qualities.replace(n_to_replace, "{}*n".format(num_n)))
        
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('alt bases', alt_bases)
        print('ref bases', ref_bases)


    if hamming_check:
        num_deletions = finalized_fixed_aligned_seq.count('*')
        hamming_distance = get_hamming_distance(finalized_fixed_aligned_seq, fixed_reference_seq) - num_deletions
        print("Hamming distance: {}".format(hamming_distance))
        assert(hamming_distance == len(alt_bases))

    # Make positions 1-based instead of 0-based
    global_positions_replaced_1based = [g+1 for g in global_positions_replaced]
    
    return alt_bases, ref_bases, qualities, global_positions_replaced_1based
    
    
def get_edit_information_wrapper(read, hamming_check=False, verbose=False):
    md_tag = read.get_tag('MD')
    cigarstring = read.cigarstring
       
    cigar_tuples = read.cigartuples
    aligned_seq = read.get_forward_sequence()
    query_qualities = read.query_qualities

    if read.is_reverse:
        aligned_seq = reverse_complement(aligned_seq)
    
    reference_seq = read.get_reference_sequence().lower()
    
    if verbose:
        print("MD tag:", md_tag)
        print("CIGAR string", cigarstring)
        print("Reference seq:", reference_seq.upper())
        print("Aligned seq:", aligned_seq)
        print("Qualities:", query_qualities)
    
    return(get_edit_information(md_tag,
                                cigar_tuples, 
                                aligned_seq, 
                                reference_seq,
                                query_qualities,
                                hamming_check=hamming_check,
                                verbose=verbose
                               ))