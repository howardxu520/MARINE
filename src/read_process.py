import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import sparse

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


def incorporate_insertions_and_deletions(aligned_sequence, cigar_tuples):
    """
    Update an aligned sequence to reflect any insertions (take away those positions) such
    that it can be better compared base-to-base to a reference sequence.
    """
    new_seq = ''
    
    current_pos = 0
    for mod, num_bases in cigar_tuples:
        if mod == 0:
            # match
            new_seq += aligned_sequence[current_pos:current_pos+num_bases]
            current_pos += num_bases
        if mod in [1, 4]:
            # insertion or soft-clipping
            current_pos += num_bases
        if mod in [2]:
            # deletion
            new_seq += ''.join(['*' for r in range(num_bases)])
        if mod in [3]:
            # N
            new_seq += ''.join(['N' for r in range(num_bases)])
            
    return new_seq


def get_contig_lengths_dict(bam_handle):
    """
    Given a bam file handle, read the header to return a dictionary
    mapping contig names to lengths.
    """
    header_lines = bam_handle.text.split("\t")
    contig_lengths = {}
    found_sn = False
    found_ln = False
    for line in header_lines:
        # Clean up line
        line = line.split("\n")[0]
        
        # Extract contig names and lengths
        if line.startswith("SN"):
            found_sn = True
            contig_name = line.split("SN:")[1]
        if found_sn == True:
            if line.startswith("LN"):
                length = int(line.split("LN:")[1])
                contig_lengths[contig_name] = length
                found_sn = False
                
    return contig_lengths


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
    md_tag = read.get_tag('MD')
    if ('G' in md_tag or 'A' in md_tag or 'T' in md_tag or 'C' in md_tag):
        # No edits present in this read, based on MD tag contents
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
    barcode = read.get_tag('CB')
    print('MD tag', md_tag)
    print("CIGAR tag", cigar_string)
    print('barcode', barcode)
    

def update_coverage_array(read, contig, contig_length, barcode_to_coverage_dict):
    # Check if we already are tracking coverage for this cell, and if not set up a new array
    barcode = read.get_tag('CB')
    position_coverage_tracker_for_contig = barcode_to_coverage_dict.get(barcode, [])
    if len(position_coverage_tracker_for_contig) == 0: 
        position_coverage_tracker_for_contig = np.zeros(contig_length)
    
    reference_positions_covered_by_read = read.get_reference_positions()
    position_coverage_tracker_for_contig[reference_positions_covered_by_read] += 1
    
    barcode_to_coverage_dict[barcode] = sparse.csr_matrix(position_coverage_tracker_for_contig)
    
    return reference_positions_covered_by_read
    
    
def get_read_information(read, contig, verbose=False):
    read_barcode = read.get_tag('CB')
    is_reverse = read.is_reverse
    reference_start = read.reference_start
    reference_end = read.reference_end
    read_id = read.query_name
    mapq = read.mapping_quality
    cigarstring = read.cigarstring
    
    if verbose:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print("Read ID:", read_id)
        print("----------------------------")
      
    # ERROR CHECKS, WITH RETURN CODE SPECIFIED
    if not has_edits(read):
        return 'no_edits', [], {}
    
    if read.is_secondary:
        return 'secondary', [], {}
       
    #if 'N' in cigarstring:
    #    return 'N', [], {}
    
    # PROCESS READ TO EXTRACT EDIT INFORMATION
    strand = '+'
    if is_reverse:
        strand = '-'
            
    alt_bases, ref_bases, qualities, positions_replaced = get_edit_information_wrapper(read, not is_reverse, verbose=verbose)

    if len(alt_bases) == 0:
        # These are reads that had deletions, and no edits.
        # They are categorized later because it is hard to tell from the MD tag if they have
        # edits at first when deletions are also indicated.
        return 'no_edits', [], {}
    
    num_edits_of_each_type = defaultdict(lambda:0)
    
    list_of_rows = []
    for alt, ref, qual, pos in zip(alt_bases, ref_bases, qualities, positions_replaced):
        assert(alt != ref)
        updated_position = pos+reference_start
        if is_reverse:
            alt = reverse_complement(alt)
            ref = reverse_complement(ref)

        distance_from_read_end = np.min([updated_position - reference_start, reference_end - updated_position])

        values_to_store = strand, distance_from_read_end, qual

        list_of_rows.append([
            read_barcode, str(contig), str(updated_position), ref, alt, read_id, strand, str(distance_from_read_end), str(qual), str(mapq)
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
        print("Failed on {}, {}".format(md_tag_parsed, e))
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

def get_edit_information(md_tag, cigar_tuples, aligned_seq, reference_seq, query_qualities, hamming_check=False, verbose=False):   
    fixed_aligned_seq = incorporate_insertions_and_deletions(aligned_seq, cigar_tuples)
    positions_replaced = get_positions_from_md_tag(md_tag, verbose=verbose)
    

    indicated_aligned_seq, alt_bases = incorporate_replaced_pos_info(fixed_aligned_seq, positions_replaced)
    
    # Account for deletions
    if '*' in indicated_aligned_seq:
        positions_deleted = find(indicated_aligned_seq, '*')
    else:
        positions_deleted = []
        
    indicated_reference_seq, ref_bases = incorporate_replaced_pos_info(reference_seq, positions_replaced, positions_deleted=positions_deleted)
    indicated_qualities, qualities = incorporate_replaced_pos_info(query_qualities, positions_replaced, qualities=True)
    
    if verbose:
        print('CIGAR tuples', cigar_tuples)
        print(fixed_aligned_seq)

        print(indicated_reference_seq)
        print(indicated_aligned_seq)
        print('alt bases', alt_bases)
        print('ref bases', ref_bases)
        
    if hamming_check:
        num_deletions = indicated_aligned_seq.count('*')
        hamming_distance = get_hamming_distance(indicated_aligned_seq, indicated_reference_seq) - num_deletions
        print("Hamming distance: {}".format(hamming_distance))
        assert(hamming_distance == len(alt_bases))
        
    return alt_bases, ref_bases, qualities, positions_replaced
    
    
def get_edit_information_wrapper(read, reverse, hamming_check=False, verbose=False):
    md_tag = read.get_tag('MD')
    cigarstring = read.cigarstring
       
    cigar_tuples = read.cigartuples
    aligned_seq = read.get_forward_sequence()
    query_qualities = read.query_qualities
    
    if not reverse:
        aligned_seq = reverse_complement(aligned_seq)
    
    reference_seq = read.get_reference_sequence().lower()
    
    if verbose:
        print("MD tag:", md_tag)
        print("CIGAR string", cigarstring)
        print("Reference seq:", reference_seq.upper())
        print("Aligned seq:", aligned_seq)
    
    return(get_edit_information(md_tag,
                                cigar_tuples, 
                                aligned_seq, 
                                reference_seq,
                                query_qualities,
                                hamming_check=hamming_check,
                                verbose=verbose
                               ))