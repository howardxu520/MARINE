import pandas as pd
import numpy as np

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


def get_positions_from_md_tag(md_tag):
    """
    Figure out which positions are replaced, from the MD tag.
    """
    md_tag_parsed = []
    for c in md_tag:
        try:
            md_tag_parsed.append(str(int(c)))
        except Exception as e:
            md_tag_parsed.append('-')

    positions = []

    try:
        position_splitters = [int(i) for i in ''.join(md_tag_parsed).split('-')]
    except Exception as e:
        print("Failed on {}, {}".format(md_tag_parsed, e))
        return None
    
    for s in position_splitters:
        if len(positions) == 0:
            positions.append(s)
        else:
            positions.append(positions[-1] + s + 1)
    return positions


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
            
    return new_seq


def incorporate_replaced_pos_info(aligned_seq, positions_replaced, qualities=False):
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
        if i in positions_replaced:
            indicated_aligned_seq.append(differences_function(a))
            if not qualities:
                bases_at_pos.append(a.upper())
            else:
                bases_at_pos.append(str(a))
        else:
            indicated_aligned_seq.append(others_function(a))
    return ''.join(indicated_aligned_seq), bases_at_pos


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


def get_edit_information(md_tag, cigar_tuples, aligned_seq, reference_seq, query_qualities, hamming_check=False):    
    fixed_aligned_seq = incorporate_insertions_and_deletions(aligned_seq, cigar_tuples)
    positions_replaced = get_positions_from_md_tag(md_tag)
    
    indicated_aligned_seq, alt_bases = incorporate_replaced_pos_info(fixed_aligned_seq, positions_replaced)
    indicated_reference_seq, ref_bases = incorporate_replaced_pos_info(reference_seq, positions_replaced)
    indicated_qualities, qualities = incorporate_replaced_pos_info(query_qualities, positions_replaced, qualities=True)
    
    if hamming_check:
        hamming_distance = get_hamming_distance(indicated_aligned_seq, indicated_reference_seq)
        print("Hamming distance: {}".format(hamming_distance))
        assert(hamming_distance == len(alt_bases))
        
    return alt_bases, ref_bases, qualities, positions_replaced
    
    
def get_edit_information_wrapper(read, reverse, hamming_check=False):
    md_tag = read.get_tag('MD')
    cigar_tuples = read.cigartuples
    aligned_seq = read.get_forward_sequence()
    query_qualities = read.query_qualities
    
    if not reverse:
        aligned_seq = reverse_complement(aligned_seq)
    
    reference_seq = read.get_reference_sequence().lower()
    return(get_edit_information(md_tag,
                                cigar_tuples, 
                                aligned_seq, 
                                reference_seq,
                                query_qualities,
                                hamming_check=hamming_check))
    
    
def has_edits(read):
    # Are there any replacements?
    md_tag = read.get_tag('MD')
    if ('G' in md_tag or 'A' in md_tag or 'T' in md_tag or 'C' in md_tag):
        # No edits present in this read, based on MD tag contents
        return True
    
    
def get_dataframe_from_barcode_dict(barcode_to_position_to_alts):
    all_rows = []
    for barcode, contig_dict in barcode_to_position_to_alts.items():
        for contig, pos_dict in contig_dict.items():
                for pos, alt_dict in pos_dict.items():
                    for alt, read_dict in alt_dict.items():
                        for read, value in read_dict.items():
                            new_row = (barcode, contig, pos, alt, read, value[0], value[1], value[2])
                            all_rows.append(new_row)

    example_dataframe = pd.DataFrame(all_rows, columns=['barcode', 'contig', 'position_ref', 'alt', 'read_id', 'strand', 'dist_from_read_end', 'quality'])
    example_dataframe['ref'] = [c.split('_')[1] for c in example_dataframe['position_ref']]
    example_dataframe['position'] = [int(c.split('_')[0]) for c in example_dataframe['position_ref']]

    return example_dataframe


def get_total_coverage_for_contig_at_position(r, coverage_dict):
    position = r.position
    contig = r.contig
    return coverage_dict.get(contig)[position]


def print_read_info(read):
    md_tag = read.get_tag('MD')
    read_id = read.query_name
    cigar_string = read.cigarstring
    barcode = read.get_tag('CR')
    print('MD tag', md_tag)
    print("CIGAR tag", cigar_string)
    print('barcode', barcode)
    
    
def update_coverage_array(read, position_coverage_tracker_for_contig):
    reference_positions_covered_by_read = read.get_reference_positions()
    position_coverage_tracker_for_contig[reference_positions_covered_by_read] += 1
    
    
def add_read_information_to_barcode_dict(read, contig, barcode_to_position_to_alts):
    read_barcode = read.get_tag('CR')
    is_reverse = read.is_reverse
    reference_start = read.reference_start
    reference_end = read.reference_end
    read_id = read.query_name
        
    # ERROR CHECKS, WITH RETURN CODE SPECIFIED
    if not has_edits(read):
        return 'no_edits'
        
    if '^' in read.get_tag('MD'):
        if '^' in read.get_tag('MD'):
            # FOR NOW SKIP DELETIONS, THEY ARE TRICKY TO PARSE...
            return 'deletion'
    
    if read.is_secondary:
        return 'secondary'
        
    # PROCESS READ TO EXTRACT EDIT INFORMATION
    strand = '+'
    if is_reverse:
        strand = '-'
            
    alt_bases, ref_bases, qualities, positions_replaced = get_edit_information_wrapper(read, not is_reverse)

    for alt, ref, qual, pos in zip(alt_bases, ref_bases, qualities, positions_replaced):
        assert(alt != ref)
        updated_position = pos+reference_start
        if is_reverse:
            alt = reverse_complement(alt)
            ref = reverse_complement(ref)

        distance_from_read_end = np.min([updated_position - reference_start, reference_end - updated_position])

        values_to_store = strand, distance_from_read_end, qual

        barcode_to_position_to_alts[read_barcode][contig]['{}_{}'.format(updated_position, ref)]\
        [alt][read_id] = values_to_store