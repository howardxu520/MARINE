
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


def incorporate_replaced_pos_info(aligned_seq, positions_replaced):
    """
    Return the aligned sequence string, with specified positions indicated as upper case
    and others as lower case. Also returns the bases at these positions themselves.
    """
    def upper(x): return x.upper()
    def lower(x): return x.lower()
    
    differences_function = upper
    others_function = lower
    
    indicated_aligned_seq = []
    bases_at_pos = []
    for i, a in enumerate(aligned_seq):
        if i in positions_replaced:
            indicated_aligned_seq.append(differences_function(a))
            bases_at_pos.append(a.upper())
        else:
            indicated_aligned_seq.append(others_function(a))
    return indicated_aligned_seq, bases_at_pos


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
