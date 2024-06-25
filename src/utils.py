import math
from glob import glob
import os
import pysam
import polars as pl
import pandas as pd
import numpy as np
import sys
from collections import OrderedDict, defaultdict

def get_contigs_that_need_bams_written(expected_contigs, split_bams_folder, barcode_tag='CB', number_of_expected_bams=4):
    bam_indices_written = [f.split('/')[-1].split('.bam')[0] for f in glob('{}/*/*.sorted.bam.bai'.format(split_bams_folder))]

    subsets_per_contig = defaultdict(lambda:0)
    for bam_index_written in bam_indices_written:
        contig_label = bam_index_written.split('_')[0]
        subsets_per_contig[contig_label] += 1

    if barcode_tag == 'CB':
        number_of_expected_bams = 4
    else:
        number_of_expected_bams = number_of_expected_bams
        
    contigs_to_write_bams_for = []
    for c in expected_contigs:
        num_written_indices = subsets_per_contig.get(c, 0)
        if num_written_indices < number_of_expected_bams:
            print("Contig {} has {}/{} bams generated".format(c, num_written_indices, number_of_expected_bams))
            contigs_to_write_bams_for.append(c)
    
    return contigs_to_write_bams_for


def make_edit_finding_jobs(bampath, output_folder, reverse_stranded=True, barcode_tag="CB", barcode_whitelist=None, contigs=[], num_intervals_per_contig=16, verbose=False, min_read_quality = 0):
    jobs = []
    
    samfile = pysam.AlignmentFile(bampath, "rb")
    contig_lengths_dict = get_contig_lengths_dict(samfile)
    print('contig_lengths_dict:{}'.format(contig_lengths_dict))
    if len(contigs) == 0:
        contigs_to_use = set(contig_lengths_dict.keys())
    else:
        contigs_to_use = set(contigs)
    for contig in contig_lengths_dict.keys():
        # Skip useless contigs
        if len(contig) > 5:
            continue
            
        if contig == 'Stamp':# or contig != '17':
            continue
            
        if contig not in contigs_to_use:
            continue

        pretty_print("\tContig {}".format(contig))
        
        contig_length = contig_lengths_dict.get(contig)
        intervals_for_contig = get_intervals(contig, contig_lengths_dict, num_intervals_per_contig)
        #print('\t\tintervals_for_contig: {}'.format(intervals_for_contig))
        # Set up for pool
        for split_index, interval in enumerate(intervals_for_contig):
            split_index = str(split_index).zfill(3)
            parameters = [bampath, contig, split_index, interval[0], interval[1], output_folder, reverse_stranded, barcode_tag, barcode_whitelist, verbose, min_read_quality]
            jobs.append(parameters)
    return jobs


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

def read_barcode_whitelist_file(barcode_whitelist_file):
    barcode_whitelist = set(pd.read_csv(barcode_whitelist_file, names=['barcodes']).barcodes.tolist())
    
    pretty_print("Barcodes in whitelist: {}".format(len(barcode_whitelist)))
    return barcode_whitelist


def pretty_print(contents, style=''):
    if type(contents) == list:
        for item in contents:
            pretty_print(item)
        sys.stdout.write("\n")
        
    else:
        to_write = '{}\n'.format(contents)
        
        before_line = None
        after_line = None
        
        styled_line = ''.join([style for i in range(len(to_write))])
        
        if style != '':
            # Line before
            pretty_print(styled_line)
            
        sys.stdout.write(to_write)
        
        if style:
            # Line after 
            pretty_print(styled_line)
    
def get_intervals(contig, contig_lengths_dict, num_intervals=4):
    """
    Splits contig coordinates into a list of {num_intervals} coordinates of equal size.
    """
    contig_length = contig_lengths_dict.get(contig)
    interval_length = math.ceil(contig_length/num_intervals)
    start = 0
    end = interval_length

    intervals = []
    while start < contig_length:
        if end > contig_length:
            end == contig_length

        interval = [start, end]
        intervals.append(interval)

        start = end
        end = start + interval_length
    return intervals


def index_bam(barcode_bam_file_path):
    index_path = '{}.bai'.format(barcode_bam_file_path)
    if not os.path.exists(index_path):
        pysam.index(barcode_bam_file_path)
    elif os.path.getmtime(index_path) < os.path.getmtime(barcode_bam_file_path):
        print("index is older than bam for {}... reindexing".format(barcode_bam_file_path))
        pysam.index(barcode_bam_file_path)
            
            
def write_rows_to_info_file(list_of_rows, f):
    for info_list in list_of_rows:
        info_line = '\t'.join(info_list) + '\n'
        f.write(info_line)
        
def write_header_to_edit_info(f):
    f.write('barcode\tcontig\tposition\tref\talt\tread_id\tstrand\tdist_from_end\tbase_quality\tmapping_quality\n')
    
def write_read_to_bam_file(read, bam_handles_for_barcodes, barcode_bam_file_path, samfile_template):
    bam_for_barcode = bam_handles_for_barcodes.get(barcode_bam_file_path)

    if not bam_for_barcode:
        bam_for_barcode = pysam.AlignmentFile(barcode_bam_file_path, "wb", template=samfile_template)
        bam_handles_for_barcodes[barcode_bam_file_path] = bam_for_barcode

    bam_for_barcode.write(read)
    
def remove_file_if_exists(file_path):
    if os.path.exists(file_path):
        #print("{} exists... deleting".format(file_path))
        os.remove(file_path)
        
def make_folder(folder_path):
    if not os.path.exists(folder_path):
        try:
            os.mkdir(folder_path)
        except Exception as e:
            pass
            #sys.stderr.write('{}_{}\n'.format(folder_path, e))

def only_keep_positions_for_region(contig, output_folder, positions_for_barcode, verbose=False):
    contig_index = str(contig.split("_")[-1]).zfill(3)
    contig_base = contig.split("_")[0]
    
    #print("Contig: {}, contig_base: {}, contig_index: {}".format(contig, contig_base, contig_index))

    edit_info_filepath_regex = "{}/edit_info/{}_{}*".format(output_folder, contig_base, contig_index)
    #print("edit_info_filepath_regex: {}".format(edit_info_filepath_regex))

    edit_info_filepath = 'Not yet set'
    try:
        edit_info_filepath = glob(edit_info_filepath_regex)[0].split('/')[-1]
        #sys.stderr.write("Contig: {}, edit_info_filepath: {}\n".format(contig, edit_info_filepath))
        
        second_half_of_filepath = edit_info_filepath.split('_{}_'.format(contig_index))[-1]
        #print("Contig: {}, second_half_of_filepath: {}".format(contig,second_half_of_filepath))
        min_pos, max_pos = second_half_of_filepath.split('_edit_info')[0].split('_')
    
        min_pos = int(min_pos)
        max_pos = int(max_pos)
        #print("Contig: {}. Min pos: {}. Max pos: {}.".format(contig, min_pos, max_pos))

        return [p for p in positions_for_barcode if p > min_pos and p < max_pos]

    except Exception as e:
        sys.stderr.write('{}, {}, {}, {}, {}, regex: {}\n'.format(e, contig, contig_base, contig_index, edit_info_filepath, edit_info_filepath_regex))
        return []
        
def check_read(read):
    return True

def get_bulk_coverage_at_pos(samfile_for_barcode, contig_bam, just_contig, pos, paired_end=False, verbose=False):
    if not paired_end:
        if verbose:
            print("~~~~~~\n!!!!SINGLE END!!!!!\n~~~~~~~`")
        # count_coverage does not work with paired end reads, because it counts both ends of the read twice...
        # For more information: https://github.com/pysam-developers/pysam/issues/744 
        coverage_at_pos = np.sum(samfile_for_barcode.count_coverage(just_contig, 
                                                                    pos-1, 
                                                                    pos, 
                                                                    quality_threshold=0, # base quality
                                                                    read_callback='all'
                                                                   ))
        return coverage_at_pos
    else:
        if verbose:
            print("~~~~~~\n!!!!PAIRED END!!!!!\npos: {}\n {}:{}~~~~~~~`".format(contig_bam, just_contig, pos))
        
        # For paired-end sequencing, it is more accurate to use the pileup function
        pileupcolumn_iter = samfile_for_barcode.pileup(just_contig, 
                                                       pos-1, 
                                                       pos, 
                                                       stepper='nofilter',
                                                       #stepper='all', 
                                                       truncate=True, 
                                                       min_base_quality=0, # base quality
                                                       max_depth=1000000)
        
        unique_read_ids = set()
        for pileupcolumn in pileupcolumn_iter:
            for pileupread in pileupcolumn.pileups:
                #if verbose:
                #    print("at pos {}:{}:".format(just_contig, pos), pileupread)
                if not pileupread.is_del and not pileupread.is_refskip:                                        
                    unique_read_ids.add(pileupread.alignment.query_name)
                    
        coverage_at_pos = len(unique_read_ids)

        if verbose:
            print('coverage_at_pos', coverage_at_pos)
            
        return coverage_at_pos

def get_coverage_df(edit_info, contig, output_folder, barcode_tag='CB', paired_end=False, 
                    verbose=False):

    if barcode_tag:
        # Single-cell, contig will include barcode ending-based suffix 
        bam_subfolder = "{}/split_bams/{}".format(output_folder, contig)
        contig_bam = '{}/{}.bam.sorted.bam'.format(bam_subfolder, contig)
    else:
        # Bulk -- no barcode ending-based suffixes, just splits
        just_contig = contig.split('_')[0]
        bam_subfolder = "{}/split_bams/{}".format(output_folder, just_contig)
        contig_bam = '{}/{}.bam.sorted.bam'.format(bam_subfolder, contig)

    #print("Contig {}. Loading {} bamfile...".format(contig, contig_bam))
    try:
        index_bam(contig_bam)
        samfile_for_barcode = pysam.AlignmentFile(contig_bam, "rb")
    except OSError:  # If MARINE cannot find a bam to index, return empty
        return pd.DataFrame(columns=['coverage', 'source'])
        
    #print("Contig {}. Loaded bamfile...".format(contig))
    unique_barcodes = sorted(list(edit_info.unique("barcode")["barcode"]))
    coverage_dict = defaultdict(lambda:{})
    
    # print("Contig {}. Iterating through barcodes...".format(contig))
    for i, barcode in enumerate(unique_barcodes):
        if i % 300 == 0:
            #print("{}/{} barcodes for {}...".format(i+1, len(unique_barcodes), contig))
            pass
            
        edit_info_for_barcode = edit_info.filter(pl.col("barcode") == barcode)
        positions_for_barcode = set(edit_info_for_barcode["position"].unique())
                    
        num_positions = len(positions_for_barcode)
        for i, pos in enumerate(positions_for_barcode):
            
            if barcode_tag:
                # For single-cell
                
                barcode_specific_contig = '{}_{}'.format(contig, barcode)
                # Alter from 3_C_AAACCCAAGAACTTCC-1, for example, to 3_AAACCCAAGAACTTCC-1'
                barcode_specific_contig_split = barcode_specific_contig.split("_")
                barcode_specific_contig_without_subdivision = "{}_{}".format(barcode_specific_contig_split[0], barcode_specific_contig_split[2])

                if verbose:
                    #print("contig_bam:", contig_bam)
                    #print("barcode_specific_contig:", barcode_specific_contig)
                    print("barcode_specific_contig_without_subdivision:", barcode_specific_contig_without_subdivision)
                    
                coverage_at_pos = np.sum(samfile_for_barcode.count_coverage(barcode_specific_contig_without_subdivision, 
                                                                            pos-1, 
                                                                            pos,  
                                                                            quality_threshold=0, # base quality
                                                                            read_callback='all'
                                                                           ))

                coverage_dict['{}:{}'.format(barcode, pos)]['coverage'] = coverage_at_pos
                coverage_dict['{}:{}'.format(barcode, pos)]['source'] = contig
                
            else:
                
                # For bulk, no barcodes, we will just have for example 19_no_barcode to convert to 19 to get coverage at that chrom
                just_contig = contig.split('_')[0]
                try:
                    coverage_at_pos = get_bulk_coverage_at_pos(samfile_for_barcode, contig_bam, just_contig, pos, 
                                                               paired_end=paired_end, verbose=verbose)
                                
                    coverage_dict['{}:{}'.format('no_barcode', pos)]['coverage'] = coverage_at_pos
                    coverage_dict['{}:{}'.format('no_barcode', pos)]['source'] = contig

                except Exception as e:
                    print("contig {}".format(contig), "Failed to get coverage", e)
                    return pd.DataFrame(columns=['coverage', 'source'])  # Should we just return empty? Or why would there be an exception here?
                
    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index')
    
    return coverage_df


def get_coverage_wrapper(parameters):
    edit_info, contig, output_folder, barcode_tag, paired_end, verbose = parameters

    output_filename = '{}/coverage/{}.tsv'.format(output_folder, contig, header=False)
    if os.path.exists(output_filename):
        return output_filename
        
    edit_info = edit_info.with_columns(
    pl.concat_str(
        [
            pl.col("barcode"),
            pl.col("position")
        ],
        separator=":",
    ).alias("barcode_position"))

    edit_info_df = edit_info.to_pandas()
    edit_info_df.index = edit_info_df['barcode_position']
    
    coverage_df = get_coverage_df(edit_info, contig, output_folder, barcode_tag=barcode_tag, 
                                  paired_end=paired_end, verbose=verbose)

    # Combine edit information with coverage information
    edit_info_and_coverage_joined = edit_info_df.join(coverage_df, how='inner')

    edit_info_and_coverage_joined['position_barcode'] = edit_info_and_coverage_joined['position'].astype(str) + '_' + edit_info_and_coverage_joined['barcode'].astype(str)
    edit_info_and_coverage_joined.to_csv(output_filename, sep='\t', header=False)
    
    return output_filename
    


def sort_bam(bam_file_name):
    output_name = bam_file_name + ".sorted.bam"
    pysam.sort("-o", output_name, bam_file_name)  
    return output_name


def rm_bam(bam_file_name):
    os.remove(bam_file_name)

def write_reads_to_file(reads, bam_file_name, header_string, barcode_tag="BC"):
    header = pysam.AlignmentHeader.from_text(header_string)
    
    header_dict = header.as_dict()
    lengths_for_sn = {}
    
    header_dict_sq = header_dict.get("SQ")
    for s in header_dict_sq:
        sn = s.get("SN")
        ln = s.get("LN")
        lengths_for_sn[sn] = ln
        
    #print("\tCurrent header length for {}: {}".format(bam_file_name, len(lengths_for_sn)))
    
    all_barcodes_for_contig = set([r.split('\t')[2] for r in reads])
    #print("\tNum barcodes for {}: {}".format(bam_file_name, len(all_barcodes_for_contig)))
        
    for new_sn in all_barcodes_for_contig:
        original_chrom = new_sn.split("_")[0]
        if not barcode_tag:
            new_sn = new_sn.split("_")[0]
            
        if new_sn not in lengths_for_sn:
            new_ln = lengths_for_sn.get(original_chrom)
            new_entry = {"SN": new_sn, "LN": new_ln}
            header_dict_sq.append(new_entry)
    
    #print("\tExample new entries: {}".format(header_dict_sq[-4:]))
    header_dict['SQ'] = header_dict_sq
    
    #print("\tNew header length: {}".format(len(header_dict.get("SQ"))))
    
    new_header = pysam.AlignmentHeader.from_dict(header_dict)
    
    num_reads = len(reads)
    
    with pysam.AlignmentFile(bam_file_name, "wb", text=str(new_header)) as bam_handle:
        for i, read_str in enumerate(reads):
            if i > 0 and i % 100000 == 0:
                sys.stdout.write('file {}: {}/{} reads\n'.format(bam_file_name.split('/')[-1], i, num_reads))
                
            try:
                read = pysam.AlignedSegment.fromstring(read_str, new_header)
                bam_handle.write(read) 
            except Exception as e:
                sys.stdout.err('{}\n\nfile {}: Failed to write read with str representation of:\n\t {}\n'.format(e,
                                                                                                      bam_file_name.split('/')[-1],
                                                                                                read_str))
                sys.exit(1)
                
            
            
    bam_handle.close()
        
        
        
def concat_and_write_bams(contig, df_dict, header_string, split_bams_folder, barcode_tag='CB', number_of_expected_bams=4, verbose=False):
    job_params = []
    
    # Sort the subcontig regions such that the reads are properly ordered 
    sorted_subcontig_names = sorted(df_dict.keys())
    sorted_subcontig_dfs = []
    for n in sorted_subcontig_names:
        subcontig = df_dict.get(n)
        if len(subcontig) > 0:
            sorted_subcontig_dfs.append(subcontig)
        
    if len(sorted_subcontig_dfs) == 0:
        if verbose:
            print("{} was empty".format(contig))
        return []

    if verbose:
        print("\t{}: num subcontigs to concat: {}".format(contig, len(sorted_subcontig_dfs)))
        # All of the reads for all of the barcodes are in this dataframe
        print("\t{}: concatting".format(contig))
        
    all_contents_df = pl.concat(sorted_subcontig_dfs)
                
    # Combine the reads (in string representation) for all rows corresponding to a barcode  
    if barcode_tag:
        suffix_options = ["A-1", "C-1", "G-1", "T-1"]
    else:
        # If we are not splitting up contigs by their barcode ending, instead let's do it by random bucket
        range_for_suffixes = number_of_expected_bams
        suffix_options = range(0, range_for_suffixes)
        
    for suffix in suffix_options:
        if barcode_tag:
            all_contents_for_suffix = all_contents_df.filter(pl.col('barcode').str.ends_with(suffix))
        else:
            all_contents_for_suffix = all_contents_df.filter(pl.col('bucket') == suffix)
            

        if verbose:
            print("\tcontig: {} suffix: {}, all_contents_df length: {}, all_contents_for_suffix length: {}".format(
                    contig,
                    suffix,
                    len(all_contents_df),
                    len(all_contents_for_suffix)
                    ))
        
        try:
            reads_deduped = list(OrderedDict.fromkeys(all_contents_for_suffix.transpose().with_columns(
                pl.concat_str(
                    [pl.col(c) for c in all_contents_for_suffix.transpose().columns],
                    separator="\n"
                     ).alias("combined_text")
            )[['combined_text']][1].item().split('\n')))
        except Exception as e:
            reads_count_for_suffix = len(all_contents_for_suffix)
            
            if reads_count_for_suffix == 0:
                print("\tWARNING: No reads found in region {}:{}... assuming this is not an issue, but perhaps worth confirming manually using samtools.".format(contig, suffix))
                reads_deduped = []
            else:
                print("\t\t### ERROR EMPTY?: {}, contig: {} suffix: {}, all_contents_df length: {}, all_contents_for_suffix length: {}".format(
                    e,
                    contig,
                    suffix,
                    len(all_contents_df),
                    reads_count_for_suffix
                    ))
                sys.exit(1)
            
        
        # Make a sub-subfolder to put the bams for this specific contig
        contig_folder = '{}/{}_{}/'.format(split_bams_folder, contig, suffix)
        if not os.path.exists(contig_folder):
            os.mkdir(contig_folder)
            
            
        bam_file_name = '{}/{}_{}.bam'.format(contig_folder, contig, suffix)
        
        # Write, sort and index bam immediately
        write_reads_to_file(reads_deduped, bam_file_name, header_string) 
        try:
            # print("\tSorting {}...".format(bam_file_name))
            sorted_bam_file_name = sort_bam(bam_file_name)
            # print("\tIndexing {}...".format(sorted_bam_file_name))
            index_bam(sorted_bam_file_name)
            rm_bam(bam_file_name)
        except Exception as e:
            print("Failed at indexing {}".format(bam_file_name))
            
    
def concat_and_write_bams_wrapper(params):
    contig, df_dict, header_string, split_bams_folder, barcode_tag, number_of_expected_bams, verbose = params
    
    # print("df_dict keys: {}".format(df_dict.keys()))
    concat_and_write_bams(contig, df_dict, header_string, split_bams_folder, barcode_tag=barcode_tag, number_of_expected_bams=number_of_expected_bams, verbose=verbose)
