import math
from glob import glob
import os
import pysam
import polars as pl
import pandas as pd
import numpy as np
import sys
from collections import OrderedDict, defaultdict


def get_contigs_that_need_bams_written(expected_contigs, split_bams_folder, barcode_tag='CB'):
    bam_indices_written = [f.split('/')[-1].split('.bam')[0] for f in glob('{}/*/*.sorted.bam.bai'.format(split_bams_folder))]

    subsets_per_contig = defaultdict(lambda:0)
    for bam_index_written in bam_indices_written:
        contig_label = bam_index_written.split('_')[0]
        subsets_per_contig[contig_label] += 1

    if barcode_tag == 'CB':
        number_of_expected_bams = 4
    else:
        number_of_expected_bams = 1
        
    contigs_to_write_bams_for = []
    for c in expected_contigs:
        num_written_indices = subsets_per_contig.get(c, 0)
        if num_written_indices < number_of_expected_bams:
            print("Contig {} has {}/{} bams generated".format(c, num_written_indices, number_of_expected_bams))
            contigs_to_write_bams_for.append(c)
    
    return contigs_to_write_bams_for


def make_edit_finding_jobs(bampath, output_folder, reverse_stranded=True, barcode_tag="CB", barcode_whitelist=None, contigs=[], num_intervals_per_contig=16, verbose=False):
    jobs = []
    
    samfile = pysam.AlignmentFile(bampath, "rb")
    contig_lengths_dict = get_contig_lengths_dict(samfile)
    
    if len(contigs) == 0:
        contigs_to_use = set(contig_lengths_dict.keys())
    else:
        contigs_to_use = set(contigs)
        
    for contig in contig_lengths_dict.keys():
        # Skip useless contigs
        if len(contig) > 5 or contig == 'Stamp':# or contig != '17':
            continue
            
        if contig not in contigs_to_use:
            continue

        pretty_print("\tContig {}".format(contig))
        
        contig_length = contig_lengths_dict.get(contig)
        intervals_for_contig = get_intervals(contig, contig_lengths_dict, num_intervals_per_contig)

        # Set up for pool
        for split_index, interval in enumerate(intervals_for_contig):
            split_index = str(split_index).zfill(3)
            parameters = [bampath, contig, split_index, interval[0], interval[1], output_folder, reverse_stranded, barcode_tag, barcode_whitelist, verbose]
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
        os.mkdir(folder_path)

def get_coverage_df(edit_info, contig, output_folder, barcode_tag='CB'):
    
    bam_subfolder = "{}/split_bams/{}".format(output_folder, contig)
    contig_bam = '{}/{}.bam.sorted.bam'.format(bam_subfolder, contig)

    print("Contig {}. Loading {} bamfile...".format(contig, contig_bam))
    index_bam(contig_bam)
    samfile_for_barcode = pysam.AlignmentFile(contig_bam, "rb")

    print("Contig {}. Loaded bamfile...".format(contig))
    unique_barcodes = sorted(list(edit_info.unique("barcode")["barcode"]))
    
    coverage_dict = {}
    
    print("Contig {}. Iterating through barcodes...".format(contig))
    for i, barcode in enumerate(unique_barcodes):
        if i % 300 == 0:
            print("{}/{} barcodes for {}...".format(i+1, len(unique_barcodes), contig))
            
        edit_info_for_barcode = edit_info.filter(pl.col("barcode") == barcode)
        positions_for_barcode = set(edit_info_for_barcode["position"].unique())
        
        num_positions = len(positions_for_barcode)
        if not barcode_tag:
            print("{}:\tTotal edit site positions: {}".format(contig, num_positions))
            
        for i, pos in enumerate(positions_for_barcode):
            
            if not barcode_tag:
                if i%10000 == 0:
                    print("\t{}:\tCoverage computed at {}/{} positions...".format(contig, i, num_positions))
            
            if barcode_tag:
                barcode_specific_contig = '{}_{}'.format(contig, barcode)
                # Alter from 3_C_AAACCCAAGAACTTCC-1, for example, to 3_AAACCCAAGAACTTCC-1'
                barcode_specific_contig_split = barcode_specific_contig.split("_")
                barcode_specific_contig_without_subdivision = "{}_{}".format(barcode_specific_contig_split[0], barcode_specific_contig_split[2])

                coverage_at_pos = np.sum(samfile_for_barcode.count_coverage(barcode_specific_contig_without_subdivision, pos-1, pos, quality_threshold=0))
                coverage_dict['{}:{}'.format(barcode, pos)] = coverage_at_pos
            else:
                # For bulk, no barcodes, we will just have for example 19_no_barcode to convert to 19
                just_contig = contig.split('_')[0]
                
                coverage_at_pos = np.sum(samfile_for_barcode.count_coverage(just_contig, pos-1, pos, quality_threshold=0))
                coverage_dict['{}:{}'.format(contig, pos)] = coverage_at_pos
            
                
    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index')
    return coverage_df


def get_edit_info_for_barcode_in_contig_wrapper(parameters):
    edit_info, contig, output_folder, barcode_tag = parameters
    coverage_df = get_coverage_df(edit_info, contig, output_folder, barcode_tag=barcode_tag)
    coverage_df.columns = ['coverage']
    
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
        
    return edit_info_df.join(coverage_df)


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
        
    print("\tCurrent header length for {}: {}".format(bam_file_name, len(lengths_for_sn)))
    
    all_barcodes_for_contig = set([r.split('\t')[2] for r in reads])
    print("\tNum barcodes for {}: {}".format(bam_file_name, len(all_barcodes_for_contig)))
        
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
    
    print("\tNew header length: {}".format(len(header_dict.get("SQ"))))
    
    new_header = pysam.AlignmentHeader.from_dict(header_dict)
    
    num_reads = len(reads)
    
    with pysam.AlignmentFile(bam_file_name, "wb", text=str(new_header)) as bam_handle:
        for i, read_str in enumerate(reads):
            if i % 100000 == 0:
                print('file {}: {}/{} reads'.format(bam_file_name.split('/')[-1], i, num_reads))
                
            try:
                read = pysam.AlignedSegment.fromstring(read_str, new_header)
                bam_handle.write(read) 
            except Exception as e:
                print('{}\n\nfile {}: Failed to write read with str representation of:\n\t {}'.format(e,
                                                                                                      bam_file_name.split('/')[-1],
                                                                                                read_str))
                sys.exit(1)
                
            
            
    bam_handle.close()
        
        
        
def concat_and_write_bams(contig, df_dict, header_string, split_bams_folder, barcode_tag='CB'):
    job_params = []
    
    # Sort the subcontig regions such that the reads are properly ordered 
    sorted_subcontig_names = sorted(df_dict.keys())
    sorted_subcontig_dfs = []
    for n in sorted_subcontig_names:
        subcontig = df_dict.get(n)
        if len(subcontig) > 0:
            sorted_subcontig_dfs.append(subcontig)
        
    if len(sorted_subcontig_dfs) == 0:
        print("Empty")
        return []
    
    print("\t{}: num subcontigs to concat: {}".format(contig, len(sorted_subcontig_dfs)))
    # All of the reads for all of the barcodes are in this dataframe
    print("\t{}: concatting".format(contig))
    all_contents_df = pl.concat(sorted_subcontig_dfs)
        
    # Combine the reads (in string representation) for all rows corresponding to a barcode  
    if barcode_tag:
        suffix_options = ["A-1", "C-1", "G-1", "T-1"]
    else:
        suffix_options = ['no_barcode']
        
    for suffix in suffix_options:
        all_contents_for_suffix = all_contents_df.filter(pl.col('barcode').str.ends_with(suffix))
        
        reads_deduped = list(OrderedDict.fromkeys(all_contents_for_suffix.transpose().with_columns(
            pl.concat_str(
                [pl.col(c) for c in all_contents_for_suffix.transpose().columns],
                separator="\n"
                 ).alias("combined_text")
        )[['combined_text']][1].item().split('\n')))
        
        # Make a sub-subfolder to put the bams for this specific contig
        contig_folder = '{}/{}_{}/'.format(split_bams_folder, contig, suffix)
        if not os.path.exists(contig_folder):
            os.mkdir(contig_folder)
            
            
        bam_file_name = '{}/{}_{}.bam'.format(contig_folder, contig, suffix)
        
        # Write, sort and index bam immediately
        write_reads_to_file(reads_deduped, bam_file_name, header_string)
        try:
            print("\tSorting {}...".format(bam_file_name))
            sorted_bam_file_name = sort_bam(bam_file_name)
            print("\tIndexing {}...".format(sorted_bam_file_name))
            index_bam(sorted_bam_file_name)
            rm_bam(bam_file_name)
        except Exception as e:
            print("Failed at indexing {}".format(bam_file_name))
            
    
def concat_and_write_bams_wrapper(params):
    contig, df_dict, header_string, split_bams_folder, barcode_tag = params
    concat_and_write_bams(contig, df_dict, header_string, split_bams_folder, barcode_tag=barcode_tag)
