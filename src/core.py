import pysam
import os
import sys
from sys import getsizeof
import time
import pandas as pd
import polars as pl
from collections import defaultdict
import multiprocessing
from tqdm import tqdm

from read_process import incorporate_replaced_pos_info,incorporate_insertions_and_deletions,\
get_positions_from_md_tag,reverse_complement,get_edit_information,get_edit_information_wrapper,\
has_edits,get_total_coverage_for_contig_at_position,\
print_read_info, update_coverage_array, get_read_information, get_hamming_distance, remove_softclipped_bases,find

from utils import get_contig_lengths_dict, get_intervals, index_bam, write_rows_to_info_file, write_header_to_edit_info, \
write_read_to_bam_file, remove_file_if_exists, make_folder, concat_and_write_bams_wrapper, make_edit_finding_jobs

import os, psutil


def run_edit_identifier(bampath, output_folder, barcode_whitelist, num_intervals_per_contig, verbose=False):
    # Make subfolder in which to information about edits
    edit_info_subfolder = '{}/edit_info'.format(output_folder)
    make_folder(edit_info_subfolder)
    
    edit_finding_jobs = make_edit_finding_jobs(bampath, output_folder, barcode_whitelist, num_intervals_per_contig, verbose)
    print("{} total jobs".format(len(edit_finding_jobs)))

    # Performance statistics
    overall_total_reads = 0
    total_seconds_for_reads = {0: 1}
    
    # Result containers
    overall_label_to_list_of_contents = defaultdict(lambda:{})
    results = []
    
    start_time = time.perf_counter()
    with multiprocessing.Pool(processes=16) as p:
        max_ = len(edit_finding_jobs)
        with tqdm(total=max_) as pbar:
            for _ in p.imap_unordered(find_edits_and_split_bams_wrapper, edit_finding_jobs):
                pbar.update()

                overall_label_to_list_of_contents[_[0]][_[1]] =  _[2]
                results.append([_[3], _[4], _[5], _[6]])

                total_reads = _[3]
                total_time = time.perf_counter() - start_time

                overall_total_reads += total_reads

                total_seconds_for_reads[overall_total_reads] = total_time

    overall_time = time.perf_counter() - start_time 
    return overall_label_to_list_of_contents, results, overall_time, overall_total_reads, total_seconds_for_reads




def find_edits(bampath, contig, split_index, start, end, output_folder, barcode_whitelist=None, verbose=False):  
    edit_info_subfolder = '{}/edit_info'.format(output_folder)
        
    time_reporting = {}
    start_time = time.perf_counter()
    
    samfile = pysam.AlignmentFile(bampath, "rb")
        
    counts = defaultdict(lambda:defaultdict(lambda:0))
    total_reads = 0
    
    bam_handles_for_barcodes = {}
    read_lists_for_barcodes = defaultdict(lambda:[])
    
    reads_for_contig = samfile.fetch(contig, start, end, multiple_iterators=True)
    
    output_file = '{}/{}_{}_{}_{}_edit_info.tsv'.format(edit_info_subfolder, contig, split_index, start, end)
    remove_file_if_exists(output_file)

    with open(output_file, 'w') as f:        
        write_header_to_edit_info(f)

        for i, read in enumerate(reads_for_contig):
            total_reads += 1
            
            if total_reads % 1000 == 0:
                time_reporting[total_reads] = time.perf_counter() - start_time

            barcode = read.get_tag("CB")
            if barcode_whitelist:
                if barcode not in barcode_whitelist:
                    counts[contig]['Barcode Filtered'] += 1
                    continue

            verbose = False
            
            try:
                error_code, list_of_rows, num_edits_of_each_type = get_read_information(read, contig, verbose=verbose)
            except Exception as e:
                print("Failed on\n{}".format(read.to_string()))
                break
                
            if error_code:
                counts[contig][error_code] += 1
            else:
                counts[contig]['edited'] += 1
                write_rows_to_info_file(list_of_rows, f)
            
            # Store each read using its string representation
            read_as_string = read.to_string()
            read_tab_separated = read_as_string.split('\t')
     
            second_new_contig_section = '{}_{}'.format(contig, barcode)
            read_tab_separated[2] = second_new_contig_section
            
            read_as_string = '\t'.join(read_tab_separated)
            
            read_lists_for_barcodes[barcode].append(read_as_string)
            
    
    # Add all reads to dictionary for contig and barcode, in their string representation
    barcode_to_concatted_reads = {}
    for barcode, read_list in read_lists_for_barcodes.items():        
        
        # Concatenate the string representations of all reads for each bam-contig combination
        all_reads_concatted = '\n'.join(read_list)
            
        # Save this concatenated block of text to dictionary
        barcode_to_concatted_reads[barcode] = all_reads_concatted
        
    time_reporting[total_reads] = time.perf_counter() - start_time
    
    samfile.close()
    
    return barcode_to_concatted_reads, total_reads, counts, time_reporting




def find_edits_and_split_bams(bampath, contig, split_index, start, end, output_folder, barcode_whitelist=None, verbose=False):
    barcode_to_concatted_reads, total_reads, counts, time_reporting = find_edits(bampath, contig, split_index,
                                                                         start, end, output_folder, barcode_whitelist=barcode_whitelist, verbose=verbose)    
    return barcode_to_concatted_reads, total_reads, counts, time_reporting
    
    
    
    
def find_edits_and_split_bams_wrapper(parameters):
    try:
        start_time = time.perf_counter()
        bampath, contig, split_index, start, end, output_folder, barcode_whitelist, verbose = parameters
        label = '{}({}):{}-{}'.format(contig, split_index, start, end)

        #print("{} ({}):{}-{}\tfind_edits_and_split_bams".format(contig, split_index, start, end))
        barcode_to_concatted_reads, total_reads, counts, time_reporting = find_edits_and_split_bams(bampath, contig, split_index, start, end,                                                                                        
                                                                                                              output_folder, 
                                                                                                              barcode_whitelist=barcode_whitelist,
                                                                                                              verbose=False)
        counts_df = pd.DataFrame.from_dict(counts)
        time_df = pd.DataFrame.from_dict(time_reporting, orient='index')
        if len(barcode_to_concatted_reads) > 0:
            barcode_to_concatted_reads_pl = pl.from_dict(barcode_to_concatted_reads).transpose(include_header=True, header_name='barcode').rename({"column_0": "contents"})
        else:
            # No transposes are allowed on empty dataframes
            barcode_to_concatted_reads_pl = pl.from_dict(barcode_to_concatted_reads)
            
        total_time = time.perf_counter() - start_time
        return contig, label, barcode_to_concatted_reads_pl, total_reads, counts_df, time_df, total_time
    except Exception as e:
        print('Contig {}: {}'.format(label, e))
        return 0, pd.DataFrame(), label, pd.DataFrame()