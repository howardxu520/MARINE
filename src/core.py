import pysam
import os
from glob import glob
import sys
from sys import getsizeof
import time
import numpy as np
import pandas as pd
import polars as pl
from collections import defaultdict
import multiprocessing
from tqdm import tqdm
from multiprocessing import get_context
import shutil
import random
import gzip

from read_process import incorporate_replaced_pos_info,incorporate_insertions_and_deletions,\
get_positions_from_md_tag,reverse_complement,get_edit_information,get_edit_information_wrapper,\
has_edits,get_total_coverage_for_contig_at_position,\
print_read_info, get_read_information, get_hamming_distance, remove_softclipped_bases,find

from utils import get_contig_lengths_dict, get_intervals, index_bam, write_rows_to_info_file, write_header_to_edit_info, \
write_read_to_bam_file, remove_file_if_exists, make_folder, concat_and_write_bams_wrapper, make_edit_finding_jobs, pretty_print, get_contigs_that_need_bams_written, split_bed_file, \
get_coverage_wrapper, write_reads_to_file, sort_bam, rm_bam, suffixes, get_broken_up_contigs, run_command, \
make_depth_command_script_single_cell, concatenate_files, generate_and_run_bash_merge, read_barcode_whitelist_file

import os, psutil


def generate_depths(output_folder, bam_filepaths, original_bam_filepath, paired_end=False, barcode_tag=None, cores=1):
    
    coverage_start_time = time.perf_counter()

    all_depth_commands = []

    combine_edit_sites_command = (
        "echo 'concatenating bed file...';"
        "for file in {}/edit_info/*edit_info.tsv; do "
        "awk 'NR > 1 {{print $2, $4-1, $4}}' OFS='\t' \"$file\"; "
        "done | sort -k1,1 -k2,2n -u > {}/combined_source_cells.bed;"
    ).format(output_folder, output_folder)

    if not os.path.exists(f'{output_folder}/combined_source_cells.bed'):
        run_command(combine_edit_sites_command)
    
    all_depth_commands.append(combine_edit_sites_command)

    output_suffix = 'source_cells'
    
    if barcode_tag:
        coverage_subfolder = '{}/coverage'.format(output_folder)
        make_folder(coverage_subfolder)

        # Single cell mode
        split_bed_file(
            f"{output_folder}/combined_{output_suffix}.bed",
            f"{output_folder}/combined_{output_suffix}_split_by_suffix",
            bam_filepaths,
            output_suffix=output_suffix
        )
        
        make_depth_command_script_single_cell(paired_end, bam_filepaths, output_folder,
                                  all_depth_commands=all_depth_commands, output_suffix='source_cells', run=True, processes=cores, barcode_tag=barcode_tag)
        
    else:
        if paired_end:
            paired_end_flag = '-s '
        else:
            paired_end_flag = ''
            
        # Bulk mode, we will not split the bed and simply use samtools depth on the combined.bed
        samtools_depth_command = f"samtools depth {paired_end_flag}-a -b {output_folder}/combined_source_cells.bed {original_bam_filepath} > {output_folder}/depths_source_cells.txt"
        run_command(samtools_depth_command)
        

    print("Concatenating edit info files...")
    concatenate_files(output_folder, "edit_info/*edit_info.tsv",
                      "{}/final_edit_info_no_coverage.tsv".format(output_folder),
                     run=True)

    print("Append the depth columns to the concatenated final_edit_info file...")

    header_columns = ['barcode', 'contig', 'position', 'ref', 'alt',
                      'read_id', 'strand', 'coverage']


    generate_and_run_bash_merge(output_folder,
                                '{}/final_edit_info_no_coverage.tsv'.format(output_folder),
                            '{}/depths_source_cells.txt'.format(output_folder), 
                            '{}/final_edit_info.tsv'.format(output_folder), 
                                header_columns, barcode_tag=barcode_tag)
    
    coverage_total_time = time.perf_counter() - coverage_start_time
    
    total_seconds_for_contig_df = pd.DataFrame({'coverage_total_time': [coverage_total_time]})
    return coverage_total_time, total_seconds_for_contig_df


# We will write the reads that contain edits to bam
def bam_processing(bam_filepath, interval_metadata_list, output_folder, barcode_tag='CB', cores=1, number_of_expected_bams=4, verbose=False):
    split_bams_folder = '{}/split_bams'.format(output_folder)
    make_folder(split_bams_folder)

    contigs_from_metadata = list(set([m['contig'] for m in interval_metadata_list if m]))
    contigs_to_generate_bams_for = get_contigs_that_need_bams_written(
        contigs_from_metadata,
        split_bams_folder, 
        barcode_tag=barcode_tag,
        number_of_expected_bams=number_of_expected_bams
    )
    
    if verbose:
        pretty_print("Will split and reconfigure the following contigs: {}".format(",".join(contigs_to_generate_bams_for)))
    
    if len(contigs_to_generate_bams_for) > 0:
        print("Processing reads from disk files for BAM generation...")
        overall_label_to_list_of_contents = process_reads_from_disk_for_bam_creation(
            interval_metadata_list, output_folder, barcode_tag, verbose
        )
        
        filtered_contents = {k: v for k, v in overall_label_to_list_of_contents.items() 
                           if k in contigs_to_generate_bams_for}
        
        if filtered_contents:
            total_bam_generation_time, total_seconds_for_bams = run_bam_reconfiguration(
                split_bams_folder, 
                bam_filepath, 
                filtered_contents, 
                contigs_to_generate_bams_for, 
                barcode_tag=barcode_tag, 
                cores=cores,
                number_of_expected_bams=number_of_expected_bams,
                verbose=verbose)
        else:
            total_bam_generation_time = 0
            total_seconds_for_bams = {}
        
        del overall_label_to_list_of_contents
        
    else:
        total_bam_generation_time = 0
        total_seconds_for_bams = {}
        
    total_seconds_for_bams_df = pd.DataFrame.from_dict(total_seconds_for_bams, orient='index')
    total_seconds_for_bams_df.columns = ['seconds']
    total_seconds_for_bams_df['contigs'] = total_seconds_for_bams_df.index
    total_seconds_for_bams_df.index = range(len(total_seconds_for_bams_df))
    
    return total_bam_generation_time, total_seconds_for_bams_df


# Cleanup for intermediate reads
def cleanup_intermediate_reads(output_folder, verbose=False):
    reads_subfolder = f"{output_folder}/intermediate_reads"
    if os.path.exists(reads_subfolder):
        if verbose:
            print(f"Cleaning up intermediate read files in {reads_subfolder}")
        try:
            shutil.rmtree(reads_subfolder)
        except Exception as e:
            if verbose:
                print(f"Error cleaning up {reads_subfolder}: {e}")
    
# We will call run_edit_identifier to find edits on the contig
def edit_finder(bam_filepath, output_folder, strandedness, barcode_tag="CB", barcode_whitelist=None, contigs=[],
                verbose=False, cores=64, min_read_quality=0, min_base_quality=0, dist_from_end=0, interval_length=2000000):
    
    pretty_print("Each contig is being split into subsets of length...".format(interval_length))
    
    interval_metadata_list, results, overall_time, overall_total_reads, \
    total_seconds_for_reads, counts_summary_df = run_edit_identifier(
        bam_filepath, 
        output_folder, 
        strandedness=strandedness,
        barcode_tag=barcode_tag,
        barcode_whitelist=barcode_whitelist,
        contigs=contigs,
        verbose=verbose,
        cores=cores,
        min_read_quality=min_read_quality,
        min_base_quality=min_base_quality,
        dist_from_end=dist_from_end,
        interval_length=interval_length
    )
    
    pretty_print(
        [
            "Reads processed:\t{}".format(overall_total_reads), 
            "Time to process reads in min:\t{}".format(round(overall_time/60, 5)),
            "Read Summary:\n{}".format(counts_summary_df)
        ],
        style="-"
    )
    
    total_seconds_for_reads_df = pd.DataFrame.from_dict(total_seconds_for_reads, orient='index')
    total_seconds_for_reads_df.columns = ['seconds']
    total_seconds_for_reads_df['reads'] = total_seconds_for_reads_df.index
    total_seconds_for_reads_df.index = range(len(total_seconds_for_reads_df))
    
    return interval_metadata_list, results, total_seconds_for_reads_df, overall_total_reads, counts_summary_df

    
# We will split up contig according to interval length, and on each split interval, we will call find_edits_and_split_bams_wrapper identify edits separately, we will merge them after all edits are finished processing.
def run_edit_identifier(bampath, output_folder, strandedness, barcode_tag="CB", barcode_whitelist=None, contigs=[], verbose=False, cores=64, min_read_quality=0, min_base_quality=0, dist_from_end=0, interval_length=2000000):
    # Make subfolder in which to information about edits
    edit_info_subfolder = '{}/edit_info'.format(output_folder)
    make_folder(edit_info_subfolder)
    
    edit_finding_jobs = make_edit_finding_jobs(bampath, output_folder, strandedness, barcode_tag, barcode_whitelist, contigs, verbose, min_read_quality, min_base_quality, dist_from_end, interval_length=interval_length)
    pretty_print("{} total jobs".format(len(edit_finding_jobs)))
    
    # Performance statistics
    overall_total_reads = 0
    total_seconds_for_reads = {0: 1}
    # Metadata point to where the reads are
    interval_metadata_list = []
    results = []
    
    start_time = time.perf_counter()
    overall_counts_summary = defaultdict(lambda: 0)
    
    with get_context("spawn").Pool(processes=cores, maxtasksperchild=4) as p:
        max_ = len(edit_finding_jobs)
        with tqdm(total=max_) as pbar:
            for result in p.imap_unordered(find_edits_and_split_bams_wrapper, edit_finding_jobs):
                pbar.update()

                contig, label, interval_metadata, total_reads, counts_summary_df, time_df, total_time = result
                
                if interval_metadata and barcode_tag:
                    interval_metadata_list.append(interval_metadata)
                if not counts_summary_df.empty:
                    for contig_key in counts_summary_df.columns:
                        for count_type in counts_summary_df.index:
                            overall_counts_summary[count_type] += counts_summary_df.loc[count_type, contig_key]
                
                overall_total_reads += total_reads
                total_seconds_for_reads[overall_total_reads] = time.perf_counter() - start_time
                
                import gc
                gc.collect()

    overall_time = time.perf_counter() - start_time 
    
    overall_count_summary_df = pd.Series(overall_counts_summary)

    return interval_metadata_list, results, overall_time, overall_total_reads, total_seconds_for_reads, overall_count_summary_df


# # Read read files from disk
# def process_reads_from_disk_for_bam_creation(interval_metadata_list, output_folder, barcode_tag, verbose=False):

#     overall_label_to_list_of_contents = defaultdict(lambda: {})
#     contig_intervals = defaultdict(list)
#     for metadata in interval_metadata_list:
#         if metadata:
#             contig_intervals[metadata['contig']].append(metadata)
    
#     for contig, intervals in contig_intervals.items():
#         if verbose:
#             print(f"Processing {len(intervals)} intervals for contig {contig}")
        
#         all_chunk_files = []
#         for interval in intervals:
#             if 'chunk_files' in interval:
#                 all_chunk_files.extend(interval['chunk_files'])
        
#         if verbose:
#             print(f"Found {len(all_chunk_files)} chunk files for {contig}")
        
#         suffix_options = suffixes.get(barcode_tag, [])
        
#         for suffix in suffix_options:
#             reads_for_suffix = defaultdict(list)

#             for chunk_file in all_chunk_files:
#                 try:
#                     with open(chunk_file, 'r') as f:
#                         for line in f:
#                             line = line.strip()
#                             if line:
#                                 parts = line.split('\t', 1)
#                                 if len(parts) == 2:
#                                     barcode, read_string = parts

#                                     if barcode.endswith(suffix):
#                                         reads_for_suffix[barcode].append(read_string)
                
#                 except Exception as e:
#                     if verbose:
#                         print(f"Error reading chunk file {chunk_file}: {e}")
#                     continue
            
#             if reads_for_suffix:
#                 df_data = []
#                 for barcode, reads in reads_for_suffix.items():
#                     df_data.append({
#                         'barcode': barcode,
#                         'contents': '\n'.join(reads),
#                         'bucket': suffix_options.index(suffix) if suffix in suffix_options else 0
#                     })
                
#                 if df_data:
#                     df = pl.DataFrame(df_data)
#                     label = f"{contig}_{suffix}"
#                     overall_label_to_list_of_contents[contig][label] = df
                    
#                     if verbose:
#                         print(f"Created DataFrame for {contig}_{suffix} with {len(df)} barcodes, total reads: {sum(len(reads) for reads in reads_for_suffix.values())}")
                    
#                     del reads_for_suffix
#                     del df_data
    
#     return overall_label_to_list_of_contents




def process_reads_from_disk_for_bam_creation(interval_metadata_list, output_folder, barcode_tag, verbose=False):
    overall_label_to_list_of_contents = defaultdict(lambda: {})
    contig_intervals = defaultdict(list)
    for metadata in interval_metadata_list:
        if metadata:
            contig_intervals[metadata['contig']].append(metadata)
    
    for contig, intervals in contig_intervals.items():
        if verbose:
            print(f"Processing {len(intervals)} intervals for contig {contig}")
        
        all_chunk_files = []
        for interval in intervals:
            if 'chunk_files' in interval:
                all_chunk_files.extend(interval['chunk_files'])
        
        if verbose:
            print(f"Found {len(all_chunk_files)} chunk files for {contig}")
        
        suffix_options = suffixes.get(barcode_tag, [])
        
        for suffix in suffix_options:
            reads_for_suffix = defaultdict(list)

            for chunk_file in all_chunk_files:
                try:
                    with open(chunk_file, 'r') as f:
                        for line in f:
                            line = line.strip()
                            if line:
                                parts = line.split('\t', 1)
                                if len(parts) == 2:
                                    barcode, read_string = parts

                                    if barcode.endswith(suffix):
                                        reads_for_suffix[barcode].append(read_string)
                
                except Exception as e:
                    if verbose:
                        print(f"Error reading chunk file {chunk_file}: {e}")
                    continue
            
            if reads_for_suffix:
                df_data = []
                for barcode, reads in reads_for_suffix.items():
                    joined_reads = '\n'.join(reads)
                    compressed_reads = gzip.compress(joined_reads.encode('utf-8'))
                    
                    df_data.append({
                        'barcode': barcode,
                        'contents_compressed': compressed_reads,
                        'read_count': len(reads), 
                        'bucket': suffix_options.index(suffix) if suffix in suffix_options else 0
                    })
                
                if df_data:
                    df = pl.DataFrame(df_data)
                    label = f"{contig}_{suffix}"
                    overall_label_to_list_of_contents[contig][label] = df
                    
                    if verbose:
                        total_reads = sum(row['read_count'] for row in df_data)
                        original_size = sum(len('\n'.join(reads).encode('utf-8')) for reads in reads_for_suffix.values())
                        compressed_size = sum(len(row['contents_compressed']) for row in df_data)
                        compression_ratio = compressed_size / original_size if original_size > 0 else 0
                        print(f"Created DataFrame for {contig}_{suffix} with {len(df)} barcodes, "
                              f"total reads: {total_reads}, compression: {compression_ratio:.2%}")
                    
                    del reads_for_suffix
                    del df_data
    
    return overall_label_to_list_of_contents


def run_bam_reconfiguration(split_bams_folder, bampath, overall_label_to_list_of_contents, contigs_to_generate_bams_for, barcode_tag='CB', cores=1, number_of_expected_bams=4, verbose=False):
    start_time = time.perf_counter()
    # Get the bam header, which will be used for each of the split bams too
    with pysam.AlignmentFile(bampath, "rb") as samfile:
        header_string = str(samfile.header)
    # num_processes = np.max([len(contigs_to_generate_bams_for), 32])
    num_processes = cores
    
    total_seconds_for_bams = {0: 1}
    total_bams = 0
    with get_context("spawn").Pool(processes=num_processes, maxtasksperchild=4) as p:
        max_ = len(contigs_to_generate_bams_for)
        with tqdm(total=max_) as pbar:
            for _ in p.imap_unordered(concat_and_write_bams_wrapper, [[i[0], i[1], header_string, split_bams_folder, 
                                                                       barcode_tag, number_of_expected_bams, verbose] for i in overall_label_to_list_of_contents.items() if i[0] in contigs_to_generate_bams_for]):
                pbar.update()
                
                total_bams += 1
                total_time = time.perf_counter() - start_time
                total_seconds_for_bams[total_bams] = total_time

    total_bam_generation_time = time.perf_counter() - start_time
    return total_bam_generation_time, total_seconds_for_bams


# This is the main function for edit finding that will be called in the marine.py
def run_edit_finding(barcode_tag,
                     barcode_whitelist_file, 
                     contigs, 
                     num_per_sublist,
                     bam_filepath, 
                     output_folder, 
                     strandedness,
                     min_read_quality,
                     min_base_quality,
                     min_dist_from_end,
                     interval_length,
                     number_of_expected_bams,
                     cores,
                     logging_folder,
                     verbose=False
                    ):
    overall_total_reads_processed = 0
    if barcode_whitelist_file:
        barcode_whitelist = read_barcode_whitelist_file(barcode_whitelist_file)
    else:
        barcode_whitelist = None

    if len(contigs) == 0:
        # Take care of the case where no contigs are specified, so that all contigs available are processed
        broken_up_contigs = [[]]
    else:
        if barcode_tag:
            # For single cell sequencing we will only process this many contigs at a time
            broken_up_contigs = get_broken_up_contigs(contigs, num_per_sublist)
        else:
            # For bulk sequencing we will just process all contigs 
            broken_up_contigs = [contigs]

    print('Contig groups to be processed:', broken_up_contigs)
    
    overall_counts_summary_df = defaultdict(lambda:0)
    overall_total_reads_processed = 0
    all_interval_metadata = []
    
    for subcontig_list in broken_up_contigs:
            
        interval_metadata_list, results, total_seconds_for_reads_df, total_reads_processed, counts_summary_df = edit_finder(
            bam_filepath, 
            output_folder, 
            strandedness,
            barcode_tag,
            barcode_whitelist,
            subcontig_list,
            verbose,
            cores=cores,
            min_read_quality=min_read_quality,
            min_base_quality=min_base_quality,
            dist_from_end=min_dist_from_end,
            interval_length=interval_length
        )

        all_interval_metadata.extend(interval_metadata_list)
        
        for k,v in counts_summary_df.items():
            overall_counts_summary_df[k] += v
            
        overall_total_reads_processed += total_reads_processed
        
        if barcode_tag:
            # Make a subfolder into which the split bams will be placed
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            pretty_print("Splitting and reconfiguring BAMs to optimize coverage calculations", style="~")
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            total_bam_generation_time, total_seconds_for_bams_df = bam_processing(
                bam_filepath, 
                all_interval_metadata,
                output_folder, 
                barcode_tag=barcode_tag, 
                cores=cores, 
                number_of_expected_bams=number_of_expected_bams, 
                verbose=verbose
            )
            
            pretty_print("Total time to concat and write bams: {} minutes".format(round(total_bam_generation_time/60, 3)))

    cleanup_intermediate_reads(output_folder, verbose)
    
    with open('{}/manifest.txt'.format(logging_folder), 'a+') as f:
        f.write(f'total_reads_processed\t{overall_total_reads_processed}\n') 
        for k, v in overall_counts_summary_df.items():
            f.write(f'{k}\t{v}\n') 

        if overall_total_reads_processed > 0:
            raw_epr = overall_counts_summary_df.get("total_edits")/overall_total_reads_processed
        else:
            raw_epr = 0
            
        f.write(f'edits per read (EPR)\t{raw_epr}\n')


# Add contig barcode info into read
def incorporate_barcode(read_as_string, contig, barcode):
    read_tab_separated = read_as_string.split('\t')
    contig_section = '{}_{}'.format(contig, barcode)
    read_tab_separated[2] = contig_section
    read_as_string = '\t'.join(read_tab_separated)
    return read_as_string


# Write reads that have edits to bam file
def write_bam_file(reads, bam_file_name, header_string):
    # Write, sort and index bam immediately
    write_reads_to_file(reads, bam_file_name, header_string)
    try:
        sorted_bam_file_name = sort_bam(bam_file_name)
        index_bam(sorted_bam_file_name)
        rm_bam(bam_file_name)
    except Exception as e:
        print("Failed at indexing {}, {}".format(bam_file_name, e))
    
    
# Look for edits and save the reads that have edits to intermediate reads files.
def find_edits(bampath, contig, split_index, start, end, output_folder, barcode_tag="CB", strandedness=0, barcode_whitelist=None, verbose=False, min_read_quality=0, min_base_quality=0, dist_from_end=0):  
    edit_info_subfolder = '{}/edit_info'.format(output_folder)
    
    reads_subfolder = '{}/intermediate_reads'.format(output_folder)
    make_folder(reads_subfolder)
        
    time_reporting = {}
    start_time = time.perf_counter()
    
    samfile = pysam.AlignmentFile(bampath, "rb")
    
    counts = defaultdict(lambda:defaultdict(lambda:0))
    counts[contig]['total_edits'] = 0
    counts[contig]['edited'] = 0
    
    total_reads = 0
    
    MAX_READS_PER_FILE = 5000000
    current_chunk = 0
    current_chunk_reads = 0
    current_chunk_file = None
    reads_written_per_barcode = defaultdict(int)
    chunk_files_created = []
    
    def get_next_chunk_file():
        nonlocal current_chunk, current_chunk_reads, current_chunk_file
        
        if current_chunk_file:
            current_chunk_file.close()
        
        chunk_filename = f"{reads_subfolder}/{contig}_{split_index}_chunk_{current_chunk:04d}.reads"
        current_chunk_file = open(chunk_filename, 'w')
        chunk_files_created.append(chunk_filename)
        current_chunk_reads = 0
        current_chunk += 1
        
        if verbose:
            print(f"Created new chunk file: {chunk_filename}")
        
        return current_chunk_file
    
    reads_for_contig = samfile.fetch(contig, start, end, multiple_iterators=True)
        
    output_file = '{}/{}_{}_{}_{}_edit_info.tsv'.format(edit_info_subfolder, contig, split_index, start, end)
    
    remove_file_if_exists(output_file)
    
    with open(output_file, 'w') as f:        
        write_header_to_edit_info(f)

        for i, read in enumerate(reads_for_contig):
            total_reads += 1
            
            if total_reads % 1000 == 0:
                time_reporting[total_reads] = time.perf_counter() - start_time

            if barcode_tag is None:
                barcode = 'no_barcode'
            elif read.has_tag(barcode_tag):
                barcode = read.get_tag(barcode_tag)
            else:
                barcode = 'missing'
                
            if barcode_whitelist and barcode != 'no_barcode':
                if barcode_tag == "IB" and read.has_tag("CB"):
                    # Single-cell long reads
                    cell_barcode = read.get_tag("CB")
                else:
                    cell_barcode = barcode
                    
                if cell_barcode not in barcode_whitelist:
                    counts[contig]['Barcode Filtered'] += 1
                    continue
            
            try:
                error_code, list_of_rows, num_edits_of_each_type = get_read_information(read, contig, strandedness=strandedness, barcode_tag=barcode_tag, verbose=verbose, min_read_quality=min_read_quality, min_base_quality=min_base_quality, dist_from_end=dist_from_end)
            except Exception as e:
                print("Failed getting read info on\n{}, {}".format(read.to_string(), e))
                break
                
            if error_code:
                counts[contig][error_code] += 1
            else:
                counts[contig]['edited'] += 1
                counts[contig]['total_edits'] += len(list_of_rows)
                write_rows_to_info_file(list_of_rows, f)
            
            # Save each read to intermediate reads file only if there was not an mapq_low error code while processing the read
            if (error_code != 'mapq_low') and (barcode_tag):
                read_to_string = read.to_string()
                # Add the contig barcode information into the read
                read_as_string = incorporate_barcode(read_to_string, contig, barcode)
                
                if current_chunk_file is None or current_chunk_reads >= MAX_READS_PER_FILE:
                    get_next_chunk_file()
                 
                # Write reads to intermediate read chunk file
                current_chunk_file.write(f"{barcode}\t{read_as_string}\n")
                current_chunk_reads += 1
                reads_written_per_barcode[barcode] += 1
                    
            elif error_code != 'mapq_low':
                split_bams_subfolder = '{}/split_bams'.format(output_folder)
                make_folder(split_bams_subfolder)
                make_folder('{}/{}'.format(split_bams_subfolder, contig))

    if current_chunk_file:
        current_chunk_file.close()
        
    time_reporting[total_reads] = time.perf_counter() - start_time
    samfile.close()
    
    if verbose:
        print(f"Created {len(chunk_files_created)} chunk files for {contig}_{split_index}")
    
    interval_metadata = {
        'contig': contig,
        'split_index': split_index,
        'start': start,
        'end': end,
        'reads_subfolder': reads_subfolder,
        'chunk_files': chunk_files_created,
        'total_chunk_files': len(chunk_files_created),
        'reads_written_per_barcode': dict(reads_written_per_barcode)
    }
    
    return interval_metadata, total_reads, counts, time_reporting
    

# This function calls find edits to look for edit through reads
def find_edits_and_split_bams_wrapper(parameters):
    try:
        start_time = time.perf_counter()
        bampath, contig, split_index, start, end, output_folder, strandedness, barcode_tag, barcode_whitelist, verbose, min_read_quality, min_base_quality, dist_from_end = parameters
        label = '{}({}):{}-{}'.format(contig, split_index, start, end)
        
        interval_metadata, total_reads, counts, time_reporting = find_edits(
            bampath, contig, split_index, start, end, output_folder, 
            barcode_tag=barcode_tag, strandedness=strandedness,
            barcode_whitelist=barcode_whitelist, verbose=verbose,
            min_read_quality=min_read_quality, min_base_quality=min_base_quality,
            dist_from_end=dist_from_end
        )
        
        counts_df = pd.DataFrame.from_dict(counts)
        
        if verbose:
            print("{}:{}, total reads: {}, reads written to files".format(contig, split_index, total_reads))
        
        time_df = pd.DataFrame.from_dict(time_reporting, orient='index')
        total_time = time.perf_counter() - start_time
        
        return contig, label, interval_metadata, total_reads, counts_df, time_df, total_time
        
    except Exception as e:
        print('Contig {}: {}'.format(label, e))
        return 0, label, None, 0, pd.DataFrame(), pd.DataFrame(), 0
    
    
    
def run_coverage_calculator(edit_info_grouped_per_contig_combined, 
                            output_folder, 
                            barcode_tag='CB',
                            paired_end=False, 
                            verbose=False,
                            processes=16
                           ):
    coverage_counting_job_params = get_job_params_for_coverage_for_edits_in_contig(
        edit_info_grouped_per_contig_combined, 
        output_folder,
        barcode_tag=barcode_tag,
        paired_end=paired_end,
        verbose=verbose
    )
    
    start_time = time.perf_counter()

    total_seconds_for_contig = {0: 1}
    total_contigs = 0
    
    results = []
    # Spawn has to be used instead of the default fork when using the polars library
    with get_context("spawn").Pool(processes=processes) as p:
        max_ = len(coverage_counting_job_params)
        with tqdm(total=max_) as pbar:
            for _ in p.imap_unordered(get_coverage_wrapper, coverage_counting_job_params):
                pbar.update()

                #print(_.columns)
                results.append(_)

                total_contigs += 1
                total_time = time.perf_counter() - start_time
                total_seconds_for_contig[total_contigs] = total_time
                
    total_time = time.perf_counter() - start_time
    return results, total_time, total_seconds_for_contig


def get_job_params_for_coverage_for_edits_in_contig(edit_info_grouped_per_contig_combined, output_folder,
                                                    barcode_tag='CB', paired_end=False, verbose=False):
    job_params = []
    
    for contig, edit_info in edit_info_grouped_per_contig_combined.items():
                    
        job_params.append([edit_info, contig, output_folder, barcode_tag, paired_end, verbose])  
        
    return job_params

    
def gather_edit_information_across_subcontigs(output_folder, barcode_tag='CB', number_of_expected_bams=4):
    
    splits = [i.split("/")[-1].split('_edit')[0] for i in sorted(glob('{}/edit_info/*'.format(output_folder)))]

    all_edit_info_for_barcodes = []

    edit_info_grouped_per_contig = defaultdict(lambda:[])
    edit_info_grouped_per_contig_combined = defaultdict(lambda:[])

    num_splits = len(splits)
    # print("Grouping edit information outputs by contig...")
    if num_splits > 500:
        interval = 400
    else:
        interval = 10
    
    for i, split in enumerate(splits):
        if i%interval == 0:
            print("\tsplit {}, {}/{}...".format(split, i, num_splits))

        contig = split.split("_")[0]

        edit_info_file = '{}/edit_info/{}_edit_info.tsv'.format(output_folder, split)
        edit_info_df = pd.read_csv(edit_info_file, sep='\t')
        edit_info_df['contig'] = edit_info_df['contig'].astype(str)
        edit_info_df['position'] = edit_info_df['position'].astype(int)

        edit_info = pl.from_pandas(edit_info_df) 
        
        if barcode_tag in ['CB', 'IS', 'IB']:
            # For single-cell or long read data we have to split the edit infos for each contig additionally by cell barcode 
            # or isoform ID in order to make the coverage calculation efficient per cell and/or per isoform
            suffix_options = suffixes.get(barcode_tag)
            for suffix in suffix_options:
                edit_info_subset = edit_info.filter(pl.col("barcode").str.ends_with(suffix))
                edit_info_grouped_per_contig["{}_{}".format(contig, suffix)].append(edit_info_subset)

        else:
            # For bulk data we don't have to do any filtering, the edits calculated for each split interval
            # can be paired directly with the split bam for that interval
            edit_info_grouped_per_contig["{}".format(split)].append(edit_info)
            
        del edit_info_df

    print("Done grouping! Concatenating ...")
    for contig, list_of_edit_info_dfs in edit_info_grouped_per_contig.items():
        edit_info_grouped_per_contig_combined[contig] = pl.concat(list_of_edit_info_dfs)
    print("Done concatenating!")
        
    return edit_info_grouped_per_contig_combined



def add_site_id(all_edit_info):
    if len(all_edit_info) == 0:
        return all_edit_info

    return all_edit_info.with_columns(
        pl.concat_str(
            [
                pl.col("barcode"),
                pl.col("contig"),
                pl.col("position"),
                pl.col("ref"),
                pl.col("alt"),
                pl.col("strand")
            ],
            separator="_",
        ).alias("site_id"))



def get_count_and_coverage_per_site(all_edit_info, skip_coverage=False):
    num_edits_df = all_edit_info.group_by("site_id").count()

    if not skip_coverage:
        coverage_df = all_edit_info.group_by("site_id").agg(pl.col("coverage").max())
        return num_edits_df.join(coverage_df, on='site_id')
    else:
        return num_edits_df



def generate_site_level_information(all_edit_info, skip_coverage=False):
    number_of_edits = len(all_edit_info)
    
    all_edit_info = add_site_id(all_edit_info)
    
    if len(all_edit_info) == 0:
        print("DataFrame is empty. Returning an empty DataFrame.")
        return all_edit_info

    identifying_values = ["site_id", "barcode", "contig", "position", "ref", "alt", "strand"]

    unique_rows_df = all_edit_info.select(identifying_values).unique(subset=identifying_values, maintain_order=True)

    count_and_coverage_at_site_df = get_count_and_coverage_per_site(all_edit_info, skip_coverage=skip_coverage)
    final_site_information_df = unique_rows_df.join(count_and_coverage_at_site_df, on='site_id')

    final_site_information_df = final_site_information_df.with_columns(
        pl.concat_str(
            [
                pl.col("ref"),
                pl.col("alt")
            ],
            separator=">",
        ).alias("conversion"))

    if not skip_coverage:
        number_of_sites = len(all_edit_info.group_by("site_id").agg(pl.col("coverage").unique()))
        assert(len(final_site_information_df) == number_of_sites)
    

    return final_site_information_df