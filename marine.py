import argparse
from multiprocessing import Pool
import multiprocessing
import numpy as np
import os
import pandas as pd
import polars as pl
import psutil
import pysam
import sys
from sys import getsizeof
import time
from tqdm import tqdm

sys.path.append('src/')

from read_process import incorporate_replaced_pos_info,incorporate_insertions_and_deletions,\
get_positions_from_md_tag,reverse_complement,get_edit_information,get_edit_information_wrapper,\
has_edits,get_total_coverage_for_contig_at_position,\
print_read_info, update_coverage_array, get_read_information, get_hamming_distance, \
remove_softclipped_bases,find

from utils import get_intervals, index_bam, write_rows_to_info_file, write_header_to_edit_info, \
write_read_to_bam_file, remove_file_if_exists, make_folder, concat_and_write_bams_wrapper, \
pretty_print, read_barcode_whitelist_file, get_contigs_that_need_bams_written

from core import run_edit_identifier, run_bam_reconfiguration, \
gather_edit_information_across_subcontigs, run_coverage_calculator

            
        
def edit_finder(bam_filepath, output_folder, barcode_whitelist, contigs=[], num_intervals_per_contig=16, 
                verbose=False):
    
    pretty_print("Each contig is being split into {} subsets...".format(num_intervals_per_contig))
    
    overall_label_to_list_of_contents, results, overall_time, overall_total_reads, \
    total_seconds_for_reads = run_edit_identifier(
        bam_filepath, 
        output_folder, 
        barcode_whitelist,
        contigs=contigs,
        num_intervals_per_contig=num_intervals_per_contig, 
        verbose=verbose
    )
    
    
    pretty_print(
        [
            "Reads processed:\t{}".format(overall_total_reads), 
            "Time to process reads in min:\t{}".format(round(overall_time/60, 5))
        ],
        style="-"
    )
    
    
    total_seconds_for_reads_df = pd.DataFrame.from_dict(total_seconds_for_reads, orient='index')
    total_seconds_for_reads_df.columns = ['seconds']
    total_seconds_for_reads_df['reads'] = total_seconds_for_reads_df.index
    total_seconds_for_reads_df.index = range(len(total_seconds_for_reads_df))
    
    
    return overall_label_to_list_of_contents, results, total_seconds_for_reads_df


def bam_processing(overall_label_to_list_of_contents, output_folder):
    split_bams_folder = '{}/split_bams'.format(output_folder)
    make_folder(split_bams_folder)
    contigs_to_generate_bams_for = get_contigs_that_need_bams_written(overall_label_to_list_of_contents, split_bams_folder)
    pretty_print("Will split and reconfigure the following contigs: {}".format(",".join(contigs_to_generate_bams_for)))
    
    
    # BAM Generation
    total_bam_generation_time = run_bam_reconfiguration(split_bams_folder, bam_filepath, overall_label_to_list_of_contents, contigs_to_generate_bams_for)
    return total_bam_generation_time
    
    
    
def coverage_processing(output_folder):
    edit_info_grouped_per_contig_combined = gather_edit_information_across_subcontigs(output_folder)
    results, total_time = run_coverage_calculator(edit_info_grouped_per_contig_combined, output_folder)
    return results, total_time

def run(bam_filepath, output_folder, contigs=[], num_intervals_per_contig=16, barcode_whitelist_file=None, verbose=False):
    logging_folder = "{}/metadata/".format(output_folder)
    make_folder(logging_folder)
    
    barcode_whitelist = read_barcode_whitelist_file(barcode_whitelist_file)
    
    # Edit identification
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pretty_print("Identifying edits", style="~")
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    overall_label_to_list_of_contents, results, total_seconds_for_reads_df = edit_finder(
        bam_filepath, 
        output_folder, 
        barcode_whitelist,
        contigs,
        num_intervals_per_contig,
        verbose
    )
    
    total_seconds_for_reads_df.to_csv("{}/edit_finder_timing.tsv".format(logging_folder), sep='\t')

    # Make a subfolder into which the split bams will be placed
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pretty_print("Contigs processed:\n\n\t{}".format(sorted(list(overall_label_to_list_of_contents.keys()))))
    pretty_print("Splitting and reconfiguring BAMs to optimize coverage calculations", style="~")
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    total_bam_generation_time = bam_processing(overall_label_to_list_of_contents, output_folder)
    
    
    # Coverage calculation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pretty_print("Total time to concat and write bams: {} minutes".format(round(total_bam_generation_time/60, 3)))
    pretty_print("Calculating coverage at edited sites", style='~')
    
    results, total_time = coverage_processing(output_folder)
    
    pretty_print("Total time to calculate coverage: {} minutes".format(round(total_time/60, 3)))
    all_edit_info = pd.concat(results)
    all_edit_info.to_csv('{}/all_edit_info.tsv'.format(output_folder), sep='\t')
    pretty_print("Done!", style="+")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MARINE')
    parser.add_argument('--bam_filepath', type=str, default='/projects/ps-yeolab3/ekofman/sailor2/data/groups_0_1_2_3_4_5_6_7_8_9_10_11_merged.bam')

    parser.add_argument('--output_folder', type=str, default='/projects/ps-yeolab3/ekofman/sailor2/scripts/script_test')
    
    parser.add_argument('--barcode_whitelist_file', type=str, default=None)
    
    parser.add_argument('--cores', type=int, default=multiprocessing.cpu_count())
    
    args = parser.parse_args()
    bam_filepath = args.bam_filepath
    output_folder = args.output_folder
    barcode_whitelist_file = args.barcode_whitelist_file
    cores = args.cores
    
    if cores is None:
        cores = 16
    pretty_print("Assuming {} cores available for multiprocessing".format(cores))
   
    # TEMP
    barcode_whitelist_file = '/projects/ps-yeolab3/ekofman/Sammi/MouseBrainEF1A_SingleCell_EPR_batch2/cellranger/results/ms_hippo_stamp_EIF4A_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    
    pretty_print(["Arguments:",
                  "\tBAM filepath:\t{}".format(bam_filepath), 
                  "\tOutput folder:\t{}".format(output_folder),
                  "\tBarcode whitelist:\t{}".format(barcode_whitelist_file),
                 ])
    
    run(bam_filepath, 
        output_folder, 
        contigs=['1'],
        barcode_whitelist_file=barcode_whitelist_file,
        num_intervals_per_contig=16
       )