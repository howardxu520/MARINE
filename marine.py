import argparse
import os
import pandas as pd
import polars as pl
import psutil
import pysam
import sys
from sys import getsizeof
import time

sys.path.append('src/')

from read_process import get_contig_lengths_dict,\
incorporate_replaced_pos_info,incorporate_insertions_and_deletions,\
get_positions_from_md_tag,reverse_complement,get_edit_information,get_edit_information_wrapper,\
has_edits,get_total_coverage_for_contig_at_position,\
print_read_info, update_coverage_array, get_read_information, get_hamming_distance, remove_softclipped_bases,find

from utils import get_intervals, index_bam, write_rows_to_info_file, write_header_to_bam, \
write_read_to_bam_file, remove_file_if_exists, make_folder, concat_and_write_bams_wrapper

def pretty_print(contents):
    if type(contents) == list:
        for item in contents:
            pretty_print(item)
    else:
        to_write = '{}\n'.format(contents)
        sys.stdout.write(to_write)
    

def run(bam_filepath):
    pretty_print(["Running algorithm", "why not"])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MARINE')
    parser.add_argument('--bam_filepath', type=str, default='/projects/ps-yeolab5/ekofman/Sammi/MouseBrainEF1A_SingleCell_EPR_batch2/filtered_possorted_ms_hippo_stamp_bam/filtered_keep_xf25_possorted_genome_with_header.bam_MD.bam')

    args = parser.parse_args()
    bam_filepath = args.bam_filepath
    
    run(bam_filepath)