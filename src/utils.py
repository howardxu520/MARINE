import math
import os
import pysam
import polars as pl
import pandas as pd
import numpy as np

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
    if not os.path.exists('{}.bai'.format(barcode_bam_file_path)):
        pysam.index(barcode_bam_file_path)
        
def write_rows_to_info_file(list_of_rows, f):
    for info_list in list_of_rows:
        info_line = '\t'.join(info_list) + '\n'
        f.write(info_line)
        
def write_header_to_bam(f):
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
        
        
def get_edit_info_for_barcode_in_contig(edit_info, contig, barcode, output_folder):
    
    bam_subfolder = "{}/split_bams/{}".format(output_folder, contig)
    barcode_bam = '{}/{}_{}.bam'.format(bam_subfolder, contig, barcode)

    samfile_for_barcode = pysam.AlignmentFile(barcode_bam, "rb")

    edit_info_for_barcode = edit_info.filter(pl.col("barcode") == barcode)
    
    positions_for_barcode = set(edit_info_for_barcode["position"].unique())

    coverage_dict = {}
    for pos in positions_for_barcode:
        coverage_at_pos = np.sum(samfile_for_barcode.count_coverage(contig, pos-1, pos, quality_threshold=0))
        coverage_dict[pos] = coverage_at_pos

    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index')
    return edit_info_for_barcode.to_pandas(), coverage_df


def get_edit_info_for_barcode_in_contig_wrapper(parameters):
    edit_info, contig, barcode, output_folder = parameters
    edit_info_for_barcode, coverage_df = get_edit_info_for_barcode_in_contig(edit_info, contig, barcode, output_folder)
    edit_info_for_barcode['contig'] = edit_info_for_barcode.contig.astype(str)
    coverage_df['contig'] = edit_info_for_barcode.contig.astype(str)
    
    return edit_info_for_barcode, coverage_df