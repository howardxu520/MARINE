import math
import os
import pysam

def get_intervals(contig, contig_lengths_dict, num_intervals=4):
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
    f.write('barcode\tcontig\tposition\tref\talt\tread_id\tstrand\tdist_from_end\tbase_quality\n')
    
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