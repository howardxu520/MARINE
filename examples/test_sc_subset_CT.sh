#! /bin/bash

python marine.py --bam_filepath examples/data/single_cell_CT.md.subset.bam --output_folder examples/sc_subset_CT --barcode_whitelist_file examples/sc_barcodes.tsv.gz --barcode_tag "CB" --num_intervals_per_contig 16
