#! /bin/bash

python marine.py --bam_filepath examples/data/sc_lr_subset.bam --output_folder examples/sc_lr_subset_CT --barcode_whitelist_file examples/data/sc_lr_barcodes.tsv.gz --barcode_tag "IB" --min_base_quality 0 --num_intervals_per_contig 5 --strandedness 2
