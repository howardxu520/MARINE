#! /bin/bash

marine.py --bam_filepath $MARINE/examples/data/sc_lr_subset.bam --output_folder $MARINE/examples/sc_lr_subset_CT --barcode_whitelist_file $MARINE/examples/data/sc_lr_barcodes.tsv.gz --barcode_tag "IB" --min_base_quality 0 --num_intervals_per_contig 5 --strandedness 2
