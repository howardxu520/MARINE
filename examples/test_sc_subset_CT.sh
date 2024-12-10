#! /bin/bash

$MARINE/marine.py --bam_filepath $MARINE/examples/data/single_cell_CT.md.subset.bam --output_folder $MARINE/examples/sc_subset_CT --barcode_whitelist_file $MARINE/examples/data/sc_barcodes.tsv.gz --barcode_tag "CB" --num_intervals_per_contig 16 --strandedness 2 --contigs "1,2,3,4,5,6" --min_base_quality 15 --all_cells_coverage
