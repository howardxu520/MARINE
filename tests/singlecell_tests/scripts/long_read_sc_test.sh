python \
/tscc/projects/ps-yeolab3/ekofman/sailor2/marine.py \
--bam_filepath \
/tscc/projects/ps-yeolab3/ekofman/sailor2/examples/data/LR_single_cell.md.subset.filtered.sorted.bam \
--output_folder \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/singlecell_tests/long_read_sc_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
16 \
--barcode_whitelist_file /tscc/projects/ps-yeolab3/ekofman/sailor2/examples/data/sc_lr_barcodes.tsv.gz \
--barcode_tag "IB" \
--contigs "6" \
--num_intervals_per_contig 4