python \
/tscc/projects/ps-yeolab3/ekofman/sailor2/marine.py \
--bam_filepath \
/tscc/projects/ps-yeolab3/ekofman/sailor2/examples/data/LR_single_cell.md.subset.filtered.sorted.bam \
--annotation_bedfile_path \
/tscc/projects/ps-yeolab3/ekofman/sailor2/annotations/cellranger-mm10-3.0.0.annotation.genes.bed \
--output_folder \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/singlecell_tests/long_read_sc_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
16 \
--barcode_whitelist_file /tscc/projects/ps-yeolab3/ekofman/sailor2/examples/data/sc_lr_barcodes.tsv.gz \
--strandedness 2 \
--barcode_tag "IB" \
--contigs "6" \
--num_intervals_per_contig 4