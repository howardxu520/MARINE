python \
/tscc/projects/ps-yeolab3/ekofman/sailor2/marine.py \
--bam_filepath \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/singlecell_tests/bams/9_3000526_only_5_cells.bam \
--annotation_bedfile_path \
/tscc/projects/ps-yeolab3/ekofman/sailor2/annotations/cellranger-GRCh38-3.0.0.annotation.genes.bed \
--output_folder \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/singlecell_tests/only_5_cells_test/ \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
1 \
--barcode_tag "CB" \
--verbose