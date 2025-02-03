mypython=$1

$mypython $MARINE/marine.py \
--bam_filepath \
$MARINE/tests/singlecell_tests/bams/9_3000526_only_5_cells.bam \
--annotation_bedfile_path \
$MARINE/annotations/cellranger-mm10-3.0.0.annotation.genes.bed \
--output_folder \
$MARINE/tests/singlecell_tests/only_5_cells_all_cells_coverage_no_tabulation_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
16 \
--barcode_tag "CB" \
--strandedness 2 \
--contigs "9" --interval_length 32000000 --all_cells_coverage --keep_intermediate_files