mypython=$1

$mypython $MARINE/marine.py \
--bam_filepath \
$MARINE/tests/strandedness_tests/bams/pair_example_18_49488551_49590000.sorted.bam \
--annotation_bedfile_path \
$MARINE/annotations/cellranger-GRCh38-3.0.0.annotation.genes.bed \
--output_folder \
$MARINE/tests/strandedness_tests/pair_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
1 \
--paired_end \
--strandedness 2 \
--contigs 18 \
--sailor \
--keep_intermediate_files \
--num_intervals_per_contig 1