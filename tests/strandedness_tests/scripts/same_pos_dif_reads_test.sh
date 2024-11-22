mypython=$1

$mypython $MARINE/marine.py \
--bam_filepath \
$MARINE/tests/strandedness_tests/bams/same_pos_dif_reads.bam \
--annotation_bedfile_path \
$MARINE/annotations/hg38_gencode.v35.annotation.genes.bed \
--output_folder \
$MARINE/tests/strandedness_tests/same_pos_dif_reads_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
16 \
--paired_end \
--strandedness 2 \
--contigs "chr17" \
--num_intervals_per_contig 16 \
--keep_intermediate_files