mypython=$1

$mypython $MARINE/marine.py \
--bam_filepath \
$MARINE/tests/strandedness_tests/bams/citrine435.bam \
--annotation_bedfile_path \
$MARINE/annotations/hg38_gencode.v35.annotation.genes.bed \
--output_folder \
$MARINE/tests/strandedness_tests/unstranded_pair_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
1 \
--paired_end \
--strandedness 0 \
--sailor \
--num_intervals_per_contig 1