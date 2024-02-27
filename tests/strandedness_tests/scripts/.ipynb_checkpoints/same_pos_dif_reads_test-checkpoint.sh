python \
/tscc/projects/ps-yeolab3/ekofman/sailor2/marine.py \
--bam_filepath \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/strandedness_tests/bams/same_pos_dif_reads.bam \
--annotation_bedfile_path \
/tscc/projects/ps-yeolab3/ekofman/sailor2/annotations/hg38_gencode.v35.annotation.genes.bed \
--output_folder \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/strandedness_tests/same_pos_dif_reads_test \
--min_dist_from_end \
0 \
--min_base_quality \
0 \
--cores \
1 \
--paired_end \
--contigs "chr17" \
--sailor \
--verbose 