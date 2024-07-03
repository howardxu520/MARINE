python \
/tscc/projects/ps-yeolab3/ekofman/sailor2/marine.py \
--bam_filepath \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/strandedness_tests/bams/pair_example_18_49488551_49590000.sorted.bam \
--annotation_bedfile_path \
/tscc/projects/ps-yeolab3/ekofman/sailor2/annotations/cellranger-GRCh38-3.0.0.annotation.genes.bed \
--output_folder \
/tscc/projects/ps-yeolab3/ekofman/sailor2/tests/strandedness_tests/pair_test \
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
--verbose \
--num_intervals_per_contig 1