#! /bin/bash

marine.py --bam_filepath $MARINE/examples/data/bulk_AI.md.subset.bam --output_folder $MARINE/examples/bulk_subset_AI --strandedness 2 --cores 16 --annotation_bedfile_path $MARINE/annotations/hg38_gencode.v35.annotation.genes.bed --contigs "1" --bedgraphs "AI" --min_base_quality 15 --min_dist_from_end 5