#! /bin/bash

python marine.py --bam_filepath examples/data/bulk_AI.md.subset.bam --output_folder examples/bulk_subset_AI --strandedness 2 --cores 16 --annotation_bedfile_path \
/tscc/projects/ps-yeolab3/ekofman/sailor2/annotations/hg38_gencode.v35.annotation.genes.bed --contigs "1"