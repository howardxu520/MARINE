# sailor2

Adding MD tag to bulk file:

samtools calmd -bAr /projects/ps-yeolab4/ekofman/Hugo/RBFOX2_bams/Hugo-A1Aligned.sortedByCoord.out.bam /projects/ps-yeolab3/ekofman/ReferenceData/hg38/cellranger-GRCh38-3.0.0/fasta/genome.fa > data/Hugo-A1Aligned.sortedByCoord.out.md.bam

Adding MD tag to Orel's file:

samtools calmd -bAr /home/omizrahi/scratch/orel_dms_old_data/bam/mrna_dms_treated_orel_mc2018_hg19_and_tb40.sorted.bam /projects/ps-yeolab3/ekofman/ReferenceData/hg38/cellranger-GRCh38-3.0.0/fasta/genome.fa > data/mrna_dms_treated_orel_mc2018_hg19_and_tb40.sorted.md.bam