# MARINE: 
Multi-core Algorithm for Rapid Identification of Nucleotide Edits
------------------------------------------------------
MARINE relies on the MD tag being present in your .bam file. To add the MD tag and then index the processed bam, use the following samtools command templates:

```
samtools calmd -bAr input.bam reference_genome.fa > input.md.bam
samtools index input.md.bam
```

Use the provided conda environment to ensure you have all required dependencies for MARINE:

```
conda env create  --file=marine_environment.yml
conda activate marine_environment
```

Notes:
* Ensure that your annotation bedfile has the same chromosome nomenclature (e.g., "9" vs "chr9") as your bam
* The more cores used, the faster MARINE will run
* The annotation file should be tab-separated and should have a non-standard bed column ordering, as follows:
```
1       29554   31109   MIR1302-2HG     +       lincRNA
1       34554   36081   FAM138A         -       lincRNA
1       65419   71585   OR4F5           +       protein_coding
1       89295   133723  AL627309.1      -       lincRNA
1       89551   91105   AL627309.3      -       lincRNA
```

### Command parameters
```
usage: marine.py [-h] [--bam_filepath BAM_FILEPATH] [--annotation_bedfile_path ANNOTATION_BEDFILE_PATH] [--output_folder OUTPUT_FOLDER]
                 [--barcode_whitelist_file BARCODE_WHITELIST_FILE] [--cores CORES] [--reverse_stranded] [--coverage] [--filtering] [--annotation]
                 [--barcode_tag BARCODE_TAG] [--min_dist_from_end MIN_DIST_FROM_END] [--min_base_quality MIN_BASE_QUALITY] [--contigs CONTIGS]
                 [--sailor] [--verbose] [--paired_end]

Run MARINE

optional arguments:
  -h, --help            show this help message and exit
  --bam_filepath BAM_FILEPATH
  --annotation_bedfile_path ANNOTATION_BEDFILE_PATH
  --output_folder OUTPUT_FOLDER
  --barcode_whitelist_file BARCODE_WHITELIST_FILE
  --cores CORES
  --reverse_stranded
  --coverage
  --filtering
  --annotation
  --barcode_tag BARCODE_TAG
  --min_dist_from_end MIN_DIST_FROM_END
  --min_base_quality MIN_BASE_QUALITY
  --contigs CONTIGS
  --sailor
  --verbose
  --paired_end
```

## Single-cell example MARINE command
```
python marine.py \
--bam_filepath /path/to/input.md.bam \
--annotation_bedfile_path /path/to/annotations/GRCh38-3.0.0.annotation.genes.bed \
--output_folder /path/to/existing/output_folder \
--barcode_tag "CB" \
--min_base_quality 15 \
--cores 32 \
--barcode_whitelist_file /path/to/filtered_feature_bc_matrix/barcodes.tsv.gz
```

## Bulk (single-end) example MARINE command
```
python marine.py \
--bam_filepath /path/to/input.md.bam \
--annotation_bedfile_path /path/to/annotations/GRCh38-3.0.0.annotation.genes.bed \
--output_folder /path/to/existing/output_folder \
--min_dist_from_end 10 \
--min_base_quality 15 \
--cores 32
```
