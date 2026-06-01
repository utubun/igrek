#!/bin/bash
set -e


echo "Creating output directories..."
mkdir -p cln/ref
mkdir -p cln/preprocessing/SH612x19
mkdir -p cln/preprocessing/SH612x20
mkdir -p cln/mapping
mkdir -p cln/variants

echo "Cleaning output subdirectories before analysis..."
find cln/preprocessing -type f -delete
find cln/mapping -type f -delete
find cln/variants -type f -delete
find cln/ref -type f -delete


echo "Running fastp for SH612x19..."
fastp \
  --in1 raw/SH612x19_251023_LH00255_A235FKJLT3_R1.fastq.gz \
  --in2 raw/SH612x19_251023_LH00255_A235FKJLT3_R2.fastq.gz \
  --unpaired1 cln/preprocessing/SH612x19/unpaired.fastq.gz \
  --unpaired2 cln/preprocessing/SH612x19/unpaired.fastq.gz \
  --failed_out cln/preprocessing/SH612x19/failed.fastq.gz \
  --dont_overwrite \
  --overrepresentation_analysis \
  --dedup \
  --trim_poly_g \
  --trim_poly_x \
  --low_complexity_filter \
  --correction \
  --merge \
  --merged_out cln/preprocessing/SH612x19/merged.fastq.gz \
  --json cln/preprocessing/SH612x19/report.json \
  --html cln/preprocessing/SH612x19/report.html


echo "Running fastp for SH612x20..."
fastp \
  --in1 raw/SH612x20_251023_LH00255_A235FKJLT3_R1.fastq.gz \
  --in2 raw/SH612x20_251023_LH00255_A235FKJLT3_R2.fastq.gz \
  --unpaired1 cln/preprocessing/SH612x20/unpaired.fastq.gz \
  --unpaired2 cln/preprocessing/SH612x20/unpaired.fastq.gz \
  --failed_out cln/preprocessing/SH612x20/failed.fastq.gz \
  --dont_overwrite \
  --overrepresentation_analysis \
  --dedup \
  --trim_poly_g \
  --trim_poly_x \
  --low_complexity_filter \
  --correction \
  --merge \
  --merged_out cln/preprocessing/SH612x20/merged.fastq.gz \
  --include_unmerged \
  --json cln/preprocessing/SH612x20/report.json \
  --html cln/preprocessing/SH612x20/report.html


echo "Downloading reference genome (PA14, GCF_000014625.1)..."
wget -q -O cln/ref/pa14.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/014/625/GCF_000014625.1_ASM1462v1/GCF_000014625.1_ASM1462v1_genomic.fna.gz" \
  && gunzip -f cln/ref/pa14.fna.gz

echo "Indexing reference genome for read alignment..."
minimap2 -d cln/ref/pa14.mmi cln/ref/pa14.fna

echo "Aligning SH612x19 reads to reference..."
minimap2 -ax sr cln/ref/pa14.mmi cln/preprocessing/SH612x19/merged.fastq.gz \
  | samtools sort -o cln/mapping/SH612x19.bam
samtools index cln/mapping/SH612x19.bam

echo "Aligning SH612x20 reads to reference..."
minimap2 -ax sr cln/ref/pa14.mmi cln/preprocessing/SH612x20/merged.fastq.gz \
  | samtools sort -o cln/mapping/SH612x20.bam
samtools index cln/mapping/SH612x20.bam

echo "Calling variants jointly for both samples..."
bcftools mpileup -a FORMAT/AD -f cln/ref/pa14.fna \
  cln/mapping/SH612x19.bam cln/mapping/SH612x20.bam \
  | bcftools call -mv -o cln/variants/joint.vcf

echo "Extracting per-sample allele counts..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
  cln/variants/joint.vcf > cln/variants/counts.tsv