.PHONY: build run up down clean package

build:
	docker compose build

run:
	docker compose up
	Rscript src/main.R
	Rscript src/update_readme.R

up:
	docker compose up --build

down:
	docker compose down

package:
	@echo "Writing data/DATA_README.txt..."
	@python3 -c "\
import textwrap, pathlib; \
pathlib.Path('data/DATA_README.txt').write_text(textwrap.dedent('''\
igrek: SH612x20 vs SH612x19 SNP Analysis - Data Archive\n\
=========================================================\n\
\n\
Input data (data/raw/)\n\
----------------------\n\
Raw paired-end FASTQ files for two Pseudomonas aeruginosa samples:\n\
  SH612x19  wildtype (wt)\n\
  SH612x20  mutant (mt)\n\
\n\
Reference data (data/cln/ref/)\n\
-------------------------------\n\
  pa14.fna        PA14 reference genome FASTA\n\
  pa14.fna.fai    samtools FASTA index\n\
  pa14.mmi        minimap2 index\n\
\n\
Preprocessing (data/cln/preprocessing/)\n\
-----------------------------------------\n\
fastp QC outputs per sample:\n\
  SH612x19/report.html    HTML QC report\n\
  SH612x19/report.json    JSON QC metrics\n\
\n\
Mapping (data/cln/mapping/)\n\
-----------------------------\n\
  SH612x19.bam      sorted BAM alignment\n\
  SH612x19.bam.bai  BAM index\n\
\n\
Variant calling (data/cln/variants/)\n\
--------------------------------------\n\
  joint.vcf              bcftools joint VCF (wt + mt)\n\
  counts.tsv             allele counts extracted from VCF\n\
  qc_all.csv             combined per-site table before filtering\n\
  final_candidates.csv   final filtered SNP candidates\n\
\n\
Analysis (src/main.R)\n\
----------------------\n\
  Fisher exact test on wt vs mt allele counts\n\
  Benjamini-Hochberg FDR adjustment\n\
  Coding consequence annotation (synonymous / nonsynonymous)\n\
  QC flags: depth, mapping quality, strand support\n\
'''))"
	@echo "Creating igrek_data.zip..."
	@zip -r igrek_data.zip data/ -x "*/.DS_Store" -x "data/.DS_Store"
	@echo "Done: igrek_data.zip"

clean:
	docker compose down --rmi local
	find data/cln -type f -delete
