# igrek: SH612x20 vs SH612x19 SNP Analysis

## Summary

This project compares two Pseudomonas aeruginosa samples against the PA14 reference genome:

- SH612x19: wildtype (wt)
- SH612x20: mutant (mt)

Pipeline steps:

1. fastp read preprocessing
2. minimap2 alignment + samtools BAM processing
3. bcftools joint variant calling
4. R-based postprocessing:
   - per-site Fisher exact test on wt vs mt allele counts
   - Benjamini-Hochberg p-value adjustment
   - coding consequence annotation (synonymous vs nonsynonymous)
   - QC flags from VCF fields (depth, mapping quality, strand support)

---

## Where the data is

Input data:

- `data/raw/`: raw FASTQ files (not tracked in Git)

Runtime/generated data:

- `data/cln/ref/`: reference FASTA and minimap2 index
- `data/cln/preprocessing/`: fastp outputs
- `data/cln/mapping/`: BAM and BAI files
- `data/cln/variants/joint.vcf`: joint VCF for wt and mt
- `data/cln/variants/counts.tsv`: extracted allele counts used by R

Analysis outputs:

- `data/cln/variants/qc_all.csv`: combined table before final filtering
- `data/cln/variants/final_candidates.csv`: final filtered candidate table

---

## Main scripts and config

- `entrypoint.sh`: end-to-end container pipeline
- `src/main.R`: statistical analysis + annotation + QC integration
- `Dockerfile`: analysis environment
- `docker-compose.yml`: container run config
- `Makefile`: convenience commands

---

## How to run

```bash
make build
make run
```

This runs the containerized workflow and writes outputs under `data/cln/`.

---

## Reproducibility Checklist

Use this checklist before sharing results:

1. `make build` completed successfully (Docker image up to date).
2. `make run` completed without errors.
3. These files exist and are fresh:
   - `data/cln/variants/joint.vcf`
   - `data/cln/variants/counts.tsv`
   - `data/cln/variants/qc_all.csv`
   - `data/cln/variants/final_candidates.csv`
4. `final_candidates.csv` matches the filtering policy currently documented in `src/main.R`.
5. README final table matches `final_candidates.csv`.
6. Raw inputs (`data/raw/`) are not staged in Git.

---

## Final candidate table

Current final table (nonsynonymous candidates after statistical filtering):

| gene_id | gene_name | position | codone | aminoacid | locus | quality | strand | depth | pval |
|---|---|---|---|---|---|---|---|---|---|
| PA14_RS18220 | puuE | 3995746 | CTG → CCG | L → P | 207 | pass | pass | pass | 2.2173514976159e-48 |
| PA14_RS19850 | PA14_RS19850 | 4345430 | ATG → ATA | M → I | 258 | pass | pass | pass | 2.1393148982346e-32 |
| PA14_RS14050 | PA14_RS14050 | 3065715 | GGC → GCC | G → A | 106 | pass | fail | pass | 0.00558944648590917 |
| PA14_RS16675 | PA14_RS16675 | 3670539 | CCT → CTT | P → L | 180 | pass | pass | pass | 0.0389268040518032 |
| PA14_RS04145 | PA14_RS04145 | 883169 | GGC → GCC | G → A | 37 | pass | pass | pass | 0.0409566020740156 |

Interpretation notes:

- `quality`: mapping-quality-based pass flag
- `strand`: strand-support pass flag (ALT supported on both strands)
- `depth`: coverage pass flag for both samples

---

## Statistical and biological caveat

This is a two-sample comparison (one wt and one mt). Results are read-supported differences between these sequenced samples, not population-level estimates.
