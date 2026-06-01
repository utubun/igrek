# dependencies
library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Paeruginosa.PA14)
library(BSgenome.Paeruginosa.NCBI.PA14)


# clean for snps only, with depth >= 10
dat <- read_tsv(
    'data/cln/variants/counts.tsv',
    col_names = c('chrom', 'pos', 'ref', 'alt', 'wt_add', 'mt_add'),
    show_col_types = FALSE,
    progress = FALSE
    ) |>
    tidyr::separate(col = wt_add, into = c('wt_ref', 'wt_alt'), sep = ',', remove = TRUE, convert = TRUE) |>
    tidyr::separate(col = mt_add, into = c('mt_ref', 'mt_alt'), sep = ',', remove = TRUE, convert = TRUE) |>
    dplyr::filter(
        (wt_ref + wt_alt >= 10) & (mt_ref + mt_alt >= 10),
        !(wt_ref == 0 & mt_ref == 0),
        (nchar(ref) == 1) & (nchar(alt) == 1) 
    ) |>
    suppressMessages()

# fisher test
pval <- mapply(
    \(wr, wa, mr, ma) {
        fisher.test(
            matrix(c(wr, wa, mr, ma), nrow = 2)
        )$p.value
    },
    dat$wt_ref,
    dat$wt_alt,
    dat$mt_ref,
    dat$mt_alt
)

# add raw, and adjusted p-values to the data
dat <- dat |>
    dplyr::mutate(
        query_id        = dplyr::row_number(),
        raw_pvalue      = pval,
        adjusted_pvalue = p.adjust(pval, method = 'BH')
    )

# map snps
snp <- VRanges(
    seqnames = dat$chrom,
    ranges   = IRanges(start = dat$pos, end = dat$pos),
    ref      = dat$ref,
    alt      = dat$alt,
    sampleNames =  Rle("mt_vs_wt")
)

# add adjusted p-values to the snp map
mcols(snp)$adjusted_pvalue <- dat$adjusted_pvalue

# extract codones and AAs
coding <- predictCoding(
    query     = snp,
    subject   = TxDb.Paeruginosa.PA14,
    seqSource = BSgenome.Paeruginosa.NCBI.PA14
)

res <- coding |>
    data.frame(stringsAsFactors =  FALSE) |>
    dplyr::select(
        query_id          = QUERYID,
        seqnames,
        gene_id           = GENEID,
        position          = start,
        strand,
        base_variant      = varAllele,
        codone_reference  = REFCODON,
        codone_variant    = VARCODON,
        aa_reference      = REFAA,
        aa_variant        = VARAA,
        protein_locus     = PROTEINLOC,
        synonymous = CONSEQUENCE,
        adjusted_pvalue
    ) |>
    dplyr::mutate(
        synonymous    = as.character(synonymous) == 'synonymous',
        codone        = sprintf("%s → %s", codone_reference, codone_variant),
        aminoacid     = sprintf("%s → %s", aa_reference, aa_variant),
    ) |>
    dplyr::select(!matches("adjusted_pvalue"))

# build QC metrics for selected variants from VCF INFO fields
vcf <- VariantAnnotation::readVcf('data/cln/variants/joint.vcf')

qc <- data.frame(rowRanges(vcf)) |>
    dplyr::mutate(
        dp  = info(vcf)$DP,
        mq  = info(vcf)$MQ,
    ) |>
    bind_cols(
        as.data.frame(as.matrix(info(vcf)$DP4)) |>
          setNames(c("ref_fwd", "ref_rev", "alt_fwd", "alt_rev"))
    ) |>
    dplyr::rename_with(tolower) |>
    dplyr::mutate(
        strand_ok = alt_fwd >= 3 & alt_rev >= 3,
        mq_ok     = mq > 40
    ) |>
    dplyr::select(
        seqnames,
        pos = start,
        dp,
        mq,
        ref_fwd,
        alt_fwd,
        ref_rev,
        alt_rev,
        strand_ok
    )

# grab dat and res metrics
qc <- dat |>
  dplyr::mutate(
    wt_depth = wt_ref + wt_alt,
    mt_depth = mt_ref + mt_alt,
    wt_af     = wt_alt / wt_depth,
    mt_af     = mt_alt / mt_depth
  ) |>
  dplyr::left_join(
    distinct(
        res,
        query_id,
        seqnames,
        gene_id,
        position,
        strand,
        codone,
        aminoacid,
        protein_locus = protein_locus,
        synonymous
    ),
    by = "query_id"
  ) |>
  dplyr::left_join(qc, by = c("seqnames", "pos")) |>
  dplyr::mutate(
    depth_ok  = wt_depth >= 10 & mt_depth >= 10
  )

# write intermediate results
readr::write_csv(qc, "data/cln/variants/qc_all.csv")

# filter by adj. p-value and synonymous
qc <- dplyr::filter(qc, adjusted_pvalue <= 0.05, !synonymous)

# write final result
readr::write_csv(qc, "data/cln/variants/qc_filt.csv")

# Get transcript metadata for genes from significant non-synonymous hits
meta <- AnnotationDbi::select(
    TxDb.Paeruginosa.PA14,
    keys = unique(as.character(qc$gene_id)),
    columns = c('GENEID', 'TXNAME', 'TXCHROM', 'TXSTART', 'TXEND', 'TXSTRAND'),
    keytype = 'GENEID'
) |>
dplyr::rename_with(tolower) |>
dplyr::select(
    gene_id   = geneid,
    gene_name = txname,
    start     = txstart,
    end       = txend
)

# final table
out <- meta |>
  dplyr::left_join(qc, by = "gene_id") |>
  dplyr::mutate(
    qc       = c("fail", "pass")[1 + (mq >= 40)],
    strand   = c("fail", "pass")[1 + strand_ok],
    depth    = c("fail", "pass")[1 + depth_ok],
    locus    = purrr::map_int(protein_locus, ~ as.integer(.x)[1]),
  ) |>
  dplyr::arrange(adjusted_pvalue, wt_af, 1 - mt_af) |>
  dplyr::select(
    gene_id,
    gene_name,
    position = pos,
    codone,
    aminoacid,
    locus,
    quality  = qc,
    strand,
    depth,
    pval     = adjusted_pvalue
  )

# write data
readr::write_csv(out, "data/cln/variants/final_candidates.csv")

source("src/update_readme.R")
