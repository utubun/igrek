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
        queryid         = dplyr::row_number(),
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
    dplyr::mutate(synonymous = as.character(synonymous) == 'synonymous') |>
    dplyr::filter(adjusted_pvalue <= 0.05, !synonymous)

# build QC metrics for selected variants from VCF INFO fields
vcf <- VariantAnnotation::readVcf('data/cln/variants/joint.vcf')
dp4 <- info(vcf)$DP4

pluck_int1 <- function(x) {
    vapply(
        as.list(x),
        function(el) {
            if (length(el) == 0) {
                return(NA_integer_)
            }
            as.integer(el[[1]])
        },
        integer(1)
    )
}

pluck_num1 <- function(x) {
    vapply(
        as.list(x),
        function(el) {
            if (length(el) == 0) {
                return(NA_real_)
            }
            as.numeric(el[[1]])
        },
        numeric(1)
    )
}

dp <- pluck_int1(info(vcf)$DP)
mq <- pluck_num1(info(vcf)$MQ)

dp4_mat <- t(vapply(
    as.list(dp4),
    function(el) {
        z <- as.integer(el)
        out <- rep(NA_integer_, 4)
        out[seq_len(min(length(z), 4))] <- z[seq_len(min(length(z), 4))]
        out
    },
    integer(4)
))

vcf_qc <- tibble::tibble(
    seqnames = as.character(GenomicRanges::seqnames(rowRanges(vcf))),
    pos = GenomicRanges::start(rowRanges(vcf)),
    dp = dp,
    mq = mq,
    ref_fwd = dp4_mat[, 1],
    ref_rev = dp4_mat[, 2],
    alt_fwd = dp4_mat[, 3],
    alt_rev = dp4_mat[, 4]
)

# Backward-compatible guards for interactive sessions with stale objects
if (!'queryid' %in% names(dat)) {
    dat <- dat |>
        dplyr::mutate(queryid = dplyr::row_number())
}

if (!'queryid' %in% names(res)) {
    if ('QUERYID' %in% names(res)) {
        res <- res |>
            dplyr::mutate(queryid = QUERYID)
    } else {
        warning('res has no queryid/QUERYID; using fallback join by seqnames+pos+alt. Rerun full script for deterministic QUERYID mapping.')
    }
}

if ('queryid' %in% names(res)) {
    qc <- res |>
        dplyr::distinct(queryid, seqnames, pos, gene, alt, pvalue) |>
        dplyr::left_join(
            dat |>
                dplyr::transmute(
                    queryid,
                    wt_ref,
                    wt_alt,
                    mt_ref,
                    mt_alt,
                    wt_depth = wt_ref + wt_alt,
                    mt_depth = mt_ref + mt_alt,
                    wt_af = wt_alt / (wt_ref + wt_alt),
                    mt_af = mt_alt / (mt_ref + mt_alt)
                ),
            by = 'queryid'
        )
} else {
    qc <- res |>
        dplyr::distinct(seqnames, pos, gene, alt, pvalue) |>
        dplyr::left_join(
            dat |>
                dplyr::transmute(
                    seqnames = chrom,
                    pos,
                    alt,
                    wt_ref,
                    wt_alt,
                    mt_ref,
                    mt_alt,
                    wt_depth = wt_ref + wt_alt,
                    mt_depth = mt_ref + mt_alt,
                    wt_af = wt_alt / (wt_ref + wt_alt),
                    mt_af = mt_alt / (mt_ref + mt_alt)
                ),
            by = c('seqnames', 'pos', 'alt')
        )
}

qc <- qc |>
    dplyr::left_join(vcf_qc, by = c('seqnames', 'pos')) |>
    dplyr::mutate(
        strand_ok = alt_fwd >= 3 & alt_rev >= 3,
        mq_ok = mq >= 40,
        depth_ok = wt_depth >= 10 & mt_depth >= 10,
        neighborhood_density = purrr::map2_int(
            seqnames,
            pos,
            \(chr, p) {
                rr <- GenomicRanges::GRanges(
                    seqnames = chr,
                    ranges = IRanges(start = max(1L, p - 20L), end = p + 20L)
                )
                sum(IRanges::overlapsAny(rowRanges(vcf), rr))
            }
        )
    ) |>
    dplyr::arrange(pvalue)

# Get transcript metadata for genes from significant non-synonymous hits
meta <- AnnotationDbi::select(
    TxDb.Paeruginosa.PA14,
    keys = unique(as.character(qc$gene)),
    columns = c('GENEID', 'TXNAME', 'TXCHROM', 'TXSTART', 'TXEND', 'TXSTRAND'),
    keytype = 'GENEID'
)
