primary_csv_path <- "data/cln/variants/final_candidates.csv"
fallback_csv_path <- "data/cln/variants/qc_filt.csv"
readme_path <- "README.md"
start_marker <- "Current final table (nonsynonymous candidates after statistical filtering):"
end_marker <- "Interpretation notes:"

csv_path <- if (file.exists(primary_csv_path)) {
  primary_csv_path
} else if (file.exists(fallback_csv_path)) {
  fallback_csv_path
} else {
  stop(
    sprintf(
      "Missing CSV file: %s (fallback also missing: %s)",
      primary_csv_path,
      fallback_csv_path
    )
  )
}

if (!file.exists(readme_path)) {
  stop(sprintf("Missing README file: %s", readme_path))
}

df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)

md_table <- c(
  paste0("| ", paste(names(df), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
)

if (nrow(df) > 0) {
  for (i in seq_len(nrow(df))) {
    row_vals <- vapply(df[i, , drop = FALSE], as.character, character(1))
    md_table <- c(md_table, paste0("| ", paste(row_vals, collapse = " | "), " |"))
  }
}

lines <- readLines(readme_path, warn = FALSE)
start_idx <- which(lines == start_marker)
end_idx <- which(lines == end_marker)

if (length(start_idx) != 1) {
  stop("README update failed: could not find unique table start marker")
}

end_idx <- end_idx[end_idx > start_idx]
if (length(end_idx) < 1) {
  stop("README update failed: could not find table end marker")
}
end_idx <- end_idx[1]

new_lines <- c(
  lines[seq_len(start_idx)],
  "",
  md_table,
  "",
  lines[end_idx:length(lines)]
)

writeLines(new_lines, readme_path)
cat(sprintf("README final table refreshed from %s\n", csv_path))
