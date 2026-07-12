#!/usr/bin/env Rscript
# extract_single_guide_results.R
#
# Extract a compact, self-describing table from a single guide TSV:
#  - choose mode: adjusted|positional
#  - choose effect: raw|shrunk
#  - keep chosen columns with their ORIGINAL names (no renaming)
#  - optionally drop PAM-overlapping variants (ranges and/or a PAM flag column)
#  - write a separate TSV of removed PAM-overlapping variants (if any removed)
#  - exclude variants whose consequence is "Others" or "Backup_gRNA"

# Supports:
#   -f <input.tsv>
#   -o <output.tsv>
#   --mode adjusted|positional
#   --effect raw|shrunk
#   --pam_regions "start-end,start-end"
#   --pam_col column_name
#   --exclude_pam TRUE|FALSE
#   --removed_out removed.tsv

# Example:
# Rscript extract_single_guide_results.R \
#   -f all_deseq2_results_condition_Day21_vs_Plasmid.tsv \
#   -o Day21_vs_Plasmid_pos_raw_filtered.tsv \
#   --mode positional \
#   --effect raw \
#   --pam_regions "100-103,210-213" \

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Argument parsing ----
args_raw <- commandArgs(trailingOnly = TRUE)

if (length(args_raw) == 0) {
  stop("No arguments provided.")
}

args <- list()
i <- 1
while (i <= length(args_raw)) {
  key <- args_raw[i]

  # Handle short flags -f and -o
  if (key %in% c("-f", "-o")) {
    if (i == length(args_raw)) stop("Missing value for ", key)
    args[[key]] <- args_raw[i + 1]
    i <- i + 2
  } else if (grepl("^--", key)) {
    if (i == length(args_raw)) stop("Missing value for ", key)
    clean_key <- sub("^--", "", key)
    args[[clean_key]] <- args_raw[i + 1]
    i <- i + 2
  } else {
    stop("Unrecognized argument: ", key)
  }
}

# ---- Required args ----
if (is.null(args[["-f"]])) {
  stop("You must provide -f <input.tsv>")
}

file <- args[["-f"]]
out  <- if (!is.null(args[["-o"]])) args[["-o"]] else "single_output.tsv"

# ---- Defaults ----
mode  <- if (!is.null(args$mode)) args$mode else "positional"
effect <- if (!is.null(args$effect)) args$effect else "raw"
drop_prefix_fields <- if (!is.null(args$drop_prefix_fields)) as.integer(args$drop_prefix_fields) else 2
pos_col <- if (!is.null(args$pos_col)) args$pos_col else "position"
pam_regions <- if (!is.null(args$pam_regions)) args$pam_regions else NULL
pam_col <- if (!is.null(args$pam_col)) args$pam_col else NULL
exclude_pam <- if (!is.null(args$exclude_pam)) as.logical(args$exclude_pam) else TRUE
removed_out <- if (!is.null(args$removed_out)) args$removed_out else NULL

if (!(mode %in% c("adjusted", "positional"))) {
  stop("--mode must be 'adjusted' or 'positional'")
}
if (!(effect %in% c("raw", "shrunk"))) {
  stop("--effect must be 'raw' or 'shrunk'")
}

# ---- Load data ----
dt <- tryCatch(
  fread(file, sep = "\t", header = TRUE, data.table = TRUE),
  error = function(e) stop("Failed reading file: ", file, "\n", conditionMessage(e))
)
setDT(dt)

# ---- Column mapping ----
get_lfc_colname <- function(mode, effect) {
  if (mode == "adjusted") paste0("adj_log2FoldChange_", effect)
  else paste0("pos_adj_log2FoldChange_", effect)
}
get_pval_colname <- function(mode, effect) {
  if (mode == "adjusted") paste0("adj_pval_", effect)
  else paste0("pos_adj_pval_", effect)
}
get_fdr_colname <- function(mode, effect) {
  if (mode == "adjusted") paste0("adj_fdr_", effect)
  else paste0("pos_adj_fdr_", effect)
}
get_stat_colname <- function(mode, effect) {
  if (mode == "adjusted") paste0("stat_adj_", effect)
  else paste0("stat_pos_", effect)
}
get_se_colname <- function(mode, effect) {
  if (mode == "positional") {
    if (effect == "raw") return("pos_total_se_raw")
    return("pos_total_se_shrunk")
  } else {
    if (effect == "raw") return("lfcSE_raw")
    return("lfcSE_shrunk")
  }
}

# If the new pos_total_se_* cols aren't present, fall back to lfcSE_* with a
# warning (mirrors merge_and_combine.R's resolve_se_col behavior, so the SE
# convention is consistent whether a variant came through the merge script or
# this single-guide script).
resolve_se_col <- function(dt, mode, effect) {
  desired <- get_se_colname(mode, effect)
  if (desired %in% names(dt)) return(desired)

  if (mode == "positional") {
    fallback <- if (effect == "raw") "lfcSE_raw" else "lfcSE_shrunk"
    if (fallback %in% names(dt)) {
      warning("Missing '", desired, "'; falling back to '", fallback,
              "'. (You may be using outputs generated before pos_total_se_* was added.)")
      return(fallback)
    }
  }
  stop("Missing SE column: ", desired)
}

lfc_col  <- get_lfc_colname(mode, effect)
se_col   <- resolve_se_col(dt, mode, effect)
pval_col <- get_pval_colname(mode, effect)
fdr_col  <- get_fdr_colname(mode, effect)
stat_col <- get_stat_colname(mode, effect)

required_cols <- c(lfc_col, pval_col, fdr_col, stat_col)

if (!("oligo_name" %in% names(dt))) stop("Input must contain 'oligo_name'.")
if (!("sequence" %in% names(dt))) dt[, sequence := NA_character_]
if (!("consequence" %in% names(dt))) dt[, consequence := NA_character_]

missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Filter unwanted consequence types ----
bad_cons <- c("Others", "Backup_gRNA")
dt <- dt[!(consequence %in% bad_cons)]

# ---- Position ----
if (!(pos_col %in% names(dt))) {
  warning("Position column not found: ", pos_col, " (creating NA column)")
  dt[, (pos_col) := NA_integer_]
}
dt[, (pos_col) := as.integer(get(pos_col))]

# ---- Variant key ----
make_variant_key_drop_prefix <- function(oligo_name, drop_n = 2L) {
  vapply(oligo_name, function(x) {
    parts <- strsplit(as.character(x), "_", fixed = TRUE)[[1]]
    if (length(parts) <= drop_n) as.character(x)
    else paste(parts[(drop_n + 1):length(parts)], collapse = "_")
  }, FUN.VALUE = character(1))
}
dt[, variant_key := make_variant_key_drop_prefix(oligo_name, drop_prefix_fields)]

# ---- PAM handling ----
parse_ranges <- function(ranges_str) {
  if (is.null(ranges_str) || ranges_str == "") return(NULL)
  pieces <- strsplit(ranges_str, ",", fixed = TRUE)[[1]]
  out <- lapply(pieces, function(p) {
    r <- strsplit(trimws(p), "-", fixed = TRUE)[[1]]
    if (length(r) != 2) stop("Invalid range: ", p)
    as.integer(r)
  })
  do.call(rbind, out)
}

mark_pam_by_ranges <- function(positions, ranges) {
  if (is.null(ranges)) return(rep(FALSE, length(positions)))
  res <- rep(FALSE, length(positions))
  for (i in seq_len(nrow(ranges))) {
    r0 <- ranges[i, ]
    res <- res | (positions >= r0[1] & positions <= r0[2])
  }
  res
}

pam_ranges_obj <- parse_ranges(pam_regions)
dt[, pam_range_flag := mark_pam_by_ranges(get(pos_col), pam_ranges_obj)]

if (!is.null(pam_col) && pam_col %in% names(dt)) {
  dt[, pam_col_flag := !is.na(get(pam_col)) &
                    get(pam_col) != "" &
                    get(pam_col) != "FALSE"]
} else {
  dt[, pam_col_flag := FALSE]
}

dt[, pam_flag := pam_range_flag | pam_col_flag]

removed_dt <- NULL
if (isTRUE(exclude_pam)) {
  # Previously this dropped PAM-overlapping variants from dt entirely.
  # Now they are retained in the main output (flagged via pam_flag) and
  # additionally written to removed_out as an audit list, but no longer removed.
  removed_dt <- dt[pam_flag %in% TRUE]
}

# ---- Output ----
out_cols <- c(
  "variant_key",
  "oligo_name",
  "consequence",
  "sequence",
  pos_col,
  "pam_flag",
  lfc_col,
  se_col,
  pval_col,
  fdr_col,
  stat_col
)

fwrite(dt[, ..out_cols], file = out, sep = "\t", quote = FALSE)
cat("Written filtered output to ", out, "\n", sep = "")

if (isTRUE(exclude_pam)) {
  if (!is.null(removed_dt) && nrow(removed_dt) > 0) {
    if (is.null(removed_out)) {
      removed_out <- sub("\\.tsv$", "_pam_flagged.tsv", out)
      if (identical(removed_out, out)) removed_out <- paste0(out, "_pam_flagged.tsv")
    }
    fwrite(removed_dt[, ..out_cols], file = removed_out, sep = "\t", quote = FALSE)
    cat("Written PAM-overlapping variants (audit copy; still retained in main output) to ", removed_out, "\n", sep = "")
    cat("Flagged ", nrow(removed_dt), " variants as PAM-overlapping (pam_flag = TRUE); they remain in the main output.\n", sep = "")
  } else {
    cat("No PAM-overlapping variants were found.\n")
  }
}
