#!/usr/bin/env Rscript
# merge_and_combine.R  (dependency-free argument parsing)
#
# Merge two guide/library result TSVs and compute combined LFC via inverse-variance weighting.
# Key features:
#  - merge key created by dropping first N underscore-separated prefixes from oligo_name (default N=2)
#  - choose mode: adjusted | positional
#  - choose effect: raw | shrunk
#  - specify PAM windows as numeric ranges per library (--pam_regions1, --pam_regions2)
#    OR supply a PAM flag column name in each file (--pam_col1/2)
#  - optionally zero weights for PAM-overlapping variants (default TRUE)
#  - exclude variants whose consequence is "Others" or "Backup_gRNA"
#  - keeps variants present in only one file (unless excluded by PAM/flags)
# Output columns include: sequences from both guides, reconciled position column named 'position',
# plus status1/status2 (selected based on --mode/--effect).
#
# Dependencies: data.table (required). ggplot2 optional (only if --plot_prefix set).
#
# Help:
# Rscript merge_and_combine.R --help
#
# Example:
# Rscript merged/merge_and_combine.R /
# -f1 sg7/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Day4.tsv /
# -f2 sg8/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Day4.tsv /
# -o merged/attempt3/DDX3X_merged_D16vD4.tsv --mode adjusted --effect raw /
# --pam_regions1 "41334188-41334210" --pam_regions2 "41334263-41334285" /
# --plot_prefix sg7_sg8_D16_adjusted_diagnostics
#
# çopy:
# Rscript merged/merge_and_combine.R -f1 sg7/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Day4.tsv -f2 sg8/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Day4.tsv -o merged/attempt3/DDX3X_merged_D16vD4.tsv --mode positional --effect raw --pam_regions1 "41334188-41334210" --pam_regions2 "41334263-41334285" --plot_prefix sg7_sg8_D16_positional_diagnostics


suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("This script requires the 'data.table' package. Install with install.packages('data.table').")
  }
})
library(data.table)

# ----------------------------
# CLI arg parsing (base R)
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat(
"merge_and_combine_guides.R

Required:
  -f1 <file1.tsv>      Guide/library 1 results TSV
  -f2 <file2.tsv>      Guide/library 2 results TSV

Optional:
  -o  <out.tsv>        Output TSV (default: combined_output.tsv)

  --mode <adjusted|positional>   (default: adjusted)
  --effect <raw|shrunk>          (default: raw)

  --drop_prefix_fields <int>     Number of underscore-delimited fields to drop from oligo_name (default: 2)
  --pos_col1 <name>              Position column name in file1 (default: position)
  --pos_col2 <name>              Position column name in file2 (default: position)

  PAM flagging (choose either ranges or columns; both allowed):
  --pam_regions1 \"start-end,start-end\"  PAM ranges for file1 (e.g. \"100-103,210-213\")
  --pam_regions2 \"start-end,start-end\"  PAM ranges for file2
  --pam_col1 <name>              Column in file1 indicating PAM overlap (truthy)
  --pam_col2 <name>              Column in file2 indicating PAM overlap (truthy)

  --exclude_pam <TRUE|FALSE>     If TRUE, weights for PAM-overlapping variants are set to 0 (default: TRUE)

  Threshold for combined_status:
  --pcut <double>     FDR cutoff (default 0.05)

  Optional diagnostics (requires ggplot2):
  --plot_prefix <prefix>   Writes diagnostics next to output by default:
                            <out_dir>/<prefix>_lfc_scatter.png
                            <out_dir>/<prefix>_combined_SE_hist.png
                            <out_dir>/<prefix>_pos_discordant.tsv
                           If prefix includes a path separator, that path is respected.

Examples:
  Rscript merge_and_combine_guides.R -f1 sg7.tsv -f2 sg8.tsv -o results/combined.tsv --mode adjusted --effect raw --pam_regions1 \"100-102\" --pam_regions2 \"145-147\" --plot_prefix diag

",
  sep = ""
  )
}

get_flag <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[1] == length(args)) stop("Flag provided without value: ", flag)
  args[idx[1] + 1]
}
has_flag <- function(flag) flag %in% args

as_logical_flag <- function(x, default = TRUE) {
  if (is.null(x)) return(default)
  x2 <- toupper(trimws(as.character(x)))
  if (x2 %in% c("TRUE","T","1","YES","Y")) return(TRUE)
  if (x2 %in% c("FALSE","F","0","NO","N")) return(FALSE)
  stop("Invalid logical value: ", x, " (expected TRUE/FALSE)")
}
as_int_flag <- function(x, default) {
  if (is.null(x)) return(default)
  suppressWarnings(v <- as.integer(x))
  if (is.na(v)) stop("Invalid integer value: ", x)
  v
}
as_num_flag <- function(x, default) {
  if (is.null(x)) return(default)
  suppressWarnings(v <- as.numeric(x))
  if (is.na(v)) stop("Invalid numeric value: ", x)
  v
}

if (has_flag("--help") || has_flag("-h")) {
  print_help()
  quit(save = "no", status = 0)
}

# Required
f1 <- get_flag("-f1")
f2 <- get_flag("-f2")
if (is.null(f1) || is.null(f2)) {
  print_help()
  stop("Missing required args: -f1 and/or -f2")
}

# Optional
out <- get_flag("-o", "combined_output.tsv")

mode <- get_flag("--mode", "adjusted")
effect <- get_flag("--effect", "raw")
drop_prefix_fields <- as_int_flag(get_flag("--drop_prefix_fields"), 2L)

pos_col1 <- get_flag("--pos_col1", "position")
pos_col2 <- get_flag("--pos_col2", "position")

pam_regions1 <- get_flag("--pam_regions1", NULL)
pam_regions2 <- get_flag("--pam_regions2", NULL)
pam_col1 <- get_flag("--pam_col1", NULL)
pam_col2 <- get_flag("--pam_col2", NULL)

exclude_pam <- as_logical_flag(get_flag("--exclude_pam", NULL), default = TRUE)

pcut <- as_num_flag(get_flag("--pcut", NULL), 0.05)

plot_prefix <- get_flag("--plot_prefix", NULL)

# validate mode/effect
mode <- match.arg(tolower(mode), c("adjusted", "positional"))
effect <- match.arg(tolower(effect), c("raw", "shrunk"))

# ----------------------------
# Helpers
# ----------------------------
read_tsv <- function(path) {
  tryCatch(
    fread(path, sep = "\t", header = TRUE, data.table = TRUE),
    error = function(e) stop("Failed reading TSV: ", path, "\n", conditionMessage(e))
  )
}
write_tsv <- function(dt, path) {
  fwrite(dt, file = path, sep = "\t", quote = FALSE)
}

make_variant_key_drop_prefix <- function(oligo_name, drop_n = 2L) {
  sapply(oligo_name, function(x) {
    parts <- unlist(strsplit(as.character(x), "_", fixed = TRUE))
    if (length(parts) <= drop_n) as.character(x) else paste(parts[(drop_n + 1):length(parts)], collapse = "_")
  }, USE.NAMES = FALSE)
}

parse_ranges <- function(ranges_str) {
  if (is.null(ranges_str) || ranges_str == "") return(NULL)
  pieces <- unlist(strsplit(ranges_str, ",", fixed = TRUE))
  out <- lapply(pieces, function(p) {
    r <- unlist(strsplit(trimws(p), "-", fixed = TRUE))
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

`%notin%` <- function(a, b) !(a %in% b)

truthy_col <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  x2 <- as.character(x)
  !is.na(x2) & x2 != "" & toupper(x2) %notin% c("FALSE", "F", "0", "NO", "N")
}

get_lfc_colname <- function(mode, effect) {
  if (mode == "adjusted") {
    paste0("adj_log2FoldChange_", ifelse(effect == "raw", "raw", "shrunk"))
  } else {
    paste0("pos_adj_log2FoldChange_", ifelse(effect == "raw", "raw", "shrunk"))
  }
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

# If the new pos_total_se_* cols aren't present, fall back to lfcSE_* with a warning
resolve_se_col <- function(dt, mode, effect, which_file = "file") {
  desired <- get_se_colname(mode, effect)
  if (desired %in% names(dt)) return(desired)

  if (mode == "positional") {
    fallback <- if (effect == "raw") "lfcSE_raw" else "lfcSE_shrunk"
    if (fallback %in% names(dt)) {
      warning(which_file, " missing '", desired, "'; falling back to '", fallback,
              "'. (You may be using outputs generated before pos_total_se_* was added.)")
      return(fallback)
    }
  }
  stop("Missing SE column in ", which_file, ": ", desired)
}
get_status_colname <- function(mode, effect) {
  # expected names from your pipeline:
  # adjusted -> stat_adj_raw / stat_adj_shrunk
  # positional -> stat_pos_raw / stat_pos_shrunk
  if (mode == "adjusted") {
    paste0("stat_adj_", ifelse(effect == "raw", "raw", "shrunk"))
  } else {
    paste0("stat_pos_", ifelse(effect == "raw", "raw", "shrunk"))
  }
}

# Diagnostics output location helper:
# If plot_prefix has path separators, respect it; otherwise put next to output file.
make_diag_prefix_path <- function(out_path, plot_prefix) {
  if (is.null(plot_prefix)) return(NULL)
  has_sep <- grepl("[/\\\\]", plot_prefix)
  if (has_sep) return(plot_prefix)
  out_dir <- dirname(out_path)
  file.path(out_dir, plot_prefix)
}

# ----------------------------
# Load & filter
# ----------------------------
dt1 <- read_tsv(f1)
dt2 <- read_tsv(f2)

bad_cons <- c("Others", "Backup_gRNA")
if ("consequence" %in% names(dt1)) dt1 <- dt1[!consequence %in% bad_cons]
if ("consequence" %in% names(dt2)) dt2 <- dt2[!consequence %in% bad_cons]

if (!("oligo_name" %in% names(dt1)) || !("oligo_name" %in% names(dt2))) {
  stop("Both inputs must contain 'oligo_name' column.")
}

# Ensure sequence exists (requested output)
if (!("sequence" %in% names(dt1))) {
  warning("File1 missing 'sequence' column; sequence1 will be NA.")
  dt1[, sequence := NA_character_]
}
if (!("sequence" %in% names(dt2))) {
  warning("File2 missing 'sequence' column; sequence2 will be NA.")
  dt2[, sequence := NA_character_]
}

# Variant keys
dt1[, variant_key := make_variant_key_drop_prefix(oligo_name, drop_prefix_fields)]
dt2[, variant_key := make_variant_key_drop_prefix(oligo_name, drop_prefix_fields)]

# Positions
if (!(pos_col1 %in% names(dt1))) {
  warning("pos_col1 not found in file1; adding NA column: ", pos_col1)
  dt1[, (pos_col1) := NA_integer_]
}
if (!(pos_col2 %in% names(dt2))) {
  warning("pos_col2 not found in file2; adding NA column: ", pos_col2)
  dt2[, (pos_col2) := NA_integer_]
}
dt1[, (pos_col1) := as.integer(get(pos_col1))]
dt2[, (pos_col2) := as.integer(get(pos_col2))]

# PAM flags
ranges1 <- parse_ranges(pam_regions1)
ranges2 <- parse_ranges(pam_regions2)

if (is.null(pam_regions1) && is.null(pam_col1)) warning("No PAM specification provided for file1 (--pam_regions1 or --pam_col1).")
if (is.null(pam_regions2) && is.null(pam_col2)) warning("No PAM specification provided for file2 (--pam_regions2 or --pam_col2).")

dt1[, pam_flag := mark_pam_by_ranges(get(pos_col1), ranges1)]
dt2[, pam_flag := mark_pam_by_ranges(get(pos_col2), ranges2)]

if (!is.null(pam_col1) && pam_col1 %in% names(dt1)) {
  dt1[, pam_flag := pam_flag | truthy_col(get(pam_col1))]
}
if (!is.null(pam_col2) && pam_col2 %in% names(dt2)) {
  dt2[, pam_flag := pam_flag | truthy_col(get(pam_col2))]
}

# ----------------------------
# Choose columns
# ----------------------------
lfc_col <- get_lfc_colname(mode, effect)

# Resolve SE column per file (allows fallback to differ per input)
se_col1 <- resolve_se_col(dt1, mode, effect, which_file = "file1")
se_col2 <- resolve_se_col(dt2, mode, effect, which_file = "file2")

status_col <- get_status_colname(mode, effect)

# Output-facing names (reflect chosen input column names)
lfc1_out <- paste0(lfc_col, "_1")
lfc2_out <- paste0(lfc_col, "_2")
se1_out  <- paste0(se_col1, "_1")
se2_out  <- paste0(se_col2, "_2")
status1_out <- paste0(status_col, "_1")
status2_out <- paste0(status_col, "_2")

# Required columns check
if (!(lfc_col %in% names(dt1))) stop("Missing column in file1: ", lfc_col)
if (!(lfc_col %in% names(dt2))) stop("Missing column in file2: ", lfc_col)
# se_col1/se_col2 existence already guaranteed by resolve_se_col()

# status columns: optional; if missing, fill with NA and warn
if (!(status_col %in% names(dt1))) {
  warning("File1 missing status column '", status_col, "'. status1 will be NA.")
  dt1[, (status_col) := NA_character_]
}
if (!(status_col %in% names(dt2))) {
  warning("File2 missing status column '", status_col, "'. status2 will be NA.")
  dt2[, (status_col) := NA_character_]
}

# Rename position columns to standard names for merge
setnames(dt1, pos_col1, "position1")
setnames(dt2, pos_col2, "position2")

# Subsets for merge (keep generic internal names: lfc1/lfc2 etc.)
df1 <- dt1[, .(variant_key,
               oligo_name1 = oligo_name,
               sequence1 = sequence,
               consequence1 = if ("consequence" %in% names(dt1)) as.character(consequence) else NA_character_,
               position1,
               pam1 = pam_flag,
               lfc1 = get(lfc_col),
               se1 = get(se_col1),
               status1 = as.character(get(status_col)))]

df2 <- dt2[, .(variant_key,
               oligo_name2 = oligo_name,
               sequence2 = sequence,
               consequence2 = if ("consequence" %in% names(dt2)) as.character(consequence) else NA_character_,
               position2,
               pam2 = pam_flag,
               lfc2 = get(lfc_col),
               se2 = get(se_col2),
               status2 = as.character(get(status_col)))]

merged <- merge(df1, df2, by = "variant_key", all = TRUE)
merged[, consequence := fifelse(!is.na(consequence1) & consequence1 != "", consequence1, consequence2)]

# ----------------------------
# Reconcile position (+ record source)
# ----------------------------
merged[, c("position", "position_source") := {
  p1 <- position1; p2 <- position2
  pos <- NA_integer_
  src <- NA_character_

  if (is.na(p1) && is.na(p2)) { pos <- NA_integer_; src <- "none" }
  else if (!is.na(p1) && is.na(p2)) { pos <- p1; src <- "only_file1" }
  else if (is.na(p1) && !is.na(p2)) { pos <- p2; src <- "only_file2" }
  else if (p1 == p2) { pos <- p1; src <- "both_equal" }
  else {
    if (!is.na(pam1) && pam1 && (is.na(pam2) || !pam2)) { pos <- p2; src <- "file2_nonPAM" }
    else if (!is.na(pam2) && pam2 && (is.na(pam1) || !pam1)) { pos <- p1; src <- "file1_nonPAM" }
    else { pos <- p2; src <- "tie_default_file2" }
  }
  list(pos, src)
}, by = seq_len(nrow(merged))]

# ----------------------------
# Weights & combine
# ----------------------------
merged[, w1 := fifelse(!is.na(lfc1) & !is.na(se1) & is.finite(se1) & se1 > 0, 1 / (se1^2), 0)]
merged[, w2 := fifelse(!is.na(lfc2) & !is.na(se2) & is.finite(se2) & se2 > 0, 1 / (se2^2), 0)]

if (exclude_pam) {
  merged[!is.na(pam1) & pam1 == TRUE, w1 := 0]
  merged[!is.na(pam2) & pam2 == TRUE, w2 := 0]
}

merged[, sum_w := w1 + w2]

merged[, combined_LFC := NA_real_]
merged[, combined_SE := NA_real_]

## NOTE: lfc1/lfc2 are coalesced to 0 here purely to avoid R's NA * 0 = NA
## behaviour. This is safe because whenever lfc1 (or lfc2) is NA, its
## corresponding weight w1 (or w2) is already guaranteed to be exactly 0
## (see the fifelse() definitions of w1/w2 above), so the coalesced
## placeholder value is always multiplied by zero and never contributes
## to the weighted average. Without this, variants present in only one
## guide's library (the other guide having no data at all, not merely a
## PAM-zeroed weight) were incorrectly assigned combined_LFC = NA instead
## of correctly inheriting the single available guide's LFC.
merged[sum_w > 0, combined_LFC := ((fcoalesce(lfc1, 0) * w1) + (fcoalesce(lfc2, 0) * w2)) / sum_w]
merged[sum_w > 0, combined_SE := 1 / sqrt(sum_w)]

merged[, combined_Z := fifelse(!is.na(combined_LFC) & !is.na(combined_SE) & combined_SE > 0,
                               combined_LFC / combined_SE, NA_real_)]
merged[, combined_p := fifelse(!is.na(combined_Z), 2 * pnorm(abs(combined_Z), lower.tail = FALSE), NA_real_)]
merged[, combined_FDR := p.adjust(combined_p, method = "BH")]

merged[, combined_status := "no impact"]
merged[!is.na(combined_FDR) & combined_FDR < pcut & !is.na(combined_LFC) & combined_LFC > 0, combined_status := "enriched"]
merged[!is.na(combined_FDR) & combined_FDR < pcut & !is.na(combined_LFC) & combined_LFC < 0, combined_status := "depleted"]
merged[, combined_status := factor(combined_status, levels = c("no impact", "depleted", "enriched"))]

# ----------------------------
# Create output-facing alias columns (do NOT rename internal lfc1/lfc2 etc.)
# ----------------------------
merged[, (lfc1_out) := lfc1]
merged[, (lfc2_out) := lfc2]
merged[, (se1_out)  := se1]
merged[, (se2_out)  := se2]
merged[, (status1_out) := status1]
merged[, (status2_out) := status2]

# ----------------------------
# Output
# ----------------------------
out_cols <- c("variant_key", "consequence",
              "oligo_name1", "sequence1",
              "oligo_name2", "sequence2",
              "position1", "position2", "position", "position_source",
              "pam1", "pam2",
              lfc1_out, se1_out, status1_out, "w1",
              lfc2_out, se2_out, status2_out, "w2",
              "combined_LFC", "combined_SE", "combined_Z", "combined_p", "combined_FDR", "combined_status")

for (c in out_cols) if (!(c %in% names(merged))) merged[, (c) := NA]

write_tsv(merged[, ..out_cols], out)
cat("Written combined output to ", out, "\n", sep = "")

# ----------------------------
# Optional diagnostics (ggplot2)
# ----------------------------
diag_prefix_path <- make_diag_prefix_path(out, plot_prefix)

if (!is.null(diag_prefix_path)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping plots. Install with install.packages('ggplot2') if needed.")
  } else {
    suppressPackageStartupMessages(library(ggplot2))
    pref <- diag_prefix_path

    # scatter LFC1 vs LFC2 styled like your SLC2A1 plot:
    # - blue points
    # - red lm line (no SE ribbon)
    # - dashed y=x reference line
    # - Pearson r and p-value annotation
    sc_dt <- merged[is.finite(lfc1) & is.finite(lfc2)]
    if (nrow(sc_dt) > 1) {

      # Pearson correlation (base R)
      cor_test <- suppressWarnings(cor.test(sc_dt$lfc1, sc_dt$lfc2, method = "pearson"))
      cor_value <- round(unname(cor_test$estimate), 3)
      p_value <- signif(cor_test$p.value, 3)

      # Place annotation bottom-right 
      x_pos <- max(sc_dt$lfc1, na.rm = TRUE)
      y_pos <- min(sc_dt$lfc2, na.rm = TRUE)
      
      cap_first <- function(x) paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))

      p1 <- ggplot(sc_dt, aes(x = lfc1, y = lfc2)) +
        geom_point(alpha = 0.6, color = "blue") +
        geom_smooth(method = "lm", color = "red", se = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        theme_minimal(base_size = 10) +
        theme(
        	plot.title = element_text(size = 11),
        	axis.title = element_text(size = 10), 
        	axis.text = element_text(size = 9)
        ) +
        labs(
          title = paste0(cap_first(mode), " ", cap_first(effect), " LFCs Between Pools"),
          x = "Pool 1 LFC",
          y = "Pool 2 LFC"
        ) +
        annotate(
          "text",
          x = x_pos, y = y_pos,
          label = paste0("r = ", cor_value, "\n p = ", p_value),
          hjust = 1, vjust = 0, size = 3.5, color = "black"
        )

      ggsave(filename = paste0(pref, "_lfc_scatter.png"), plot = p1, width = 6, height = 5)
    } else {
      warning("Not enough paired finite LFC values to compute correlation / draw scatter.")
    }

    # combined SE histogram
    if (any(!is.na(merged$combined_SE))) {
      p2 <- ggplot(merged[!is.na(combined_SE)], aes(x = combined_SE)) +
        geom_histogram(bins = 60) +
        labs(title = "Combined SE distribution", x = "combined_SE")
      ggsave(filename = paste0(pref, "_combined_SE_hist.png"), plot = p2, width = 6, height = 4)
    }

    # position discordance table (include selected position + source)
    pd <- merged[!is.na(position1) & !is.na(position2) & position1 != position2,
                 .(variant_key, position1, position2, position, position_source, pam1, pam2)]
    if (nrow(pd) > 0) fwrite(pd, paste0(pref, "_pos_discordant.tsv"), sep = "\t", quote = FALSE)

    cat("Diagnostics written with prefix ", pref, "\n", sep = "")
  }
}