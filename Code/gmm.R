#!/usr/bin/env Rscript
# gmm.R
#
# Fits a Gaussian Mixture Model (mclust) on an LFC column from either:
#  - merged output (combined_LFC + combined_status)
#  - single-guide output (one of: adj_log2FoldChange_raw/shrunk, pos_adj_log2FoldChange_raw/shrunk
#    with matching status column: stat_adj_raw/shrunk or stat_pos_raw/shrunk)
#
# Key behavior:
#  - G is fixed to 2 (no CLI option)
#  - ALWAYS writes a diagnostic plot (requires ggplot2 + ragg)
#  - Fits GMM only on rows with status in {depleted, no impact}
#  - Leaves enriched variants alone (not fit; status unchanged)
#  - Output = original columns + new GMM columns appended
#
# Dependencies: data.table, mclust, ggplot2, ragg
#
# Usage:
#   Rscript gmm.R -f in.tsv [-o out.tsv] [--lfc_col COL] [--status_col COL]
#     [--seed 1] [--plot_prefix prefix]
#
# Notes:
#  - If --plot_prefix is not provided, plot is written next to output as:
#      <out_base>_GMM_plot.png
#
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Requires data.table")
  if (!requireNamespace("mclust", quietly = TRUE)) stop("Requires mclust")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Requires ggplot2")
  if (!requireNamespace("ragg", quietly = TRUE)) stop("Requires ragg")
})

library(data.table)
library(mclust)

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat(
"gmm_on_sge_output.R

Required:
  -f <input.tsv>

Optional:
  -o <output.tsv>        Default: <input_basename>_with_GMM.tsv (in same directory)
  --lfc_col <name>       LFC column (auto-detect if omitted)
  --status_col <name>    Status column (auto-detect if omitted; can be derived from lfc_col for single-guide)
  --seed <int>           RNG seed (default: 1)
  --plot_prefix <path>   Plot prefix (default: next to output: <out_base>_GMM)

Behavior:
  - G is fixed to 2
  - GMM fit set = rows with status in {depleted, no impact}
  - enriched rows are not fit and left unchanged

Examples:
  # merged output (auto-detect combined_LFC + combined_status)
  Rscript gmm_on_sge_output.R -f merged.tsv

  # single-guide output (explicit)
  Rscript gmm_on_sge_output.R -f single.tsv --lfc_col pos_adj_log2FoldChange_raw --status_col stat_pos_raw

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

as_int_flag <- function(x, default) {
  if (is.null(x)) return(default)
  suppressWarnings(v <- as.integer(x))
  if (is.na(v)) stop("Invalid integer value: ", x)
  v
}

if (has_flag("--help") || has_flag("-h")) {
  print_help()
  quit(save = "no", status = 0)
}

# ---- Required
f_in <- get_flag("-f")
if (is.null(f_in)) {
  print_help()
  stop("Missing required arg: -f <input.tsv>")
}
if (!file.exists(f_in)) stop("Input file not found: ", f_in)

# ---- Optional
out <- get_flag("-o", NULL)
if (is.null(out)) {
  base <- sub("\\.tsv$", "", basename(f_in), ignore.case = TRUE)
  out <- file.path(dirname(f_in), paste0(base, "_with_GMM.tsv"))
}

lfc_col_user <- get_flag("--lfc_col", NULL)
status_col_user <- get_flag("--status_col", NULL)

seed <- as_int_flag(get_flag("--seed", NULL), 1L)
plot_prefix_user <- get_flag("--plot_prefix", NULL)

# ---- Load
dt <- tryCatch(
  fread(f_in, sep = "\t", header = TRUE, data.table = TRUE),
  error = function(e) stop("Failed reading TSV: ", f_in, "\n", conditionMessage(e))
)
setDT(dt)

first_present <- function(cands, nms) {
  hit <- cands[cands %in% nms]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

detect_lfc_col <- function(nms) {
  cands <- c(
    "combined_LFC",
    "pos_adj_log2FoldChange_raw",
    "pos_adj_log2FoldChange_shrunk",
    "adj_log2FoldChange_raw",
    "adj_log2FoldChange_shrunk"
  )
  first_present(cands, nms)
}

detect_status_col <- function(nms) {
  cands <- c(
    "combined_status",
    "stat_pos_raw",
    "stat_pos_shrunk",
    "stat_adj_raw",
    "stat_adj_shrunk"
  )
  first_present(cands, nms)
}

derive_status_from_lfc <- function(lfc_col) {
  # Only for the 4 single-guide possibilities
  map <- list(
    "pos_adj_log2FoldChange_raw"    = "stat_pos_raw",
    "pos_adj_log2FoldChange_shrunk" = "stat_pos_shrunk",
    "adj_log2FoldChange_raw"        = "stat_adj_raw",
    "adj_log2FoldChange_shrunk"     = "stat_adj_shrunk"
  )
  if (!is.null(map[[lfc_col]])) return(map[[lfc_col]])
  NULL
}

lfc_col <- if (!is.null(lfc_col_user)) lfc_col_user else detect_lfc_col(names(dt))
if (is.null(lfc_col) || !(lfc_col %in% names(dt))) {
  stop("Could not determine LFC column. Provide --lfc_col <colname>.")
}

# status col: prefer user, else:
#  - merged: combined_status if present
#  - single: derive from lfc_col, else detect from file
status_col <- NULL
if (!is.null(status_col_user)) {
  status_col <- status_col_user
} else {
  if ("combined_LFC" %in% names(dt) && lfc_col == "combined_LFC") {
    status_col <- if ("combined_status" %in% names(dt)) "combined_status" else detect_status_col(names(dt))
  } else {
    derived <- derive_status_from_lfc(lfc_col)
    if (!is.null(derived) && derived %in% names(dt)) {
      status_col <- derived
    } else {
      status_col <- detect_status_col(names(dt))
    }
  }
}
if (is.null(status_col) || !(status_col %in% names(dt))) {
  stop("Could not determine status column. Provide --status_col <colname>.")
}

# ---- Fit settings
fit_status_norm <- c("depleted", "no impact")  # fixed by your requirement

status_raw <- as.character(dt[[status_col]])
status_norm <- tolower(trimws(status_raw))

fit_set <- status_norm %in% fit_status_norm
# enriched are excluded from fit set and left unchanged
fit_set <- fit_set & (status_norm != "enriched")

x_all <- suppressWarnings(as.numeric(dt[[lfc_col]]))
fit_idx <- which(fit_set & !is.na(x_all) & is.finite(x_all))

if (length(fit_idx) < 20) {
  stop("Too few rows to fit GMM (need >= 20). Rows eligible for fit: ", length(fit_idx))
}

# ---- GMM (G fixed to 2)
G_fixed <- 2L
x <- x_all[fit_idx]

set.seed(seed)
gmm <- mclust::Mclust(x, G = G_fixed)

K <- ncol(gmm$z)
prob_cols <- paste0("GMM_prob_cluster", seq_len(K))

means <- as.numeric(gmm$parameters$mean)

# Label components (only depleted vs no impact)
label_map <- rep("no impact", K)
depleted_comp <- which.min(means)
label_map[depleted_comp] <- "depleted"

# ---- Add output cols
dt[, GMM_fit_set := FALSE]
dt[, GMM_cluster := as.integer(NA)]
dt[, GMM_label := as.character(NA)]
for (pc in prob_cols) dt[, (pc) := as.numeric(NA)]

dt[fit_idx, GMM_fit_set := TRUE]
dt[fit_idx, GMM_cluster := gmm$classification]
dt[fit_idx, GMM_label := label_map[gmm$classification]]
for (k in seq_len(K)) {
  dt[fit_idx, (prob_cols[k]) := gmm$z[, k]]
}

# ---- Sign correction ----
# Mclust(x, G = 2) is not constrained to equal variances, so it can select an
# unequal-variance model where the "depleted" (lower-mean) component is wide
# enough to absorb a strongly positive LFC outlier, purely because that point
# scores a higher posterior probability there than under the tight "no impact"
# component. This is not biologically valid: a positive LFC cannot represent
# depletion. Any fitted row labeled "depleted" with LFC > 0 is reassigned to
# "no impact" here, deterministically and identically for every row.
# GMM_cluster and GMM_prob_cluster* are left untouched (raw model output,
# kept for auditability); only GMM_label (and, downstream, GMM_status) reflect
# the override. GMM_sign_override flags exactly which rows were reassigned.
dt[, GMM_sign_override := FALSE]
sign_override_idx <- fit_idx[dt$GMM_label[fit_idx] == "depleted" & x_all[fit_idx] > 0]
if (length(sign_override_idx) > 0) {
  dt[sign_override_idx, GMM_sign_override := TRUE]
  dt[sign_override_idx, GMM_label := "no impact"]
}

# Post-GMM status:
# - default to original
# - for fitted rows, override with GMM label
dt[, GMM_status := status_raw]
dt[fit_idx, GMM_status := GMM_label]

# Status changed flag
dt[, GMM_status_changed := FALSE]
orig_norm <- tolower(trimws(as.character(status_raw)))
new_norm  <- tolower(trimws(as.character(dt$GMM_status)))
dt[!is.na(orig_norm) & !is.na(new_norm) & (orig_norm != new_norm), GMM_status_changed := TRUE]

# ---- Write output (all original cols + appended new cols)
fwrite(dt, file = out, sep = "\t", quote = FALSE)
cat("Wrote: ", out, "\n", sep = "")
cat("Used LFC column: ", lfc_col, "\n", sep = "")
cat("Used status column: ", status_col, "\n", sep = "")
cat("Fit rows: ", length(fit_idx), " / ", nrow(dt), "\n", sep = "")
cat("Sign-corrected rows (depleted -> no impact, LFC > 0): ", length(sign_override_idx),
    " / ", length(fit_idx), " fit rows\n", sep = "")
cat("G (fixed): ", gmm$G, "   model: ", gmm$modelName, "\n", sep = "")
cat("Component means (fit set only):\n")
print(data.frame(component = seq_len(K), mean_LFC = means, label = label_map)[order(means), ], row.names = FALSE)

# ---- ALWAYS write plot
# If plot_prefix not provided, default next to output with base name
out_base <- sub("\\.tsv$", "", basename(out), ignore.case = TRUE)
default_pref <- file.path(dirname(out), paste0(out_base, "_GMM"))
plot_prefix <- if (!is.null(plot_prefix_user)) plot_prefix_user else default_pref
plot_file <- paste0(plot_prefix, "_plot.png")

plot_dt <- copy(dt)
plot_dt[, lfc_val := suppressWarnings(as.numeric(get(lfc_col)))]
plot_dt <- plot_dt[!is.na(lfc_val) & !is.na(GMM_cluster) & GMM_fit_set == TRUE]

if (nrow(plot_dt) > 0) {
  plot_dt <- plot_dt[order(lfc_val)]
  plot_dt[, variant_index := seq_len(.N)]

  # Cluster mean lines reflect the raw (uncorrected) GMM components, since
  # those are the model's actual fitted components - the sign correction is
  # a labeling override on top, not a re-fit.
  plot_dt[, cluster_label := paste0("Cluster ", GMM_cluster)]

  # Points get a third category for sign-corrected rows, so they're visually
  # distinguishable from both raw clusters even though GMM_cluster (used for
  # the mean lines above) still reflects their original raw assignment.
  override_label <- "Sign-corrected (was depleted, LFC > 0)"
  plot_dt[, plot_group := fifelse(GMM_sign_override == TRUE, override_label, cluster_label)]

  all_levels <- c(sort(unique(as.character(plot_dt$cluster_label))), override_label)
  plot_dt[, cluster_label := factor(cluster_label, levels = all_levels)]
  plot_dt[, plot_group := factor(plot_group, levels = all_levels)]

  means_dt <- plot_dt[, .(cluster_mean = mean(lfc_val, na.rm = TRUE)), by = cluster_label]

  cluster_levels <- setdiff(all_levels, override_label)
  point_colors <- setNames(c("#1B9E77", "#7570B3")[seq_along(cluster_levels)], cluster_levels)
  point_colors[override_label] <- "#E41A1C"

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = variant_index, y = lfc_val, color = plot_group)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::geom_hline(
      data = means_dt,
      ggplot2::aes(yintercept = cluster_mean, color = cluster_label),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::scale_color_manual(values = point_colors, drop = !any(plot_dt$GMM_sign_override)) +
    ggplot2::labs(
      title = paste0("GMM (G=2) on ", lfc_col, " | fit: depleted + no impact"),
      x = "Row (sorted by LFC within fit set)",
      y = lfc_col,
      color = "GMM cluster"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(plot_file, p, width = 10, height = 5, dpi = 300, device = ragg::agg_png)
  cat("Wrote plot: ", plot_file, "\n", sep = "")
} else {
  warning("No fitted rows available for plotting after filtering (unexpected). Plot not written.")
}