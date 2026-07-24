#!/usr/bin/env Rscript

# signal_to_noise.R
#
# Signal-to-noise QC metrics across analysis conditions, for the 5'UTR pilot
# screen (CTCF, DDX3X, SON, and 1 other gene; 3 of which have 2 guides merged
# via merge_and_combine.R, 1 has a single guide via extract_single_guide.R).
# Each gene's guide(s) were then run through gmm.R, producing *_with_GMM.tsv
# files.
#
# This is intentionally simple: it reports the same per-file metrics as your
# original signal_to_noise.R (median |Z| overall/controls/functional, %
# controls significant, Cohen's d) for every gene x condition combination in
# ONE combined output table, so you can eyeball comparisons across all 4
# genes x up to 3 conditions (Day4_positional, Day4_adjusted,
# Plasmid_positional) yourself. There is NO automatic composite score or
# ranking - deliberately, since a hierarchical rank across adjusted vs
# positional isn't a fair like-for-like comparison (residual positional bias
# is structurally biased toward whichever candidate already had positional
# correction applied), and with only 4 genes manual inspection is more
# transparent than an automated decision layer.
#
# Column selection prefers GMM_status (from gmm.R) for significance, and
# combined_LFC/combined_SE (merged, two-guide genes) falling back to
# pos_adj_log2FoldChange_raw / adj_log2FoldChange_raw + pos_total_se_raw /
# lfcSE_raw (single-guide gene) - i.e. "GMM status and combined LFC, or
# pos_adj LFCs for the single-guide case".
#
# Gene / reference / metric are inferred per file (not passed via a single
# global --metric flag), since one run compares all 3 conditions together:
#   gene   = filename prefix before "_merged" or "_extracted" (e.g. "CTCF", "SON")
#   ref    = "Day4" or "Plasmid", detected from those tokens anywhere in the
#            full file path (case-insensitive)
#   metric = "positional" or "adjusted", detected as a token in the filename
# Override any of these with --manifest (cols: file, gene, ref, metric) if
# auto-detection gets a specific file wrong.
#
# File discovery is RECURSIVE: everything before the first wildcard in
# --input is treated as a base directory, and the last path segment is a
# filename pattern matched at ANY depth beneath it. So "*_with_GMM.tsv" (no
# directory wildcards at all) searches the whole current directory tree, and
# "CTCF/*_with_GMM.tsv" searches recursively under CTCF/ only. You do not
# need to match your actual folder depth with repeated "*/*/*" segments, and
# different genes/conditions can sit at different depths without breaking
# anything.
#
# --dry_run prints the inferred (gene, ref, metric, condition) per file and
# exits without computing any metrics, so you can check parsing (AND check
# it didn't pick up unrelated files elsewhere in the tree) before running the
# full thing.
#
# Usage example:
#   Rscript signal_to_noise.R \
#     --input "*_with_GMM.tsv" \
#     --out_prefix out/utr_signal_noise \
#     --control "SNV" --signal "LOF"
#
#   # Check parsing first:
#   Rscript signal_to_noise.R --input "*_with_GMM.tsv" --dry_run
#
#   # CTCF has no LOF-labeled variants; its closest functional analog is
#   # START_Insertion. Override just CTCF's signal_set with a manifest TSV
#   # (a "gene"-keyed row applies to every one of that gene's files at once,
#   # so you don't need to list all 3 of CTCF's files individually):
#   #   gene<TAB>signal_set
#   #   CTCF<TAB>START_Insertion
#   Rscript signal_to_noise.R --input "*_with_GMM.tsv" --manifest ctcf_override.tsv \
#     --out_prefix out/utr_signal_noise --control "SNV" --signal "LOF"
#
#   # Further restrict CTCF's signal set to only strong-Kozak-context start
#   # insertions (consequence alone doesn't distinguish Kozak strength - that
#   # lives in variant_key, e.g. "...StrongKozak..."), by adding a
#   # signal_key_pattern column to the same manifest:
#   #   gene<TAB>signal_set<TAB>signal_key_pattern
#   #   CTCF<TAB>START_Insertion<TAB>StrongKozak
#   Rscript signal_to_noise.R --input "*_with_GMM.tsv" --manifest ctcf_override.tsv \
#     --out_prefix out/utr_signal_noise --control "SNV" --signal "LOF"
#
#   # If the current directory also has unrelated files/folders (e.g. an
#   # SLC2A1 folder from a different project), check with --dry_run first;
#   # if it picks up files you don't want, point --input at a narrower base
#   # directory instead (glob2rx-based matching does not support brace
#   # expansion like "{CTCF,DDX3X}", so list the genes in the path or run per
#   # gene and combine the resulting .metrics.tsv files afterward).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# ----------------------------
# DEFAULT SETS (UTR pilot - override with --control/--signal, or per-gene via
# --manifest, if a gene's consequence labels differ - e.g. CTCF has no "LOF"
# variants, so its closest functional analog is "START_Insertion". See
# --manifest below for how to set this per-gene.)
# ----------------------------

# ----------------------------
# HELPERS
# ----------------------------

parse_csv_arg <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

# Sys.glob() requires an exact directory-depth match (no recursive "**"), which
# is fragile when different genes/conditions sit at different folder depths
# (e.g. merged_adjusted/ vs sg3/ vs an extracted single-guide folder). This
# instead treats --input as: everything before the first wildcard is a base
# directory, and the last path segment is a filename pattern matched
# RECURSIVELY at any depth under that base directory.
find_input_files <- function(input_pattern) {
  parts <- strsplit(input_pattern, "/", fixed = TRUE)[[1]]
  wildcard_idx <- which(grepl("[*?]", parts))[1]

  if (is.na(wildcard_idx)) {
    # No wildcard at all - treat as a literal path/glob as before.
    return(Sys.glob(input_pattern))
  }

  base_dir <- if (wildcard_idx == 1) "." else paste(parts[seq_len(wildcard_idx - 1)], collapse = "/")
  filename_glob <- parts[length(parts)]
  filename_regex <- utils::glob2rx(filename_glob)

  list.files(path = base_dir, pattern = filename_regex, recursive = TRUE, full.names = TRUE)
}

pick_first_existing <- function(dt, candidates) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_median <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_real_)
  stats::median(v, na.rm = TRUE)
}

# ----------------------------
# METADATA INFERENCE (gene / ref / metric)
# ----------------------------

infer_gene <- function(path, gene_regex = NA_character_) {
  bn <- sub("\\.tsv$", "", basename(path), ignore.case = TRUE)
  bn <- sub("_with_GMM$", "", bn, ignore.case = TRUE)

  if (!is.na(gene_regex) && nzchar(gene_regex)) {
    m <- regmatches(bn, regexec(gene_regex, bn, ignore.case = TRUE))[[1]]
    if (length(m) >= 2 && nzchar(m[2])) return(m[2])
  }

  sub("_(merged|extracted|single).*$", "", bn, ignore.case = TRUE)
}

infer_ref <- function(path, day4_regex, plasmid_regex) {
  if (grepl(plasmid_regex, path, ignore.case = TRUE, perl = TRUE)) return("Plasmid")
  if (grepl(day4_regex,    path, ignore.case = TRUE, perl = TRUE)) return("Day4")
  "Other"
}

infer_metric <- function(path, positional_regex, adjusted_regex) {
  bn <- basename(path)
  is_pos <- grepl(positional_regex, bn, ignore.case = TRUE, perl = TRUE)
  is_adj <- grepl(adjusted_regex,    bn, ignore.case = TRUE, perl = TRUE)
  if (is_pos && !is_adj) return("positional")
  if (is_adj && !is_pos) return("adjusted")
  NA_character_
}

build_file_metadata <- function(files, gene_regex, day4_regex, plasmid_regex,
                                positional_regex, adjusted_regex, manifest_path,
                                default_control_str, default_signal_str,
                                default_signal_key_pattern = "") {
  meta <- data.table(
    file          = files,
    gene          = vapply(files, infer_gene, character(1), gene_regex = gene_regex),
    ref           = vapply(files, infer_ref, character(1),
                           day4_regex = day4_regex, plasmid_regex = plasmid_regex),
    metric        = vapply(files, infer_metric, character(1),
                           positional_regex = positional_regex, adjusted_regex = adjusted_regex),
    control_set_str      = default_control_str,
    signal_set_str        = default_signal_str,
    signal_key_pattern    = default_signal_key_pattern
  )

  # Manifest rows can target a specific file (col "file") OR an entire gene
  # (col "gene", with "file" left blank) - the latter is the convenient way
  # to set e.g. CTCF's signal_set once and have it apply to all of CTCF's
  # conditions, rather than repeating it per file.
  if (!is.null(manifest_path) && nzchar(manifest_path)) {
    if (!file.exists(manifest_path)) stop("Manifest file not found: ", manifest_path)
    man <- fread(manifest_path, colClasses = "character")
    if (!("file" %in% names(man)) && !("gene" %in% names(man))) {
      stop("Manifest must contain a 'file' column and/or a 'gene' column.")
    }

    for (i in seq_len(nrow(man))) {
      key_file <- if ("file" %in% names(man)) man$file[i] else NA_character_
      key_gene <- if ("gene" %in% names(man)) man$gene[i] else NA_character_

      if (!is.na(key_file) && nzchar(key_file)) {
        idx <- which(meta$file == key_file | basename(meta$file) == basename(key_file))
        if (length(idx) == 0) {
          warning("Manifest row not matched to any input file: ", key_file)
          next
        }
      } else if (!is.na(key_gene) && nzchar(key_gene)) {
        idx <- which(meta$gene == key_gene)
        if (length(idx) == 0) {
          warning("Manifest row not matched to any file for gene: ", key_gene)
          next
        }
      } else {
        warning("Manifest row ", i, " has neither 'file' nor 'gene' set; skipped.")
        next
      }

      if ("gene"        %in% names(man) && nzchar(man$gene[i]))          meta$gene[idx]            <- man$gene[i]
      if ("ref"         %in% names(man) && nzchar(man$ref[i]))           meta$ref[idx]             <- man$ref[i]
      if ("metric"      %in% names(man) && nzchar(man$metric[i]))        meta$metric[idx]          <- man$metric[i]
      if ("control_set" %in% names(man) && nzchar(man$control_set[i]))   meta$control_set_str[idx] <- man$control_set[i]
      if ("signal_set"  %in% names(man) && nzchar(man$signal_set[i]))    meta$signal_set_str[idx]  <- man$signal_set[i]
      if ("signal_key_pattern" %in% names(man) && nzchar(man$signal_key_pattern[i]))
        meta$signal_key_pattern[idx] <- man$signal_key_pattern[i]
    }
  }

  meta[, condition := paste0(ref, "_", metric)]
  meta[]
}

# ----------------------------
# COLUMN SELECTION
# ----------------------------

choose_cols <- function(dt, metric) {
  metric <- tolower(metric)

  cfg <- switch(metric,
    pos = , positional = list(
      pfx      = "pos",
      score_fb = "pos_adj_log2FoldChange_raw",
      se_fb    = "pos_total_se_raw",
      fdr_fb   = "pos_adj_fdr_raw",
      stat_fb  = "stat_pos_raw",
      label    = "positional"
    ),
    adj = , adjusted = list(
      pfx      = "adj",
      score_fb = "adj_log2FoldChange_raw",
      se_fb    = "lfcSE_raw",
      fdr_fb   = "adj_fdr_raw",
      stat_fb  = "stat_adj_raw",
      label    = "adjusted"
    ),
    stop("Unknown metric '", metric, "'. Must be 'positional' or 'adjusted'.")
  )

  p <- cfg$pfx
  score_col <- pick_first_existing(dt, c(
    paste0(p, "_combined_LFC"), paste0("combined_LFC_", p), "combined_LFC"
  ))
  se_col <- pick_first_existing(dt, c(
    paste0(p, "_combined_SE"), paste0("combined_SE_", p), "combined_SE"
  ))
  fdr_col <- pick_first_existing(dt, c(
    paste0(p, "_combined_FDR"), paste0("combined_FDR_", p), "combined_FDR",
    cfg$fdr_fb, "padj_raw"
  ))
  stat_col <- pick_first_existing(dt, c(
    "GMM_status",
    paste0(p, "_combined_status"), paste0("combined_status_", p), "combined_status",
    cfg$stat_fb
  ))

  if (is.na(score_col)) score_col <- cfg$score_fb
  if (is.na(se_col))    se_col    <- cfg$se_fb

  if (!all(c(score_col, se_col) %in% names(dt))) {
    stop(
      "Metric '", cfg$label, "' requires score/se columns. Expected either ",
      "combined LFC/SE columns or ", cfg$score_fb, " and ", cfg$se_fb, "."
    )
  }
  list(score = score_col, se = se_col, fdr = fdr_col, stat = stat_col)
}

# ----------------------------
# COHEN'S D
# ----------------------------

cohens_d <- function(x_signal, x_control) {
  x_signal  <- x_signal[is.finite(x_signal)]
  x_control <- x_control[is.finite(x_control)]
  if (length(x_signal) < 2 || length(x_control) < 2) return(NA_real_)

  s1 <- stats::var(x_signal)
  s0 <- stats::var(x_control)
  sp <- sqrt(
    ((length(x_signal)  - 1) * s1 +
     (length(x_control) - 1) * s0) /
    (length(x_signal) + length(x_control) - 2)
  )
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x_signal) - mean(x_control)) / sp
}

# ----------------------------
# PER-FILE SUMMARY
# ----------------------------

summarise_file <- function(dt, gene, ref, metric, condition,
                           control_set, signal_set,
                           signal_key_pattern = "",
                           fdr_cutoff  = 0.05,
                           min_n_for_d = 10) {

  if (!("consequence" %in% names(dt))) {
    stop("Input must contain a 'consequence' column.")
  }

  cols <- choose_cols(dt, metric)
  score_col <- cols$score
  se_col    <- cols$se
  fdr_col   <- cols$fdr
  stat_col  <- cols$stat

  x <- as.data.table(copy(dt))
  x <- x[!is.na(get(score_col)) & !is.na(get(se_col)) & get(se_col) > 0]
  x[, z     := get(score_col) / get(se_col)]
  x[, abs_z := abs(z)]

  median_abs_z <- safe_median(x$abs_z)

  xc <- x[consequence %in% control_set]
  xs <- x[consequence %in% signal_set]

  # Optional additional filter on the signal set, matched against variant_key
  # (e.g. restricting CTCF's START_Insertion signal set to only the
  # strong-Kozak-context variants, whose variant_key contains "StrongKozak",
  # since the consequence label alone does not distinguish Kozak strength).
  if (nzchar(signal_key_pattern)) {
    if (!("variant_key" %in% names(xs))) {
      stop("signal_key_pattern was supplied but input has no 'variant_key' column.")
    }
    xs <- xs[grepl(signal_key_pattern, variant_key, ignore.case = TRUE)]
  }

  n_controls              <- nrow(xc)
  n_signal                <- nrow(xs)
  median_abs_z_controls   <- safe_median(xc$abs_z)
  median_abs_z_functional <- safe_median(xs$abs_z)

  pct_controls_significant <- NA_real_
  control_sig_basis        <- NA_character_

  if (n_controls > 0) {
    if (!is.na(stat_col) && stat_col %in% names(xc)) {
      pct_controls_significant <-
        100 * mean(tolower(trimws(as.character(xc[[stat_col]]))) != "no impact", na.rm = TRUE)
      control_sig_basis <- paste0("stat:", stat_col)
    } else if (!is.na(fdr_col) && fdr_col %in% names(xc)) {
      pct_controls_significant <-
        100 * mean(xc[[fdr_col]] < fdr_cutoff, na.rm = TRUE)
      control_sig_basis <- paste0("fdr:", fdr_col, "<", fdr_cutoff)
    }
  }

  d <- NA_real_
  if (n_signal >= min_n_for_d && n_controls >= min_n_for_d) {
    d <- cohens_d(xs$z, xc$z)
  }

  data.table(
    targeton        = gene,
    condition       = condition,
    ref             = ref,
    metric          = metric,
    input_n         = nrow(dt),
    n_variants_used = nrow(x),
    score_col       = score_col,
    se_col          = se_col,
    fdr_col         = ifelse(is.na(fdr_col),  NA_character_, fdr_col),
    stat_col        = ifelse(is.na(stat_col), NA_character_, stat_col),

    median_abs_z = median_abs_z,

    control_set              = paste(control_set, collapse = ","),
    n_controls               = n_controls,
    median_abs_z_controls    = median_abs_z_controls,
    pct_controls_significant = pct_controls_significant,
    control_sig_basis        = control_sig_basis,

    signal_set                    = paste(signal_set, collapse = ","),
    signal_key_pattern             = ifelse(nzchar(signal_key_pattern), signal_key_pattern, NA_character_),
    n_signal                      = n_signal,
    median_abs_z_functional       = median_abs_z_functional,
    cohens_d_z_signal_vs_control  = d,

    input_file = NA_character_  # filled in by caller
  )
}

# ----------------------------
# CLI
# ----------------------------

opt_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input TSV(s). Glob, e.g. \"*/*/*_with_GMM.tsv\""),
  make_option(c("-o", "--out_prefix"), type = "character", default = "signal_noise",
              help = "Output prefix [default %default]"),

  make_option(c("--manifest"), type = "character", default = "",
              help = paste0(
                "Optional TSV overriding auto-detected metadata. Columns: ",
                "'file' (match a specific file, by exact path or basename) and/or ",
                "'gene' (match ALL files already inferred as this gene - convenient ",
                "for setting one gene's signal_set/control_set at once). Optional ",
                "override columns: gene, ref, metric, control_set, signal_set, ",
                "signal_key_pattern (comma-separated consequence labels for ",
                "control_set/signal_set, e.g. for CTCF's lack of LOF variants: ",
                "gene=CTCF, signal_set=START_Insertion; signal_key_pattern is a ",
                "regex/substring matched against variant_key, applied on top of ",
                "signal_set, e.g. signal_key_pattern=StrongKozak)"
              )),
  make_option(c("--dry_run"), action = "store_true", default = FALSE,
              help = "Print inferred (gene, ref, metric, condition) per file and exit, without computing metrics"),

  make_option(c("--gene_regex"), type = "character", default = "",
              help = "Optional regex with 1 capture group to extract gene from basename (overrides default prefix logic)"),
  make_option(c("--ref_day4_regex"), type = "character", default = "day4|d4",
              help = "Regex (case-insensitive) matched against full path to detect Day4 reference [default %default]"),
  make_option(c("--ref_plasmid_regex"), type = "character", default = "plasmid",
              help = "Regex (case-insensitive) matched against full path to detect Plasmid reference [default %default]"),
  make_option(c("--metric_positional_regex"), type = "character", default = "positional",
              help = "Regex (case-insensitive) matched against filename to detect positional metric [default %default]"),
  make_option(c("--metric_adjusted_regex"), type = "character", default = "adjusted",
              help = "Regex (case-insensitive) matched against filename to detect adjusted metric [default %default]"),

  make_option(c("--control"), type = "character", default = "SNV",
              help = "Comma-separated neutral/control consequence labels; global default, override per-gene via --manifest [default %default]"),
  make_option(c("--signal"), type = "character", default = "LOF",
              help = "Comma-separated signal consequence labels used for Cohen's d; global default, override per-gene via --manifest [default %default]"),
  make_option(c("--signal_key_pattern"), type = "character", default = "",
              help = paste0(
                "Optional regex/substring matched against variant_key, applied as an ",
                "ADDITIONAL filter on top of --signal (consequence-based). Global ",
                "default is empty (no extra filtering); override per-gene via ",
                "--manifest's signal_key_pattern column, e.g. for CTCF: ",
                "gene=CTCF, signal_key_pattern=StrongKozak, to restrict its ",
                "START_Insertion signal set to strong-Kozak-context variants only ",
                "[default none]"
              )),
  make_option(c("--fdr_cutoff"), type = "double", default = 0.05,
              help = "FDR cutoff for fallback significance calling [default %default]"),
  make_option(c("--min_n_for_d"), type = "integer", default = 10,
              help = "Minimum n per group for Cohen's d [default %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$input) || !nzchar(opt$input)) stop("Please provide --input")

files <- find_input_files(opt$input)
if (length(files) == 0) stop(paste0("No files matched: ", opt$input))

gene_regex_arg <- if (nzchar(opt$gene_regex)) opt$gene_regex else NA_character_

meta <- build_file_metadata(
  files               = files,
  gene_regex          = gene_regex_arg,
  day4_regex          = opt$ref_day4_regex,
  plasmid_regex       = opt$ref_plasmid_regex,
  positional_regex    = opt$metric_positional_regex,
  adjusted_regex      = opt$metric_adjusted_regex,
  manifest_path       = opt$manifest,
  default_control_str = opt$control,
  default_signal_str  = opt$signal,
  default_signal_key_pattern = opt$signal_key_pattern
)

message("Inferred metadata per file:")
print(
  meta[, .(file = basename(file), gene, ref, metric, condition,
           control_set = control_set_str, signal_set = signal_set_str,
           signal_key_pattern)],
  row.names = FALSE
)

bad_metric <- meta[is.na(metric)]
if (nrow(bad_metric) > 0) {
  stop(
    "Could not determine metric (positional/adjusted) for the following file(s). ",
    "Fix the filename, adjust --metric_positional_regex/--metric_adjusted_regex, ",
    "or supply --manifest to override:\n",
    paste0("  ", bad_metric$file, collapse = "\n")
  )
}

bad_ref <- meta[ref == "Other"]
if (nrow(bad_ref) > 0) {
  warning(
    "Could not determine reference (Day4/Plasmid) for the following file(s); ",
    "they are labeled 'Other' and will still be processed but won't appear ",
    "in the standard Day4/Plasmid comparisons:\n",
    paste0("  ", bad_ref$file, collapse = "\n"),
    call. = FALSE
  )
}

if (isTRUE(opt$dry_run)) {
  message("\n--dry_run set: exiting without computing metrics.")
  quit(save = "no", status = 0)
}

all_summ <- vector("list", nrow(meta))

for (i in seq_len(nrow(meta))) {
  f  <- meta$file[i]
  dt <- fread(f)

  s <- summarise_file(
    dt          = dt,
    gene        = meta$gene[i],
    ref         = meta$ref[i],
    metric      = meta$metric[i],
    condition   = meta$condition[i],
    control_set = parse_csv_arg(meta$control_set_str[i]),
    signal_set  = parse_csv_arg(meta$signal_set_str[i]),
    signal_key_pattern = meta$signal_key_pattern[i],
    fdr_cutoff  = opt$fdr_cutoff,
    min_n_for_d = opt$min_n_for_d
  )

  s[, input_file := f]
  all_summ[[i]] <- s
}

summary_dt <- rbindlist(all_summ, fill = TRUE)
setorder(summary_dt, targeton, ref, metric)

out_tsv <- paste0(opt$out_prefix, ".metrics.tsv")
fwrite(summary_dt, out_tsv, sep = "\t")
message("Wrote: ", out_tsv)

# ----------------------------
# PLOT: median |Z| by condition, one panel per gene
# ----------------------------

plot_dt <- summary_dt[is.finite(median_abs_z)]

if (nrow(plot_dt) > 0) {
  p <- ggplot(plot_dt, aes(x = condition, y = median_abs_z, color = metric, shape = ref)) +
    geom_point(size = 2.5) +
    facet_wrap(~ targeton, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Signal-to-noise summary: median |Z| by condition",
      x     = "Condition",
      y     = "median(|Z|)",
      color = "Metric",
      shape = "Reference"
    )

  out_png <- paste0(opt$out_prefix, ".median_abs_z.png")
  ggsave(
    filename = out_png, plot = p,
    width  = max(8, 2.2 * length(unique(plot_dt$condition))),
    height = max(5, 2.5 * ceiling(length(unique(plot_dt$targeton)) / 2)),
    units  = "in", dpi = 300
  )
  message("Wrote: ", out_png)
}
