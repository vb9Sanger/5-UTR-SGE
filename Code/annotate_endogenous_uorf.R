#!/usr/bin/env Rscript

# Run example:
# Rscript annotate_endogenous_uorf_minimal.R \
#   CTCF_annotated.tsv \
#   CTCF_with_uorf_status.tsv \
#   CTCF
#
# Note: DDX3X was removed from TARGET_CONFIGS. Its endogenous start codon is a
# non-ATG (GTG) with no real evidence of translation, so there's no reliable
# endogenous uORF to check disruption against for that gene.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# ============================================================
# CONFIG
# ============================================================

TARGET_CONFIGS <- list(

  CTCF = list(
    has_endogenous_uorf = TRUE,

    wt_utr_seq1 = "AATGATTACGGACCTGAAGTCAAAAAATAAGATGCGCTAGTGGACAGATTGCTGACCAGGGGCTTGAGAGCTGGGTTCTATTTTCCCTCCTCAAACTGACTTTGCAGCCACGGAGAG",
    wt_utr_seq2 = "AATGATTACGGACCTGAAGCCAAAGAACAAGATGCGCTAGTGGACAGATTGCTGACCAGGGGCTTGAGAGCTGGGTTCTATTTTCCCTCCTCAAACTGACTTTGCAGCCACGGAGAG",

    upstream_anchor_seq = "TAAAATGATTTGTGCTTTTCTTTTAG",

    downstream_anchor_seq1 = "GTAAGTGCTATTTCCATTTTTTCTATCTTAAAAATGGTAGTAATATAACCCTCCTGGAGTTGGAGGATTT",
    downstream_anchor_seq2 = "GTAAGTGCTATTTCCATTTTTTCTATCTTAAAAATGGTAGTAATATAACCCTCCCGGAATTAGAGGATTT",

    # Known endogenous uORF start (1-based in WT UTR)
    endogenous_uorf_start_local1 = 2L,
    endogenous_uorf_start_local2 = 2L,

    # Expected start codon at the endogenous uORF start
    endogenous_uorf_start_codon1 = "ATG",
    endogenous_uorf_start_codon2 = "ATG"
  )

)

# ============================================================
# Helpers
# ============================================================

choose_final_status <- function(df) {
  status_cols <- c("GMM_status", "combined_status", "stat_pos_raw", "stat_adj_raw")
  present_cols <- status_cols[status_cols %in% names(df)]

  if (length(present_cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }

  out <- as.character(df[[present_cols[1]]])

  if (length(present_cols) > 1) {
    for (col in present_cols[-1]) {
      out <- dplyr::coalesce(out, as.character(df[[col]]))
    }
  }

  out
}

is_missing_string <- function(x) {
  is.na(x) || x == ""
}

find_fixed_once <- function(seq, pattern) {
  if (is_missing_string(seq) || is_missing_string(pattern)) return(NA_integer_)
  hit <- regexpr(pattern, seq, fixed = TRUE)[1]
  if (hit < 0) return(NA_integer_)
  as.integer(hit)
}

find_first_inframe_stop <- function(seq, start_idx, max_end) {
  if (is_missing_string(seq) || is.na(start_idx)) return(NA_integer_)

  stop_codons <- c("TAA", "TAG", "TGA")

  for (i in seq(start_idx + 3L, max_end, by = 3L)) {
    codon <- substr(seq, i, i + 2L)
    if (codon %in% stop_codons) return(i)
  }

  NA_integer_
}

extract_utr_region <- function(seq,
                               upstream_anchor_seq,
                               downstream_anchor_seq) {

  if (is_missing_string(seq)) {
    return(NA_character_)
  }

  left <- 1L
  if (!is_missing_string(upstream_anchor_seq)) {
    hit <- find_fixed_once(seq, upstream_anchor_seq)
    if (!is.na(hit)) {
      left <- hit + nchar(upstream_anchor_seq)
    }
  }

  right <- nchar(seq)
  if (!is_missing_string(downstream_anchor_seq)) {
    hit <- find_fixed_once(seq, downstream_anchor_seq)
    if (!is.na(hit)) {
      right <- hit - 1L
    }
  }

  if (right < left) return(NA_character_)

  substr(seq, left, right)
}

get_context <- function(df, i, cfg) {

  src <- df$annotation_source[i] %||% NA_character_

  if (!is.na(src) && src == "sequence1") {
    return(list(
      sequence = df$sequence1[i] %||% df$sequence[i],
      wt = cfg$wt_utr_seq1,
      down = cfg$downstream_anchor_seq1,
      start = cfg$endogenous_uorf_start_local1,
      start_codon = cfg$endogenous_uorf_start_codon1 %||% "ATG"
    ))
  }

  if (!is.na(src) && src == "sequence2") {
    return(list(
      sequence = df$sequence2[i] %||% df$sequence[i],
      wt = cfg$wt_utr_seq2,
      down = cfg$downstream_anchor_seq2,
      start = cfg$endogenous_uorf_start_local2,
      start_codon = cfg$endogenous_uorf_start_codon2 %||% "ATG"
    ))
  }

  # fallback
  return(list(
    sequence = df$sequence[i],
    wt = cfg$wt_utr_seq2 %||% cfg$wt_utr_seq1,
    down = cfg$downstream_anchor_seq2 %||% cfg$downstream_anchor_seq1,
    start = cfg$endogenous_uorf_start_local2 %||% cfg$endogenous_uorf_start_local1,
    start_codon = cfg$endogenous_uorf_start_codon2 %||% cfg$endogenous_uorf_start_codon1 %||% "ATG"
  ))
}

annotate_row <- function(ctx, cfg) {

  seq <- ctx$sequence
  wt_utr <- ctx$wt
  start <- ctx$start
  start_codon <- ctx$start_codon %||% "ATG"

  if (is_missing_string(seq) || is_missing_string(wt_utr) || is.na(start)) {
    return(tibble(
      endogenous_uorf_status = NA_character_,
      uorf_disrupting = NA,
      endogenous_uorf_wt_seq = NA_character_,
      endogenous_uorf_alt_seq = NA_character_
    ))
  }

  alt_utr <- extract_utr_region(
    seq,
    cfg$upstream_anchor_seq,
    ctx$down
  )

  if (is_missing_string(alt_utr)) {
    return(tibble(
      endogenous_uorf_status = NA_character_,
      uorf_disrupting = NA,
      endogenous_uorf_wt_seq = NA_character_,
      endogenous_uorf_alt_seq = NA_character_
    ))
  }

  # Sanity check: expected WT start codon should be present
  wt_start_ok <- substr(wt_utr, start, start + 2L) == start_codon

  # WT uORF
  wt_stop <- if (wt_start_ok) {
    find_first_inframe_stop(wt_utr, start, nchar(wt_utr) - 2L)
  } else {
    NA_integer_
  }

  wt_seq <- if (!is.na(wt_stop)) {
    substr(wt_utr, start, wt_stop + 2L)
  } else {
    NA_character_
  }

  wt_len <- if (!is.na(wt_stop)) {
    (wt_stop + 2L) - start + 1L
  } else {
    NA_integer_
  }

  # ALT uORF
  alt_start_ok <- substr(alt_utr, start, start + 2L) == start_codon

  alt_stop <- if (alt_start_ok) {
    find_first_inframe_stop(alt_utr, start, nchar(alt_utr) - 2L)
  } else {
    NA_integer_
  }

  alt_seq <- if (!is.na(alt_stop)) {
    substr(alt_utr, start, alt_stop + 2L)
  } else {
    NA_character_
  }

  alt_len <- if (!is.na(alt_stop)) {
    (alt_stop + 2L) - start + 1L
  } else {
    NA_integer_
  }

  status <- case_when(
    !wt_start_ok ~ "wt_start_not_found",
    !alt_start_ok ~ "start_lost",
    is.na(alt_stop) ~ "no_inframe_stop",
    alt_len < wt_len ~ "shortened",
    alt_len > wt_len ~ "lengthened",
    alt_len == wt_len ~ "intact",
    TRUE ~ NA_character_
  )

  disrupting <- status %in% c("start_lost", "no_inframe_stop", "shortened", "lengthened")

  tibble(
    endogenous_uorf_status = status,
    uorf_disrupting = disrupting,
    endogenous_uorf_wt_seq = wt_seq,
    endogenous_uorf_alt_seq = alt_seq
  )
}

# ============================================================
# Main
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript annotate_endogenous_uorf_minimal.R input.tsv output.tsv <target>")
}

input <- args[1]
output <- args[2]
target <- args[3]

if (!target %in% names(TARGET_CONFIGS)) {
  stop(
    "Unknown target: ", target,
    ". Available targets: ", paste(names(TARGET_CONFIGS), collapse = ", ")
  )
}

cfg <- TARGET_CONFIGS[[target]]
df <- read_tsv(input, show_col_types = FALSE)

out <- bind_cols(
  df,
  bind_rows(lapply(seq_len(nrow(df)), function(i) {
    ctx <- get_context(df, i, cfg)
    annotate_row(ctx, cfg)
  }))
)

# Ensure final_status exists and is last column
if (!"final_status" %in% names(out)) {
  out <- out %>%
    mutate(final_status = choose_final_status(.))
}

out <- out %>%
  select(-any_of("final_status"), everything(), final_status)

write_tsv(out, output)
message("Wrote: ", output)