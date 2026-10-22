#!/usr/bin/env Rscript

# Run example:
# Rscript annotate_endogenous_uorf_minimal.R \
#   CTCF_annotated.tsv \
#   CTCF_with_uorf_status.tsv \
#   CTCF
#
# Also supports SON:
# Rscript annotate_endogenous_uorf_minimal.R \
#   SON_annotated.tsv \
#   SON_with_uorf_status.tsv \
#   SON
#
# Note: DDX3X was removed from TARGET_CONFIGS. Its endogenous start codon is a
# non-ATG (GTG) with no real evidence of translation, so there's no reliable
# endogenous uORF to check disruption against for that gene.
#
# Note on SON: this isn't a uORF fully contained in the 5' UTR like CTCF's -
# it's a putative overlapping ORF (oORF) that starts at the very first base
# of exon 1 and reads out of frame straight through the canonical CDS start
# codon, terminating on a TAG a little further into the CDS. The same
# start/stop disruption logic applies either way (start_lost, shortened,
# lengthened, no_inframe_stop/stop_lost, intact) - the only difference is
# that "wt_utr_seq2" here is set to span from the exon-1 start through past
# that in-CDS stop codon, rather than stopping at the UTR/CDS boundary.

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
  ),

  SON = list(
    has_endogenous_uorf = TRUE,

    # SON only has a populated "sequence2" reference context (see notes in
    # annotate_UTR.R) - there is no genuine seq1 side, so those fields are
    # left NA and the script's existing %||% fallbacks handle that cleanly.
    wt_utr_seq1 = NA_character_,

    # Spans from the first base of exon 1 (the oORF start) through the
    # canonical CDS start codon and on to the oORF's own in-frame TAG stop,
    # i.e. wt_utr_seq2 + cds_anchor_seq2 from annotate_UTR.R's SON config:
    #   ATGCTGGGAGCCTGGAGGACTAGCGAGGAGGAGTTGAGAGAACGGAGCGGACGCC (55nt UTR)
    # + ATGGCGACCAACATCGAGCAGATTTTTAGGTCTTTCGTGGTCAGTAAGTTTCGAGAAATTCAACAGGAGCTTTCCAG (start of CDS)
    wt_utr_seq2 = paste0(
      "ATGCTGGGAGCCTGGAGGACTAGCGAGGAGGAGTTGAGAGAACGGAGCGGACGCC",
      "ATGGCGACCAACATCGAGCAGATTTTTAGGTCTTTCGTGGTCAGTAAGTTTCGAGAAATTCAACAGGAGCTTTCCAG"
    ),

    # Same upstream anchor used in annotate_UTR.R for SON - marks the end of
    # exon 1's upstream flank / the start of the oORF.
    upstream_anchor_seq = "TTGGCCGCTGCGCCTCCTCCCGAGGC",

    # No downstream anchor needed: with none supplied, extract_utr_region()
    # simply extends the right boundary to the end of the row's own
    # sequence, which is exactly what we want here since the oORF's stop
    # codon sits inside the CDS, past where a UTR-only downstream anchor
    # would normally cut off.
    downstream_anchor_seq1 = NA_character_,
    downstream_anchor_seq2 = NA_character_,

    # oORF starts at the very first base of exon 1 (local position 1 in
    # wt_utr_seq2 above), not position 2 like CTCF's endogenous uORF.
    endogenous_uorf_start_local1 = NA_integer_,
    endogenous_uorf_start_local2 = 1L,

    endogenous_uorf_start_codon1 = NA_character_,
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

# Shortest anchor fragment accepted when the full anchor cannot be matched.
MIN_ANCHOR_BP <- 10L

# Locate a region boundary from its anchor sequence.
#
# A variant falling inside the anchor itself breaks an exact match. Previously
# the boundary then defaulted to the start (or end) of the oligo, which
# silently returned a window beginning in the wrong place and produced
# spurious start_lost calls for variants that lie outside the region
# altogether. Instead, progressively shorter fragments of the anchor are tried
# from its region-facing end inwards, so a boundary whose outer bases were
# altered can still be found; shortened fragments are accepted only where they
# match in exactly one place. Where no fragment matches, NA is returned and the
# variant is reported as unassessable rather than misclassified.
find_anchor <- function(seq, anchor, from_3prime) {
  n <- nchar(anchor)

  for (k in seq(n, min(MIN_ANCHOR_BP, n), by = -1L)) {
    pat <- if (from_3prime) substr(anchor, n - k + 1L, n) else substr(anchor, 1L, k)
    hits <- gregexpr(pat, seq, fixed = TRUE)[[1]]
    hits <- hits[hits > 0L]

    # the full-length anchor is distinctive enough to take the first match,
    # as before; shortened fragments must be unique
    if (length(hits) >= 1L && (k == n || length(hits) == 1L)) {
      return(if (from_3prime) hits[1] + k else hits[1])
    }
  }

  NA_integer_
}

extract_utr_region <- function(seq,
                               upstream_anchor_seq,
                               downstream_anchor_seq,
                               region_len = NA_integer_) {

  if (is_missing_string(seq)) {
    return(NA_character_)
  }

  left <- if (is_missing_string(upstream_anchor_seq)) 1L else
    find_anchor(seq, upstream_anchor_seq, from_3prime = TRUE)

  right <- if (is_missing_string(downstream_anchor_seq)) nchar(seq) else {
    hit <- find_anchor(seq, downstream_anchor_seq, from_3prime = FALSE)
    if (is.na(hit)) NA_integer_ else hit - 1L
  }

  # A variant lying wholly within one anchor cannot be matched even by its
  # shortest fragment, but by the same token it falls outside the region, which
  # is therefore unchanged. Where the opposite boundary was found, the missing
  # one is placed at the wild-type region length rather than defaulting to the
  # end of the oligo.
  # This only applies where the opposite boundary came from an anchor match of
  # its own; a boundary that merely defaulted to the end of the oligo (as for a
  # target with no downstream anchor) carries no positional information.
  if (!is.na(region_len)) {
    if (is.na(left) && !is.na(right) && !is_missing_string(downstream_anchor_seq)) {
      left <- right - region_len + 1L
    } else if (is.na(right) && !is.na(left) && !is_missing_string(upstream_anchor_seq)) {
      right <- left + region_len - 1L
    }
  }

  if (is.na(left) || is.na(right)) return(NA_character_)
  if (left < 1L || right > nchar(seq) || right < left) return(NA_character_)

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
    ctx$down,
    region_len = nchar(wt_utr)
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