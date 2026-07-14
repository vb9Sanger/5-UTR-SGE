#!/usr/bin/env Rscript

# Ex:
# Rscript annotate_sge_utr_simple.R input.tsv output.tsv EXAMPLE_CODING_TARGET

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# ============================================================
# CONFIG
# Define one entry per target/gene.
#
# Sequence contexts:
#   - cds_anchor_seq1 / cds_anchor_seq2 define where the CDS begins
#   - wt_utr_seq1 / wt_utr_seq2 define the corresponding WT UTR sequence
#
# For merged files:
#   - sequence2 is preferred
#   - cds_anchor_seq2 and wt_utr_seq2 are used when sequence2 is present
#   - cds_anchor_seq1 and wt_utr_seq1 are used only if sequence2 is missing
#
# For single-column inputs (sequence + oligo_name):
#   - the script still prefers seq2-style references if available
#   - if seq2 config values exist, cds_anchor_seq2 / wt_utr_seq2 are used
#   - otherwise cds_anchor_seq1 / wt_utr_seq1 are used
#
# Note on single-column input:
#   - If a target's config has BOTH seq1 and seq2 genuinely populated (as is
#     the case for CTCF, DDX3X, and SLC2A1 - only SON has just one side
#     populated), single-column "sequence" input is ambiguous and raises an
#     explicit error, since there's no way to tell which guide's reference
#     context a given row belongs to. Use merged sequence1/sequence2 input
#     for these targets instead. (An older version of this script
#     disambiguated single-column SLC2A1 input specifically via "sg5"/"sg6"
#     in oligo_name - that workflow is no longer used, and the same ambiguity
#     risk existed for CTCF/DDX3X too, just unguarded.)
#
# Notes:
#   - wt_utr_seq1 / wt_utr_seq2 are only needed for targets with endogenous
#     UTR ATGs/uORFs (used to detect gained UTR ATGs)
#
#   - canonical_kozak_15bp should be the reference 15 bp surrounding
#     the canonical CDS ATG:
#       9 nt upstream + ATG + 3 nt downstream
#
#   - if needed, this can be made guide-specific using
#       canonical_kozak_15bp1 / canonical_kozak_15bp2
#
#   - for targets without a captured CDS anchor, frame can be derived from
#     a downstream anchor found in the observed oligo:
#       * downstream_anchor_seq1 / downstream_anchor_seq2 locate a real
#         downstream position in the oligo
#       * downstream_frame_offset1 / downstream_frame_offset2 specify the
#         number of nt from the start of that anchor to the true frame
#         reference position
#
#     Example:
#       if the first base of downstream_anchor_seq is already in the desired
#       CDS frame, use downstream_frame_offset = 0L
#       if the frame reference is 2 nt after the anchor start, use 2L
#
#   - the script automatically adds a 'final_status' column at output,
#     prioritising the following columns if present:
#       GMM_status -> combined_status -> stat_pos_raw -> stat_adj_raw
# ============================================================

# ============================================================
# CONFIG
# ============================================================

TARGET_CONFIGS <- list(

  DDX3X = list(
    wt_utr_seq1 = NA_character_,
    wt_utr_seq2 = NA_character_,
    wt_utr_seq = NA_character_,  # backwards compatibility fallback
    has_endogenous_utr_atg = FALSE,

    cds_anchor_seq1 = "ATGAGTCATGTGGCAGTGGAAAATGCGCTCGGGCTGGACCAGCAG",
    cds_anchor_seq2 = "ATGAGTCATGTGGCAGTGGAAAACGCGCTCGGACTGGACCAGCAG",

    upstream_anchor_seq = "GGCTTTCCAGCGGGTATATTAGATCCGTGGCCGCGCGGTGCGCTCC",
    downstream_anchor_seq1 = NA_character_,
    downstream_anchor_seq2 = NA_character_,
    downstream_anchor_seq = NA_character_,
    downstream_frame_offset1 = NA_integer_,
    downstream_frame_offset2 = NA_integer_,
    downstream_frame_offset = NA_integer_,

    canonical_kozak_15bp = "TCTTCAGGGATGAGT",
    canonical_kozak_15bp1 = NA_character_,
    canonical_kozak_15bp2 = NA_character_
  ),

  SLC2A1 = list(
    wt_utr_seq1 = "AAGAGGCAAGAGGTAGCAACCGCGAGCGTGCCGGTCGCTAGTCGCGGGTCCCCGAGTGAGCACGCCAGGGAGCAGGAGAGCAAACAACGGAGGTCGGAGTCAGAGTCGCAGTGGGAGTCCCCGGACCGGAGCACGAGCCTGAGCGGGAGAGCGCCGCTCGCACGCCCGTCGCCACCCGCGTACCCGGCGCAGCCAGAGCCACCAGCGCAGCGCTGCC",
    wt_utr_seq2 = "AAGAGGCAAGAGGTAGCAACCGCGAGCGTGCCGGTCGCTAGTCGCGGGTCCCCGAGTGAGCACGCCAGGGAGCAGGAGACCAAACGACGGGGGTCGGAGTCAGAGTCGCAGTGGGAGTCCCCGGACCGGAGCACGAGCCTGAGCGGGAGAGCGCCGCTCGCACGCCCGTCGCCACCCGCGTACCCGGCGCAGCCAGAGCCACCAGCGCAGCGCTGCC",
    wt_utr_seq = NA_character_,  # backwards compatibility fallback
    has_endogenous_utr_atg = FALSE,

    cds_anchor_seq1 = "ATGGAGCCCAGCAGCAAG",
    cds_anchor_seq2 = "ATGGAGCCTAGTAGTAAG",

    upstream_anchor_seq = NA_character_,
    downstream_anchor_seq1 = NA_character_,
    downstream_anchor_seq2 = NA_character_,
    downstream_anchor_seq = NA_character_,
    downstream_frame_offset1 = NA_integer_,
    downstream_frame_offset2 = NA_integer_,
    downstream_frame_offset = NA_integer_,

    canonical_kozak_15bp = "AGCGCTGCCATGGAG",
    canonical_kozak_15bp1 = NA_character_,
    canonical_kozak_15bp2 = NA_character_
  ),

  SON = list(
    wt_utr_seq1 = NA_character_,
    wt_utr_seq2 = "ATGCTGGGAGCCTGGAGGACTAGCGAGGAGGAGTTGAGAGAACGGAGCGGACGCC",
    wt_utr_seq = NA_character_,  # backwards compatibility fallback
    has_endogenous_utr_atg = TRUE,

    cds_anchor_seq1 = NA_character_,
    cds_anchor_seq2 = "ATGGCGACCAACATCGAGCAGATTTTTAGGTCTTTCGTGGTCAGTAAGTTTCGAGAAATTCAACAGGAGCTTTCCAG",

    upstream_anchor_seq = "TTGGCCGCTGCGCCTCCTCCCGAGGC",
    downstream_anchor_seq1 = NA_character_,
    downstream_anchor_seq2 = NA_character_,
    downstream_anchor_seq = NA_character_,
    downstream_frame_offset1 = NA_integer_,
    downstream_frame_offset2 = NA_integer_,
    downstream_frame_offset = NA_integer_,

    canonical_kozak_15bp = "GCGGACGCCATGGCG",
    canonical_kozak_15bp1 = NA_character_,
    canonical_kozak_15bp2 = NA_character_
  ),

  CTCF = list(
    wt_utr_seq1 = "AATGATTACGGACCTGAAGTCAAAAAATAAGATGCGCTAGTGGACAGATTGCTGACCAGGGGCTTGAGAGCTGGGTTCTATTTTCCCTCCTCAAACTGACTTTGCAGCCACGGAGAG",
    wt_utr_seq2 = "AATGATTACGGACCTGAAGCCAAAGAACAAGATGCGCTAGTGGACAGATTGCTGACCAGGGGCTTGAGAGCTGGGTTCTATTTTCCCTCCTCAAACTGACTTTGCAGCCACGGAGAG",
    wt_utr_seq = NA_character_,
    has_endogenous_utr_atg = TRUE,

    cds_anchor_seq1 = NA_character_,
    cds_anchor_seq2 = NA_character_,

    upstream_anchor_seq = "TAAAATGATTTGTGCTTTTCTTTTAG",
    downstream_anchor_seq1 = "GTAAGTGCTATTTCCATTTTTTCTATCTTAAAAATGGTAGTAATATAACCCTCCTGGAGTTGGAGGATTT",
    downstream_anchor_seq2 = "GTAAGTGCTATTTCCATTTTTTCTATCTTAAAAATGGTAGTAATATAACCCTCCCGGAATTAGAGGATTT",
    downstream_anchor_seq = NA_character_,
    downstream_frame_offset1 = 0L,
    downstream_frame_offset2 = 0L,
    downstream_frame_offset = NA_integer_,

    canonical_kozak_15bp = NA_character_,
    canonical_kozak_15bp1 = NA_character_,
    canonical_kozak_15bp2 = NA_character_
  )
)

# ============================================================
# Helpers
# ============================================================

is_missing_string <- function(x) {
  is.na(x) || x == ""
}

find_fixed_once <- function(seq, pattern) {
  if (is_missing_string(seq) || is_missing_string(pattern)) return(NA_integer_)
  hit <- regexpr(pattern, seq, fixed = TRUE)[1]
  if (hit < 0) return(NA_integer_)
  as.integer(hit)
}

find_all_codons <- function(seq, codon) {
  if (is_missing_string(seq) || nchar(seq) < 3) return(integer(0))
  hits <- gregexpr(codon, seq, fixed = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) return(integer(0))
  as.integer(hits)
}

get_subseq_safe <- function(seq, start, end) {
  if (is_missing_string(seq) || is.na(start) || is.na(end)) return(NA_character_)
  if (start < 1 || end > nchar(seq) || start > end) return(NA_character_)
  substr(seq, start, end)
}

get_kozak15 <- function(seq, start_pos) {
  # 9 upstream + start codon (3 nt) + 3 downstream = 15 bp
  get_subseq_safe(seq, start_pos - 9L, start_pos + 5L)
}

classify_kozak_from_oligo <- function(oligo_name) {
  x <- oligo_name %||% ""
  case_when(
    str_detect(x, regex("StrongKozak", ignore_case = TRUE)) ~ "Strong",
    str_detect(x, regex("WeakKozak", ignore_case = TRUE)) ~ "Weak",
    TRUE ~ NA_character_
  )
}

frame_relative_to_reference <- function(start_idx, reference_idx) {
  if (is.na(start_idx) || is.na(reference_idx)) return(NA_integer_)
  ((start_idx - reference_idx) %% 3L) + 1L
}

find_first_inframe_stop <- function(seq, start_idx, max_end = NA_integer_) {
  if (is_missing_string(seq) || is.na(start_idx)) return(NA_integer_)

  stop_codons <- c("TAA", "TAG", "TGA")
  seq_len <- nchar(seq)

  if (is.na(max_end)) {
    max_end <- seq_len - 2L
  } else {
    max_end <- min(max_end, seq_len - 2L)
  }

  first_codon <- start_idx + 3L
  if (first_codon > max_end) return(NA_integer_)

  for (i in seq(first_codon, max_end, by = 3L)) {
    codon <- substr(seq, i, i + 2L)
    if (codon %in% stop_codons) return(i)
  }

  NA_integer_
}

has_valid_orf <- function(seq, start_idx, max_end = NA_integer_) {
  if (is_missing_string(seq) || is.na(start_idx)) return(NA)

  stop_idx <- find_first_inframe_stop(
    seq = seq,
    start_idx = start_idx,
    max_end = max_end
  )

  if (is.na(stop_idx)) return(FALSE)

  # Require at least one sense codon between start and stop.
  # Excludes immediate start-stop cases:
  #   ATGTGA, ATGTAA, ATGTAG
  (stop_idx - start_idx) >= 6L
}

extract_orf_sequence <- function(seq, start_idx, max_end = NA_integer_) {
  stop_idx <- find_first_inframe_stop(seq, start_idx, max_end = max_end)

  if (is.na(stop_idx)) return(NA_character_)

  # Exclude start-stop mini-ORFs from returned ORF sequence as well.
  if ((stop_idx - start_idx) < 6L) return(NA_character_)

  substr(seq, start_idx, stop_idx + 2L)
}

is_startstop_at <- function(seq, start_idx) {
  if (is_missing_string(seq) || is.na(start_idx)) return(NA)

  stop_codons <- c("TAA", "TAG", "TGA")
  seq_len <- nchar(seq)

  if ((start_idx + 5L) > seq_len) return(FALSE)

  next_codon <- substr(seq, start_idx + 3L, start_idx + 5L)
  next_codon %in% stop_codons
}

score_kozak <- function(kozak15) {
  if (length(kozak15) == 0) return(numeric(0))

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    return(rep(NA_real_, length(kozak15)))
  }

  # ORFik human PFM
  pfm <- t(matrix(
    as.integer(c(
      20,20,21,21,19,24,46,29,19,22,28,16,
      27,33,32,23,32,38,10,38,45,15,39,26,
      35,29,28,39,30,26,37,20,28,49,18,37,
      18,18,19,17,19,12, 7,13, 8,14,15,21
    )),
    ncol = 4
  ))
  rownames(pfm) <- c("A", "C", "G", "T")

  pwm <- Biostrings::PWM(pfm)

  out <- rep(NA_real_, length(kozak15))

  ok <- !is.na(kozak15) & nchar(kozak15) == 15L
  if (!any(ok)) return(out)

  # Remove ATG (positions 10–12), keep 12 flanking bases
  flank12 <- paste0(substr(kozak15[ok], 1, 9),
                    substr(kozak15[ok], 13, 15))

  raw_scores <- vapply(
    flank12,
    function(s) {
      if (nchar(s) != 12L) return(NA_real_)
      as.numeric(Biostrings::PWMscoreStartingAt(pwm, s, starting.at = 1))
    },
    numeric(1)
  )

  out[ok] <- raw_scores
  out
}

extract_utr_region <- function(seq,
                               upstream_anchor_seq = NA_character_,
                               cds_anchor_seq = NA_character_,
                               downstream_anchor_seq = NA_character_,
                               downstream_frame_offset = NA_integer_) {
  if (is_missing_string(seq)) {
    return(list(
      utr_seq = NA_character_,
      utr_start_in_oligo = NA_integer_,
      utr_end_in_oligo = NA_integer_,
      cds_atg_idx = NA_integer_,
      cds_end_idx = NA_integer_,
      downstream_anchor_idx = NA_integer_,
      downstream_frame_ref_idx = NA_integer_
    ))
  }

  left_bound <- 1L
  if (!is_missing_string(upstream_anchor_seq)) {
    up_hit <- find_fixed_once(seq, upstream_anchor_seq)
    if (!is.na(up_hit)) {
      left_bound <- up_hit + nchar(upstream_anchor_seq)
    }
  }

  cds_atg_idx <- NA_integer_
  cds_end_idx <- NA_integer_
  downstream_anchor_idx <- NA_integer_
  downstream_frame_ref_idx <- NA_integer_
  right_bound <- nchar(seq)

  if (!is_missing_string(cds_anchor_seq)) {
    cds_hit <- find_fixed_once(seq, cds_anchor_seq)
    if (!is.na(cds_hit)) {
      cds_atg_idx <- cds_hit
      cds_end_idx <- cds_hit + nchar(cds_anchor_seq) - 1L
      right_bound <- cds_hit - 1L
    }
  } else if (!is_missing_string(downstream_anchor_seq)) {
    down_hit <- find_fixed_once(seq, downstream_anchor_seq)
    if (!is.na(down_hit)) {
      downstream_anchor_idx <- down_hit

      if (!is.na(downstream_frame_offset)) {
        downstream_frame_ref_idx <- down_hit + downstream_frame_offset
      }

      right_bound <- down_hit - 1L
    }
  }

  if (right_bound < left_bound) {
    return(list(
      utr_seq = NA_character_,
      utr_start_in_oligo = NA_integer_,
      utr_end_in_oligo = NA_integer_,
      cds_atg_idx = cds_atg_idx,
      cds_end_idx = cds_end_idx,
      downstream_anchor_idx = downstream_anchor_idx,
      downstream_frame_ref_idx = downstream_frame_ref_idx
    ))
  }

  list(
    utr_seq = substr(seq, left_bound, right_bound),
    utr_start_in_oligo = left_bound,
    utr_end_in_oligo = right_bound,
    cds_atg_idx = cds_atg_idx,
    cds_end_idx = cds_end_idx,
    downstream_anchor_idx = downstream_anchor_idx,
    downstream_frame_ref_idx = downstream_frame_ref_idx
  )
}

pick_gained_start <- function(utr_seq, consequence, has_endogenous_utr_atg, wt_utr_seq = NA_character_) {
  alt_atg <- find_all_codons(utr_seq, "ATG")

  gained_atg <- integer(0)

  if (!isTRUE(has_endogenous_utr_atg)) {
    gained_atg <- alt_atg
  } else {
    wt_atg <- find_all_codons(wt_utr_seq, "ATG")
    gained_atg <- setdiff(alt_atg, wt_atg)
  }

  if (identical(consequence, "START_Insertion") && length(gained_atg) == 0L) {
    gained_atg <- alt_atg
  }

  if (length(gained_atg) > 0L) {
    return(list(pos_local = min(gained_atg)))
  }

  list(pos_local = NA_integer_)
}

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

pick_single_reference_context <- function(target_name, oligo_name, cfg) {
  oligo_name <- oligo_name %||% ""

  use_seq2 <- !is_missing_string(cfg$cds_anchor_seq2 %||% NA_character_) ||
    !is_missing_string(cfg$wt_utr_seq2 %||% NA_character_) ||
    !is_missing_string(cfg$downstream_anchor_seq2 %||% NA_character_) ||
    !is.na(cfg$downstream_frame_offset2 %||% NA_integer_) ||
    !is_missing_string(cfg$canonical_kozak_15bp2 %||% NA_character_)

  use_seq1 <- !is_missing_string(cfg$cds_anchor_seq1 %||% NA_character_) ||
    !is_missing_string(cfg$wt_utr_seq1 %||% NA_character_) ||
    !is_missing_string(cfg$downstream_anchor_seq1 %||% NA_character_) ||
    !is.na(cfg$downstream_frame_offset1 %||% NA_integer_) ||
    !is_missing_string(cfg$canonical_kozak_15bp1 %||% NA_character_)

  # If BOTH seq1 and seq2 are genuinely populated for this target, a
  # single-column "sequence" input is inherently ambiguous - there's no way
  # to tell, from the config alone, which guide's reference context a given
  # row belongs to. This is not unique to any one gene: CTCF (wt_utr_seq1 !=
  # wt_utr_seq2, downstream_anchor_seq1 != downstream_anchor_seq2), DDX3X
  # (cds_anchor_seq1 != cds_anchor_seq2), and SLC2A1 (cds_anchor_seq1 !=
  # cds_anchor_seq2) all have this structure - only SON does not (only seq2
  # is populated there). Silently defaulting to seq2 (as the old code below
  # would do) risks using the wrong guide's reference for a given row without
  # any warning, so this fails loudly and generically instead.
  if (use_seq1 && use_seq2) {
    stop(
      "Target '", target_name, "' has two genuinely different reference ",
      "contexts (both seq1 and seq2 config values are populated), so ",
      "single-column 'sequence' input is ambiguous and not supported. ",
      "Provide merged 'sequence1'/'sequence2' columns instead."
    )
  }

  if (use_seq2) {
    return(list(
      ref_label = "sequence2",
      wt_utr_seq = cfg$wt_utr_seq2 %||% cfg$wt_utr_seq %||% NA_character_,
      cds_anchor_seq = cfg$cds_anchor_seq2 %||% cfg$cds_anchor_seq %||% NA_character_,
      downstream_anchor_seq = cfg$downstream_anchor_seq2 %||% cfg$downstream_anchor_seq %||% NA_character_,
      downstream_frame_offset = cfg$downstream_frame_offset2 %||% cfg$downstream_frame_offset %||% NA_integer_,
      canonical_kozak_15bp = cfg$canonical_kozak_15bp2 %||% cfg$canonical_kozak_15bp %||% NA_character_
    ))
  }

  if (use_seq1) {
    return(list(
      ref_label = "sequence1",
      wt_utr_seq = cfg$wt_utr_seq1 %||% cfg$wt_utr_seq %||% NA_character_,
      cds_anchor_seq = cfg$cds_anchor_seq1 %||% cfg$cds_anchor_seq %||% NA_character_,
      downstream_anchor_seq = cfg$downstream_anchor_seq1 %||% cfg$downstream_anchor_seq %||% NA_character_,
      downstream_frame_offset = cfg$downstream_frame_offset1 %||% cfg$downstream_frame_offset %||% NA_integer_,
      canonical_kozak_15bp = cfg$canonical_kozak_15bp1 %||% cfg$canonical_kozak_15bp %||% NA_character_
    ))
  }

  return(list(
    ref_label = "sequence",
    wt_utr_seq = cfg$wt_utr_seq %||% NA_character_,
    cds_anchor_seq = cfg$cds_anchor_seq %||% NA_character_,
    downstream_anchor_seq = cfg$downstream_anchor_seq %||% NA_character_,
    downstream_frame_offset = cfg$downstream_frame_offset %||% NA_integer_,
    canonical_kozak_15bp = cfg$canonical_kozak_15bp %||% NA_character_
  ))
}

annotate_one <- function(variant_key,
                         consequence,
                         oligo_name,
                         sequence,
                         cfg,
                         wt_utr_seq = NA_character_,
                         cds_anchor_seq = NA_character_,
                         downstream_anchor_seq = NA_character_,
                         downstream_frame_offset = NA_integer_,
                         canonical_kozak_15bp = NA_character_,
                         annotation_source = NA_character_) {
  is_utr_variant <- str_detect(variant_key %||% "", "UTR")

  blank <- tibble(
    is_utr_variant = is_utr_variant,
    is_start_insertion = FALSE,
    inserted_kozak_strength = NA_character_,
    gained_start_pos = NA_integer_,
    gained_start_frame = NA_integer_,
    gained_start_kozak15 = NA_character_,
    gained_start_kozak_score = NA_real_,
    gained_atg = FALSE,
    creates_startstop = NA,
    creates_uorf = NA,
    uorf_seq = NA_character_,
    creates_oorf = NA,
    oorf_seq = NA_character_,
    disrupts_endogenous_kozak = NA,
    endogenous_kozak_alt15 = NA_character_,
    endogenous_kozak_alt_score = NA_real_,
    annotation_source = annotation_source
  )

  if (!isTRUE(is_utr_variant)) return(blank)

  region <- extract_utr_region(
    seq = sequence,
    upstream_anchor_seq = cfg$upstream_anchor_seq %||% NA_character_,
    cds_anchor_seq = cds_anchor_seq,
    downstream_anchor_seq = downstream_anchor_seq,
    downstream_frame_offset = downstream_frame_offset
  )

  utr_seq <- region$utr_seq
  utr_offset <- region$utr_start_in_oligo
  cds_atg_idx <- region$cds_atg_idx
  cds_end_idx <- region$cds_end_idx
  downstream_frame_ref_idx <- region$downstream_frame_ref_idx

  is_start_insertion <- identical(consequence, "START_Insertion")
  inserted_kozak_strength <- if (is_start_insertion) classify_kozak_from_oligo(oligo_name) else NA_character_

  gained <- pick_gained_start(
    utr_seq = utr_seq,
    consequence = consequence,
    has_endogenous_utr_atg = cfg$has_endogenous_utr_atg,
    wt_utr_seq = wt_utr_seq
  )

  gained_start_pos <- if (!is.na(gained$pos_local) && !is.na(utr_offset)) {
    utr_offset + gained$pos_local - 1L
  } else {
    NA_integer_
  }

  gained_start_frame <- if (!is.na(cds_atg_idx) && !is.na(gained_start_pos)) {
    frame_relative_to_reference(gained_start_pos, cds_atg_idx)
  } else if (!is.na(downstream_frame_ref_idx) && !is.na(gained_start_pos)) {
    frame_relative_to_reference(gained_start_pos, downstream_frame_ref_idx)
  } else {
    NA_integer_
  }

  gained_start_kozak15 <- if (!is.na(gained_start_pos)) {
    get_kozak15(sequence, gained_start_pos)
  } else {
    NA_character_
  }

  gained_start_kozak_score <- score_kozak(gained_start_kozak15)

  creates_startstop <- if (!is.na(gained_start_pos)) {
    is_startstop_at(sequence, gained_start_pos)
  } else {
    NA
  }

  uorf_search_end <- if (!is.na(cds_atg_idx)) {
    cds_atg_idx - 1L
  } else {
    region$utr_end_in_oligo
  }

  upstream_stop_idx <- if (!is.na(gained_start_pos) && !is.na(uorf_search_end)) {
    find_first_inframe_stop(
      seq = sequence,
      start_idx = gained_start_pos,
      max_end = uorf_search_end
    )
  } else {
    NA_integer_
  }

  has_upstream_inframe_stop <- !is.na(upstream_stop_idx)

  has_valid_uorf <- if (!is.na(gained_start_pos) && !is.na(uorf_search_end)) {
    has_valid_orf(
      seq = sequence,
      start_idx = gained_start_pos,
      max_end = uorf_search_end
    )
  } else {
    NA
  }

  creates_uorf <- if (isTRUE(creates_startstop)) {
    FALSE
  } else {
    has_valid_uorf
  }

  uorf_seq <- if (isTRUE(creates_uorf) && !is.na(gained_start_pos) && !is.na(uorf_search_end)) {
    extract_orf_sequence(
      sequence,
      gained_start_pos,
      max_end = uorf_search_end
    )
  } else {
    NA_character_
  }

  creates_oorf <- if (!is.na(gained_start_pos) &&
                      !is.na(gained_start_frame) &&
                      (!is.na(cds_atg_idx) || !is.na(downstream_frame_ref_idx))) {
    if (isTRUE(creates_startstop)) {
      FALSE
    } else if (gained_start_frame == 1L) {
      FALSE
    } else {
      !has_upstream_inframe_stop
    }
  } else {
    NA
  }

  oorf_seq <- if (isTRUE(creates_oorf) &&
                  !is.na(gained_start_pos) &&
                  !is.na(cds_end_idx)) {
    extract_orf_sequence(sequence, gained_start_pos, max_end = cds_end_idx)
  } else {
    NA_character_
  }

  endogenous_kozak_alt15 <- if (!is.na(cds_atg_idx)) {
    get_kozak15(sequence, cds_atg_idx)
  } else {
    NA_character_
  }

  endogenous_kozak_alt_score <- score_kozak(endogenous_kozak_alt15)

  disrupts_endogenous_kozak <- if (!is.na(cds_atg_idx) &&
                                   !is_missing_string(canonical_kozak_15bp) &&
                                   !is.na(endogenous_kozak_alt15)) {
    endogenous_kozak_alt15 != canonical_kozak_15bp
  } else {
    NA
  }

  tibble(
    is_utr_variant = TRUE,
    is_start_insertion = is_start_insertion,
    inserted_kozak_strength = inserted_kozak_strength,
    gained_start_pos = gained_start_pos,
    gained_start_frame = gained_start_frame,
    gained_start_kozak15 = gained_start_kozak15,
    gained_start_kozak_score = gained_start_kozak_score,
    gained_atg = !is.na(gained_start_pos),
    creates_startstop = creates_startstop,
    creates_uorf = creates_uorf,
    uorf_seq = uorf_seq,
    creates_oorf = creates_oorf,
    oorf_seq = oorf_seq,
    disrupts_endogenous_kozak = disrupts_endogenous_kozak,
    endogenous_kozak_alt15 = endogenous_kozak_alt15,
    endogenous_kozak_alt_score = endogenous_kozak_alt_score,
    annotation_source = annotation_source
  )
}

annotate_single <- function(df, cfg, target_name) {
  out_list <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    ctx <- pick_single_reference_context(
      target_name = target_name,
      oligo_name = df$oligo_name[i],
      cfg = cfg
    )

    ann <- annotate_one(
      variant_key = df$variant_key[i],
      consequence = df$consequence[i],
      oligo_name = df$oligo_name[i],
      sequence = df$sequence[i],
      cfg = cfg,
      wt_utr_seq = ctx$wt_utr_seq,
      cds_anchor_seq = ctx$cds_anchor_seq,
      downstream_anchor_seq = ctx$downstream_anchor_seq,
      downstream_frame_offset = ctx$downstream_frame_offset,
      canonical_kozak_15bp = ctx$canonical_kozak_15bp,
      annotation_source = ctx$ref_label
    )
    out_list[[i]] <- ann
  }

  bind_cols(df, bind_rows(out_list))
}

annotate_merged_prefer_seq2 <- function(df, cfg) {
  out_list <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    seq2 <- df$sequence2[i] %||% NA_character_
    oligo2 <- df$oligo_name2[i] %||% NA_character_

    seq1 <- df$sequence1[i] %||% NA_character_
    oligo1 <- df$oligo_name1[i] %||% NA_character_

    use_seq2 <- !is_missing_string(seq2)

    if (use_seq2) {
      ann <- annotate_one(
        variant_key = df$variant_key[i],
        consequence = df$consequence[i],
        oligo_name = oligo2,
        sequence = seq2,
        cfg = cfg,
        wt_utr_seq = cfg$wt_utr_seq2 %||% cfg$wt_utr_seq %||% NA_character_,
        cds_anchor_seq = cfg$cds_anchor_seq2 %||% cfg$cds_anchor_seq %||% NA_character_,
        downstream_anchor_seq = cfg$downstream_anchor_seq2 %||% cfg$downstream_anchor_seq %||% NA_character_,
        downstream_frame_offset = cfg$downstream_frame_offset2 %||% cfg$downstream_frame_offset %||% NA_integer_,
        canonical_kozak_15bp = cfg$canonical_kozak_15bp2 %||% cfg$canonical_kozak_15bp %||% NA_character_,
        annotation_source = "sequence2"
      )
    } else {
      ann <- annotate_one(
        variant_key = df$variant_key[i],
        consequence = df$consequence[i],
        oligo_name = oligo1,
        sequence = seq1,
        cfg = cfg,
        wt_utr_seq = cfg$wt_utr_seq1 %||% cfg$wt_utr_seq %||% NA_character_,
        cds_anchor_seq = cfg$cds_anchor_seq1 %||% cfg$cds_anchor_seq %||% NA_character_,
        downstream_anchor_seq = cfg$downstream_anchor_seq1 %||% cfg$downstream_anchor_seq %||% NA_character_,
        downstream_frame_offset = cfg$downstream_frame_offset1 %||% cfg$downstream_frame_offset %||% NA_integer_,
        canonical_kozak_15bp = cfg$canonical_kozak_15bp1 %||% cfg$canonical_kozak_15bp %||% NA_character_,
        annotation_source = "sequence1"
      )
    }

    out_list[[i]] <- ann
  }

  bind_cols(df, bind_rows(out_list))
}

# ============================================================
# Main
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript annotate_sge_utr_simple.R <input.tsv> <output.tsv> <target_name>",
      "\nExample:",
      "Rscript annotate_sge_utr_simple.R input.tsv output.tsv SLC2A1"
    )
  )
}

input_file <- args[1]
output_file <- args[2]
target_name <- args[3]

if (!target_name %in% names(TARGET_CONFIGS)) {
  stop("target_name must be one of: ", paste(names(TARGET_CONFIGS), collapse = ", "))
}

cfg <- TARGET_CONFIGS[[target_name]]
df <- read_tsv(input_file, show_col_types = FALSE)

required_common <- c("variant_key", "consequence")
if (!all(required_common %in% names(df))) {
  stop("Input file must contain columns: ", paste(required_common, collapse = ", "))
}

is_single <- all(c("oligo_name", "sequence") %in% names(df))
is_merged <- all(c("oligo_name1", "sequence1", "oligo_name2", "sequence2") %in% names(df))

if (is_single) {
  out <- annotate_single(df, cfg, target_name)
} else if (is_merged) {
  out <- annotate_merged_prefer_seq2(df, cfg)
} else {
  stop(
    paste(
      "Input format not recognised.",
      "Need either:",
      "(variant_key, consequence, oligo_name, sequence)",
      "or",
      "(variant_key, consequence, oligo_name1, sequence1, oligo_name2, sequence2)"
    )
  )
}

out <- out %>%
  mutate(final_status = choose_final_status(.)) %>%
  select(-any_of("final_status"), everything(), final_status)

write_tsv(out, output_file)
message("Wrote: ", output_file)
