"""
extract_vep_annotations.py
--------------------------
Extract SpliceAI (max delta score) and UTRAnnotator annotations from a pair of
VEP-annotated VCF files (two guides per screen) and merge into an existing
annotation TSV.

Usage:
    python extract_vep_annotations.py \
        --vcf1  path/to/gene_sg7_vep.vcf \
        --vcf2  path/to/gene_sg8_vep.vcf \
        --tsv   path/to/merged_annotated.tsv \
        --out   path/to/output_annotated.tsv

    Single-guide genes (e.g. SON): omit --vcf2 entirely.
        python extract_vep_annotations.py \
            --vcf1 path/to/SON_sgX_vep.vcf \
            --tsv  path/to/SON_extracted_annotated.tsv \
            --out  path/to/SON_output_annotated.tsv

    Genes with no SpliceAI data but with Pangolin scores instead (e.g. SON):
    pass a flat, already-split-VEP TSV (one row per variant, columns already
    split out rather than a packed VCF CSQ field) via --pangolin_tsv in place
    of --vcf1/--vcf2. UTRAnnotator columns are read from this same file.
        python extract_vep_annotations.py \
            --pangolin_tsv path/to/SON_sg10_splitvep_pangolin.tsv \
            --tsv          path/to/SON_extracted_annotated.tsv \
            --out          path/to/SON_output_annotated.tsv

Key decisions:
  - SpliceAI delta score: max of the 4 delta scores (DS_AG, DS_AL, DS_DG, DS_DL)
    across all transcripts in the CSQ field. This gives the single worst-case
    predicted splicing impact. Range 0–1; ≥0.2 is typically considered potentially
    pathogenic, ≥0.5 is high-confidence.
  - UTRAnnotator: consequence and annotation string are taken from the first
    transcript that has them populated (identical across protein-coding transcripts).
  - Discrepancy handling between sg7 and sg8 VCFs:
    * SpliceAI: take the MAX delta score across both VCFs (conservative / worst-case).
      In practice, scores are identical across guides (they depend only on the variant
      position, not the guide sequence), so this will never inflate values.
    * UTRAnnotator: consequence types are always identical. Annotation strings differ
      only in key ordering (e.g. "DistanceToCDS=...:type=..." vs "type=...:DistanceToCDS=...")
      — both are semantically identical. The script reports which VCF a discrepancy
      came from so you can audit if needed. The sg8 annotation is preferred when both
      are present; if only one guide has the annotation, that value is used.
  - Pangolin (--pangolin_tsv mode, used for genes like SON with no SpliceAI
    coverage): Pangolin reports a gain score (>=0) and a loss score (<=0)
    per variant, analogous to SpliceAI's four directional deltas. We report
    'pangolin_delta' = max(abs(gain_score), abs(loss_score)) as the single
    worst-case splicing-disruption magnitude, on the same 0-1 scale as
    spliceai_delta, so the two can be reasoned about the same way even
    though they live in separate output columns (spliceai_delta is left
    null for genes run in this mode, and pangolin_delta is left null for
    genes run with --vcf1/--vcf2 - the two scores are never both computed
    for the same gene, so there's no discrepancy-resolution step here).
    The signed gain_score/gain_pos/loss_score/loss_pos are also carried
    through unchanged for auditing.
"""

import re
import argparse
import pandas as pd
from pathlib import Path


# ── Helpers ────────────────────────────────────────────────────────────────────

def extract_variant_key(vcf_id: str) -> str:
    """Strip transcript + guide prefix from VCF ID to get the bare variant key.
    
    E.g. 'ENST00000644876.2_sg7_CDS_114_C_SNV'  →  'CDS_114_C_SNV'
    """
    match = re.search(r'_sg\d+_(.+)$', vcf_id)
    return match.group(1) if match else vcf_id


def parse_csq_header(vcf_path: str) -> list[str]:
    """Return the ordered list of CSQ sub-field names from the VCF header."""
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                fmt = re.search(r'Format: ([^"]+)"', line)
                if fmt:
                    return fmt.group(1).strip().split('|')
    raise ValueError(f"No CSQ FORMAT definition found in {vcf_path}")


def parse_vcf_annotations(vcf_path: str) -> dict:
    """
    Parse a VEP-annotated VCF and return per-variant annotations.

    Returns
    -------
    dict keyed by variant_key:
        {
            'spliceai_delta': float | None,   # max(DS_AG, DS_AL, DS_DG, DS_DL)
            'utr_consequence': str,            # e.g. '5_prime_UTR_premature_start_codon_gain_variant'
            'utr_annotation':  str,            # e.g. 'type=uORF:KozakStrength=Strong:...'
        }
    """
    csq_fields = parse_csq_header(vcf_path)

    # SpliceAI delta score columns
    spliceai_ds_cols = ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL',
                        'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']
    ds_indices = {col: csq_fields.index(col)
                  for col in spliceai_ds_cols if col in csq_fields}

    # UTRAnnotator columns
    utr_cons_idx = csq_fields.index('5UTR_consequence')   if '5UTR_consequence'  in csq_fields else None
    utr_ann_idx  = csq_fields.index('5UTR_annotation')    if '5UTR_annotation'   in csq_fields else None

    if not ds_indices:
        print(f"  WARNING: No SpliceAI DS columns found in {Path(vcf_path).name}")
    if utr_cons_idx is None:
        print(f"  WARNING: No UTRAnnotator columns found in {Path(vcf_path).name}")

    results = {}
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue

            vkey = extract_variant_key(cols[2])   # ID column
            info = cols[7]

            # Find the CSQ field in INFO
            csq_raw = next((x[4:] for x in info.split(';') if x.startswith('CSQ=')), None)
            if not csq_raw:
                continue

            max_delta     = None
            utr_consequence = ''
            utr_annotation  = ''

            for transcript in csq_raw.split(','):
                vals = transcript.split('|')

                # SpliceAI: max delta score across all 4 splice directions
                scores = []
                for col, idx in ds_indices.items():
                    if idx < len(vals) and vals[idx] not in ('', '.'):
                        try:
                            scores.append(float(vals[idx]))
                        except ValueError:
                            pass
                if scores:
                    tx_max = max(scores)
                    if max_delta is None or tx_max > max_delta:
                        max_delta = tx_max

                # UTRAnnotator: take the first non-empty value
                # (values are identical across transcripts; key order may vary)
                if utr_cons_idx is not None and not utr_consequence:
                    v = vals[utr_cons_idx] if utr_cons_idx < len(vals) else ''
                    if v:
                        utr_consequence = v
                if utr_ann_idx is not None and not utr_annotation:
                    v = vals[utr_ann_idx] if utr_ann_idx < len(vals) else ''
                    if v:
                        utr_annotation = v

            results[vkey] = {
                'spliceai_delta': max_delta,
                'utr_consequence': utr_consequence,
                'utr_annotation':  utr_annotation,
            }

    return results


def parse_pangolin_tsv(tsv_path: str) -> dict:
    """
    Parse a flat, already-split-VEP TSV (one row per variant/transcript, with
    CSQ sub-fields already split into their own columns - e.g. produced by
    bcftools +split-vep) and return per-variant Pangolin + UTRAnnotator
    annotations. Used for genes that have Pangolin scores instead of SpliceAI
    (e.g. SON).

    Returns
    -------
    dict keyed by variant_key:
        {
            'pangolin_delta':      float | None,  # max(abs(gain), abs(loss))
            'pangolin_gain_score': float | None,
            'pangolin_gain_pos':   float | None,
            'pangolin_loss_score': float | None,
            'pangolin_loss_pos':   float | None,
            'utr_consequence':     str,
            'utr_annotation':      str,
        }
    """
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)
    # Headers look like '#[1]CHROM', '[64]pangolin_gain_score', etc. - strip
    # the leading '#' and the bracketed column-number prefix.
    df.columns = [re.sub(r'^#?\[\d+\]', '', c) for c in df.columns]

    required = ['ID']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{tsv_path} is missing expected column(s): {missing}")

    pangolin_cols = ['pangolin_gain_score', 'pangolin_gain_pos',
                     'pangolin_loss_score', 'pangolin_loss_pos']
    for col in pangolin_cols:
        if col not in df.columns:
            print(f"  WARNING: column '{col}' not found in {Path(tsv_path).name}")

    utr_cons_col = '_5UTR_consequence' if '_5UTR_consequence' in df.columns else None
    utr_ann_col  = '_5UTR_annotation'  if '_5UTR_annotation'  in df.columns else None
    if utr_cons_col is None:
        print(f"  WARNING: No UTRAnnotator columns found in {Path(tsv_path).name}")

    def to_float(v):
        if v is None or (isinstance(v, float) and pd.isna(v)) or v in ('', '.'):
            return None
        try:
            return float(v)
        except ValueError:
            return None

    df['variant_key'] = df['ID'].apply(extract_variant_key)

    results = {}
    for vkey, group in df.groupby('variant_key'):
        gain_score = loss_score = gain_pos = loss_pos = None
        utr_consequence = ''
        utr_annotation  = ''

        for _, row in group.iterrows():
            gs = to_float(row.get('pangolin_gain_score')) if 'pangolin_gain_score' in df.columns else None
            ls = to_float(row.get('pangolin_loss_score')) if 'pangolin_loss_score' in df.columns else None
            gp = to_float(row.get('pangolin_gain_pos'))   if 'pangolin_gain_pos'   in df.columns else None
            lp = to_float(row.get('pangolin_loss_pos'))   if 'pangolin_loss_pos'   in df.columns else None

            # Keep the values from whichever row has the largest-magnitude
            # score seen so far (mirrors the "worst-case across transcripts"
            # logic used for SpliceAI).
            if gs is not None and (gain_score is None or abs(gs) > abs(gain_score)):
                gain_score, gain_pos = gs, gp
            if ls is not None and (loss_score is None or abs(ls) > abs(loss_score)):
                loss_score, loss_pos = ls, lp

            if utr_cons_col and not utr_consequence:
                v = row.get(utr_cons_col, '')
                if v and v != '.' and not pd.isna(v):
                    utr_consequence = v
            if utr_ann_col and not utr_annotation:
                v = row.get(utr_ann_col, '')
                if v and v != '.' and not pd.isna(v):
                    utr_annotation = v

        magnitudes = [abs(s) for s in (gain_score, loss_score) if s is not None]
        pangolin_delta = max(magnitudes) if magnitudes else None

        results[vkey] = {
            'pangolin_delta':      pangolin_delta,
            'pangolin_gain_score': gain_score,
            'pangolin_gain_pos':   gain_pos,
            'pangolin_loss_score': loss_score,
            'pangolin_loss_pos':   loss_pos,
            'utr_consequence':     utr_consequence,
            'utr_annotation':      utr_annotation,
        }

    return results


def merge_annotations(ann1: dict, ann2: dict) -> dict:
    """
    Merge annotations from one or two VCFs (one or two guides for the same
    variants). For single-guide genes, pass ann2={} - every field then
    resolves to ann1's values and no discrepancies are raised.

    Strategy:
      - spliceai_delta: take the MAX of the two (conservative; scores are
        position-based so will always agree — max is a safe fallback).
      - utr_consequence: should always agree; use whichever is non-empty.
      - utr_annotation: key ordering may differ between guides (same data);
        prefer ann2 (sg8) when both present.

    Returns a merged dict keyed by all variant_keys from either VCF.
    Also returns a list of discrepancy notes for auditing.
    """
    all_keys = set(ann1) | set(ann2)
    merged   = {}
    discrepancies = []

    for k in all_keys:
        a1 = ann1.get(k, {'spliceai_delta': None, 'utr_consequence': '', 'utr_annotation': ''})
        a2 = ann2.get(k, {'spliceai_delta': None, 'utr_consequence': '', 'utr_annotation': ''})

        # SpliceAI
        d1, d2 = a1['spliceai_delta'], a2['spliceai_delta']
        if d1 is not None and d2 is not None:
            delta = max(d1, d2)
        elif d1 is not None:
            delta = d1
        else:
            delta = d2   # may be None

        # UTR consequence
        uc1, uc2 = a1['utr_consequence'], a2['utr_consequence']
        if uc1 and uc2 and uc1 != uc2:
            discrepancies.append(
                f"  DISCREPANCY utr_consequence [{k}]: vcf1={uc1!r}  vcf2={uc2!r}  → using vcf2"
            )
        utr_consequence = uc2 or uc1   # prefer sg8

        # UTR annotation (key-order variation is OK)
        ua1, ua2 = a1['utr_annotation'], a2['utr_annotation']
        utr_annotation = ua2 or ua1   # prefer sg8

        merged[k] = {
            'spliceai_delta': delta,
            'utr_consequence': utr_consequence,
            'utr_annotation':  utr_annotation,
        }

    return merged, discrepancies


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf1', required=False, default=None,
                        help='First VEP VCF (e.g. sg7). Not used with --pangolin_tsv.')
    parser.add_argument('--vcf2', required=False, default=None,
                        help='Second VEP VCF (e.g. sg8). Omit for single-guide genes (e.g. SON). '
                             'Not used with --pangolin_tsv.')
    parser.add_argument('--pangolin_tsv', required=False, default=None,
                        help='Flat split-VEP TSV with Pangolin + UTRAnnotator columns, for genes '
                             'with no SpliceAI coverage (e.g. SON). Mutually exclusive with '
                             '--vcf1/--vcf2.')
    parser.add_argument('--tsv',  required=True, help='Existing annotation TSV to append to')
    parser.add_argument('--out',  required=True, help='Output TSV path')
    parser.add_argument('--key_col', default='variant_key',
                        help='Column in TSV that holds the variant identifier (default: variant_key)')
    args = parser.parse_args()

    if args.pangolin_tsv and (args.vcf1 or args.vcf2):
        parser.error("--pangolin_tsv cannot be combined with --vcf1/--vcf2.")
    if not args.pangolin_tsv and not args.vcf1:
        parser.error("Must provide either --vcf1 (optionally with --vcf2) or --pangolin_tsv.")

    pangolin_mode = args.pangolin_tsv is not None

    if pangolin_mode:
        print(f"Parsing {Path(args.pangolin_tsv).name} (Pangolin mode) …")
        merged = parse_pangolin_tsv(args.pangolin_tsv)
        print(f"  → {len(merged)} variants parsed")
    else:
        print(f"Parsing {Path(args.vcf1).name} …")
        ann1 = parse_vcf_annotations(args.vcf1)
        print(f"  → {len(ann1)} variants parsed")

        single_guide = args.vcf2 is None
        if single_guide:
            print("No --vcf2 given: running in single-guide mode (e.g. SON).")
            ann2 = {}
        else:
            print(f"Parsing {Path(args.vcf2).name} …")
            ann2 = parse_vcf_annotations(args.vcf2)
            print(f"  → {len(ann2)} variants parsed")

        print("Merging annotations …" if not single_guide else "Building annotation table …")
        merged, discrepancies = merge_annotations(ann1, ann2)

        if single_guide:
            # No second guide to disagree with, so discrepancy checking is moot.
            pass
        elif discrepancies:
            print(f"\n{len(discrepancies)} UTR consequence discrepancy/ies between VCFs (using vcf2):")
            for d in discrepancies:
                print(d)
        else:
            print("  No UTR consequence discrepancies between the two guides ✓")

    print(f"\nLoading TSV: {Path(args.tsv).name} …")
    tsv = pd.read_csv(args.tsv, sep='\t')
    print(f"  → {len(tsv)} rows, {len(tsv.columns)} columns")

    # Build lookup dataframe from merged annotations
    ann_df = pd.DataFrame.from_dict(merged, orient='index')
    ann_df.index.name = args.key_col
    ann_df = ann_df.reset_index()

    # Merge
    out = tsv.merge(ann_df, on=args.key_col, how='left')

    # Coverage report
    n_utr = (out['utr_consequence'] != '').sum()
    print(f"\nCoverage in output:")
    if pangolin_mode:
        n_pangolin = out['pangolin_delta'].notna().sum()
        print(f"  pangolin_delta non-null : {n_pangolin}/{len(out)} ({100*n_pangolin/len(out):.1f}%)")
    else:
        n_spliceai = out['spliceai_delta'].notna().sum()
        print(f"  spliceai_delta non-null : {n_spliceai}/{len(out)} ({100*n_spliceai/len(out):.1f}%)")
    print(f"  utr_consequence present : {n_utr}/{len(out)} ({100*n_utr/len(out):.1f}%)")

    out.to_csv(args.out, sep='\t', index=False)
    print(f"\nSaved to {args.out}")


if __name__ == '__main__':
    main()
