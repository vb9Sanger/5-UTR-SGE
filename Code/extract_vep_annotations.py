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
    parser.add_argument('--vcf1', required=True, help='First VEP VCF (e.g. sg7)')
    parser.add_argument('--vcf2', required=False, default=None,
                        help='Second VEP VCF (e.g. sg8). Omit for single-guide genes (e.g. SON).')
    parser.add_argument('--tsv',  required=True, help='Existing annotation TSV to append to')
    parser.add_argument('--out',  required=True, help='Output TSV path')
    parser.add_argument('--key_col', default='variant_key',
                        help='Column in TSV that holds the variant identifier (default: variant_key)')
    args = parser.parse_args()

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
    n_spliceai = out['spliceai_delta'].notna().sum()
    n_utr      = (out['utr_consequence'] != '').sum()
    print(f"\nCoverage in output:")
    print(f"  spliceai_delta non-null : {n_spliceai}/{len(out)} ({100*n_spliceai/len(out):.1f}%)")
    print(f"  utr_consequence present : {n_utr}/{len(out)} ({100*n_utr/len(out):.1f}%)")

    out.to_csv(args.out, sep='\t', index=False)
    print(f"\nSaved to {args.out}")


if __name__ == '__main__':
    main()
