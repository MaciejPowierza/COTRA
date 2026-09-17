#!/usr/bin/env python3
"""
00_data_prep.py — Convert raw xlsx input to TSV for the R pipeline.

The xlsx file (FINAL_RESULTS_1732_4894.functional_annotation.FINAL.with_HEK293T_WT_ATAC.xlsx)
contains 6,208 loci × 67 columns. This script reads it and saves as TSV,
which the R pipeline (r_pipeline/main.R) expects.

Usage:
    python 00_data_prep.py
"""
import pandas as pd
import os
import sys

# ============ CONFIGURATION ============
# EDIT THIS PATH: point to your raw xlsx file
XLSX_PATH = "/path/to/FINAL_RESULTS_1732_4894.functional_annotation.FINAL.with_HEK293T_WT_ATAC.xlsx"

# EDIT THIS PATH: where to save the TSV (should be r_pipeline/data/ relative to this script)
# Default: saves to r_pipeline/data/ next to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TSV_DIR = os.path.join(SCRIPT_DIR, "r_pipeline", "data")
TSV_FILENAME = "FINAL_RESULTS_1732_4894_with_ATAC.tsv"
# ========================================

def main():
    tsv_path = os.path.join(TSV_DIR, TSV_FILENAME)

    print(f"Reading xlsx: {XLSX_PATH}", flush=True)
    if not os.path.exists(XLSX_PATH):
        print(f"ERROR: File not found: {XLSX_PATH}", file=sys.stderr)
        print("Please edit XLSX_PATH at the top of this script.", file=sys.stderr)
        sys.exit(1)

    df = pd.read_excel(XLSX_PATH)
    print(f"  Loaded {len(df)} rows × {len(df.columns)} columns", flush=True)

    # Create output directory if needed
    os.makedirs(TSV_DIR, exist_ok=True)

    # Save as TSV
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"  Saved TSV: {tsv_path}", flush=True)

    # Quick sanity check
    df2 = pd.read_csv(tsv_path, sep="\t", nrows=5)
    assert len(df2.columns) == len(df.columns), "Column count mismatch after round-trip!"
    print(f"  Sanity check passed: {len(df2.columns)} columns verified", flush=True)
    print("Done.", flush=True)

if __name__ == "__main__":
    main()
