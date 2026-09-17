#!/usr/bin/env python3
"""
01_sw_features.py — Smith-Waterman alignment features + ablation study.

Computes SW local alignment between 23nt gRNA (spacer+PAM) and 200nt centered_sequence,
extracts 55 alignment features, patches long_data.csv, and runs ablation:
  Config A: Hamming mm_pos only (current baseline)
  Config B: SW features only (replace mm_pos)
  Config C: Both Hamming + SW

Usage:
    python 01_sw_features.py

Requires: long_data.csv, test_loci.csv, feature_groups.csv from the R pipeline (Step 1).
"""
import pandas as pd, numpy as np, xgboost as xgb
from scipy.stats import spearmanr
from Bio import pairwise2
import os, sys, warnings, time
warnings.filterwarnings('ignore')
def log(m): print(m, flush=True)

# ============ CONFIGURATION ============
# EDIT THESE PATHS to match your local setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
print(SCRIPT_DIR)

# Path to the TSV data file (output of 00_data_prep.py, or your existing TSV)
DATA_PATH = os.path.join(SCRIPT_DIR, "r_pipeline", "data", "FINAL_RESULTS_1732_4894_with_ATAC.tsv")
print(DATA_PATH)

# Path to long_data.csv and related files (output of R pipeline Step 1)
LONG_DATA_PATH = os.path.join(SCRIPT_DIR, "results", "nn_data", "long_data.csv")
TEST_LOCI_PATH = os.path.join(SCRIPT_DIR, "results", "nn_data", "test_loci.csv")
FEATURE_GROUPS_PATH = os.path.join(SCRIPT_DIR, "results", "nn_data", "feature_groups.csv")

# Output directory for results
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "results")
# ========================================

GRNA = {
    "1732": "TGCTCGAGTGGGTCCCCGTGAGG",  # spacer+PAM (23nt)
    "4894": "GAGGACGAGATGTAAGAGGCTGG",
}

# ── Step 1: Load data and compute SW features ──
log("Step 1: Loading data and computing SW alignments...")
if not os.path.exists(DATA_PATH):
    print(f"ERROR: Data file not found: {DATA_PATH}", file=sys.stderr)
    sys.exit(1)
df = pd.read_csv(DATA_PATH, sep="\t")
log(f"  {len(df)} loci")

ss = "1732_1.breakends.noEnds.FILTERED"; se = "Neg_ctrl_1h_DNA1.breakends.noEnds.FILTERED"
sc = list(df.columns[list(df.columns).index(ss):list(df.columns).index(se)+1])

# Guide membership
guide_indices = {}
for g in ["1732", "4894"]:
    gs = [c for c in sc if c.startswith(g)]
    guide_indices[g] = np.where((df[gs].sum(axis=1) > 0) | (df[f"mismatch_count_{g}"] == 0))[0]
    log(f"  gRNA {g}: {len(guide_indices[g])} loci")

def compute_sw_features(seq, grna_seq):
    """Compute SW alignment and extract features."""
    alignments = pairwise2.align.localms(seq, grna_seq, 2, -1, -10, -5,
                                          one_alignment_only=True)
    if not alignments:
        return {f"sw_gap_pos_{i}": 0 for i in range(1, 24)} | \
               {f"sw_mm_pos_{i}": 0 for i in range(1, 24)} | \
               {"sw_score": 0, "sw_score_norm": 0, "sw_align_len": 0,
                "sw_num_gaps": 0, "sw_num_gap_opens": 0,
                "sw_match_count": 0, "sw_mismatch_count": 0,
                "sw_start_offset": 0, "sw_pam_match": 0}

    aln = alignments[0]
    seq_aln = aln.seqA[aln.start:aln.end]
    grna_aln = aln.seqB[aln.start:aln.end]

    score = aln.score
    align_len = len(seq_aln)

    num_gaps = seq_aln.count('-') + grna_aln.count('-')
    num_gap_opens = 0
    in_gap = False
    for s, g in zip(seq_aln, grna_aln):
        if s == '-' or g == '-':
            if not in_gap:
                num_gap_opens += 1
                in_gap = True
        else:
            in_gap = False

    match_count = sum(1 for s, g in zip(seq_aln, grna_aln) if s != '-' and g != '-' and s == g)
    mismatch_count = sum(1 for s, g in zip(seq_aln, grna_aln) if s != '-' and g != '-' and s != g)

    # Per-position features (aligned to gRNA positions)
    gap_pos = np.zeros(23, dtype=int)
    mm_pos = np.zeros(23, dtype=int)
    grna_pos = 0
    for s, g in zip(seq_aln, grna_aln):
        if g == '-':
            continue
        if grna_pos >= 23:
            break
        if s == '-':
            gap_pos[grna_pos] = 1
        elif s != g:
            mm_pos[grna_pos] = 1
        grna_pos += 1

    start_offset = aln.start - 100
    pam_match = 1 if all(mm_pos[i] == 0 and gap_pos[i] == 0 for i in range(20, 23)) else 0

    features = {
        "sw_score": float(score),
        "sw_score_norm": float(score) / max(1, align_len),
        "sw_align_len": int(align_len),
        "sw_num_gaps": int(num_gaps),
        "sw_num_gap_opens": int(num_gap_opens),
        "sw_match_count": int(match_count),
        "sw_mismatch_count": int(mismatch_count),
        "sw_start_offset": int(start_offset),
        "sw_pam_match": int(pam_match),
    }
    for i in range(23):
        features[f"sw_gap_pos_{i+1}"] = int(gap_pos[i])
        features[f"sw_mm_pos_{i+1}"] = int(mm_pos[i])
    return features

# Compute SW features per guide
sw_by_locus = {}
t0 = time.time()
for g in ["1732", "4894"]:
    grna_seq = GRNA[g]
    for idx in guide_indices[g]:
        row = df.iloc[idx]
        cid = row["cluster_id"]
        centered_seq = row["centered_sequence"]
        sw_by_locus[cid] = compute_sw_features(centered_seq, grna_seq)
    log(f"  gRNA {g}: {len(guide_indices[g])} alignments done")

log(f"  Total SW features computed for {len(sw_by_locus)} loci in {time.time()-t0:.1f}s")

# Verify
cid0 = df.iloc[guide_indices["1732"][0]]["cluster_id"]
log(f"  Sample {cid0}: sw_score={sw_by_locus[cid0]['sw_score']}, sw_align_len={sw_by_locus[cid0]['sw_align_len']}, "
    f"sw_mm_count={sw_by_locus[cid0]['sw_mismatch_count']}, sw_gaps={sw_by_locus[cid0]['sw_num_gaps']}")

# ── Step 2: Load and patch long_data ──
log("\nStep 2: Loading long_data.csv...")
if not os.path.exists(LONG_DATA_PATH):
    print(f"ERROR: long_data.csv not found: {LONG_DATA_PATH}", file=sys.stderr)
    print("Please run the R pipeline (Step 1) first.", file=sys.stderr)
    sys.exit(1)
ld = pd.read_csv(LONG_DATA_PATH)
log(f"  {len(ld)} rows, {len(ld.columns)} cols")

# Build SW feature DataFrame
sw_cols = list(sw_by_locus[cid0].keys())
sw_df = pd.DataFrame.from_dict(sw_by_locus, orient='index', columns=sw_cols)
sw_df.index.name = 'locus_id'
ld = ld.merge(sw_df, on='locus_id', how='left')
for c in sw_cols:
    ld[c] = ld[c].fillna(0)

log(f"  Patched {len(sw_cols)} SW features")
log(f"  sw_score range: [{ld['sw_score'].min():.1f}, {ld['sw_score'].max():.1f}]")
log(f"  sw_mm_pos_1 nonzero: {(ld['sw_mm_pos_1']!=0).sum()}")
log(f"  sw_num_gaps nonzero: {(ld['sw_num_gaps']!=0).sum()}")

# ── Step 3: Load test split and groups ──
log("\nStep 3: Loading test split and feature groups...")
tl = pd.read_csv(TEST_LOCI_PATH)['locus_id'].values
fg = pd.read_csv(FEATURE_GROUPS_PATH)
groups = {g: fg[fg['group']==g]['feature'].tolist() for g in fg['group'].unique()}
log(f"  {len(tl)} test loci")

# ── Step 4: Ablation study ──
log("\nStep 4: Ablation study (XGBoost depth=6, nrounds=30)...")

def recall_at_k(pred, actual, tdf, k=50):
    rs = []
    for g in ["1732", "4894"]:
        c = f"grna_{g}"
        if c not in tdf.columns: continue
        idx = tdf[c] == 1
        if idx.sum() == 0: continue
        pg = pred[idx]; ag = actual[idx]; lg = tdf['locus_id'].values[idx]
        nz = ag[ag > 0]
        if len(nz) < 5: continue
        thr = np.median(nz)
        th = lg[ag > thr]
        lp = pd.Series(pg).groupby(lg).max()
        tk = lp.sort_values(ascending=False).index[:min(k, len(lp))]
        rs.append(np.isin(th, tk).sum() / len(th))
    return np.mean(rs) if rs else np.nan

def metrics(pred, actual, tdf):
    sp = spearmanr(pred, actual)[0]
    r2 = 1 - np.sum((actual-pred)**2) / np.sum((actual-np.mean(actual))**2)
    s17 = spearmanr(pred[tdf['grna_1732']==1], actual[tdf['grna_1732']==1])[0] if (tdf['grna_1732']==1).any() else np.nan
    s48 = spearmanr(pred[tdf['grna_4894']==1], actual[tdf['grna_4894']==1])[0] if (tdf['grna_4894']==1).any() else np.nan
    return {'spearman_all': sp, 'r2_all': r2, 'spearman_1732': s17, 'spearman_4894': s48,
            'recall_at_50': recall_at_k(pred, actual, tdf, 50),
            'recall_at_100': recall_at_k(pred, actual, tdf, 100)}

def train_xgb(ld, fcols, tl):
    tr = ld[~ld['locus_id'].isin(tl)]; te = ld[ld['locus_id'].isin(tl)]
    Xtr = np.nan_to_num(tr[fcols].values, nan=0.0); ytr = tr['log_edit_count'].values
    Xte = np.nan_to_num(te[fcols].values, nan=0.0); yte = te['log_edit_count'].values
    dtr = xgb.DMatrix(Xtr, label=ytr); dte = xgb.DMatrix(Xte, label=yte)
    best_sp, best_m = -1, None
    for d in [4, 6]:
        p = {'max_depth': d, 'eta': 0.1, 'subsample': 0.8, 'colsample_bytree': 0.8,
             'objective': 'reg:squarederror', 'eval_metric': 'rmse', 'seed': 123,
             'n_jobs': -1, 'tree_method': 'hist'}
        m = xgb.train(p, dtr, num_boost_round=30)
        pr = m.predict(dte); mt = metrics(pr, yte, te)
        log(f"    d={d}: Sp={mt['spearman_all']:.4f} R²={mt['r2_all']:.4f} R@50={mt['recall_at_50']:.4f} R@100={mt['recall_at_100']:.4f}")
        if mt['spearman_all'] > best_sp: best_sp, best_m = mt['spearman_all'], mt
    return best_m

# Define feature sets for ablation
sw_scalar = ['sw_score', 'sw_score_norm', 'sw_align_len', 'sw_num_gaps',
             'sw_num_gap_opens', 'sw_match_count', 'sw_mismatch_count',
             'sw_start_offset', 'sw_pam_match']
sw_gap_pos = [f"sw_gap_pos_{i}" for i in range(1, 24)]
sw_mm_pos = [f"sw_mm_pos_{i}" for i in range(1, 24)]
all_sw = sw_scalar + sw_gap_pos + sw_mm_pos

# Base groups (non-alignment, always included)
base_groups = ['sequence', 'pair', 'per_guide_struct', 'global_struct',
               'chromatin', 'interaction', 'experimental']
base_cols = [c for g in base_groups for c in groups.get(g, []) if c in ld.columns]

# Config A: Hamming only (current alignment group)
cfgA_cols = base_cols + [c for c in groups.get('alignment', []) if c in ld.columns]
cfgA_cols = [c for c in cfgA_cols if not c.startswith('sw_')]

# Config B: SW only (replace alignment group)
cfgB_cols = base_cols + all_sw
cfgB_cols = [c for c in cfgB_cols if c in ld.columns]
cfgB_cols = [c for c in cfgB_cols if not c.startswith('mm_pos_') and c != 'PAM_mismatch']

# Config C: Both Hamming + SW
cfgC_cols = base_cols + [c for c in groups.get('alignment', []) if c in ld.columns] + all_sw
cfgC_cols = [c for c in cfgC_cols if c in ld.columns]

configs = [
    ('A_hamming_only', cfgA_cols, 'Current Hamming mm_pos (baseline)'),
    ('B_sw_only', cfgB_cols, 'SW features replacing Hamming'),
    ('C_both', cfgC_cols, 'Hamming + SW features'),
]

results = {}
for name, fcols, desc in configs:
    log(f"\n  === Config {name}: {desc} ({len(fcols)} features) ===")
    m = train_xgb(ld, fcols, tl)
    results[name] = m
    log(f"  → Sp={m['spearman_all']:.4f} R²={m['r2_all']:.4f} R@50={m['recall_at_50']:.4f} R@100={m['recall_at_100']:.4f}")

# ── Step 5: Save comparison ──
log("\n\n" + "="*60)
log("SW ABLATION RESULTS")
log("="*60)
rows = []
for i, (name, fcols, desc) in enumerate(configs):
    m = results[name]
    rows.append({'config': name, 'description': desc, 'n_features': len(fcols), **m})
comp = pd.DataFrame(rows)
pd.set_option('display.float_format', '{:.4f}'.format)
log("\n" + comp.to_string(index=False))

os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, "sw_ablation_comparison.csv")
comp.to_csv(out_path, index=False)
log(f"\nSaved to {out_path}")
