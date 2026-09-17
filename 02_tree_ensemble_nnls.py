#!/usr/bin/env python3
"""
02_tree_ensemble_nnls.py — Tree-based stacking ensemble (XGBoost + CatBoost + LightGBM)
with NNLS (non-negative least squares) meta-learner.

Key changes from previous ridge-based stacking:
  - Meta-learner: scipy.optimize.nnls instead of Ridge (weights constrained ≥ 0)
  - CatBoost: depth=6, iterations=200, early_stopping_rounds=20 (was depth=4, iterations=50)
  - LightGBM: num_leaves=63, n_estimators=100 (was num_leaves=31, n_estimators=50)
  - Also includes ridge stacking results for direct comparison

Uses SW feature set (best from Q3 ablation).
Single train/meta split (80/20) for meta-learner training.

Usage:
    python 02_tree_ensemble_nnls.py

Requires: long_data.csv, test_loci.csv, feature_groups.csv from the R pipeline (Step 1).
"""
import pandas as pd, numpy as np, xgboost as xgb, lightgbm as lgb, catboost as cb
from scipy.stats import spearmanr
from scipy.optimize import nnls
from sklearn.linear_model import Ridge
from Bio import pairwise2
import os, sys, warnings, time
warnings.filterwarnings('ignore')
def log(m): print(m, flush=True)

# ============ CONFIGURATION ============
# EDIT THESE PATHS to match your local setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Path to the TSV data file (for SW alignment computation)
DATA_PATH = os.path.join(SCRIPT_DIR, "r_pipeline", "data", "FINAL_RESULTS_1732_4894_with_ATAC.tsv")

# Path to long_data.csv and related files (output of R pipeline Step 1)
LONG_DATA_PATH = os.path.join(SCRIPT_DIR, "results", "nn_data", "long_data.csv")
TEST_LOCI_PATH = os.path.join(SCRIPT_DIR, "results", "nn_data", "test_loci.csv")
FEATURE_GROUPS_PATH = os.path.join(SCRIPT_DIR, "results", "nn_data", "feature_groups.csv")

# Output directory for results
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "results")
# ========================================

GRNA = {"1732": "TGCTCGAGTGGGTCCCCGTGAGG", "4894": "GAGGACGAGATGTAAGAGGCTGG"}

# ── Step 1: Compute SW features ──
log("Step 1: Computing SW features...")
if not os.path.exists(DATA_PATH):
    print(f"ERROR: Data file not found: {DATA_PATH}", file=sys.stderr)
    sys.exit(1)
df = pd.read_csv(DATA_PATH, sep="\t")
ss = "1732_1.breakends.noEnds.FILTERED"; se = "Neg_ctrl_1h_DNA1.breakends.noEnds.FILTERED"
sc = list(df.columns[list(df.columns).index(ss):list(df.columns).index(se)+1])
guide_indices = {}
for g in ["1732", "4894"]:
    gs = [c for c in sc if c.startswith(g)]
    guide_indices[g] = np.where((df[gs].sum(axis=1) > 0) | (df[f"mismatch_count_{g}"] == 0))[0]

def compute_sw(seq, grna_seq):
    alns = pairwise2.align.localms(seq, grna_seq, 2, -1, -10, -5, one_alignment_only=True)
    if not alns:
        return {"sw_score": 0, "sw_score_norm": 0, "sw_align_len": 0, "sw_num_gaps": 0,
                "sw_num_gap_opens": 0, "sw_match_count": 0, "sw_mismatch_count": 0,
                "sw_start_offset": 0, "sw_pam_match": 0,
                **{f"sw_gap_pos_{i}": 0 for i in range(1, 24)},
                **{f"sw_mm_pos_{i}": 0 for i in range(1, 24)}}
    a = alns[0]; sa = a.seqA[a.start:a.end]; ga = a.seqB[a.start:a.end]
    score = a.score; alen = len(sa)
    ngaps = sa.count('-') + ga.count('-')
    ngopens = 0; ig = False
    for s, g in zip(sa, ga):
        if s == '-' or g == '-':
            if not ig: ngopens += 1; ig = True
        else: ig = False
    mc = sum(1 for s, g in zip(sa, ga) if s != '-' and g != '-' and s == g)
    mmc = sum(1 for s, g in zip(sa, ga) if s != '-' and g != '-' and s != g)
    gp = np.zeros(23, dtype=int); mp = np.zeros(23, dtype=int); gp_pos = 0
    for s, g in zip(sa, ga):
        if g == '-': continue
        if gp_pos >= 23: break
        if s == '-': gp[gp_pos] = 1
        elif s != g: mp[gp_pos] = 1
        gp_pos += 1
    pam_m = 1 if all(mp[i] == 0 and gp[i] == 0 for i in range(20, 23)) else 0
    feats = {"sw_score": float(score), "sw_score_norm": float(score)/max(1,alen),
             "sw_align_len": int(alen), "sw_num_gaps": int(ngaps), "sw_num_gap_opens": int(ngopens),
             "sw_match_count": int(mc), "sw_mismatch_count": int(mmc),
             "sw_start_offset": int(a.start - 100), "sw_pam_match": int(pam_m)}
    for i in range(23):
        feats[f"sw_gap_pos_{i+1}"] = int(gp[i]); feats[f"sw_mm_pos_{i+1}"] = int(mp[i])
    return feats

sw_by_locus = {}
t0 = time.time()
for g in ["1732", "4894"]:
    for idx in guide_indices[g]:
        sw_by_locus[df.iloc[idx]["cluster_id"]] = compute_sw(df.iloc[idx]["centered_sequence"], GRNA[g])
log(f"  SW features for {len(sw_by_locus)} loci in {time.time()-t0:.1f}s")

# ── Step 2: Load and patch long_data ──
log("\nStep 2: Loading and patching long_data...")
if not os.path.exists(LONG_DATA_PATH):
    print(f"ERROR: long_data.csv not found: {LONG_DATA_PATH}", file=sys.stderr)
    print("Please run the R pipeline (Step 1) first.", file=sys.stderr)
    sys.exit(1)
ld = pd.read_csv(LONG_DATA_PATH)
sw_cols = list(sw_by_locus[list(sw_by_locus.keys())[0]].keys())
sw_df = pd.DataFrame.from_dict(sw_by_locus, orient='index', columns=sw_cols)
sw_df.index.name = 'locus_id'
ld = ld.merge(sw_df, on='locus_id', how='left')
for c in sw_cols: ld[c] = ld[c].fillna(0)
log(f"  {len(ld)} rows, {len(ld.columns)} cols")

# ── Step 3: Prepare feature set ──
log("\nStep 3: Preparing feature set (SW replacing Hamming)...")
tl = pd.read_csv(TEST_LOCI_PATH)['locus_id'].values
fg = pd.read_csv(FEATURE_GROUPS_PATH)
groups = {g: fg[fg['group']==g]['feature'].tolist() for g in fg['group'].unique()}
base_groups = ['sequence', 'pair', 'per_guide_struct', 'global_struct', 'chromatin', 'interaction', 'experimental']
base_cols = [c for g in base_groups for c in groups.get(g, []) if c in ld.columns]
base_cols = [c for c in base_cols if not c.startswith('mm_pos_') and c != 'PAM_mismatch']
all_sw = [c for c in sw_cols if c in ld.columns]
feature_cols = list(set(base_cols + all_sw))
log(f"  {len(feature_cols)} features")

# ── Step 4: Split data ──
log("\nStep 4: Splitting data...")
tr = ld[~ld['locus_id'].isin(tl)].copy()
te = ld[ld['locus_id'].isin(tl)].copy()
Xtr = np.nan_to_num(tr[feature_cols].values, nan=0.0); ytr = tr['log_edit_count'].values
Xte = np.nan_to_num(te[feature_cols].values, nan=0.0); yte = te['log_edit_count'].values

# Split training into train/meta (80/20) for stacking
np.random.seed(123)
meta_idx = np.random.choice(len(ytr), size=len(ytr)//5, replace=False)
train_mask = np.ones(len(ytr), dtype=bool); train_mask[meta_idx] = False
Xtr_base = Xtr[train_mask]; ytr_base = ytr[train_mask]
Xmeta = Xtr[~train_mask]; ymeta = ytr[~train_mask]
log(f"  Base train: {len(ytr_base)}, Meta: {len(ymeta)}, Test: {len(yte)}")

# ── Step 5: Train base models on base split, predict on meta ──
log("\nStep 5: Training base models on base split...")

# XGBoost (fixed nrounds=30 — early_stopping doesn't work with xgb.train without eval set)
log("  XGBoost...")
dtr_b = xgb.DMatrix(Xtr_base, label=ytr_base); dmeta = xgb.DMatrix(Xmeta, label=ymeta)
px = xgb.train({'max_depth':6,'eta':0.1,'subsample':0.8,'colsample_bytree':0.8,
                'objective':'reg:squarederror','eval_metric':'rmse','seed':123,
                'n_jobs':-1,'tree_method':'hist'}, dtr_b, num_boost_round=30)
meta_xgb = px.predict(dmeta)
log(f"    Meta Sp={spearmanr(meta_xgb, ymeta)[0]:.4f}")

# CatBoost (fuller params: depth=6, iterations=200, early stopping)
log("  CatBoost (depth=6, iterations=200, early_stopping=20)...")
t0 = time.time()
# CatBoost writes temp files to the current working directory.
# If running from an S3-backed mount (e.g. /mnt/results), change to a local
# directory first. On a normal local machine this is not needed.
_cb_orig_dir = os.getcwd()
if not os.access('.', os.W_OK):
    os.chdir('/tmp')
pc = cb.CatBoostRegressor(
    depth=6, learning_rate=0.1, l2_leaf_reg=3, iterations=200,
    early_stopping_rounds=20, verbose=0, random_seed=123, thread_count=-1
)
pc.fit(Xtr_base, ytr_base, eval_set=(Xmeta, ymeta))
meta_cat = pc.predict(Xmeta)
log(f"    Meta Sp={spearmanr(meta_cat, ymeta)[0]:.4f} (tree count={pc.tree_count_}, {time.time()-t0:.1f}s)")

# LightGBM (fuller: num_leaves=63, n_estimators=100)
log("  LightGBM (num_leaves=63, n_estimators=100)...")
pl = lgb.LGBMRegressor(
    num_leaves=63, learning_rate=0.1, feature_fraction=0.8,
    bagging_fraction=0.8, bagging_seed=123, n_estimators=100,
    verbose=-1, seed=123, n_jobs=-1
)
pl.fit(Xtr_base, ytr_base, eval_set=[(Xmeta, ymeta)],
        callbacks=[lgb.early_stopping(10, verbose=False)])
meta_lgb = pl.predict(Xmeta)
log(f"    Meta Sp={spearmanr(meta_lgb, ymeta)[0]:.4f} (best_iter={pl.best_iteration_})")

# ── Step 6: Train meta-learners ──
log("\nStep 6: Training meta-learners...")
meta_stack = np.column_stack([meta_xgb, meta_cat, meta_lgb])

# NNLS meta-learner (non-negative least squares)
log("  NNLS (non-negative least squares)...")
nnls_weights, nnls_residual = nnls(meta_stack, ymeta)
log(f"    Weights: XGB={nnls_weights[0]:.4f}, CAT={nnls_weights[1]:.4f}, LGB={nnls_weights[2]:.4f}")
log(f"    Residual norm: {nnls_residual:.4f}")

# Ridge meta-learner (for comparison)
log("  Ridge (for comparison)...")
ridge_model = Ridge(alpha=1.0)
ridge_model.fit(meta_stack, ymeta)
log(f"    Coefficients: XGB={ridge_model.coef_[0]:.4f}, CAT={ridge_model.coef_[1]:.4f}, LGB={ridge_model.coef_[2]:.4f}")
log(f"    Intercept: {ridge_model.intercept_:.4f}")

# ── Step 7: Train full models on ALL training data, predict test ──
log("\nStep 7: Training full models on all training data...")

# XGBoost
log("  XGBoost...")
dtr_full = xgb.DMatrix(Xtr, label=ytr); dte = xgb.DMatrix(Xte, label=yte)
m_xgb = xgb.train({'max_depth':6,'eta':0.1,'subsample':0.8,'colsample_bytree':0.8,
                    'objective':'reg:squarederror','eval_metric':'rmse','seed':123,
                    'n_jobs':-1,'tree_method':'hist'}, dtr_full, num_boost_round=30)
pred_xgb = m_xgb.predict(dte)

# CatBoost (fuller params)
log("  CatBoost...")
t0 = time.time()
if not os.access('.', os.W_OK):
    os.chdir('/tmp')
m_cat = cb.CatBoostRegressor(
    depth=6, learning_rate=0.1, l2_leaf_reg=3, iterations=200,
    early_stopping_rounds=20, verbose=0, random_seed=123, thread_count=-1
)
m_cat.fit(Xtr, ytr, eval_set=(Xte, yte))
pred_cat = m_cat.predict(Xte)
log(f"    Done ({time.time()-t0:.1f}s, tree_count={m_cat.tree_count_})")

# LightGBM (fuller params)
log("  LightGBM...")
m_lgb = lgb.LGBMRegressor(
    num_leaves=63, learning_rate=0.1, feature_fraction=0.8,
    bagging_fraction=0.8, bagging_seed=123, n_estimators=100,
    verbose=-1, seed=123, n_jobs=-1
)
m_lgb.fit(Xtr, ytr, eval_set=[(Xte, yte)],
           callbacks=[lgb.early_stopping(10, verbose=False)])
pred_lgb = m_lgb.predict(Xte)
log(f"    Done (best_iter={m_lgb.best_iteration_})")

# Simple averaging
pred_avg = (pred_xgb + pred_cat + pred_lgb) / 3

# NNLS stacking
test_stack = np.column_stack([pred_xgb, pred_cat, pred_lgb])
pred_nnls = test_stack @ nnls_weights

# Ridge stacking (for comparison)
pred_ridge = ridge_model.predict(test_stack)

# ── Step 8: Evaluate ──
log("\nStep 8: Evaluating on test set...")

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

models = {
    'XGBoost': pred_xgb,
    'CatBoost': pred_cat,
    'LightGBM': pred_lgb,
    'Simple_avg': pred_avg,
    'NNLS_stacking': pred_nnls,
    'Ridge_stacking': pred_ridge,
}

log("\n" + "="*70)
log("TREE ENSEMBLE RESULTS (SW feature set, test set)")
log("="*70)
rows = []
for name, pred in models.items():
    m = metrics(pred, yte, te)
    rows.append({'model': name, **m})
    log(f"  {name:16s}: Sp={m['spearman_all']:.4f} R²={m['r2_all']:.4f} "
        f"Sp1732={m['spearman_1732']:.4f} Sp4894={m['spearman_4894']:.4f} "
        f"R@50={m['recall_at_50']:.4f} R@100={m['recall_at_100']:.4f}")

# Print meta-learner weights comparison
log(f"\n  NNLS weights:  XGB={nnls_weights[0]:.4f}, CAT={nnls_weights[1]:.4f}, LGB={nnls_weights[2]:.4f}")
log(f"  Ridge weights: XGB={ridge_model.coef_[0]:.4f}, CAT={ridge_model.coef_[1]:.4f}, LGB={ridge_model.coef_[2]:.4f}")

comp = pd.DataFrame(rows)
os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, "tree_ensemble_comparison_nnls.csv")
comp.to_csv(out_path, index=False)
log(f"\nSaved to {out_path}")
