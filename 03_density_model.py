#!/usr/bin/env python3
"""
03_density_model.py — Zero-inflated mixture model for density estimation.

Component 1 (Classifier): P(edit_count > 0 | features)
Component 2 (Regressor): E[log_edit_count | edit_count > 0, features]
Combined: P(edited) × E[count|edited]

Also: cross-guide transfer analysis (train on 1732, test on 4894, and vice versa).
Uses SW feature set (best from Q3 ablation).

Usage:
    python 03_density_model.py

Requires: long_data.csv, test_loci.csv, feature_groups.csv from the R pipeline (Step 1).
"""
import pandas as pd, numpy as np, xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
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

# ── SW feature computation ──
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
        return {"sw_score":0,"sw_score_norm":0,"sw_align_len":0,"sw_num_gaps":0,
                "sw_num_gap_opens":0,"sw_match_count":0,"sw_mismatch_count":0,
                "sw_start_offset":0,"sw_pam_match":0,
                **{f"sw_gap_pos_{i}":0 for i in range(1,24)},
                **{f"sw_mm_pos_{i}":0 for i in range(1,24)}}
    a=alns[0]; sa=a.seqA[a.start:a.end]; ga=a.seqB[a.start:a.end]
    score=a.score; alen=len(sa)
    ngaps=sa.count('-')+ga.count('-')
    ngopens=0; ig=False
    for s,g in zip(sa,ga):
        if s=='-' or g=='-':
            if not ig: ngopens+=1; ig=True
        else: ig=False
    mc=sum(1 for s,g in zip(sa,ga) if s!='-' and g!='-' and s==g)
    mmc=sum(1 for s,g in zip(sa,ga) if s!='-' and g!='-' and s!=g)
    gp=np.zeros(23,dtype=int); mp=np.zeros(23,dtype=int); gpp=0
    for s,g in zip(sa,ga):
        if g=='-': continue
        if gpp>=23: break
        if s=='-': gp[gpp]=1
        elif s!=g: mp[gpp]=1
        gpp+=1
    pam=1 if all(mp[i]==0 and gp[i]==0 for i in range(20,23)) else 0
    f={"sw_score":float(score),"sw_score_norm":float(score)/max(1,alen),
       "sw_align_len":int(alen),"sw_num_gaps":int(ngaps),"sw_num_gap_opens":int(ngopens),
       "sw_match_count":int(mc),"sw_mismatch_count":int(mmc),
       "sw_start_offset":int(a.start-100),"sw_pam_match":int(pam)}
    for i in range(23): f[f"sw_gap_pos_{i+1}"]=int(gp[i]); f[f"sw_mm_pos_{i+1}"]=int(mp[i])
    return f

sw_by_locus = {}
t0=time.time()
for g in ["1732","4894"]:
    for idx in guide_indices[g]:
        sw_by_locus[df.iloc[idx]["cluster_id"]]=compute_sw(df.iloc[idx]["centered_sequence"],GRNA[g])
log(f"  SW features for {len(sw_by_locus)} loci in {time.time()-t0:.1f}s")

# ── Load and patch long_data ──
log("\nStep 2: Loading long_data...")
if not os.path.exists(LONG_DATA_PATH):
    print(f"ERROR: long_data.csv not found: {LONG_DATA_PATH}", file=sys.stderr)
    print("Please run the R pipeline (Step 1) first.", file=sys.stderr)
    sys.exit(1)
ld=pd.read_csv(LONG_DATA_PATH)
sw_cols=list(sw_by_locus[list(sw_by_locus.keys())[0]].keys())
sw_df=pd.DataFrame.from_dict(sw_by_locus,orient='index',columns=sw_cols);sw_df.index.name='locus_id'
ld=ld.merge(sw_df,on='locus_id',how='left')
for c in sw_cols: ld[c]=ld[c].fillna(0)

# ── Feature set ──
log("\nStep 3: Preparing features...")
tl=pd.read_csv(TEST_LOCI_PATH)['locus_id'].values
fg=pd.read_csv(FEATURE_GROUPS_PATH)
groups={g:fg[fg['group']==g]['feature'].tolist() for g in fg['group'].unique()}
base_groups=['sequence','pair','per_guide_struct','global_struct','chromatin','interaction','experimental']
base_cols=[c for g in base_groups for c in groups.get(g,[]) if c in ld.columns]
base_cols=[c for c in base_cols if not c.startswith('mm_pos_') and c!='PAM_mismatch']
all_sw=[c for c in sw_cols if c in ld.columns]
feature_cols=list(set(base_cols+all_sw))
log(f"  {len(feature_cols)} features")

# ── Helpers ──
def recall_at_k(pred,actual,tdf,k=50):
    rs=[]
    for g in["1732","4894"]:
        c=f"grna_{g}"
        if c not in tdf.columns: continue
        idx=tdf[c]==1
        if idx.sum()==0: continue
        pg=pred[idx];ag=actual[idx];lg=tdf['locus_id'].values[idx]
        nz=ag[ag>0]
        if len(nz)<5: continue
        thr=np.median(nz);th=lg[ag>thr]
        lp=pd.Series(pg).groupby(lg).max()
        tk=lp.sort_values(ascending=False).index[:min(k,len(lp))]
        rs.append(np.isin(th,tk).sum()/len(th))
    return np.mean(rs) if rs else np.nan

def metrics(pred,actual,tdf):
    sp=spearmanr(pred,actual)[0]
    r2=1-np.sum((actual-pred)**2)/np.sum((actual-np.mean(actual))**2)
    return{'spearman_all':sp,'r2_all':r2,
           'recall_at_50':recall_at_k(pred,actual,tdf,50),
           'recall_at_100':recall_at_k(pred,actual,tdf,100)}

# ── Step 4: Train mixture model ──
log("\nStep 4: Training zero-inflated mixture model...")
tr=ld[~ld['locus_id'].isin(tl)].copy(); te=ld[ld['locus_id'].isin(tl)].copy()
Xtr=np.nan_to_num(tr[feature_cols].values,nan=0.0)
ytr_raw=tr['edit_count'].values; ytr_log=tr['log_edit_count'].values
Xte=np.nan_to_num(te[feature_cols].values,nan=0.0)
yte_raw=te['edit_count'].values; yte_log=te['log_edit_count'].values
ytr_bin=(ytr_raw>0).astype(int); yte_bin=(yte_raw>0).astype(int)
log(f"  Train: {len(ytr_bin)} pairs, {ytr_bin.sum()} edited ({ytr_bin.mean()*100:.1f}%)")
log(f"  Test:  {len(yte_bin)} pairs, {yte_bin.sum()} edited ({yte_bin.mean()*100:.1f}%)")

# Component 1: Classifier
log("\n  Component 1: XGBoost classifier P(edited>0)...")
dtr_c=xgb.DMatrix(Xtr,label=ytr_bin); dte_c=xgb.DMatrix(Xte,label=yte_bin)
best_auc,best_cp,best_cr=0,None,0
for d in[4,6]:
    p={'max_depth':d,'eta':0.1,'subsample':0.8,'colsample_bytree':0.8,
       'objective':'binary:logistic','eval_metric':'auc','seed':123,'n_jobs':-1,'tree_method':'hist'}
    cv=xgb.cv(params=p,dtrain=dtr_c,num_boost_round=50,nfold=5,early_stopping_rounds=10,verbose_eval=False)
    bi=cv['test-auc-mean'].idxmax(); auc=cv['test-auc-mean'].max()
    log(f"    depth={d}: AUC={auc:.4f} @ round {bi}")
    if auc>best_auc: best_auc,best_cp,best_cr=auc,p,bi
clf=xgb.train(best_cp,dtr_c,num_boost_round=max(best_cr+1,1))
pred_edited=clf.predict(dte_c)
auc_test=roc_auc_score(yte_bin,pred_edited)
log(f"  Test AUC-ROC: {auc_test:.4f}")

# Component 2: Regressor (nonzero only)
log("\n  Component 2: XGBoost regressor E[count|edited>0]...")
nz_tr=ytr_raw>0; nz_te=yte_raw>0
dtr_r=xgb.DMatrix(Xtr[nz_tr],label=ytr_log[nz_tr])
dte_r=xgb.DMatrix(Xte[nz_te],label=yte_log[nz_te])
dte_full=xgb.DMatrix(Xte)
best_rmse,best_rp,best_rr=float('inf'),None,0
for d in[4,6]:
    p={'max_depth':d,'eta':0.1,'subsample':0.8,'colsample_bytree':0.8,
       'objective':'reg:squarederror','eval_metric':'rmse','seed':123,'n_jobs':-1,'tree_method':'hist'}
    cv=xgb.cv(params=p,dtrain=dtr_r,num_boost_round=50,nfold=5,early_stopping_rounds=10,verbose_eval=False)
    bi=cv['test-rmse-mean'].idxmin(); rmse=cv['test-rmse-mean'].min()
    log(f"    depth={d}: RMSE={rmse:.4f} @ round {bi}")
    if rmse<best_rmse: best_rmse,best_rp,best_rr=rmse,p,bi
reg=xgb.train(best_rp,dtr_r,num_boost_round=max(best_rr+1,1))
pred_mag=reg.predict(dte_full)  # E[count|edited] for ALL test pairs

# Regressor on nonzero only
m_nz=metrics(pred_mag[nz_te],yte_log[nz_te],te[nz_te])
log(f"  Regressor on nonzero: Sp={m_nz['spearman_all']:.4f} R²={m_nz['r2_all']:.4f}")

# ── Step 5: Evaluate ranking strategies ──
log("\n\n"+"="*70)
log("MIXTURE MODEL: Ranking strategy comparison (test set)")
log("="*70)
log(f"  Classifier AUC-ROC: {auc_test:.4f}")

pred_density=pred_edited
pred_magnitude=pred_mag
pred_combined=pred_edited*pred_mag

strategies={'density (P_edited)':pred_density,'magnitude (E[count|ed])':pred_magnitude,
            'combined (P×E)':pred_combined}
for name,pred in strategies.items():
    m=metrics(pred,yte_log,te)
    log(f"  {name:30s}: Sp={m['spearman_all']:.4f} R²={m['r2_all']:.4f} "
        f"R@50={m['recall_at_50']:.4f} R@100={m['recall_at_100']:.4f}")

# ── Step 6: Cross-guide transfer ──
log("\n\n"+"="*70)
log("CROSS-GUIDE TRANSFER ANALYSIS")
log("="*70)

cross_guide_results = []
for tg, gg in [("1732","4894"),("4894","1732")]:
    log(f"\n  Train on gRNA {tg} → Test on gRNA {gg}:")
    tr_m=tr[f'grna_{tg}']==1; te_m=te[f'grna_{gg}']==1
    Xg_tr=np.nan_to_num(tr[tr_m][feature_cols].values,nan=0.0)
    yg_tr_raw=tr[tr_m]['edit_count'].values; yg_tr_log=tr[tr_m]['log_edit_count'].values
    yg_tr_bin=(yg_tr_raw>0).astype(int)
    Xg_te=np.nan_to_num(te[te_m][feature_cols].values,nan=0.0)
    yg_te_raw=te[te_m]['edit_count'].values; yg_te_log=te[te_m]['log_edit_count'].values
    yg_te_bin=(yg_te_raw>0).astype(int)
    log(f"    Train: {len(yg_tr_bin)} ({yg_tr_bin.sum()} edited), Test: {len(yg_te_bin)} ({yg_te_bin.sum()} edited)")

    # Classifier
    dg_c=xgb.DMatrix(Xg_tr,label=yg_tr_bin); dg_ct=xgb.DMatrix(Xg_te)
    cg=xgb.train({'max_depth':6,'eta':0.1,'subsample':0.8,'colsample_bytree':0.8,
                  'objective':'binary:logistic','eval_metric':'auc','seed':123,
                  'n_jobs':-1,'tree_method':'hist'},dg_c,num_boost_round=30)
    pg_ed=cg.predict(dg_ct)
    auc_g=roc_auc_score(yg_te_bin,pg_ed) if yg_te_bin.sum()>0 else np.nan

    # Regressor
    nz_g=yg_tr_raw>0; nzg_te=yg_te_raw>0
    if nz_g.sum()>10 and nzg_te.sum()>10:
        dg_r=xgb.DMatrix(Xg_tr[nz_g],label=yg_tr_log[nz_g])
        rg=xgb.train({'max_depth':6,'eta':0.1,'subsample':0.8,'colsample_bytree':0.8,
                      'objective':'reg:squarederror','eval_metric':'rmse','seed':123,
                      'n_jobs':-1,'tree_method':'hist'},dg_r,num_boost_round=30)
        pg_mag=rg.predict(dg_ct)
        sp_g=spearmanr(pg_mag[nzg_te],yg_te_log[nzg_te])[0]
    else:
        pg_mag=np.zeros(len(yg_te_log)); sp_g=np.nan

    sp_c=spearmanr(pg_ed*pg_mag,yg_te_log)[0]
    log(f"    Classifier AUC: {auc_g:.4f}")
    log(f"    Regressor Sp (nonzero): {sp_g:.4f}")
    log(f"    Combined Sp (all): {sp_c:.4f}")
    cross_guide_results.append({
        'train_guide': tg, 'test_guide': gg,
        'classifier_auc': auc_g, 'regressor_sp_nz': sp_g, 'combined_sp': sp_c
    })

# ── Save ──
log("\n\nSaving...")
res={'classifier_auc':auc_test,
     'density_sp':metrics(pred_density,yte_log,te)['spearman_all'],
     'magnitude_sp':metrics(pred_magnitude,yte_log,te)['spearman_all'],
     'combined_sp':metrics(pred_combined,yte_log,te)['spearman_all'],
     'density_r50':metrics(pred_density,yte_log,te)['recall_at_50'],
     'magnitude_r50':metrics(pred_magnitude,yte_log,te)['recall_at_50'],
     'combined_r50':metrics(pred_combined,yte_log,te)['recall_at_50'],
     'density_r100':metrics(pred_density,yte_log,te)['recall_at_100'],
     'magnitude_r100':metrics(pred_magnitude,yte_log,te)['recall_at_100'],
     'combined_r100':metrics(pred_combined,yte_log,te)['recall_at_100'],
     'regressor_nz_sp':m_nz['spearman_all'],'regressor_nz_r2':m_nz['r2_all']}
pd.DataFrame([res]).to_csv(os.path.join(OUTPUT_DIR,"density_model_comparison.csv"),index=False)
pd.DataFrame(cross_guide_results).to_csv(os.path.join(OUTPUT_DIR,"cross_guide_transfer.csv"),index=False)
log(f"Saved to {OUTPUT_DIR}/density_model_comparison.csv")
log(f"Saved to {OUTPUT_DIR}/cross_guide_transfer.csv")
log("\nDone.")
