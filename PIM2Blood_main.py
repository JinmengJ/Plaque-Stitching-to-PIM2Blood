#!/usr/bin/env python3
"""
This script is intentionally self-contained for manuscript code release.

"""


from pathlib import Path
import json, re, warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import NMF
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score

# -----------------------------
# Paths and constants
# -----------------------------
BASE = Path('/data2/jiajm/Atherosclerosis_Manuscript/new_panel/final/RADAR_ModulesNet_142')
SAVE = BASE / 'ProportionalSplit_DevTrainExternal_FinalHoldout'
SAVE.mkdir(parents=True, exist_ok=True)

DEV_RAW_FP = BASE / 'train_merged_raw.csv'
EXT_RAW_FP = BASE / 'external_merged_raw.csv'
PIMS_FP = BASE / 'pims_markers_intersect.csv'
UKB_COHORT_FP = Path('/data2/jiajm/Atherosclerosis_Manuscript/new_panel/incident_filtered_cohort_with_full_olink.csv')

META = ['SampleID', 'SubjectID', 'Group', 'Time']
RANDOM_STATE = 20260615
DEV_HOLDOUT_SIZE = 0.30

NMF_K_LIST = [10, 12, 14, 16]
NMF_REPEATS = 5
NMF_ALPHA_H = 0.3
NMF_L1_RATIO = 0.9
NMF_MAX_ITER = 800
MODULE_TARGET_SIZE = 30
MODULE_MIN_SIZE = 10
MODULE_MAX_SIZE = 80
MAX_GENE_MODULE_MEMBERSHIP = 2
N_CORE_MODULES = 4

GENE_GATE_TRAIN_AUC_MIN = 0.60
GENE_GATE_EXTERNAL_AUC_MIN = 0.60
GENE_GATE_TRAIN_COHEN_MIN = 0.40
GENE_GATE_EXTERNAL_COHEN_MIN = 0.40
GENE_GATE_TEMPORAL_D_MIN = 0.30
K_MAX = 40
MIN_GENES_PER_MODULE_FOR_SELECTION = 5

# -----------------------------
# Utilities
# -----------------------------
def normalize_group_time(df):
    x = df.copy()
    x['Group'] = x['Group'].astype(str).str.strip()
    g = x['Group'].str.lower()
    x['Group'] = np.where(g.eq('ami'), 'AMI', np.where(g.eq('stable'), 'Stable', x['Group']))
    x['Time'] = x['Time'].astype(str).str.strip().str.lower()
    tmap = {
        'admission':'adm','admit':'adm','baseline':'adm','0':'adm','0d':'adm','day0':'adm',
        'discharge':'dis','dis':'dis','4-6d':'dis','day5':'dis','5d':'dis',
        '1month':'1m','1_month':'1m','1-mo':'1m','1 mo':'1m',
        '6month':'6m','6_month':'6m','6-mo':'6m','6 mo':'6m',
    }
    x['Time'] = x['Time'].map(lambda v: tmap.get(v, v))
    x.loc[x['Group'].eq('Stable'), 'Time'] = 'adm'
    return x

def zscore_series(s):
    s = pd.to_numeric(s, errors='coerce')
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s*0.0
    return (s - s.mean()) / sd

def zscore_df_by_columns(df):
    return df.apply(zscore_series, axis=0)

def cohen_d(x1, x0):
    x1 = pd.to_numeric(pd.Series(x1), errors='coerce').dropna().values
    x0 = pd.to_numeric(pd.Series(x0), errors='coerce').dropna().values
    if len(x1) < 2 or len(x0) < 2:
        return np.nan
    sp = np.sqrt(((len(x1)-1)*np.var(x1, ddof=1) + (len(x0)-1)*np.var(x0, ddof=1)) / max(len(x1)+len(x0)-2, 1))
    return np.nan if sp == 0 else float((np.mean(x1)-np.mean(x0))/sp)

def safe_ttest(x1, x0):
    x1 = pd.to_numeric(pd.Series(x1), errors='coerce').dropna().values
    x0 = pd.to_numeric(pd.Series(x0), errors='coerce').dropna().values
    if len(x1) < 2 or len(x0) < 2:
        return np.nan
    return float(stats.ttest_ind(x1, x0, equal_var=False).pvalue)

def safe_auc(y, s):
    y = np.asarray(y).astype(int)
    s = pd.to_numeric(pd.Series(s), errors='coerce').values
    ok = np.isfinite(s)
    y, s = y[ok], s[ok]
    if len(np.unique(y)) < 2:
        return np.nan
    return float(roc_auc_score(y, s))

def safe_ap(y, s):
    y = np.asarray(y).astype(int)
    s = pd.to_numeric(pd.Series(s), errors='coerce').values
    ok = np.isfinite(s)
    y, s = y[ok], s[ok]
    if len(np.unique(y)) < 2:
        return np.nan
    return float(average_precision_score(y, s))

def make_easy_mask(df):
    g = df['Group'].astype(str).str.lower()
    t = df['Time'].astype(str).str.lower()
    return (g.eq('ami') & t.eq('adm')) | (g.eq('stable') & t.eq('adm'))

def make_easy_y(df):
    g = df['Group'].astype(str).str.lower()
    t = df['Time'].astype(str).str.lower()
    return (g.eq('ami') & t.eq('adm')).astype(int)

def paired_adm_later_stats(df, score_col, later_time='dis'):
    d = df.copy()
    d = d[d['Group'].astype(str).str.lower().eq('ami')].copy()
    d = d[d['Time'].isin(['adm', later_time])].copy()
    if d.empty:
        return {'n_pair':0, f'adm_minus_{later_time}_mean':np.nan, f'adm_minus_{later_time}_d':np.nan, f'adm_minus_{later_time}_p':np.nan}
    piv = d.pivot_table(index='SubjectID', columns='Time', values=score_col, aggfunc='mean')
    if 'adm' not in piv.columns or later_time not in piv.columns:
        return {'n_pair':0, f'adm_minus_{later_time}_mean':np.nan, f'adm_minus_{later_time}_d':np.nan, f'adm_minus_{later_time}_p':np.nan}
    piv = piv[['adm', later_time]].dropna()
    if len(piv) < 2:
        return {'n_pair':int(len(piv)), f'adm_minus_{later_time}_mean':np.nan, f'adm_minus_{later_time}_d':np.nan, f'adm_minus_{later_time}_p':np.nan}
    diff = piv['adm'] - piv[later_time]
    sd = diff.std(ddof=1)
    dval = diff.mean()/sd if sd and np.isfinite(sd) else np.nan
    pval = stats.ttest_rel(piv['adm'], piv[later_time]).pvalue
    return {'n_pair':int(len(piv)), f'adm_minus_{later_time}_mean':float(diff.mean()), f'adm_minus_{later_time}_d':float(dval) if np.isfinite(dval) else np.nan, f'adm_minus_{later_time}_p':float(pval)}

def assoc_stats_for_score(df, score_col, prefix):
    mask = make_easy_mask(df)
    y = make_easy_y(df)
    s = pd.to_numeric(df[score_col], errors='coerce')
    x1 = s[mask & (y == 1)]
    x0 = s[mask & (y == 0)]
    auc = safe_auc(y[mask], s[mask])
    direction = 'AMI_higher' if (x1.mean() - x0.mean()) >= 0 else 'Stable_higher'
    out = {
        f'{prefix}_n_easy': int(mask.sum()),
        f'{prefix}_n_AMI_adm': int((mask & (y == 1)).sum()),
        f'{prefix}_n_Stable': int((mask & (y == 0)).sum()),
        f'{prefix}_mean_AMIadm': float(x1.mean()) if len(x1) else np.nan,
        f'{prefix}_mean_Stable': float(x0.mean()) if len(x0) else np.nan,
        f'{prefix}_diff_AMI_minus_Stable': float(x1.mean() - x0.mean()) if len(x1) and len(x0) else np.nan,
        f'{prefix}_cohen_d': cohen_d(x1, x0),
        f'{prefix}_p_AMI_vs_Stable': safe_ttest(x1, x0),
        f'{prefix}_AUC': auc,
        f'{prefix}_AP': safe_ap(y[mask], s[mask]),
        f'{prefix}_baseline_AP': float(np.mean(y[mask])) if mask.sum() else np.nan,
        f'{prefix}_direction': direction,
    }
    pair = paired_adm_later_stats(df.assign(__score__=s.values), '__score__', 'dis')
    out[f'{prefix}_adm_minus_dis_mean'] = pair['adm_minus_dis_mean']
    out[f'{prefix}_adm_minus_dis_d'] = pair['adm_minus_dis_d']
    out[f'{prefix}_adm_minus_dis_p'] = pair['adm_minus_dis_p']
    out[f'{prefix}_n_pair_dis'] = pair['n_pair']
    return out

def vector_z(v, positive_only=False):
    s = pd.to_numeric(pd.Series(v), errors='coerce').fillna(0.0)
    if positive_only:
        s = s.clip(lower=0)
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean()) / sd

def canonical(x):
    return re.sub(r'[^a-z0-9]+', '', str(x).lower())

def find_olink_col(df, gene):
    candidates = [f'olink_{gene.lower()}', gene.lower(), gene.upper(), gene]
    for c in candidates:
        if c in df.columns:
            return c
    targets = {canonical(f'olink_{gene}'), canonical(gene)}
    for c in df.columns:
        if canonical(c) in targets:
            return c
    return None

# -----------------------------
# 1. Load and split Dev first
# -----------------------------
dev_full = normalize_group_time(pd.read_csv(DEV_RAW_FP, low_memory=False))
ext_raw = normalize_group_time(pd.read_csv(EXT_RAW_FP, low_memory=False))
pim_candidates = pd.read_csv(PIMS_FP)['gene_in_data'].astype(str).str.strip().tolist()

gene_cols_dev = [c for c in dev_full.columns if c not in META]
gene_cols_ext = [c for c in ext_raw.columns if c not in META]
common_microarray_genes = set(gene_cols_dev) & set(gene_cols_ext)

if UKB_COHORT_FP.exists():
    ukb_header = pd.read_csv(UKB_COHORT_FP, nrows=0, low_memory=False)
    ukb_header_df = pd.DataFrame(columns=ukb_header.columns)
    pim_genes = [
        g for g in pim_candidates
        if g in common_microarray_genes and find_olink_col(ukb_header_df, g) is not None
    ]
else:
    pim_genes = [g for g in pim_candidates if g in common_microarray_genes]

pd.DataFrame({'gene': pim_genes}).to_csv(SAVE / 'pims_olink_microarray_intersection.csv', index=False)
if len(pim_genes) == 0:
    raise ValueError('No PIM/Olink/microarray intersecting markers were found.')

subject_table = dev_full.groupby('SubjectID').agg(
    label=('Group', lambda x: 'AMI' if (x.astype(str).str.lower() == 'ami').any() else 'Stable'),
    n_samples=('SampleID', 'size')
).reset_index()
train_subjects, holdout_subjects = train_test_split(
    subject_table['SubjectID'],
    test_size=DEV_HOLDOUT_SIZE,
    random_state=RANDOM_STATE,
    stratify=subject_table['label']
)
train_subjects = set(train_subjects)
holdout_subjects = set(holdout_subjects)
dev_train = dev_full[dev_full['SubjectID'].isin(train_subjects)].copy().reset_index(drop=True)
dev_holdout = dev_full[dev_full['SubjectID'].isin(holdout_subjects)].copy().reset_index(drop=True)

split_summary = pd.DataFrame([
    {'dataset':'Dev_full', 'n_samples':len(dev_full), 'n_subjects':dev_full['SubjectID'].nunique(), 'n_easy':int(make_easy_mask(dev_full).sum()), 'n_AMI_adm':int(make_easy_y(dev_full[make_easy_mask(dev_full)]).sum()), 'n_Stable':int((1-make_easy_y(dev_full[make_easy_mask(dev_full)])).sum())},
    {'dataset':'Dev_train', 'n_samples':len(dev_train), 'n_subjects':dev_train['SubjectID'].nunique(), 'n_easy':int(make_easy_mask(dev_train).sum()), 'n_AMI_adm':int(make_easy_y(dev_train[make_easy_mask(dev_train)]).sum()), 'n_Stable':int((1-make_easy_y(dev_train[make_easy_mask(dev_train)])).sum())},
    {'dataset':'Dev_holdout', 'n_samples':len(dev_holdout), 'n_subjects':dev_holdout['SubjectID'].nunique(), 'n_easy':int(make_easy_mask(dev_holdout).sum()), 'n_AMI_adm':int(make_easy_y(dev_holdout[make_easy_mask(dev_holdout)]).sum()), 'n_Stable':int((1-make_easy_y(dev_holdout[make_easy_mask(dev_holdout)])).sum())},
    {'dataset':'External_constraint', 'n_samples':len(ext_raw), 'n_subjects':ext_raw['SubjectID'].nunique(), 'n_easy':int(make_easy_mask(ext_raw).sum()), 'n_AMI_adm':int(make_easy_y(ext_raw[make_easy_mask(ext_raw)]).sum()), 'n_Stable':int((1-make_easy_y(ext_raw[make_easy_mask(ext_raw)])).sum())},
])
split_summary.to_csv(SAVE / 'dev_subject_split_summary.csv', index=False)
dev_train.to_csv(SAVE / 'dev_train_raw.csv', index=False)
dev_holdout.to_csv(SAVE / 'dev_holdout_raw.csv', index=False)
ext_raw.to_csv(SAVE / 'external_constraint_raw.csv', index=False)
print('[Split]')
print(split_summary.to_string(index=False))
print('PIM genes used:', len(pim_genes))

# -----------------------------
# 2. NMF on Dev_train only
# -----------------------------
X_raw = dev_train[pim_genes].apply(pd.to_numeric, errors='coerce')
X_raw = X_raw.fillna(X_raw.mean(axis=0))
X = (X_raw - X_raw.min(axis=0)) + 1e-9

def hoyer_sparsity(A):
    A = np.asarray(A); n = A.size
    if n == 0: return 0.0
    l1 = np.sum(np.abs(A)); l2 = np.sqrt(np.sum(A*A)) + 1e-12
    return float((np.sqrt(n) - (l1/l2)) / (np.sqrt(n) - 1 + 1e-12))

nmf_rows = []
best = None
for K in NMF_K_LIST:
    for r in range(NMF_REPEATS):
        seed = RANDOM_STATE + 100*K + r
        nmf = NMF(n_components=K, init='nndsvda', random_state=seed,
                  alpha_W=0.0, alpha_H=NMF_ALPHA_H, l1_ratio=NMF_L1_RATIO,
                  beta_loss='frobenius', solver='cd', max_iter=NMF_MAX_ITER)
        W = nmf.fit_transform(X.values)
        H = nmf.components_
        recon = float(nmf.reconstruction_err_)
        sparsity = hoyer_sparsity(H)
        score = sparsity / (recon + 1e-9)
        row = {'K':K, 'repeat':r, 'seed':seed, 'reconstruction_err':recon, 'sparsity_H':sparsity, 'selection_score':score}
        nmf_rows.append(row)
        cand = {**row, 'W':W, 'H':H}
        if best is None or cand['selection_score'] > best['selection_score']:
            best = cand

K_NMF = int(best['K'])
mods = [f'Mod{i:02d}' for i in range(K_NMF)]
H = best['H']
W = best['W']
H_mod_gene = pd.DataFrame(H, index=mods, columns=pim_genes)
H_gene_mod = H_mod_gene.T.copy()

pd.DataFrame(nmf_rows).to_csv(SAVE / 'nmf_K_selection_devtrain_only.csv', index=False)
H_mod_gene.to_csv(SAVE / 'nmf_H_gene_loadings_devtrain_only.csv')
pd.DataFrame(W, index=dev_train['SampleID'], columns=mods).to_csv(SAVE / 'nmf_W_sample_scores_devtrain_only.csv')

assigned = {g:0 for g in pim_genes}
module_map = {}
for m in mods:
    vals = H_gene_mod[m].sort_values(ascending=False)
    glist = []
    for g, val in vals.head(MODULE_MAX_SIZE).items():
        if val <= 0:
            continue
        if assigned[g] < MAX_GENE_MODULE_MEMBERSHIP:
            glist.append(g)
            assigned[g] += 1
        if len(glist) >= MODULE_TARGET_SIZE:
            break
    module_map[m] = glist

# Keep modules with enough positive selected genes.
valid_modules = [m for m, gl in module_map.items() if len(gl) >= MODULE_MIN_SIZE]
H_sel = pd.DataFrame(0.0, index=pim_genes, columns=valid_modules)
for m in valid_modules:
    H_sel.loc[module_map[m], m] = H_gene_mod.loc[module_map[m], m]
H_w = H_sel.divide(H_sel.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
H_w.to_csv(SAVE / 'module_membership_weighted_devtrainNMF.csv')
with open(SAVE / 'module_map_devtrainNMF.json', 'w') as f:
    json.dump({m: module_map[m] for m in valid_modules}, f, indent=2)
pd.DataFrame({'module':valid_modules, 'n_genes':[len(module_map[m]) for m in valid_modules]}).to_csv(SAVE / 'module_sizes_devtrainNMF.csv', index=False)
print('[NMF] chosen K=', K_NMF, 'valid modules=', valid_modules)

# -----------------------------
# 3. Module scores and module priority using Dev_train + External only
# -----------------------------
def score_modules_for_dataset(df, weights):
    out = df[META].copy()
    genes = list(weights.index)
    Xd = df[genes].apply(pd.to_numeric, errors='coerce')
    Zd = zscore_df_by_columns(Xd)
    for m in weights.columns:
        w = weights[m]
        use = w[w > 0]
        if use.empty:
            out[m] = np.nan
        else:
            out[m] = Zd[use.index].multiply(use, axis=1).sum(axis=1) / use.sum()
    return out

mod_score_devtrain = score_modules_for_dataset(dev_train, H_w)
mod_score_external = score_modules_for_dataset(ext_raw, H_w)
mod_score_holdout = score_modules_for_dataset(dev_holdout, H_w)
mod_score_devtrain.to_csv(SAVE / 'module_scores_devtrain_weighted.csv', index=False)
mod_score_external.to_csv(SAVE / 'module_scores_external_constraint_weighted.csv', index=False)
mod_score_holdout.to_csv(SAVE / 'module_scores_devholdout_weighted_final_only.csv', index=False)

module_rows = []
for m in valid_modules:
    row = {'module':m, 'n_module_genes':int((H_w[m] > 0).sum())}
    row.update(assoc_stats_for_score(mod_score_devtrain, m, 'devtrain'))
    row.update(assoc_stats_for_score(mod_score_external, m, 'external'))
    row.update(assoc_stats_for_score(mod_score_holdout, m, 'devholdout_shown_not_used'))
    module_rows.append(row)
module_rank = pd.DataFrame(module_rows)
module_rank['module_priority_score'] = (
    0.35*vector_z(module_rank['devtrain_cohen_d'], True).values +
    0.25*vector_z(module_rank['devtrain_adm_minus_dis_d'], True).values +
    0.25*vector_z(module_rank['external_cohen_d'], True).values +
    0.15*vector_z(module_rank['external_AUC'].fillna(0.5), False).values
)
module_rank['eligible_core_module'] = (
    module_rank['devtrain_direction'].eq('AMI_higher') &
    module_rank['external_direction'].eq('AMI_higher') &
    (module_rank['devtrain_cohen_d'] > 0) &
    (module_rank['external_cohen_d'] > 0)
)
module_rank = module_rank.sort_values(['eligible_core_module','module_priority_score'], ascending=[False, False]).reset_index(drop=True)
module_rank['module_rank'] = np.arange(1, len(module_rank)+1)
module_rank.to_csv(SAVE / 'module_rank_devtrain_external.csv', index=False)
core_modules = module_rank.head(N_CORE_MODULES)['module'].tolist()
pd.DataFrame({'core_module_order':np.arange(1,len(core_modules)+1), 'module':core_modules}).to_csv(SAVE / 'core_modules_selected_devtrain_external.csv', index=False)
print('[Core modules]', core_modules)
print(module_rank[['module','eligible_core_module','module_priority_score','devtrain_AUC','external_AUC','devtrain_cohen_d','external_cohen_d']].head(8).to_string(index=False))

# -----------------------------
# 4. Gene stats, gate, and ranking using Dev_train + External only
# -----------------------------
def gene_stats_one_dataset(df, genes, prefix):
    rows = []
    mask = make_easy_mask(df)
    y = make_easy_y(df)
    for g in genes:
        s = pd.to_numeric(df[g], errors='coerce')
        x1 = s[mask & (y == 1)]
        x0 = s[mask & (y == 0)]
        diff = x1.mean() - x0.mean()
        row = {
            'gene':g,
            f'{prefix}_gene_mean_AMIadm':float(x1.mean()) if len(x1) else np.nan,
            f'{prefix}_gene_mean_Stable':float(x0.mean()) if len(x0) else np.nan,
            f'{prefix}_gene_diff_AMI_minus_Stable':float(diff) if np.isfinite(diff) else np.nan,
            f'{prefix}_gene_cohen_d':cohen_d(x1,x0),
            f'{prefix}_gene_p':safe_ttest(x1,x0),
            f'{prefix}_gene_auc':safe_auc(y[mask], s[mask]),
            f'{prefix}_gene_direction':'AMI_higher' if diff >= 0 else 'Stable_higher',
        }
        pair = paired_adm_later_stats(df.assign(__score__=s.values), '__score__', 'dis')
        row[f'{prefix}_gene_n_pair'] = pair['n_pair']
        row[f'{prefix}_gene_adm_minus_dis_mean'] = pair['adm_minus_dis_mean']
        row[f'{prefix}_gene_adm_minus_dis_d'] = pair['adm_minus_dis_d']
        row[f'{prefix}_gene_adm_minus_dis_p'] = pair['adm_minus_dis_p']
        rows.append(row)
    return pd.DataFrame(rows)

candidate_genes = sorted(set(H_w.index[H_w[core_modules].sum(axis=1) > 0]))
gs_train = gene_stats_one_dataset(dev_train, candidate_genes, 'devtrain')
gs_ext = gene_stats_one_dataset(ext_raw, candidate_genes, 'external')
gene_stats = gs_train.merge(gs_ext, on='gene', how='outer')

records = []
for m_i, m in enumerate(core_modules, start=1):
    w = H_w[m]
    genes_m = w[w > 0].sort_values(ascending=False).index.tolist()
    d = pd.DataFrame({'gene':genes_m, 'module_use':m, 'module_order':m_i, 'weighted_membership':w.loc[genes_m].values})
    d = d.merge(gene_stats, on='gene', how='left')
    d['gene_direction_use'] = np.where(d['devtrain_gene_direction'].eq('AMI_higher'), 'AMI_higher', d['external_gene_direction'])
    d = d[d['gene_direction_use'].eq('AMI_higher')].copy()
    tr_auc = pd.to_numeric(d['devtrain_gene_auc'], errors='coerce')
    ex_auc = pd.to_numeric(d['external_gene_auc'], errors='coerce')
    tr_d = pd.to_numeric(d['devtrain_gene_cohen_d'], errors='coerce')
    ex_d = pd.to_numeric(d['external_gene_cohen_d'], errors='coerce')
    temporal = pd.to_numeric(d['devtrain_gene_adm_minus_dis_d'], errors='coerce')
    d['passed_gene_gate'] = (
        (tr_auc >= GENE_GATE_TRAIN_AUC_MIN) |
        (ex_auc >= GENE_GATE_EXTERNAL_AUC_MIN) |
        (tr_d >= GENE_GATE_TRAIN_COHEN_MIN) |
        (ex_d >= GENE_GATE_EXTERNAL_COHEN_MIN) |
        (temporal >= GENE_GATE_TEMPORAL_D_MIN)
    ).fillna(False)
    reasons = []
    for _, r in d.iterrows():
        rr = []
        if pd.notna(r.get('devtrain_gene_auc')) and r['devtrain_gene_auc'] >= GENE_GATE_TRAIN_AUC_MIN: rr.append('devtrain_auc')
        if pd.notna(r.get('external_gene_auc')) and r['external_gene_auc'] >= GENE_GATE_EXTERNAL_AUC_MIN: rr.append('external_auc')
        if pd.notna(r.get('devtrain_gene_cohen_d')) and r['devtrain_gene_cohen_d'] >= GENE_GATE_TRAIN_COHEN_MIN: rr.append('devtrain_cohen_d')
        if pd.notna(r.get('external_gene_cohen_d')) and r['external_gene_cohen_d'] >= GENE_GATE_EXTERNAL_COHEN_MIN: rr.append('external_cohen_d')
        if pd.notna(r.get('devtrain_gene_adm_minus_dis_d')) and r['devtrain_gene_adm_minus_dis_d'] >= GENE_GATE_TEMPORAL_D_MIN: rr.append('temporal_d')
        reasons.append(';'.join(rr) if rr else 'none')
    d['gene_gate_reason'] = reasons
    d['gene_score'] = (
        0.30*vector_z(d['weighted_membership']).values +
        0.20*vector_z(d['devtrain_gene_auc'].fillna(0.5)).values +
        0.15*vector_z(d['external_gene_auc'].fillna(0.5)).values +
        0.15*vector_z(d['devtrain_gene_cohen_d'], True).values +
        0.10*vector_z(d['external_gene_cohen_d'], True).values +
        0.10*vector_z(d['devtrain_gene_adm_minus_dis_d'], True).values
    )
    d = d.sort_values(['passed_gene_gate','gene_score','weighted_membership'], ascending=[False,False,False]).reset_index(drop=True)
    d['rank_within_module'] = np.arange(1, len(d)+1)
    records.append(d)

gene_pool = pd.concat(records, ignore_index=True)
gene_pool.to_csv(SAVE / 'candidate_gene_pool_devtrain_external.csv', index=False)
print('[Gene pool]', gene_pool.shape)
print(gene_pool.groupby(['module_use','passed_gene_gate']).size().reset_index(name='n').to_string(index=False))

# -----------------------------
# 5. K scan and panel selection using Dev_train + External only
# -----------------------------
def select_panel_for_k(gene_pool, core_modules, K):
    base_n = K // len(core_modules)
    extra = K % len(core_modules)
    quotas = {m: base_n + (1 if i < extra else 0) for i, m in enumerate(core_modules)}
    rows = []
    used = set()
    for m in core_modules:
        dm = gene_pool[(gene_pool['module_use'] == m) & gene_pool['passed_gene_gate']].sort_values('rank_within_module')
        take = []
        for _, r in dm.iterrows():
            if r['gene'] in used:
                continue
            take.append(r)
            used.add(r['gene'])
            if len(take) >= quotas[m]:
                break
        rows.extend(take)
    # Backfill from passed-gate genes if duplicated membership creates shortfall.
    if len(rows) < K:
        allg = gene_pool[gene_pool['passed_gene_gate']].sort_values(['module_order','rank_within_module'])
        for _, r in allg.iterrows():
            if r['gene'] in used:
                continue
            rows.append(r)
            used.add(r['gene'])
            if len(rows) >= K:
                break
    panel = pd.DataFrame(rows).copy()
    if panel.empty or len(panel) != K:
        raise ValueError(f'Could not select K={K}, selected {len(panel)}')
    panel = panel.sort_values(['module_order','rank_within_module']).reset_index(drop=True)
    panel['panel_rank'] = np.arange(1, len(panel)+1)
    panel['selection_depth'] = panel.groupby('module_use').cumcount() + 1
    panel['K'] = K
    panel['equal_weight'] = 1.0 / K
    counts = panel.groupby('module_use').size().to_dict()
    panel['final_module_count'] = panel['module_use'].map(counts)
    return panel, counts

def fit_panel_params(df_fit, genes):
    Xp = df_fit[genes].apply(pd.to_numeric, errors='coerce')
    mu = Xp.mean(axis=0)
    sd = Xp.std(axis=0, ddof=0).replace(0, np.nan).fillna(1.0)
    return mu, sd

def score_panel(df, genes, mu, sd, score_col='PanelScore'):
    Z = (df[genes].apply(pd.to_numeric, errors='coerce') - mu[genes]) / sd[genes]
    out = df[META].copy()
    out[score_col] = Z.mean(axis=1, skipna=True)
    out[f'{score_col}_n_genes_used'] = Z.notna().sum(axis=1)
    return out

def evaluate_panel_score(scored, score_col='PanelScore'):
    m = assoc_stats_for_score(scored, score_col, 'tmp')
    return {k.replace('tmp_', ''):v for k,v in m.items()}

K_values = list(range(len(core_modules), K_MAX + 1, len(core_modules)))
all_panels = []
k_rows = []
for K in K_values:
    if K / len(core_modules) < MIN_GENES_PER_MODULE_FOR_SELECTION:
        eligible = False
    else:
        eligible = True
    try:
        panel, counts = select_panel_for_k(gene_pool, core_modules, K)
    except Exception as e:
        print('[K skip]', K, e)
        continue
    genes = panel['gene'].astype(str).tolist()
    dev_fit_for_scaling = pd.concat([dev_train, ext_raw], ignore_index=True)
    mu, sd = fit_panel_params(dev_fit_for_scaling, genes)
    sc_dev = score_panel(dev_train, genes, mu, sd, 'PanelScore')
    sc_ext = score_panel(ext_raw, genes, mu, sd, 'PanelScore')
    sc_hold = score_panel(dev_holdout, genes, mu, sd, 'PanelScore')
    md = evaluate_panel_score(sc_dev, 'PanelScore')
    me = evaluate_panel_score(sc_ext, 'PanelScore')
    mh = evaluate_panel_score(sc_hold, 'PanelScore')
    row = {
        'K':K,
        'genes_per_module_nominal':K/len(core_modules),
        'eligible_for_K_selection':eligible,
        'core_modules':';'.join(core_modules),
        'panel_genes':';'.join(genes),
        'module_counts':';'.join([f'{m}:{counts.get(m,0)}' for m in core_modules]),
        'n_passed_gene_gate':int(panel['passed_gene_gate'].sum()),
    }
    for prefix, mm in [('devtrain', md), ('external', me), ('devholdout', mh)]:
        for k, v in mm.items():
            row[f'{prefix}_{k}'] = v
    k_rows.append(row)
    all_panels.append(panel)

kdf = pd.DataFrame(k_rows)
kdf['selection_score_devtrain_external'] = (
    0.30*vector_z(kdf['external_AUC'].fillna(0.5)).values +
    0.25*vector_z(kdf['external_AP'].fillna(0.0)).values +
    0.20*vector_z(kdf['external_cohen_d'], True).values +
    0.20*vector_z(kdf['devtrain_adm_minus_dis_d'], True).values +
    0.05*vector_z(kdf['devtrain_AUC'].fillna(0.5)).values
)
# Compact plateau: choose the smallest eligible K within a small margin of the best development score.
dev_only = kdf[kdf['eligible_for_K_selection']].copy()
max_score = dev_only['selection_score_devtrain_external'].max()
plateau = dev_only[dev_only['selection_score_devtrain_external'] >= max_score - 0.15].sort_values('K')
if plateau.empty:
    best = dev_only.sort_values(['selection_score_devtrain_external','K'], ascending=[False, True]).iloc[0]
    selection_rule = 'fallback_best_devtrain_external_score'
else:
    best = plateau.iloc[0]
    selection_rule = 'devtrain_external_compact_plateau'
best_k = int(best['K'])
kdf['selected_K'] = kdf['K'].eq(best_k)
kdf['selection_rule'] = selection_rule
kdf.to_csv(SAVE / 'Kscan_devtrain_external_with_holdout_shown.csv', index=False)
plateau.to_csv(SAVE / 'Kscan_plateau_candidates_devtrain_external.csv', index=False)
pd.DataFrame([best]).assign(selection_rule=selection_rule, n_plateau_candidates=len(plateau)).to_csv(SAVE / 'selected_K_by_devtrain_external_compactPlateau.csv', index=False)

all_panel_df = pd.concat(all_panels, ignore_index=True)
all_panel_df.to_csv(SAVE / 'all_candidate_panels_devtrain_external.csv', index=False)
frozen_panel = all_panel_df[all_panel_df['K'] == best_k].sort_values('panel_rank').reset_index(drop=True)
frozen_panel['selected_K'] = best_k
frozen_panel['frozen_panel_source'] = 'Subject-level proportional split discovery; Dev_holdout unused for selection'
frozen_panel.to_csv(SAVE / 'selected_panel_by_devtrain_external.csv', index=False)
(SAVE / 'frozen_panel_selected_genes.txt').write_text('\n'.join(frozen_panel['gene'].astype(str).tolist()) + '\n')
print('[Selected K]', best_k, selection_rule)
print(frozen_panel[['panel_rank','gene','module_use','selection_depth','gene_score','passed_gene_gate','gene_gate_reason']].to_string(index=False))

# -----------------------------
# 6. Final microarray scores and metrics for frozen panel
# -----------------------------
LOCKED = frozen_panel['gene'].astype(str).tolist()
mu, sd = fit_panel_params(pd.concat([dev_train, ext_raw], ignore_index=True), LOCKED)
params = pd.DataFrame({'gene':LOCKED, 'mean_fit_devtrain_external':mu[LOCKED].values, 'sd_fit_devtrain_external':sd[LOCKED].values})
params.to_csv(SAVE / 'frozen_panel_score_parameters_fit_on_devtrain_plus_external.csv', index=False)

score_rows = []
metric_rows = []
for role, name, df in [
    ('development', 'Dev_train', dev_train),
    ('development_constraint', 'External_constraint', ext_raw),
    ('final_test', 'Dev_holdout', dev_holdout),
]:
    sc = score_panel(df, LOCKED, mu, sd, 'FrozenPanelScore')
    sc['dataset_role'] = role
    sc['dataset_name'] = name
    score_rows.append(sc)
    metric_rows.append({'dataset_role':role, 'dataset_name':name, **evaluate_panel_score(sc, 'FrozenPanelScore')})

scores_micro = pd.concat(score_rows, ignore_index=True)
metrics_micro = pd.DataFrame(metric_rows)
scores_micro.to_csv(SAVE / 'panel_scores_microarray.csv', index=False)
metrics_micro.to_csv(SAVE / 'panel_metrics_microarray.csv', index=False)
print('[Microarray metrics]')
print(metrics_micro[['dataset_name','n_easy','n_AMI_adm','n_Stable','AUC','AP','cohen_d','adm_minus_dis_d']].to_string(index=False))

# -----------------------------
# 7. Method manifest
# -----------------------------
manifest = {
    'pipeline': 'subject-level proportional split discovery',
    'input_files': {
        'dev_raw': str(DEV_RAW_FP),
        'external_raw': str(EXT_RAW_FP),
        'pim_markers': str(PIMS_FP),
    },
    'data_roles': {
        'Dev_train': 'NMF/module discovery, module/gene/K development, score-scaling fit',
        'External_constraint': 'development stability constraint and score-scaling fit; not independent validation',
        'Dev_holdout': 'final microarray test only; not used for NMF, gene selection, K selection, or score scaling',
    },
    'random_state': RANDOM_STATE,
    'dev_holdout_size': DEV_HOLDOUT_SIZE,
    'selected_K': best_k,
    'selected_genes': LOCKED,
    'core_modules': core_modules,
}
with open(SAVE / 'pipeline_manifest.json', 'w') as f:
    json.dump(manifest, f, indent=2)
print('[DONE] results saved to', SAVE)
