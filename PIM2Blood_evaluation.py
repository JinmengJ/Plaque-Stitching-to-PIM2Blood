#!/usr/bin/env python3
"""Evaluate a frozen 20-marker panel in UKB Olink.

This script is separate from marker discovery. It reads a fixed panel generated
by proportional_split_panel_model.py and performs panel-only UKB/Olink outcome
evaluation without any marker reselection or clinical modeling.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

BASE = Path("/data2/jiajm/Atherosclerosis_Manuscript/new_panel/final/RADAR_ModulesNet_142")
DEFAULT_PANEL_FILE = BASE / "ProportionalSplit_DevTrainExternal_FinalHoldout" / "selected_panel_by_devtrain_external.csv"
DEFAULT_UKB_FILE = Path("/data2/jiajm/Atherosclerosis_Manuscript/new_panel/incident_filtered_cohort_with_full_olink.csv")
DEFAULT_OUT_DIR = BASE / "ProportionalSplit_DevTrainExternal_FinalHoldout" / "UKB_Olink_PanelEvaluation"

REQUIRED_COLUMNS = ["EID_STD", "T0", "followup_days", "t_event_days"]
SPECIFICITY_TARGETS = [0.80, 0.85, 0.90]


def canonical(x: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def find_olink_col(df: pd.DataFrame, gene: str) -> str | None:
    candidates = [f"olink_{gene.lower()}", gene.lower(), gene.upper(), gene]
    for c in candidates:
        if c in df.columns:
            return c
    targets = {canonical(f"olink_{gene}"), canonical(gene)}
    for c in df.columns:
        if canonical(c) in targets:
            return c
    return None


def read_panel(panel_file: Path) -> pd.DataFrame:
    panel = pd.read_csv(panel_file)
    if "gene" not in panel.columns:
        raise ValueError(f"Panel file must contain a 'gene' column: {panel_file}")
    panel = panel.copy()
    panel["gene"] = panel["gene"].astype(str).str.strip()
    panel = panel[panel["gene"].ne("")].drop_duplicates("gene").reset_index(drop=True)
    if "panel_rank" not in panel.columns:
        panel["panel_rank"] = np.arange(1, len(panel) + 1)
    if len(panel) != 20:
        print(f"[WARN] Panel has {len(panel)} markers, not 20. Evaluation will use all listed markers.")
    return panel


def build_labels(person_df: pd.DataFrame, horizon_days: int, gap_days: int = 0) -> tuple[np.ndarray, np.ndarray]:
    t = pd.to_numeric(person_df["t_event_days"], errors="coerce")
    follow = pd.to_numeric(person_df["followup_days"], errors="coerce")
    has_event = t.notna()
    t2 = t.fillna(np.inf)
    y_all = ((has_event) & (t2 >= gap_days) & (t2 <= horizon_days)).astype(int).to_numpy()
    eligible = ((follow >= horizon_days) | (y_all == 1))
    eligible = pd.Series(eligible).fillna(False).to_numpy(dtype=bool)
    return y_all[eligible], eligible


def compute_basic_metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    y = np.asarray(y).astype(int)
    score = pd.to_numeric(pd.Series(score), errors="coerce").to_numpy()
    ok = np.isfinite(score)
    y, score = y[ok], score[ok]
    out = {"n": int(len(y)), "events": int(y.sum()), "base_rate": float(y.mean()) if len(y) else np.nan}
    if len(np.unique(y)) < 2:
        out.update({"AUC": np.nan, "AP": np.nan})
    else:
        out.update({"AUC": float(roc_auc_score(y, score)), "AP": float(average_precision_score(y, score))})
    return out


def fixed_specificity_metrics(y: np.ndarray, score: np.ndarray, specificity_target: float) -> dict[str, float]:
    y = np.asarray(y).astype(int)
    score = pd.to_numeric(pd.Series(score), errors="coerce").to_numpy()
    ok = np.isfinite(score)
    y, score = y[ok], score[ok]
    if len(np.unique(y)) < 2:
        return {
            "selected_n": np.nan,
            "tp": np.nan,
            "fp": np.nan,
            "sensitivity": np.nan,
            "specificity": np.nan,
            "ppv": np.nan,
            "nns": np.nan,
        }

    neg_scores = score[y == 0]
    threshold = np.quantile(neg_scores, specificity_target)
    pred = score >= threshold
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    sensitivity = tp / (tp + fn) if (tp + fn) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    nns = 1 / ppv if ppv and ppv > 0 else np.inf
    return {
        "threshold": float(threshold),
        "selected_n": int(pred.sum()),
        "tp": tp,
        "fp": fp,
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "ppv": float(ppv),
        "nns": float(nns),
    }


def evaluate_ukb(panel_file: Path, ukb_file: Path, out_dir: Path, horizon_days: int) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    panel = read_panel(panel_file)
    cohort = pd.read_csv(ukb_file, low_memory=False).reset_index(drop=True)

    panel_map = panel[["panel_rank", "gene"]].copy()
    if "module_use" in panel.columns:
        panel_map["module_use"] = panel["module_use"]
    panel_map["olink_column"] = panel_map["gene"].map(lambda g: find_olink_col(cohort, g))
    panel_map.to_csv(out_dir / "panel_olink_mapping.csv", index=False)

    missing_genes = panel_map.loc[panel_map["olink_column"].isna(), "gene"].tolist()
    missing_required = [c for c in REQUIRED_COLUMNS if c not in cohort.columns]
    if missing_genes or missing_required:
        status = {
            "ran": False,
            "missing_panel_genes": missing_genes,
            "missing_required_columns": missing_required,
        }
        (out_dir / "ukb_olink_status.json").write_text(json.dumps(status, indent=2))
        return status

    genes = panel_map["gene"].tolist()
    olink_cols = panel_map["olink_column"].tolist()
    protein = cohort[olink_cols].apply(pd.to_numeric, errors="coerce")
    protein.columns = genes

    mu = protein.mean(axis=0, skipna=True)
    sd = protein.std(axis=0, skipna=True, ddof=0).replace(0, 1.0).fillna(1.0)
    z = (protein - mu) / sd
    cohort["PanelScore_equal_weight"] = z.mean(axis=1, skipna=True)
    cohort["n_panel_used"] = z.notna().sum(axis=1)
    for gene in genes:
        cohort[f"z_{gene}"] = z[gene].to_numpy()

    pd.DataFrame({
        "gene": genes,
        "olink_column": olink_cols,
        "mean_used_for_zscore": mu[genes].to_numpy(),
        "sd_used_for_zscore": sd[genes].to_numpy(),
    }).to_csv(out_dir / "panel_ukb_zscore_parameters.csv", index=False)

    y, eligible = build_labels(cohort, horizon_days=horizon_days)
    cohort_eval = cohort.loc[eligible].copy().reset_index(drop=True)
    cohort_eval[f"AMI_{horizon_days}d"] = y.astype(int)

    y = cohort_eval[f"AMI_{horizon_days}d"].to_numpy(int)
    scores = {
        "Panel-only": cohort_eval["PanelScore_equal_weight"].to_numpy(float),
    }

    score_cols = [
        "EID_STD",
        "T0",
        f"AMI_{horizon_days}d",
        "age_recruit",
        "sex",
        "t_event_days",
        "followup_days",
        "PanelScore_equal_weight",
        "n_panel_used",
    ] + [f"z_{g}" for g in genes]
    score_cols = [c for c in score_cols if c in cohort_eval.columns]
    score_df = cohort_eval[score_cols].copy()
    for model_name, score in scores.items():
        score_df[f"score_{model_name}"] = score
    score_df.to_csv(out_dir / f"scores_{horizon_days}d_panel_ukb_olink.csv", index=False)

    subgroups = {"overall": np.ones(len(cohort_eval), dtype=bool)}
    if "sex" in cohort_eval.columns:
        sex = pd.to_numeric(cohort_eval["sex"], errors="coerce")
        subgroups["male"] = sex.eq(1).to_numpy()
        subgroups["female"] = sex.eq(0).to_numpy()
    if "age_recruit" in cohort_eval.columns:
        age = pd.to_numeric(cohort_eval["age_recruit"], errors="coerce")
        subgroups["age_ge_60"] = age.ge(60).to_numpy()
    if "sex" in cohort_eval.columns and "age_recruit" in cohort_eval.columns:
        subgroups["male_age_ge_60"] = (sex.eq(1) & age.ge(60)).to_numpy()

    auc_rows = []
    specificity_rows = []
    for group_name, mask in subgroups.items():
        yy = y[mask]
        for model_name, score_all in scores.items():
            score = score_all[mask]
            auc_rows.append({
                "group": group_name,
                "model": model_name,
                "horizon_days": horizon_days,
                "panel_n": len(genes),
                "panel_genes": ";".join(genes),
                **compute_basic_metrics(yy, score),
            })
            for specificity in SPECIFICITY_TARGETS:
                specificity_rows.append({
                    "group": group_name,
                    "model": model_name,
                    "horizon_days": horizon_days,
                    "specificity_target": specificity,
                    **fixed_specificity_metrics(yy, score, specificity),
                })

    pd.DataFrame(auc_rows).to_csv(out_dir / f"AUC_AP_{horizon_days}d_panel_ukb_olink.csv", index=False)
    pd.DataFrame(specificity_rows).to_csv(out_dir / f"specificity80_85_90_{horizon_days}d_panel_ukb_olink.csv", index=False)

    status = {
        "ran": True,
        "model": "Panel-only",
        "panel_file": str(panel_file),
        "ukb_file": str(ukb_file),
        "out_dir": str(out_dir),
        "horizon_days": horizon_days,
        "n_panel_genes": len(genes),
        "n_eligible": int(len(cohort_eval)),
        "events": int(y.sum()),
    }
    (out_dir / "ukb_olink_status.json").write_text(json.dumps(status, indent=2))
    return status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate a fixed 20-marker panel-only score in UKB Olink.")
    parser.add_argument("--panel-file", type=Path, default=DEFAULT_PANEL_FILE)
    parser.add_argument("--ukb-file", type=Path, default=DEFAULT_UKB_FILE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--horizon-days", type=int, default=365)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    status = evaluate_ukb(
        panel_file=args.panel_file,
        ukb_file=args.ukb_file,
        out_dir=args.out_dir,
        horizon_days=args.horizon_days,
    )
    print("[UKB/Olink status]", status)


if __name__ == "__main__":
    main()
