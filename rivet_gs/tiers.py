from __future__ import annotations
import json
import pandas as pd

P95_FLAG = "p95_flag"
BF1_FLAG = "bestF1_flag"
TIER_COL = "tier"
PROBA_COL = "y_prob_max"
P95_THR_COL = "P95_thr"
BF1_THR_COL = "BestF1_thr"

def ensure_flags_and_tier(df: pd.DataFrame, thresholds_json: str | None) -> pd.DataFrame:
    has_p95 = P95_FLAG in df.columns
    has_bf1 = BF1_FLAG in df.columns
    has_tier = TIER_COL in df.columns

    thr_p95 = thr_bf1 = None
    if thresholds_json:
        try:
            with open(thresholds_json, "r") as f:
                data = json.load(f)
            for k in ["P95_threshold", "P95_thr", "p95", "P95"]:
                if k in data:
                    thr_p95 = float(data[k]); break
            for k in ["BestF1_threshold", "BestF1_thr", "best_f1", "F1"]:
                if k in data:
                    thr_bf1 = float(data[k]); break
        except Exception:
            pass

    if not has_p95:
        if P95_THR_COL in df.columns:
            df[P95_FLAG] = (df[PROBA_COL] >= df[P95_THR_COL]).astype(int)
        elif thr_p95 is not None and PROBA_COL in df.columns:
            df[P95_FLAG] = (df[PROBA_COL] >= thr_p95).astype(int)

    if not has_bf1:
        if BF1_THR_COL in df.columns:
            df[BF1_FLAG] = (df[PROBA_COL] >= df[BF1_THR_COL]).astype(int)
        elif thr_bf1 is not None and PROBA_COL in df.columns:
            df[BF1_FLAG] = (df[PROBA_COL] >= thr_bf1).astype(int)

    if not has_tier:
        def _tier_row(r):
            if r.get(P95_FLAG, 0) == 1: return "T1"
            if r.get(BF1_FLAG, 0) == 1: return "T2"
            return ""
        df[TIER_COL] = df.apply(_tier_row, axis=1)

    return df
