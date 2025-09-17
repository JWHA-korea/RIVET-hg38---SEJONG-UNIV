from __future__ import annotations
from typing import Union, Dict
import numpy as np
import pandas as pd

GENE_COL  = "gene"
PROBA_COL = "y_prob_max"

def _norm01(x: pd.Series) -> pd.Series:
    s = pd.Series(x, dtype=float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    mn, mx = float(s.min()), float(s.max())
    return (s - mn) / (mx - mn) if mx > mn else s*0.0

def _winsor_top(x: pd.Series, q: float = 0.99) -> pd.Series:
    s = pd.Series(x, dtype=float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    cap = float(s.quantile(q))
    return s.clip(upper=cap)

def _hpo_df_from_any(hpo_w: Union[None, Dict[str, float], pd.DataFrame, pd.Series]) -> pd.DataFrame:
    if hpo_w is None:
        return pd.DataFrame(columns=[GENE_COL, "hpo_weight"])
    if isinstance(hpo_w, dict):
        df = pd.DataFrame({GENE_COL: list(hpo_w.keys()),
                           "hpo_weight": [float(hpo_w[g]) for g in hpo_w]})
    elif isinstance(hpo_w, pd.Series):
        df = pd.DataFrame({GENE_COL: [str(k) for k in hpo_w.index],
                           "hpo_weight": [float(v) if pd.notna(v) else 0.0 for v in hpo_w.values]})
    elif isinstance(hpo_w, pd.DataFrame):
        df = hpo_w.copy()
        if GENE_COL not in df.columns:
            df = df.rename(columns={df.columns[0]: GENE_COL})
        wc = next((c for c in ["hpo_weight", "HPO", "weight"] if c in df.columns), None)
        if wc is None:
            wc = df.columns[1] if len(df.columns) >= 2 else None
        if wc is None or GENE_COL not in df.columns:
            return pd.DataFrame(columns=[GENE_COL, "hpo_weight"])
        df = df[[GENE_COL, wc]].rename(columns={wc: "hpo_weight"})
    else:
        return pd.DataFrame(columns=[GENE_COL, "hpo_weight"])
    df[GENE_COL] = df[GENE_COL].astype(str).str.upper()
    df["hpo_weight"] = pd.to_numeric(df["hpo_weight"], errors="coerce").fillna(0.0).clip(0, 1)
    return df

def combine(df_gene: pd.DataFrame,
            hpo_w: Union[None, Dict[str, float], pd.DataFrame, pd.Series],
            weights: Dict[str, float] | None,
            winsor_q: float = 0.99,
            score_lo: float = 0.10,
            score_hi: float = 0.99,
            jitter: float = 1e-6) -> pd.DataFrame:
    """
    Rerank fusion (rerank_auto.py와 동일한 스케일 전략):
        - FUNC winsorize(q) → minmax, NET/PATH/NOVEL/HPO minmax
        - FINAL_raw = Σ w_i * channel_i
        - score = lo + (hi-lo) * minmax(FINAL_raw) + jitter
    """
    df = df_gene.copy()

    for col in ["FUNC", "NET", "PATH", "NOVEL"]:
        if col not in df.columns:
            df[col] = 0.0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    H = _hpo_df_from_any(hpo_w)
    if not H.empty:
        df[GENE_COL] = df[GENE_COL].astype(str).str.upper()
        df = df.merge(H, on=GENE_COL, how="left")
    if "hpo_weight" not in df.columns:
        df["hpo_weight"] = 0.0
    df["hpo_weight"] = pd.to_numeric(df["hpo_weight"], errors="coerce").fillna(0.0).clip(0, 1)

    FUNCn  = _norm01(_winsor_top(df["FUNC"], q=winsor_q))
    NETn   = _norm01(df["NET"])
    PATHn  = _norm01(df["PATH"])
    NOVELn = _norm01(df["NOVEL"])
    HPOn   = _norm01(df["hpo_weight"])

    w = {"FUNC": 0.35, "NET": 0.28, "PATH": 0.18, "NOVEL": 0.12, "HPO": 0.07}
    if isinstance(weights, dict):
        w.update({k.upper(): float(v) for k, v in weights.items()})

    FINAL_raw = (w["FUNC"]*FUNCn + w["NET"]*NETn + w["PATH"]*PATHn +
                 w["NOVEL"]*NOVELn + w["HPO"]*HPOn)

    score01 = _norm01(FINAL_raw)
    if jitter and jitter > 0:
        rng = np.random.default_rng(42)
        score01 = score01 + rng.normal(0.0, jitter, size=len(score01))
    score01 = score01.clip(0, 1)
    score = score_lo + (score_hi - score_lo) * score01

    df["FINAL_score"] = FINAL_raw
    df["score"] = score

    if PROBA_COL not in df.columns:
        df[PROBA_COL] = df["FUNC"]

    df = df.sort_values(["score", PROBA_COL], ascending=[False, False]).reset_index(drop=True)
    return df

