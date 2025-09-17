from __future__ import annotations
import os, datetime as dt
import pandas as pd
from .utils.io_utils import ensure_dir, write_json

REPORT_COLS = [
    "gene", "y_prob_max", "p95_flag", "bestF1_flag",
    "rank", "clinvar_plp_flag", "tier", "score", "disease"
]

def _coerce_columns(df: pd.DataFrame, disease_label: str) -> pd.DataFrame:
    out = df.copy()

    # score: prefer FINAL_score -> score
    if "score" not in out.columns:
        if "FINAL_score" in out.columns:
            out["score"] = out["FINAL_score"]
        else:
            out["score"] = 0.0

    # rank: prefer existing rank -> rank_final -> 0..N-1
    if "rank" not in out.columns:
        if "rank_final" in out.columns:
            out["rank"] = out["rank_final"]
        else:
            out["rank"] = range(len(out))

    # ensure mandatory columns exist with sensible defaults
    defaults = {
        "gene": "",
        "y_prob_max": 0.0,
        "p95_flag": 0,
        "bestF1_flag": 0,
        "clinvar_plp_flag": 0,
        "tier": "",
    }
    for k, v in defaults.items():
        if k not in out.columns:
            out[k] = v

    # disease label column
    out["disease"] = disease_label

    # reorder & subset
    missing = [c for c in REPORT_COLS if c not in out.columns]
    if missing:
        # all missing should be filled above; this is just safety
        for c in missing:
            out[c] = "" if c in ("gene", "tier", "disease") else 0
    out = out[REPORT_COLS]

    # stable sort for top view
    out = out.sort_values(["score", "y_prob_max"], ascending=[False, False], kind="mergesort")
    return out

def write_outputs(df: pd.DataFrame, outdir: str, disease: str, meta: dict, disease_label: str | None = None) -> dict:
    ensure_dir(outdir)
    slug = "".join(c if c.isalnum() else "_" for c in (disease_label or disease).lower()).strip("_") or "disease"
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(outdir, f"{slug}_{ts}")
    ensure_dir(run_dir)

    # curated report
    final_df = _coerce_columns(df, disease_label or disease)
    final_path = os.path.join(run_dir, f"final_for_report_{slug}.tsv")
    final_df.to_csv(final_path, sep="\t", index=False)

    # topN with the same schema
    top_n = int(meta.get("top_n", 1000))
    top_path = os.path.join(run_dir, f"top{top_n}_{slug}.tsv")
    final_df.head(top_n).to_csv(top_path, sep="\t", index=False)

    # meta
    meta_path  = os.path.join(run_dir, f"run_{slug}.json")
    write_json(meta_path, meta)

    return {"final": final_path, "top": top_path, "meta": meta_path, "run_dir": run_dir}
