#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RIVET Gene Study CLI (HPO mapping + rerank-style fusion)

- 최소 실행 인자: --disease, --outdir (원하면 --disease-label)
- extras(STRING 네트워크/경로/문헌)는 config/paths.yaml 의 extras 블록이 존재하면 자동 병합
- 점수 스케일은 rerank_auto.py 와 유사( winsorize→minmax, 채널 가중합, 0.10~0.99 스케일 )
- 'disease' 컬럼은 **HPO 매핑된 유전자에만** 채워짐(HPO-only)

Author: JIWOOK HA (Sejong Univ. Academic Festival)
"""
from __future__ import annotations

import os, sys, json, argparse, datetime, re
from typing import Dict, Any, Set

import pandas as pd

from .auto_discovery import discover_paths, discover_defaults
from .hpo_resolver import disease_to_hpo_ids
from .hpo_mapper import map_hpo_to_gene_weights
from .tiers import ensure_flags_and_tier
from .combine_scores import combine, _hpo_df_from_any
from .extras import build_net_scores, load_path_scores, load_novel_scores

GENE_COL  = "gene"
PROBA_COL = "y_prob_max"

REPORT_COLS = [
    "gene", "y_prob_max", "p95_flag", "bestF1_flag", "rank",
    "clinvar_plp_flag", "tier", "score", "disease"
]

# ---------------------------- helpers ----------------------------

def _slugify(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "run"

def _ts() -> str:
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

def _parse_weights(s: str | None, defaults: Dict[str, Any]) -> Dict[str, float]:
    base = {k.upper(): float(v) for k, v in (defaults.get("weights") or {}).items()}
    if not s:
        return base
    out = dict(base)
    for tok in s.split(","):
        tok = tok.strip()
        if not tok or ":" not in tok:
            continue
        k, v = tok.split(":", 1)
        try:
            out[k.strip().upper()] = float(v)
        except Exception:
            pass
    return out

def _ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def _to_num(val, default=0.0, clip01=False):
    """Robust numeric cast for Series/scalar/list with optional [0,1] clip."""
    if isinstance(val, pd.Series):
        out = pd.to_numeric(val, errors="coerce").fillna(default)
        return out.clip(0, 1) if clip01 else out
    if isinstance(val, (list, tuple)):
        s = pd.Series(val)
        out = pd.to_numeric(s, errors="coerce").fillna(default)
        return out.clip(0, 1) if clip01 else out
    try:
        num = pd.to_numeric(val, errors="coerce")
    except Exception:
        num = None
    if pd.isna(num):
        num = default
    x = float(num)
    if clip01:
        if x < 0.0: x = 0.0
        if x > 1.0: x = 1.0
    return x

def _seed_genes_from_hpo(hpo_w) -> Set[str]:
    df = _hpo_df_from_any(hpo_w)
    if df.empty:
        return set()
    s = _to_num(df["hpo_weight"], default=0.0, clip01=True)
    return set(df.loc[s > 0.0, GENE_COL].astype(str).str.upper().tolist())

def _load_seed_file(path: str) -> Set[str]:
    s: Set[str] = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g:
                s.add(g.upper())
    return s

# ---------------------------- main ----------------------------

def main():
    ap = argparse.ArgumentParser(
        prog="rivet-gs",
        description="RIVET Gene Study CLI (HPO mapping → rerank-style fusion)"
    )
    # 필수
    ap.add_argument("--disease", required=True, help="Exact HPO-registered disease name (phenotype.hpoa col2)")
    ap.add_argument("--outdir", required=True, help="Output directory (a timestamped subdir will be created)")
    # 선택(표기)
    ap.add_argument("--disease-label", default=None, help="Label for the report (default: same as --disease)")
    # 선택(가중치/출력)
    ap.add_argument("--weights", default=None, help='Weights "FUNC:0.35,NET:0.28,PATH:0.18,NOVEL:0.12,HPO:0.07"')
    ap.add_argument("--top", type=int, default=None, help="Top-N (default: defaults.yaml)")
    ap.add_argument("--min-score", type=float, default=None, help="Filter by FINAL_score (raw) after combine")
    # 선택(extras: config/paths.yaml에 있으면 자동 사용, CLI는 override 용)
    ap.add_argument("--string-net", default=None, help="STRING-derived network (geneA\\tgeneB\\tweight)")
    ap.add_argument("--path-scores", default=None, help="Pathway summary scores per gene")
    ap.add_argument("--literature", default=None, help="Literature counts/PMIDs table")
    ap.add_argument("--seeds", default=None, help="Manual seed gene list (one gene per line, overrides HPO seeds)")
    # 선택(rerank knobs)
    ap.add_argument("--winsor-q", type=float, default=0.99)
    ap.add_argument("--gamma",     type=float, default=0.60, help="Personalized PageRank gamma")
    ap.add_argument("--novel-log", action="store_true",  default=True,  help="Use log1p on literature counts")
    ap.add_argument("--no-novel-log", dest="novel_log", action="store_false")
    ap.add_argument("--score-lo",  type=float, default=0.10)
    ap.add_argument("--score-hi",  type=float, default=0.99)
    ap.add_argument("--jitter",    type=float, default=1e-6)

    args = ap.parse_args()

    paths = discover_paths()
    defaults = discover_defaults() or {}
    disease_label = args.disease_label or args.disease

    gene_scores = paths.get("gene_scores")
    thresholds  = paths.get("thresholds")
    hpo_genes   = paths.get("hpo_genes")
    pheno_hpoa  = paths.get("phenotype_hpoa")

    if not (gene_scores and os.path.exists(gene_scores)):
        print("[ERROR] Gene score table not found. Set in config/paths.yaml or RIVET_GENE_SCORES.", file=sys.stderr); sys.exit(1)
    if not (thresholds and os.path.exists(thresholds)):
        print("[ERROR] thresholds.json not found. Set in config/paths.yaml or RIVET_THRESHOLDS.", file=sys.stderr); sys.exit(1)
    if not (hpo_genes and os.path.exists(hpo_genes)):
        print("[ERROR] genes_to_phenotype.txt not found. Set in config/paths.yaml or RIVET_HPO_GENES.", file=sys.stderr); sys.exit(1)
    if not (pheno_hpoa and os.path.exists(pheno_hpoa)):
        print("[ERROR] phenotype.hpoa not found. Set in config/paths.yaml or RIVET_PHENOTYPE_HPOA.", file=sys.stderr); sys.exit(1)

    print(f"[INFO] gene_scores    = {gene_scores}")
    print(f"[INFO] thresholds     = {thresholds}")
    print(f"[INFO] hpo_genes      = {hpo_genes}")
    print(f"[INFO] phenotype_hpoa = {pheno_hpoa}")

    weights = _parse_weights(args.weights, defaults)
    top_n   = int(args.top or (defaults.get("top_n") or 1000))

    # 1) HPO 매핑
    hpo_ids     = disease_to_hpo_ids(args.disease, pheno_hpoa)  # exact match, 없으면 raise
    hpo_weights = map_hpo_to_gene_weights(hpo_genes, hpo_ids)

    # 2) gene base table
    df_gene = pd.read_csv(gene_scores, sep="\t")

    # FUNC fallback (없으면 y_prob_max 사용)
    if "FUNC" not in df_gene.columns and "y_prob_max" in df_gene.columns:
        df_gene["FUNC"] = _to_num(df_gene["y_prob_max"], default=0.0, clip01=True)

    # flag/tier 정리
    df_gene = ensure_flags_and_tier(df_gene, thresholds)

    # 3) EXTRAS: config/paths.yaml → extras 블록 자동 탐지, CLI가 있으면 override
    cfg_ex = (paths.get("extras") or {})
    net_path  = args.string_net or cfg_ex.get("string_net")
    path_path = args.path_scores or cfg_ex.get("path_scores")
    lit_path  = args.literature  or cfg_ex.get("literature")

    # seeds: HPO → seed set, --seeds 파일이 있으면 override
    seed_genes = _seed_genes_from_hpo(hpo_weights)
    if args.seeds:
        if os.path.exists(args.seeds):
            seed_genes = _load_seed_file(args.seeds)
        else:
            print(f"[WARN] --seeds not found: {args.seeds} (fallback to HPO seeds)")

    # NET
    net_df = build_net_scores(net_path, seed_genes, gamma=args.gamma) if net_path else pd.DataFrame(columns=[GENE_COL, "NET_ex"])
    if not net_df.empty:
        df_gene[GENE_COL] = df_gene[GENE_COL].astype(str).str.upper()
        df_gene = df_gene.merge(net_df, on=GENE_COL, how="left")
        df_gene["NET"] = df_gene["NET"] if "NET" in df_gene.columns else df_gene["NET_ex"]
        df_gene["NET"] = _to_num(df_gene["NET"], default=0.0, clip01=True)
        df_gene.drop(columns=["NET_ex"], inplace=True, errors="ignore")
    else:
        if "NET" not in df_gene.columns: df_gene["NET"] = 0.0
        df_gene["NET"] = _to_num(df_gene["NET"], default=0.0, clip01=True)

    # PATH
    path_df = load_path_scores(path_path)
    if path_df is not None and not path_df.empty:
        df_gene = df_gene.merge(path_df, on=GENE_COL, how="left")
    if "PATH" not in df_gene.columns: df_gene["PATH"] = 0.0
    df_gene["PATH"] = _to_num(df_gene["PATH"], default=0.0, clip01=True)

    # NOVEL
    nov_df = load_novel_scores(lit_path, log_scale=args.novel_log, invert=True)
    if nov_df is not None and not nov_df.empty:
        df_gene = df_gene.merge(nov_df, on=GENE_COL, how="left")
        df_gene["NOVEL"] = df_gene["NOVEL"] if "NOVEL" in df_gene.columns else df_gene.get("NOVEL_ex", 0.0)
        if "NOVEL" not in df_gene.columns: df_gene["NOVEL"] = 0.0
        df_gene["NOVEL"] = _to_num(df_gene["NOVEL"], default=0.0, clip01=True)
        df_gene.drop(columns=["NOVEL_ex"], inplace=True, errors="ignore")
    else:
        if "NOVEL" not in df_gene.columns: df_gene["NOVEL"] = 0.0
        df_gene["NOVEL"] = _to_num(df_gene["NOVEL"], default=0.0, clip01=True)

    # 4) combine → FINAL_score & score (rerank scale)
    df_final = combine(
        df_gene, hpo_weights, weights,
        winsor_q=args.winsor_q,
        score_lo=args.score_lo, score_hi=args.score_hi, jitter=args.jitter
    )

    # optional filter by FINAL_score (raw)
    if args.min_score is not None:
        df_final = df_final[df_final["FINAL_score"] >= float(args.min_score)].copy()

    # 5) HPO-only disease labeling (라벨 모드 없음)
    try:
        df_hw = _hpo_df_from_any(hpo_weights)
        df_hw["gene"] = df_hw["gene"].astype(str).str.upper()
        df_hw["hpo_weight"] = pd.to_numeric(df_hw["hpo_weight"], errors="coerce").fillna(0.0)
        hpo_set = set(df_hw.loc[df_hw["hpo_weight"] > 0.0, "gene"].tolist())
    except Exception:
        hpo_set = set()

    genes_up = df_final["gene"].astype(str).str.upper()
    df_final["disease"] = ""  # 기본은 공란
    df_final.loc[genes_up.isin(hpo_set), "disease"] = disease_label  # HPO에만 라벨

    # Ensure report columns exist
    for c in ["y_prob_max", "p95_flag", "bestF1_flag", "rank", "clinvar_plp_flag", "tier", "score", "disease"]:
        if c not in df_final.columns:
            df_final[c] = 0 if c.endswith("_flag") else ""

    # minimal report view (정렬: score desc, y_prob_max desc)
    df_final_sorted = df_final.sort_values(["score", PROBA_COL], ascending=[False, False]).reset_index(drop=True)
    cols = [c for c in REPORT_COLS if c in df_final_sorted.columns]
    report_df = df_final_sorted[cols].copy()

    # output
    slug   = _slugify(disease_label)
    outdir = os.path.join(args.outdir, f"{slug}_{_ts()}")
    _ensure_dir(outdir)

    final_path   = os.path.join(outdir, f"final_for_report_{slug}.tsv")
    top_n        = int(args.top or (defaults.get("top_n") or 1000))
    top_path_min = os.path.join(outdir, f"top{top_n}_{slug}.tsv")
    top_path_all = os.path.join(outdir, f"top{top_n}_full_{slug}.tsv")

    report_df.to_csv(final_path, sep="\t", index=False)
    report_df.head(top_n).to_csv(top_path_min, sep="\t", index=False)
    df_final_sorted.head(top_n).to_csv(top_path_all, sep="\t", index=False)

    run_info = {
        "disease": args.disease,
        "disease_label": disease_label,
        "paths": {
            "gene_scores": gene_scores,
            "thresholds": thresholds,
            "hpo_genes": hpo_genes,
            "phenotype_hpoa": pheno_hpoa,
            "extras": {
                "string_net": net_path,
                "path_scores": path_path,
                "literature":  lit_path
            }
        },
        "weights": _parse_weights(args.weights, defaults),
        "rerank": {
            "winsor_q": args.winsor_q,
            "gamma": args.gamma,
            "novel_log": args.novel_log,
            "score_lo": args.score_lo,
            "score_hi": args.score_hi,
            "jitter": args.jitter
        },
        "top_n": top_n,
        "rows_out": int(len(df_final_sorted)),
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds")
    }
    with open(os.path.join(outdir, "run.json"), "w", encoding="utf-8") as f:
        json.dump(run_info, f, ensure_ascii=False, indent=2)

    print(f"[OK] final (report): {final_path}")
    print(f"[OK] topN (report) : {top_path_min}")
    print(f"[OK] topN (full)   : {top_path_all}")
    print(f"[OK] meta          : {os.path.join(outdir, 'run.json')}")

if __name__ == "__main__":
    main()

