# rivet_gs/extras.py
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Dict, Set, Tuple, Optional

import numpy as np
import pandas as pd


# ---------------------------
# PATH scores loader
# ---------------------------
def load_path_scores(path: str | None) -> Optional[pd.DataFrame]:
    """
    Load pathway summary scores per gene.

    Expected format (TSV, header):
        gene    PATH
    - gene: gene symbol (string)
    - PATH: float in [0,1] (will be coerced & clipped)

    Returns a DataFrame with columns ["gene", "PATH"], or None if path is missing.
    """
    if not path or not Path(path).exists():
        return None
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "gene" not in df.columns:
        # try first col as gene if header is missing
        df.columns = [c if i else "gene" for i, c in enumerate(df.columns)]
    df["gene"] = df["gene"].astype(str).str.upper()

    # find PATH column (exact or case-insensitive)
    path_col = None
    for c in df.columns:
        if c.strip().lower() == "path":
            path_col = c
            break
    if path_col is None:
        # try second column as PATH
        if df.shape[1] >= 2:
            path_col = df.columns[1]
        else:
            return None

    path_vals = pd.to_numeric(df[path_col], errors="coerce").fillna(0.0)
    path_vals = path_vals.clip(0.0, 1.0)
    out = pd.DataFrame({"gene": df["gene"], "PATH": path_vals})
    return out


# ---------------------------
# NOVEL (literature) loader
# ---------------------------
def load_novel_scores(
    path: str | None,
    log_scale: bool = True,
    invert: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Load literature-based novelty and convert to a 0..1 score (higher = more novel).

    Expected input (TSV; compression inferred):
      - Must include 'gene' column.
      - Either:
          * 'count' column with integer counts, or
          * 'pmids' column with comma/space-separated PMIDs (will be counted).

    Steps:
      1) Uppercase genes
      2) Count = numeric or len(pmids)
      3) Optional log1p(count)
      4) Invert (max - x) so more papers -> lower novelty
      5) Min-max to [0,1] as NOVEL_ex
    """
    if not path or not Path(path).exists():
        return None

    df = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    if "gene" not in df.columns:
        return None
    df["gene"] = df["gene"].astype(str).str.upper()

    if "count" in df.columns:
        cnt = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    elif "pmids" in df.columns:
        # split on comma/space; count tokens
        def _pmid_len(s: str) -> int:
            s = (s or "").strip()
            if not s:
                return 0
            toks = [x for x in re.split(r"[,\s]+", s) if x]
            return len(toks)

        cnt = df["pmids"].fillna("").map(_pmid_len).astype(float)
    else:
        # try a common alias
        for alias in ["n", "n_mentions", "papers", "hits"]:
            if alias in df.columns:
                cnt = pd.to_numeric(df[alias], errors="coerce").fillna(0)
                break
        else:
            return None

    cnt = cnt.astype(float)

    if log_scale:
        cnt = np.log1p(cnt)

    if invert:
        cnt = (cnt.max() - cnt)

    # normalize to [0,1]
    cmin, cmax = float(cnt.min()), float(cnt.max())
    if cmax > cmin:
        nov = (cnt - cmin) / (cmax - cmin)
    else:
        nov = cnt * 0.0

    return pd.DataFrame({"gene": df["gene"], "NOVEL_ex": nov.astype(float)})


# ---------------------------
# NET (STRING-derived) builder
# ---------------------------
def build_net_scores(
    net_path: str | None,
    seed_genes: Set[str],
    gamma: float = 0.60,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> pd.DataFrame:
    """
    Personalized PageRank over a simple gene-gene weighted network.

    Input:
      - net_path: TSV with columns (geneA, geneB, weight) in [0,1]; header required.
      - seed_genes: set of uppercase gene symbols used for personalization.
      - gamma: damping on graph (typical 0.60 ~ 0.85)
      - max_iter, tol: power-iteration params

    Output:
      DataFrame ["gene","NET_ex"] with 0..1 min-max scaled scores.
      If net_path is missing/empty or seeds empty → returns empty DataFrame.
    """
    if not net_path or not Path(net_path).exists():
        return pd.DataFrame(columns=["gene", "NET_ex"])

    edges = pd.read_csv(net_path, sep="\t", dtype=str)
    # sanitize/rename
    cols_lower = [c.lower() for c in edges.columns]
    col_map = {}
    for c in edges.columns:
        cl = c.lower()
        if cl in ("genea", "source", "a", "node1", "protein1"):
            col_map[c] = "geneA"
        elif cl in ("geneb", "target", "b", "node2", "protein2"):
            col_map[c] = "geneB"
        elif cl in ("weight", "w", "score", "combined_score", "s"):
            col_map[c] = "weight"
    edges = edges.rename(columns=col_map)

    # require columns
    for need in ["geneA", "geneB"]:
        if need not in edges.columns:
            return pd.DataFrame(columns=["gene", "NET_ex"])

    # weight handling
    if "weight" not in edges.columns:
        edges["weight"] = 1.0
    edges["weight"] = pd.to_numeric(edges["weight"], errors="coerce").fillna(0.0)
    # clip to [0,1]
    edges["weight"] = edges["weight"].clip(0.0, 1.0)

    # uppercase nodes
    edges["geneA"] = edges["geneA"].astype(str).str.upper()
    edges["geneB"] = edges["geneB"].astype(str).str.upper()

    # build node index
    nodes = pd.unique(pd.concat([edges["geneA"], edges["geneB"]], ignore_index=True))
    if len(nodes) == 0 or len(seed_genes) == 0:
        return pd.DataFrame(columns=["gene", "NET_ex"])

    idx: Dict[str, int] = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)

    # adjacency (row-normalized transition)
    # out-degree by A -> B edges
    out_w = np.zeros(n, dtype=float)
    # using adjacency as COO triplets: i->j
    rows, cols, vals = [], [], []
    for a, b, w in edges[["geneA", "geneB", "weight"]].itertuples(index=False):
        ia = idx.get(a)
        ib = idx.get(b)
        if ia is None or ib is None:
            continue
        if w <= 0.0:
            continue
        rows.append(ia)
        cols.append(ib)
        vals.append(float(w))
        out_w[ia] += float(w)
        # treat as undirected if not already present (STRING is usually undirected-like)
        rows.append(ib)
        cols.append(ia)
        vals.append(float(w))
        out_w[ib] += float(w)

    if len(vals) == 0:
        return pd.DataFrame(columns=["gene", "NET_ex"])

    # build row-normalized matrix in dict-of-lists form to avoid heavy deps
    # T[j, i] means prob from i->j when multiplying r_next = gamma*T @ r + (1-gamma)*p
    # We'll store forward i->j edges and use transpose multiply manually.
    neighbors: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    # First, accumulate lists per source i
    tmp: Dict[int, list] = {}
    for i, j, w in zip(rows, cols, vals):
        tmp.setdefault(i, []).append((j, w))
    for i, lst in tmp.items():
        if out_w[i] > 0:
            js = np.array([j for j, _ in lst], dtype=int)
            ws = np.array([w for _, w in lst], dtype=float)
            ws = ws / out_w[i]
            neighbors[i] = (js, ws)

    # personalization vector p
    p = np.zeros(n, dtype=float)
    seed_idx = [idx[g] for g in seed_genes if g in idx]
    if len(seed_idx) == 0:
        # if HPO seeds not in graph, bail out
        return pd.DataFrame(columns=["gene", "NET_ex"])
    p[seed_idx] = 1.0
    p = p / p.sum()

    # init
    r = p.copy()

    # power iteration: r_{t+1} = gamma * W^T r_t + (1-gamma) * p
    # W is row-normalized (sum_j W[i,j] = 1), so W^T multiply:
    for _ in range(max_iter):
        r_next = np.zeros_like(r)
        # for each source i, distribute r[i] to neighbors j
        for i, (js, ws) in neighbors.items():
            ri = r[i]
            if ri == 0:
                continue
            r_next[js] += gamma * ri * ws
        r_next += (1.0 - gamma) * p

        # check convergence (L1)
        if np.abs(r_next - r).sum() < tol:
            r = r_next
            break
        r = r_next

    # min-max to 0..1
    rmin, rmax = float(r.min()), float(r.max())
    if rmax > rmin:
        s = (r - rmin) / (rmax - rmin)
    else:
        s = r * 0.0

    out = pd.DataFrame({
        "gene": nodes.astype(str),
        "NET_ex": s.astype(float),
    })
    return out

