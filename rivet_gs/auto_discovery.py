from __future__ import annotations
import os, yaml
from typing import Dict, Optional

ENV = {
    "GENE_SCORES": "RIVET_GENE_SCORES",
    "THRESHOLDS":  "RIVET_THRESHOLDS",
    "HPO_OBO":     "RIVET_HPO_OBO",
    "HPO_GENES":   "RIVET_HPO_GENES",
    "PHENO_HPOA":  "RIVET_PHENOTYPE_HPOA",
    "VAR_SCORES":  "RIVET_VARIANT_SCORES",
}

def _find_repo_root(start: str) -> str:
    cur = os.path.abspath(start)
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "config", "paths.yaml")):
            return cur
        nxt = os.path.dirname(cur)
        if nxt == cur:
            break
        cur = nxt
    return os.path.abspath(os.path.join(start, ".."))

_THIS_DIR = os.path.dirname(__file__)
_REPO_ROOT = _find_repo_root(os.path.join(_THIS_DIR, ".."))

_DEF = {
    "paths": os.path.join(_REPO_ROOT, "config", "paths.yaml"),
    "defaults": os.path.join(_REPO_ROOT, "config", "defaults.yaml"),
}

def _read_yaml(path: str) -> Dict:
    try:
        if not os.path.exists(path):
            return {}
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return data or {}
    except Exception:
        return {}

def _norm_path(p: Optional[str], base_dir: str) -> Optional[str]:
    if not p:
        return None
    p = os.path.expandvars(os.path.expanduser(str(p)))
    if not os.path.isabs(p):
        p = os.path.abspath(os.path.join(base_dir, p))
    return p

def _pick_env(env_key: str) -> Optional[str]:
    var = ENV.get(env_key)
    if not var: return None
    val = os.environ.get(var, "").strip()
    if not val: return None
    val = os.path.expandvars(os.path.expanduser(val))
    return val if os.path.exists(val) else None

def discover_paths() -> Dict[str, Optional[str]]:
    paths_yaml = _DEF["paths"]
    cfg = _read_yaml(paths_yaml)

    # ⬇️ 핵심: 상대경로를 항상 "repo 루트" 기준으로 정규화
    base = _REPO_ROOT

    def pick(key: str, env_key: str) -> Optional[str]:
        ev = _pick_env(env_key)
        if ev:
            return ev
        val = _norm_path(cfg.get(key), base) if isinstance(cfg, dict) else None
        if val and os.path.exists(val):
            return val
        return val

    return {
        "gene_scores":    pick("gene_scores",    "GENE_SCORES"),
        "thresholds":     pick("thresholds",     "THRESHOLDS"),
        "hpo_obo":        pick("hpo_obo",        "HPO_OBO"),
        "hpo_genes":      pick("hpo_genes",      "HPO_GENES"),
        "phenotype_hpoa": pick("phenotype_hpoa", "PHENO_HPOA"),
        "variant_scores": pick("variant_scores", "VAR_SCORES"),
    }

def discover_defaults() -> Dict:
    return _read_yaml(_DEF["defaults"]) or {}
