from __future__ import annotations
from typing import Set, List, Tuple

# phenotype.hpoa columns (tab-separated):
# 0: database_id, 1: disease_name, 2: qualifier, 3: hpo_id, ...

def _norm(s: str) -> str:
    # 소문자 + 공백 정규화 + 기호 제거
    return "".join(ch for ch in s.lower().strip() if ch.isalnum() or ch.isspace())

def _scan_pheno(path: str):
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            dname = parts[1]   # disease_name
            hpo_id = parts[3]  # hpo_id
            if not (hpo_id and hpo_id.startswith('HP:')):
                continue
            yield dname, hpo_id

def disease_to_hpo_ids(disease_name: str, phenotype_hpoa_path: str) -> Set[str]:
    """정확 일치 -> 정규화 일치 -> '부분 포함' 단일 후보 -> 다중 후보면 제시."""
    target_norm = _norm(disease_name)
    exact: Set[str] = set()
    fuzzy: Set[str] = set()
    contains: dict[str, Set[str]] = {}  # disease_name -> {HP:...}

    for dname, hpo_id in _scan_pheno(phenotype_hpoa_path):
        if dname == disease_name:
            exact.add(hpo_id)
        if _norm(dname) == target_norm:
            fuzzy.add(hpo_id)
        if target_norm and target_norm in _norm(dname):
            contains.setdefault(dname, set()).add(hpo_id)

    if exact:
        return exact
    if fuzzy:
        return fuzzy
    if contains:
        # 부분일치가 여러 개면 사용자에게 후보 제시
        if len(contains) == 1:
            return next(iter(contains.values()))
        suggestions = ", ".join(list(contains.keys())[:10])
        raise ValueError(
            f"Ambiguous disease name '{disease_name}'. Did you mean one of: {suggestions} ..."
        )
    raise ValueError(
        f"Disease name '{disease_name}' not found in phenotype.hpoa. "
        f"Try a canonical label from HPO (see docs/HPO_USAGE.md)."
    )

def list_disease_names(phenotype_hpoa_path: str, limit: int | None = None) -> List[str]:
    names: List[str] = []
    seen = set()
    for dname, _ in _scan_pheno(phenotype_hpoa_path):
        if dname not in seen:
            seen.add(dname)
            names.append(dname)
            if limit and len(names) >= limit:
                break
    return names
