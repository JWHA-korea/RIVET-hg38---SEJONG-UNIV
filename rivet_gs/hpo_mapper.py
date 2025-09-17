from __future__ import annotations
from typing import Set
import pandas as pd

GENE_COL = "gene"

def map_hpo_to_gene_weights(genes_to_phenotype_path: str, hpo_ids: Set[str]) -> pd.DataFrame:
    """
    genes_to_phenotype.txt에서 주어진 HPO ID 집합과 매칭되는 유전자에 가중치를 부여해 반환.
    반환 컬럼: gene, hpo_weight (0..1 스케일, 최소 양수 0.2 보장)
    """
    # 1) 우선 헤더 유무 판단
    df = pd.read_csv(genes_to_phenotype_path, sep='\t', header=None, dtype=str, engine='c')
    first = df.iloc[0].astype(str).str.lower().tolist()
    header_like = any('hpo' in x or 'gene' in x for x in first)

    if header_like:
        df = pd.read_csv(genes_to_phenotype_path, sep='\t', header=0, dtype=str, engine='c')
        gene_cands = [c for c in df.columns if str(c).lower() in ('gene', 'gene_symbol', 'entrez_gene_symbol')]
        hpo_cands  = [c for c in df.columns if 'hpo' in str(c).lower() and 'id' in str(c).lower()]
        if gene_cands and hpo_cands:
            gcol, hcol = gene_cands[0], hpo_cands[0]
            df = df[[gcol, hcol]].copy()
            df.columns = ['gene_symbol', 'hpo_id']
        else:
            # 헤더가 애매하면 인덱스로 폴백
            df = pd.read_csv(genes_to_phenotype_path, sep='\t', header=None, dtype=str, engine='c')
            df = df[[1, 2]].copy()
            df.columns = ['gene_symbol', 'hpo_id']
    else:
        # 헤더 없음 가정: 1=gene_symbol, 2=hpo_id
        df = df[[1, 2]].copy()
        df.columns = ['gene_symbol', 'hpo_id']

    # 2) 주어진 HPO ID들만 필터
    df = df[df['hpo_id'].isin(hpo_ids)].copy()
    if df.empty:
        return pd.DataFrame({GENE_COL: [], 'hpo_weight': []})

    # 3) gene별 매칭 횟수 -> 0..1 스케일
    g_counts = df.groupby('gene_symbol').size().rename('hpo_hits').reset_index()
    mn = g_counts['hpo_hits'].min()
    mx = g_counts['hpo_hits'].max()
    if mx == mn:
        g_counts['hpo_weight'] = 1.0
    else:
        g_counts['hpo_weight'] = (g_counts['hpo_hits'] - mn) / (mx - mn)
        # 0으로 눌리는 항목에 최소 양수 가중치 부여(탐색성 확보)
        g_counts.loc[g_counts['hpo_weight'] == 0, 'hpo_weight'] = 0.2

    out = g_counts[['gene_symbol', 'hpo_weight']].copy()
    out.rename(columns={'gene_symbol': GENE_COL}, inplace=True)
    return out
