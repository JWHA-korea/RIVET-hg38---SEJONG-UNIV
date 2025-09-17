# RIVET Gene Study (HPO mapping)

Author: JIWOOK HA 
mail : wldnr0728@naver.com
Affiliation: Sejong University, Department of Smart Bio-Industry Convergence  
Purpose: Gene prioritization pipeline (module) prepared for the Sejong University Academic Festival.

## What is this?
- You pass an HPO-registered disease name; the tool resolves its HPO terms / converts them into gene-level weights / combines with your existing gene scores / produces a final, disease-specific priority table.
- Thresholding uses both BestF1 and Precision@95 to derive flags and T1/T2 tiers.
- Report schema (final TSV):  
  `gene, y_prob_max, p95_flag, bestF1_flag, rank, clinvar_plp_flag, tier, score, disease`

## Quickstart
```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .
rivet-gs --disease "Hereditary breast and/or ovarian cancer syndrome" \
         --disease-label "breast cancer" \
         --outdir runs
CLI
rivet-gs --disease "<HPO name>" --outdir runs \
         [--disease-label "<label name>"] \
         [--weights "FUNC:...,NET:...,PATH:...,NOVEL:...,HPO:..."] \
         [--top N] [--min-score X] [--threads K]
The disease name must exactly match column 2 (disease_name) in data/hpo/phenotype.hpoa.

How to find the exact HPO disease name¬
grep -P "^[^#]" data/hpo/phenotype.hpoa | grep -i "<keyword>" | cut -f2 | sort -u | head -50
Data bundle & License
This repo bundles the original HPO files (hp.obo, genes_to_phenotype.txt, phenotype.hpoa; 2023-10-09).

See NOTICE-HPO.md for source/version/terms.

## Data & Licenses
This project bundles unmodified HPO files. Gene scores are derived upstream using third-party resources
(e.g., snpEff, dbNSFP, gnomAD, ClinVar, CCR, GENCODE; optionally STRING/Reactome/GO/PubTator).
See `docs/DATA_SOURCES.md` and `NOTICE-HPO.md`. Please cite the original sources when publishing.

Citation
JIWOOK HA. RIVET Gene Study (HPO mapping)
Sejong Univ. Academic Festival Edition, 2025.

License
MIT Â© JIWOOK HA
