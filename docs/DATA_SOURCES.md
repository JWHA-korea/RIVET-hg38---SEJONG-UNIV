# Data & Tools Provenance

The `rivet-gs` module does **not** fetch external data at runtime.
It consumes precomputed gene-level scores and thresholds produced upstream.

## Bundled in this repository
- **HPO**: `hp.obo`, `genes_to_phenotype.txt`, `phenotype.hpoa`  
  Version: 2023-10-09 (unmodified). See `NOTICE-HPO.md`.

## Used upstream to generate `data/prioritized.gene.with_thresholds.tsv`
> These resources/tools were used in our upstream RIVET pipeline.  
> Please cite the original sources when publishing results.

- Variant annotation / prediction: snpEff; dbNSFP (incl. CADD, REVEL, PolyPhen-2, SIFT, MutationTaster, etc.); GERP++; phyloP/phastCons
- Population / constraint: gnomAD (exomes v4.1); pLI, LOEUF, missense Z; pext (GTEx-based)
- Clinical labels: ClinVar; ACMG gene list
- Constraint regions: CCR (Constrained Coding Regions)
- Reference & annotation: GRCh38; GENCODE v44

### Notes
- This repo ships derived tables only, not the original datasets above.
- Users must follow each sources license/terms if they download original files.
- Network/pathway features (STRING/Reactome/GO) are **not used at runtime** here; they were applied upstream (if at all) or can be added via optional local files.

