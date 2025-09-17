# HPO Usage

- Use column 2 (`disease_name`) from `data/hpo/phenotype.hpoa` exactly.
- Find candidates:
  ```bash
  grep -P "^[^#]" data/hpo/phenotype.hpoa | grep -i "<keyword>" | cut -f2 | sort -u | head -50
Example:

css
rivet-gs --disease "Hereditary breast and/or ovarian cancer syndrome" \
         --disease-label "breast cancer" \
         --outdir runs
