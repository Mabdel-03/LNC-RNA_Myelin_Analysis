# Tests

Tests are synthetic and self-contained. They do not read UK Biobank data.

Primary coverage:

- additive OLS recovery and HC3 sensitivity behavior;
- composite sign conventions and missing-component rules;
- ROI-PCA reproducibility and PC2 retention;
- REGENIE/BOLT result ingestion parsers;
- REGENIE/BOLT phenotype file formatting;
- raw-p-only annotations with no FDR/Bonferroni columns;
- AA/AG/GG pairwise genotype contrast behavior.

Run from the repo root:

```bash
python -m pytest tests -q
```

The batch pipeline can also run tests after a job with
`bash scripts/run_all.sh --with-tests`, but normal full analyses should be
submitted via the Slurm wrappers in `scripts/sbatch/`.
