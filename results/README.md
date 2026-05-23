# Results Directory

This directory stores aggregate outputs that document the rs2546890 analysis.
Only reviewable, non-individual-level summaries should be committed.

## Tracked Summary Outputs

- `genotype_qc_summary.csv`, `sample_counts.csv`, and
  `covariate_missingness.csv` document the analysis inputs and sample funnel.
- `phenotype_manifest.csv` and `phenotype_tiers.csv` document how UKB fields
  were classified into primary, secondary, replication, control, and
  exploratory families.
- `association_results_primary.csv`, `association_results_composites.csv`,
  `association_results_roi_pca.csv`, `association_results_genotypic.csv`,
  `association_results_pairwise_genotype.csv`, and
  `association_results_model_matrix.csv` are OLS-derived aggregate outputs.
- `association_results_*_lmm.csv` files are ingested REGENIE/BOLT aggregate
  outputs.
- `association_ms_risk.csv` and `association_ms_risk_lmm.csv` summarize the MS
  disease-risk analyses.
- `report.md`, `report_hierarchical.md`, and `variant_story_report.md` are the
  human-readable reports.

## Ignored Outputs

`intermediate/`, `variant/`, `lmm/`, `lmm_inputs/`, `lmm_inputs_ms/`,
`composite_phenotypes.csv`, `roi_pca_phenotypes.csv`, and
`ms_phenotype.tsv` are ignored. They are either large, per-eid, raw engine
outputs, or Slurm/runtime artifacts and should not be pushed.

## P-Value Policy

All result summaries use raw p values only. Files named
`multiple_testing_summary*.csv` are compatibility copies of raw-p summaries;
they do not contain adjusted p values.
