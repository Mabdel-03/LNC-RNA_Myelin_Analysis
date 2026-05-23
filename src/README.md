# Source Pipeline Stages

Every script accepts `--config config.yaml` and is normally invoked by a Slurm
wrapper through `scripts/run_all.sh`, except the MS disease-risk scripts, which
are invoked by `scripts/sbatch/pipeline_ms_risk_lmm.sbatch.sh`.

| Stage | Script | Purpose | Main outputs |
|---|---|---|---|
| 00 | `00_inspect_inputs.py` | Validate configured files, tools, and package versions. | `logs/inspect_report.txt`, `logs/versions.txt` |
| 01 | `01_extract_variant.py` | Extract rs2546890 A-dosage from pgen/BGEN/pre-extracted input and assign hardcall genotypes. | `results/intermediate/variant_dosage.tsv`, `results/genotype_qc_summary.csv` |
| 02 | `02_build_phenotype_matrix.py` | Build UKB MRI phenotype matrix, parse field IDs, classify IDP families, and derive RD from L2/L3. | `results/intermediate/phenotypes_wide.tsv`, `results/phenotype_manifest.csv` |
| 03 | `03_build_covariates.py` | Assemble age, sex, site, array, PCs, imaging QC, WMH/vascular covariates where available. | `results/intermediate/covariates.tsv`, `results/covariate_missingness.csv` |
| 04 | `04_merge_qc_sample.py` | Join genotype/phenotype/covariates; apply ancestry, covariate, relatedness, and genotype QC. | `results/intermediate/analysis_ready.tsv`, `keep_ols.txt`, `keep_lmm.txt`, `sample_counts.csv` |
| 04b | `04b_derive_composites_and_pca.py` | Build demyelination/free-water/axonal/geometry composites and tract ROI-PC scores. | ignored extended phenotype files, `roi_pca_loadings.csv`, `phenotype_tiers.csv` |
| 04c | `04c_prepare_lmm_inputs.py` | Materialize REGENIE/BOLT phenotype, covariate, keep, and target-variant files. | `results/lmm_inputs/` |
| 05a | `05a_run_ols.py` | Fit additive OLS, HC3 sensitivity, genotype model, pairwise AA/AG/GG, dominant/recessive tests. | `association_results_*.csv`, `raw_p_testing_summary.csv` |
| 05b | `05b_run_regenie.py` | Run REGENIE quantitative-trait step 1 and rs2546890 step 2. | raw files under ignored `results/lmm/` |
| 05c | `05c_run_bolt.py` | Optional BOLT-LMM engine for selected phenotypes. | raw files under ignored `results/lmm/` |
| 05d | `05d_ingest_lmm_results.py` | Convert REGENIE/BOLT outputs into report-ready aggregate CSVs. | `association_results_*_lmm.csv`, `raw_p_testing_summary_lmm.csv` |
| 06 | `06_plots_and_report.py` | Produce summary plots and technical Markdown reports. | `figures/*.png`, `report.md`, `report_hierarchical.md` |
| Targeted | `targeted_diffusion_check.py` | Filter OLS outputs to FA, L1, and RD targeted tract checks. | `targeted_diffusion_check.csv`, `targeted_diffusion_forest.png` |
| 07 | `07_run_ms_risk.py` | Derive MS cases and fit logistic disease-risk sensitivities. | `association_ms_risk.csv`, `ms_phenotype_audit.csv` |
| 07b | `07b_run_ms_risk_lmm.py` | Fit REGENIE binary-trait/Firth MS-risk model. | `association_ms_risk_lmm.csv` |

Per-eid tables and raw LMM engine outputs are ignored because they may contain
individual-level UKB data or large derived intermediates. Aggregate CSVs and
figures are tracked when they are needed to audit the final report.
