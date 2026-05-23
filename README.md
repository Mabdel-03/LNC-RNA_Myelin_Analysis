# LNC-RNA Myelin Analysis — rs2546890 × UKBB MRI white-matter / myelin IDPs

Reproducible local pipeline that tests whether **rs2546890** (chr5:159332892, GRCh38) allele dosage associates with UK Biobank MRI-derived white-matter / myelin phenotypes (FA, MD, L1–L3, RD, NODDI ICVF/ISOVF/ODI, WMH volume).

The default engine is a **mixed-model (REGENIE)** run in parallel with OLS for sensitivity. The pipeline organizes phenotypes into focused primary, secondary, replication, control, and exploratory families and reports **raw p values only**. It does not compute or report adjusted p values.

> **Interpretation guardrail (enforced in `report_hierarchical.md`):** Ordinary diffusion MRI signals (FA, MD, L1-L3, ICVF, OD, ISOVF) are *indirect* and should not be labeled as myelination effects unless supported by myelin-sensitive MRI (MTR, MTsat, MWF, qT1) or orthogonal validation. The current expected interpretation is that rs2546890-A appears more consistent with a free-water or tract-geometry signal than a canonical demyelination signal unless the composite analysis below shows otherwise.

---

## Layout

```
lnc_rna_mri/
├── config.yaml                # all user-editable paths and settings
├── requirements.txt
├── figures/                   # tracked final figures; see figures/README.md
├── results/                   # tracked aggregate outputs; see results/README.md
├── scripts/
│   ├── README.md              # script-level organization
│   ├── run_all.sh             # batch-only orchestrator (00 → 06 + targeted report)
│   └── sbatch/                # canonical Slurm wrappers
├── tests/                     # synthetic/unit tests; see tests/README.md
└── src/
    ├── README.md              # step-by-step source responsibilities
    ├── utils.py               # shared helpers
    ├── 00_inspect_inputs.py   # validate paths, tools, packages
    ├── 01_extract_variant.py  # rs2546890 dosage (plink2 / bgen / pre-extracted)
    ├── 02_build_phenotype_matrix.py
    ├── 03_build_covariates.py
    ├── 04_merge_qc_sample.py  # ancestry / relatedness / QC funnel
    ├── 05a_run_ols.py        # OLS additive + model variants + genotype pairwise tests
    ├── 05b_run_regenie.py    # REGENIE focused mixed-model run
    ├── 05d_ingest_lmm_results.py
    ├── 06_plots_and_report.py # QQ, manhattan-by-IDP, forest, report.md
    ├── 07_run_ms_risk.py      # MS G35 logistic disease-risk checks
    ├── 07b_run_ms_risk_lmm.py # REGENIE-BT/Firth MS risk model
    └── targeted_diffusion_check.py
```

Outputs land in `results/`, `figures/`, `logs/`. Aggregate CSV/Markdown/PNG outputs are tracked when they are useful for review; per-eid UKB tables, raw LMM internals, LOCO files, and Slurm stdout/stderr are ignored.

The current narrative report is [`results/variant_story_report.md`](results/variant_story_report.md). It ties the script provenance, phenotype construction, statistical models, figures, and MS-risk validation into the rs2546890-A story.

---

## Inputs the pipeline expects

Set the paths inside `config.yaml`. Defaults are wired to the local Kellis-lab UKB tree (`/net/bmc-lab4/...`, `/net/bmc-lab5/...`, `/home/mabdel03/...`).

| Section | Key | What it points to | Required? |
|---|---|---|---|
| `inputs` | `pgen_prefix` | plink2 pgen prefix for imputed genotypes (GRCh38) | primary path |
| `inputs` | `variant_metadata_file` | UKB MFI v3 TSV (zst ok) — for INFO / EAF | optional |
| `inputs` | `chr5_bgen` + `sample_file` | BGEN fallback if pgen absent | fallback |
| `inputs` | `preextracted_variant_dosage_file` | CSV with `eid, dosage_A` to skip extraction | optional |
| `inputs` | `basket_tab` | UKB basket TSV (used for IDPs + covariates) | primary path |
| `inputs` | `phenotype_file` | one wide CSV of phenotypes — alt to basket | optional |
| `inputs` | `data_dictionary_file` | TSV/CSV with `field_id, description` | optional (better labels) |
| `inputs` | `mri_pheno_tsv_dirs` | already-extracted per-field TSVs to enrich | optional |
| `inputs` | `genetic_covariates_file` | pre-built TSV with eid + age/sex/site/PCs | optional |
| `inputs` | `imaging_confounds_file` | official UKB imaging confounds TSV | optional |
| `inputs` | `ancestry_relatedness_file` | TSV with `eid, ancestry_label, kin_keep` | optional |

`00_inspect_inputs.py` prints which of these resolve on this machine before any heavy work runs.

---

## How to run

### Interactive checks

```bash
# Lightweight inspection is still safe interactively.
python src/00_inspect_inputs.py --config config.yaml

# Development-only local smoke checks require an explicit override.
ALLOW_LOCAL_RUN=1 bash scripts/run_all.sh --config config.yaml --dry-run
```

### SLURM on Luria — kellis partition

Full pipelines are batch-only. Submit the wrapper sbatch scripts; they
activate the right conda env, set tool paths, cap BLAS threads, and pass the
Slurm CPU allocation into OLS/REGENIE/BOLT:

```bash
# OLS only: phenotype-level OLS workers across the Slurm CPU allocation
sbatch scripts/sbatch/pipeline_ols.sbatch.sh

# Mixed-model headline (REGENIE step1+step2 + OLS sensitivity)
# 16 cpus / 300G / up to 48h — scales linearly with lmm.max_phenotypes
sbatch scripts/sbatch/pipeline_regenie.sbatch.sh

# Optional BOLT-LMM secondary engine
# 32 cpus / 300G / up to 48h
sbatch scripts/sbatch/pipeline_bolt.sbatch.sh

# MS disease-risk validation (logistic sensitivity + REGENIE-BT/Firth)
sbatch scripts/sbatch/pipeline_ms_risk_lmm.sbatch.sh
```

See [`scripts/sbatch/README.md`](scripts/sbatch/README.md) for wall-clock
budgeting (REGENIE step 1 is ~15-25 min/phenotype) and how to override the
engine at submit time.

Each script accepts `--config` and `--dry-run`. `--dry-run` overrides the value in config for that step.

---

## Variant extraction priority

`src/01_extract_variant.py` tries the following in order:

1. **Pre-extracted dosage CSV** (`inputs.preextracted_variant_dosage_file`). Must contain `eid` and `dosage_A` (a `genotype` column with `AA/AG/GG` is also consumed).
2. **plink2 on pgen** (default on this machine). The extractor probes multiple ID forms in order and falls back to a coordinate range:
   - `rs<id>` (UKB-imputed pgen often lacks rsIDs)
   - `<chr>:<pos>:<ref>:<alt>` and the swapped variant (UKB-style colon ID)
   - `chr<chr>:<pos>:<ref>:<alt>` (Ensembl-style)
   - `--chr X --from-bp/--to-bp` coordinate range as last resort

   The position used is `variant.pos_grch37` or `variant.pos_grch38` depending on `variant.pgen_build` (the imputed UKB v3 pgen here is GRCh37). The plink2 `.raw` header is parsed to detect which allele is being counted; if it's the **other** allele (G), the dosage is flipped (`dosage_A = 2 - dosage_G`).
3. **BGEN fallback** via plink2 (or qctool/bgenix if you wire one up).
4. If none of the above apply, the step writes templated extraction commands to `logs/variant_extraction_commands.sh` and exits with a clear error.

INFO and EAF are looked up from `inputs.variant_metadata_file` (UKB MFI v3) when available.

For rs2546890 on this machine, the variant matches as `5:158759900:A:G` (GRCh37); the effect allele A is the alphabetically-first allele so plink2 counts it directly, no dosage flip needed.

---

## Statistical design (one-line summary)

- **Primary:** `INRT(IDP) ~ dosage_A + age + age² + sex + age×sex + C(site) + head_size + dMRI_motion + C(array) + PC1..PCk`. `aa_vs_gg_additive = 2 × β`.
- **Model variants:** additive OLS with nonrobust and HC3 SEs, categorical genotype 2-df model, AA/AG/GG pairwise contrasts, dominant A, recessive A, and focused REGENIE mixed-model results.
- **Sensitivity (genotypic):** `INRT(IDP) ~ C(genotype, ref=GG) + same covariates`. Reports AG-vs-GG, AA-vs-GG, 2-df Wald raw p.
- **Disease risk:** `MS_G35 ~ dosage_A + covariates` via logistic regression in the unrelated cohort, plus REGENIE binary-trait/Firth on the larger kinship-tolerant set.
- **Outliers** trimmed at `|z| > 6` on the raw scale before transform.
- **WMH** uses `log1p` then INRT.
- **Multiplicity handling:** raw p values only. Summary tables count `p < 0.05`, `p < 0.01`, and `p < 0.001`; no correction is applied.

**Focused phenotype tiers:** the confirmatory family is FA/MD/RD in pre-registered TBSS skeleton callosal and projection-tract ROIs. Directional AD/L1, NODDI ICVF/ISOVF, weighted-tract replication, WMH/pathology controls, and exploratory QSM/T2*/GWC/structural outcomes are tagged separately in the phenotype manifest. Composite phenotypes and ROI-level PCA scores are still available as secondary summaries.

**Mixed-model engine:** by default (`models.engine: regenie+ols`) the pipeline runs REGENIE step 1 (LOCO null model, multi-phenotype) + step 2 (rs2546890 only) using the prior SI-loneliness REGENIE infrastructure (`/home/mabdel03/data/software/regenie/regenie`, HapMap3 BED at `…/ukb_genoHM3/ukb_genoHM3_bed`, model SNPs at `…/ukb_genoHM3_modelSNPs.txt`, imputed pgen at `/net/bmc-lab5/…/ukb_imp`). BOLT-LMM is available as a secondary engine (`engine: bolt` or `bolt+ols`) and by default runs only the 15-or-so composite + ROI-PC phenotypes.

Engine choices in `models.engine`:
- `ols`              — OLS only (back-compat)
- `regenie`          — REGENIE LMM only
- `bolt`             — BOLT-LMM only (composites + ROI PCs by default)
- `regenie+ols`      — both, REGENIE is the report headline (default)
- `bolt+ols`         — both, BOLT is the report headline
- `all`              — all three engines

---

## Outputs

| File | What |
|---|---|
| `results/genotype_qc_summary.csv` | rsID, alleles, EAF, missingness, hardcalls, HWE p, INFO, method |
| `results/phenotype_manifest.csv` | one row per IDP (panel, family, transform, include flag) |
| `results/sample_counts.csv` | sample funnel (n at every QC step) |
| `results/covariate_missingness.csv` | per-covariate missing n/% |
| `results/association_results_primary.csv` | β, SE, t, raw p, n, raw-p ranks/flags, `aa_vs_gg_additive_2beta` (OLS, single IDPs) |
| `results/association_results_genotypic.csv` | AG-vs-GG, AA-vs-GG, 2-df Wald raw p (OLS) |
| `results/association_results_model_matrix.csv` | additive, HC3, dominant, recessive, and genotypic raw-p model comparisons |
| `results/association_results_pairwise_genotype.csv` | AA-vs-GG, AG-vs-GG, AA-vs-AG covariate-conditioned raw-p contrasts |
| `results/association_results_composites.csv` | OLS on composite phenotypes (NEW) |
| `results/association_results_roi_pca.csv` | OLS on ROI-PCA phenotypes (NEW) |
| `results/association_results_*_lmm.csv` | parallel LMM (REGENIE/BOLT) results (NEW) |
| `results/association_ms_risk.csv` | logistic MS-risk sensitivity models: additive, dominant, recessive, AA/AG/GG pairwise |
| `results/association_ms_risk_lmm.csv` | REGENIE binary-trait/Firth MS-risk headline |
| `results/raw_p_testing_summary.csv` | per-family raw-p test counts + nominal hit counts (OLS) |
| `results/raw_p_testing_summary_lmm.csv` | per-family raw-p test counts + nominal hit counts (LMM) |
| `results/multiple_testing_summary*.csv` | compatibility copies of the raw-p summaries; they do not contain adjusted p values |
| `results/model_diagnostics_summary.csv` | per-phenotype r², cond no, BP/JB p (NEW) |
| `results/composite_phenotypes.csv` | per-eid composite scores (NEW) |
| `results/roi_pca_phenotypes.csv` | per-eid ROI PC scores (NEW) |
| `results/roi_pca_loadings.csv` | per-tract PC1/PC2 metric loadings (NEW) |
| `results/phenotype_tiers.csv` | phenotype → family / panel / metric / region (NEW) |
| `results/report.md` | back-compat OLS run summary |
| `results/report_hierarchical.md` | LMM headline + OLS sensitivity + composite-driven interpretation (NEW) |
| `results/variant_story_report.md` | prose scientific report explaining phenotype construction, models, results, and biological interpretation |
| `figures/qqplot_pvalues.png` | QQ of primary p-values |
| `figures/manhattan_like_idp_results.png` | per-IDP −log10(p) grouped by panel/modality |
| `figures/effect_size_forest_top_hits.png` | top-20 β±CI |
| `figures/genotype_violin_top3.png` | optional: raw distributions by AA/AG/GG for top 3 |
| `figures/effect_heatmap_by_tract_metric.png` | NEW: tracts × metrics, color=signed −log10(p) |
| `figures/composite_effects_forest.png` | NEW: composite × tract β±CI |
| `figures/top_roi_pca_loadings.png` | NEW: PC1 loadings for top-ranked tracts |
| `figures/targeted_diffusion_forest.png` | FA/L1/RD targeted tract forest plot; generated by `targeted_diffusion_check.py` |
| `logs/inspect_report.txt`, `logs/run_log.txt`, `logs/versions.txt`, `logs/model_warnings.txt` | diagnostics |
| `logs/regenie_step{1,2}.log` | LMM driver logs |

---

## Reproducibility

- `project.random_seed` seeded in 05.
- `logs/versions.txt` records Python, package, and external-tool versions per run.
- `logs/run_log.txt` captures every step's stdout/stderr with timestamps.
- Set `project.overwrite: false` to write each run into `results/YYYYmmdd_HHMMSS/`.

---

## Sample QC defaults (local)

`config.yaml` wires the standard UKB Sample-QC file as the single source for genetic covariates AND ancestry / unrelated flags:

- `inputs.genetic_covariates_file` → `sqc.20220316.tsv.gz` (488K samples, 80+ cols incl. population, sex, age, genotyping_array, used_in_pca_calculation, sr_WB, UKB_PC1..PC40).
- `inputs.relatedness_file` → `ukb_rel_a21942_s488172.dat.gz` (KING pairs).
- Step 03 maps `UKB_PC1..PC{N}` → `PC1..PC{N}`, `genotyping_array` → `array`, plus pulls imaging covariates (site, head_size, motion_dmri = field 25746 dMRI outlier slices, acq_date) from the basket.
- Step 04 applies (in order):
  1. **Ancestry**: `sample.ancestry_label="White British"` → SQC `population_MM == "WB_MM"` (also accepts `population == "WB"` or `sr_WB == 1`).
  2. **Unrelated**: if `inputs.relatedness_file` is set → greedy graph prune of KING pairs with `Kinship >= sample.kinship_threshold` (default 0.0884 = 3rd-degree). Else falls back to SQC's `used_in_pca_calculation == TRUE` flag.

## Known limitations

- BGEN code path is templated but **not exercised locally** (only pgen is available here).
- UKB imaging confounds file (Smith et al. 2020 Resource 1977) is not on disk; only core imaging covariates (age, sex, site, head size, dMRI outlier slices) are included. Drop the official file into `inputs.imaging_confounds_file` to merge it in.
- LMM engines (REGENIE primary, BOLT-LMM secondary) are now first-class run targets via `models.engine`. The legacy `models.use_lmm` flag is preserved but superseded.
- Hard-call genotypes are assigned only where dosage falls within `1 - hardcall_probability_threshold` of an integer; samples below that are NA for the genotypic model but kept for the additive one.
- The local imputed pgen is **GRCh37** despite the path naming suggesting v3. If you swap in a GRCh38 source, flip `variant.pgen_build: "GRCh38"` so the extractor queries with `pos_grch38` instead.
