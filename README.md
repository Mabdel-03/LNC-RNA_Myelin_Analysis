# LNC-RNA Myelin Analysis — rs2546890 × UKBB MRI white-matter / myelin IDPs

Reproducible local pipeline that tests whether **rs2546890** (chr5:159332892, GRCh38) allele dosage associates with UK Biobank MRI-derived white-matter / myelin phenotypes (FA, MD, L1–L3, RD, NODDI ICVF/ISOVF/ODI, WMH volume).

The primary model is **OLS additive** on the inverse-rank-normalized IDP with full UKB covariate adjustment. A **genotypic** (AA/AG/GG vs GG) model is fit as a sensitivity check so the AA-vs-GG contrast can be read directly.

---

## Layout

```
lnc_rna_mri/
├── config.yaml                # all user-editable paths and settings
├── requirements.txt
├── scripts/run_all.sh         # orchestrator (00 → 06)
└── src/
    ├── utils.py               # shared helpers
    ├── 00_inspect_inputs.py   # validate paths, tools, packages
    ├── 01_extract_variant.py  # rs2546890 dosage (plink2 / bgen / pre-extracted)
    ├── 02_build_phenotype_matrix.py
    ├── 03_build_covariates.py
    ├── 04_merge_qc_sample.py  # ancestry / relatedness / QC funnel
    ├── 05_run_association.py  # OLS + genotypic + multiple testing
    └── 06_plots_and_report.py # QQ, manhattan-by-IDP, forest, report.md
```

Outputs land in `results/`, `figures/`, `logs/`. With `project.overwrite=false` each full run is written into a timestamped subfolder of `results/`.

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

```bash
# 1) Inspect: validates paths/tools, writes logs/inspect_report.txt
python src/00_inspect_inputs.py --config config.yaml

# 2) Full dry run (no plink2/no model fits, just plan + manifests)
bash scripts/run_all.sh --config config.yaml --dry-run

# 3) Test mode — caps to first 3 phenotypes, runs the real models end-to-end
#    (also requires project.dry_run: false OR drop --dry-run + flip config)
bash scripts/run_all.sh --config config.yaml

# 4) Full run: edit config.yaml → project.dry_run: false, project.test_mode: false
bash scripts/run_all.sh --config config.yaml
```

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
- **Sensitivity (genotypic):** `INRT(IDP) ~ C(genotype, ref=GG) + same covariates`. Reports AG-vs-GG, AA-vs-GG, 2-df Wald p.
- **Outliers** trimmed at `|z| > 6` on the raw scale before transform.
- **WMH** uses `log1p` then INRT.
- **Multiple testing:** Bonferroni + BH-FDR per panel (primary / secondary).

LMM (REGENIE/BOLT) is **not** the default. Setting `models.use_lmm: true` halts with a command template rather than silently falling back.

---

## Outputs

| File | What |
|---|---|
| `results/genotype_qc_summary.csv` | rsID, alleles, EAF, missingness, hardcalls, HWE p, INFO, method |
| `results/phenotype_manifest.csv` | one row per IDP (panel, transform, include flag) |
| `results/sample_counts.csv` | sample funnel (n at every QC step) |
| `results/covariate_missingness.csv` | per-covariate missing n/% |
| `results/association_results_primary.csv` | β, SE, t, p, n, FDR, Bonferroni, `aa_vs_gg_additive_2beta` |
| `results/association_results_genotypic.csv` | AG-vs-GG, AA-vs-GG, 2-df Wald p |
| `results/report.md` | run summary, top hits, sensitivity comparison, caveats |
| `figures/qqplot_pvalues.png` | QQ of primary p-values |
| `figures/manhattan_like_idp_results.png` | per-IDP −log10(p) grouped by panel/modality |
| `figures/effect_size_forest_top_hits.png` | top-20 β±CI |
| `figures/genotype_violin_top3.png` | optional: raw distributions by AA/AG/GG for top 3 |
| `logs/inspect_report.txt`, `logs/run_log.txt`, `logs/versions.txt`, `logs/model_warnings.txt` | diagnostics |

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
- OLS only — set `use_lmm: true` to emit a REGENIE command template; no silent LMM fallback.
- Hard-call genotypes are assigned only where dosage falls within `1 - hardcall_probability_threshold` of an integer; samples below that are NA for the genotypic model but kept for the additive one.
- The local imputed pgen is **GRCh37** despite the path naming suggesting v3. If you swap in a GRCh38 source, flip `variant.pgen_build: "GRCh38"` so the extractor queries with `pos_grch38` instead.
