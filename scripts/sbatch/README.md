# sbatch wrappers for the rs2546890 × UKBB MRI pipeline

These scripts submit the pipeline to the **kellis** partition on Luria.
Each wrapper activates the right conda env (`Python_Analysis`), prepends
the required tool paths (`GWAS_env/bin` for plink2, `bolt_lmm/bin` for
BOLT-LMM), and patches the SLURM thread allocation into `config.yaml`
before invoking `bash scripts/run_all.sh`.

| Script | Engine | Threads | Mem | Wall | Use when |
|---|---|---:|---:|---:|---|
| `pipeline_ols.sbatch.sh` | `ols` | 4 | 32G | 8h | Fast OLS-only PheWAS + composites + ROI PCA. ~30 min actual. |
| `pipeline_regenie.sbatch.sh` | `regenie+ols` (default) | 16 | 128G | 48h | Mixed-model headline run. Step 1 dominates wall-clock (~15-25 min/pheno). Tune `lmm.max_phenotypes` in config. |
| `pipeline_bolt.sbatch.sh` | `bolt+ols` | 32 | 150G | 48h | Optional secondary engine. Per-phenotype BOLT runs in serial. |

## Submitting

From the repo root (`lnc_rna_mri/`):

```bash
sbatch scripts/sbatch/pipeline_ols.sbatch.sh
sbatch scripts/sbatch/pipeline_regenie.sbatch.sh
sbatch scripts/sbatch/pipeline_bolt.sbatch.sh
```

Override the engine without editing config:

```bash
sbatch --export=ENGINE_OVERRIDE=regenie     scripts/sbatch/pipeline_regenie.sbatch.sh    # REGENIE only
sbatch --export=ENGINE_OVERRIDE=regenie+ols scripts/sbatch/pipeline_ols.sbatch.sh         # promote ols job
```

## Outputs

Each job writes:
- stdout / stderr → `logs/sbatch/lncrna_<engine>_<jobid>.{out,err}`
- pipeline stdout → `logs/run_log.txt` (also gets per-step REGENIE/BOLT
  driver logs under `logs/regenie_step{1,2}.log`, `logs/bolt_<pheno>.log`)
- results CSVs + figures → `results/`, `figures/` (see top-level README)
- LMM headline results → `results/association_results_*_lmm.csv` +
  `results/report_hierarchical.md`

## Wall-clock scaling for REGENIE step 1

REGENIE step 1 cost is roughly linear in `n_phenotypes` (each gets its own
block-iterative ridge regression). Rough budget on this HM3 BED (444K SNPs)
× 39K imaging subjects × 16 threads:

| Phenotypes | Step 1 wall | Notes |
|---:|---:|---|
| 5 | ~90 min | Demo / fast check |
| 50 | ~6 hr | Bench-test scope |
| ~316 (secondary only) | ~24-36 hr | Composites + ROI PCs only |
| ~1000 (full) | ~3-5 days | Full panel — chunk via SLURM array if needed |

Set the scope in `config.yaml`:
```yaml
lmm:
  scope_families: ["primary", "secondary", "controls"]   # exploratory off by default
  scope_include_exploratory: false                        # set true to add ~700 single IDPs
  max_phenotypes: 0                                       # 0 = no cap; >0 caps for demos
```

If you need to chunk further, the simplest pattern is to split
`results/lmm_inputs/pheno_col_list.txt` into chunks and submit one step 1
job per chunk, then iterate step 2 over the resulting `_pred.list` files.
The chunking helper is not committed — add as a follow-up if needed.

## Resuming a failed job

REGENIE step 1 outputs are cached: as long as `results/lmm/regenie_step1_pred.list`
and `results/lmm/regenie_step1.pheno_sha256` exist and the phenotype list
hasn't changed, `05b_run_regenie.py` will skip step 1 and jump to step 2.
Set `lmm.regenie.skip_step1: true` to force-skip even without the cache.

## Cluster notes

- `Python_Analysis` env: pandas, numpy, scipy, statsmodels, scikit-learn, pyyaml.
- `GWAS_env` env: plink2 binary (needed by step 01 variant extraction).
- `bolt_lmm` env: BOLT-LMM v2.5 binary + plink2.
- `regenie_env` env: REGENIE v3.4.1 binary — the wrapper hard-codes this path
  via `lmm.regenie.binary` in `config.yaml`.
- HapMap3 BED (148 GB) and KING relatedness pairs live under
  `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/geno/`
  and `/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/sqc/`.
  Read paths are baked into `config.yaml` — update if you re-locate them.
