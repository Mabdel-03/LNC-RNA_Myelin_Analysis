# Pipeline Scripts

This directory contains the user-facing execution entrypoints. The pipeline is
intended to run from Slurm batch jobs, not from long interactive sessions.

## Canonical Entry Points

- `run_all.sh` is the shared orchestrator used by the MRI batch wrappers. It
  runs input inspection, variant extraction, phenotype/covariate construction,
  sample QC, composites/PCA, selected model engines, plots/reports, and the
  targeted diffusion summary when OLS outputs are produced.
- `sbatch/pipeline_ols.sbatch.sh` runs the OLS-only analysis.
- `sbatch/pipeline_regenie.sbatch.sh` runs the REGENIE+OLS headline analysis.
- `sbatch/pipeline_bolt.sbatch.sh` runs the optional BOLT+OLS engine.
- `sbatch/pipeline_ms_risk_lmm.sbatch.sh` runs the MS disease-risk validation
  scripts, including REGENIE binary-trait/Firth.

`run_all.sh` refuses non-Slurm execution unless `ALLOW_LOCAL_RUN=1` is set.
That override is only for short development smoke checks.

## Local Recovery Artifacts

Files named `pipeline_regenie_resume*.sbatch.sh` are ignored by git. They were
used for local recovery/debugging of failed REGENIE attempts and should not be
treated as public pipeline entrypoints.
