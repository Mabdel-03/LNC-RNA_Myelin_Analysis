# rs2546890-A, White-Matter Vulnerability, and MS Risk

Generated from the completed `lnc_rna_mri` pipeline outputs on 2026-05-23.

This report summarizes the evidence that **rs2546890-A** is associated with a
modest increase in multiple sclerosis risk and with tract-specific UK Biobank
MRI diffusion phenotypes that are relevant to demyelination vulnerability. All
p values are **raw p values only**. No FDR, Bonferroni, q values, or adjusted
p values are used.

## Executive Summary

The core finding is convergent rather than singular. The A allele has a small
MS-risk signal, and the same allele is associated with diffusion MRI differences
in long white-matter pathways that are plausible targets or substrates for MS
vulnerability: **left cerebral peduncle**, **left corticospinal tract**,
**posterior thalamic radiation**, **uncinate fasciculus**, **forceps/callosal
tracts**, and **superior cerebellar peduncles**.

The single-IDP pattern is strongest in the **left cerebral peduncle**, where
rs2546890-A is associated with higher FA and lower RD/MD. This is not a simple
active-demyelination pattern, because active demyelination is often expected to
increase RD. A better interpretation is that the variant tags a baseline
white-matter microstructure axis, or an injury/repair susceptibility axis, that
is visible in diffusion MRI and overlaps anatomically with MS-relevant tracts.

The disease-risk analysis anchors the imaging signal. In the REGENIE
binary-trait/Firth model, rs2546890-A was associated with MS case status at
`OR = 1.043`, 95% CI `1.000-1.088`, raw `p = 0.0492` in `n = 430,354`.
The unrelated-cohort logistic sensitivity model was similar:
`OR = 1.048`, 95% CI `1.003-1.095`, raw `p = 0.0347`.

The bottom line is that **rs2546890-A is best framed as a candidate
white-matter/MS vulnerability allele**, not as proof of direct MRI-measured
demyelination. UKB diffusion IDPs are indirect microstructure phenotypes.
They can support a demyelination-relevant story only when interpreted with
tract anatomy, NODDI/free-water context, genotype contrasts, WMH controls, and
disease-risk validation.

## Analysis Provenance

The table below maps the scientific steps in this report to the scripts that
produced them.

| step | script | purpose | key outputs used here |
|---|---|---|---|
| Input audit | [`../src/00_inspect_inputs.py`](../src/00_inspect_inputs.py) | Checks configured paths, tools, and package versions before heavy work. | `../logs/inspect_report.txt`, `../logs/versions.txt` |
| Variant extraction | [`../src/01_extract_variant.py`](../src/01_extract_variant.py) | Extracts rs2546890 dosage from UKB pgen/BGEN/pre-extracted input, flips allele coding if needed, and assigns hardcall genotypes. | `genotype_qc_summary.csv`, ignored per-eid `intermediate/variant_dosage.tsv` |
| MRI phenotype construction | [`../src/02_build_phenotype_matrix.py`](../src/02_build_phenotype_matrix.py) | Parses UKB field IDs, creates the phenotype manifest, classifies IDP families, and derives RD from L2/L3. | `phenotype_manifest.csv`, ignored per-eid `intermediate/phenotypes_wide.tsv` |
| Covariates | [`../src/03_build_covariates.py`](../src/03_build_covariates.py) | Builds age, sex, site, array, PCs, imaging-QC, pathology, and vascular covariates when available. | `covariate_missingness.csv`, ignored per-eid `intermediate/covariates.tsv` |
| Sample QC | [`../src/04_merge_qc_sample.py`](../src/04_merge_qc_sample.py) | Joins genotype/phenotype/covariates; applies ancestry, covariate, relatedness, and genotype QC. | `sample_counts.csv`, ignored keep lists |
| Composites and ROI PCA | [`../src/04b_derive_composites_and_pca.py`](../src/04b_derive_composites_and_pca.py) | Builds tract composites and ROI-PC scores from matched diffusion metrics. | `phenotype_tiers.csv`, `roi_pca_loadings.csv`, ignored per-eid composite/PCA tables |
| LMM input writing | [`../src/04c_prepare_lmm_inputs.py`](../src/04c_prepare_lmm_inputs.py) | Writes REGENIE/BOLT phenotype, covariate, keep, and target-variant files. | ignored `lmm_inputs/` |
| OLS and model variants | [`../src/05a_run_ols.py`](../src/05a_run_ols.py) | Fits additive OLS, HC3 sensitivity, genotype model, pairwise genotype contrasts, dominant/recessive models. | `association_results_primary.csv`, `association_results_model_matrix.csv`, `association_results_pairwise_genotype.csv` |
| REGENIE QT LMM | [`../src/05b_run_regenie.py`](../src/05b_run_regenie.py), [`../src/05d_ingest_lmm_results.py`](../src/05d_ingest_lmm_results.py) | Runs kinship-aware quantitative-trait LMM and ingests rs2546890 results. | `association_results_composites_lmm.csv`, `association_results_roi_pca_lmm.csv`, `association_results_controls_lmm.csv` |
| Plot/report generation | [`../src/06_plots_and_report.py`](../src/06_plots_and_report.py) | Produces technical reports and main figures. | `report.md`, `report_hierarchical.md`, `../figures/*.png` |
| Targeted diffusion check | [`../src/targeted_diffusion_check.py`](../src/targeted_diffusion_check.py) | Filters OLS results to FA, L1, and RD in targeted tracts. | `targeted_diffusion_check.csv`, `../figures/targeted_diffusion_forest.png` |
| MS logistic risk | [`../src/07_run_ms_risk.py`](../src/07_run_ms_risk.py) | Derives MS case status and fits additive/dominant/recessive/pairwise logistic models. | `association_ms_risk.csv`, `ms_phenotype_audit.csv` |
| MS REGENIE-BT risk | [`../src/07b_run_ms_risk_lmm.py`](../src/07b_run_ms_risk_lmm.py) | Runs REGENIE binary-trait/Firth MS-risk model. | `association_ms_risk_lmm.csv` |

Canonical execution is through the Slurm wrappers in
[`../scripts/sbatch/`](../scripts/sbatch/). The MRI pipeline is run by
`pipeline_ols.sbatch.sh`, `pipeline_regenie.sbatch.sh`, or
`pipeline_bolt.sbatch.sh`. The disease-risk validation is run by
`pipeline_ms_risk_lmm.sbatch.sh`.

## Variant and Sample Definition

The variant was modeled as A-allele dosage.

| item | value |
|---|---:|
| Variant | rs2546890 |
| GRCh38 position | chr5:159332892 |
| GRCh37 position used in UKB pgen | chr5:158759900 |
| Effect allele | A |
| Other allele | G |
| Estimated A-allele frequency | 0.516 |
| INFO / MFI | 1.0 |
| Hardcall AA / AG / GG | 130,414 / 242,210 / 114,785 |
| HWE p | 0.000361 |

The imaging sample was restricted to the configured white-British analysis set.
OLS used an unrelated subset; LMM used the larger kinship-tolerant keep set.

| sample funnel step | n |
|---|---:|
| Genotype input | 487,409 |
| Phenotype input | 502,357 |
| Covariate input | 488,221 |
| Merged genotype/phenotype/covariate data | 487,150 |
| White British via SQC | 430,518 |
| Required covariates non-missing | 430,518 |
| OLS unrelated set | 399,161 |
| LMM keep set | 430,354 |

## Phenotype Construction

The MRI phenotypes are UK Biobank imaging-derived phenotypes (IDPs). They are
ROI/tract summaries rather than voxelwise maps. The pipeline therefore uses a
tract/ROI design rather than a voxelwise design.

### Diffusion Metrics

The diffusion families were parsed from UKB field IDs in
[`../src/02_build_phenotype_matrix.py`](../src/02_build_phenotype_matrix.py).
The main metrics are:

| metric | construction and interpretation |
|---|---|
| FA | Fractional anisotropy from TBSS skeleton or tractography-weighted tracts. Higher FA often reflects more directional diffusion, but can also reflect tract geometry, fiber coherence, crossing fibers, or axonal packing. |
| MD | Mean diffusivity. Higher MD can reflect less restricted diffusion, free water, edema, or tissue loss. |
| L1 | Principal tensor eigenvalue, used here as axial diffusivity context. |
| L2/L3 | Secondary tensor eigenvalues, used to derive RD. |
| RD | Derived phenotype: `(L2 + L3) / 2`. RD is demyelination-relevant but not myelin-specific. |
| ICVF | Intracellular volume fraction from NODDI-style UKB IDPs, used as axon-packing context. |
| ISOVF | Isotropic volume fraction, used as free-water/extracellular-water context. |
| OD/MO | Orientation/anisotropy context used in tract-geometry composites where available. |

### Phenotype Families

The analysis was organized into interpretable families rather than a flat
all-IDP screen.

| family | construction | role in interpretation |
|---|---|---|
| Primary skeleton tensor | FA, MD, and derived RD in pre-registered TBSS skeleton callosal/projection tracts from `focused_analysis.tbss_primary_fa_base_ids`. | Main tract-level diffusion family. |
| Directional tensor context | L1 in the same primary TBSS tracts. | Helps distinguish axis-specific diffusion from RD/MD shifts. |
| Secondary mechanistic | ICVF and ISOVF in the same primary TBSS tracts. | Adds axon-packing and free-water context. |
| Replication weighted tensor | FA, MD, and RD in tractography-weighted tracts from `focused_analysis.weighted_replication_fa_base_ids`. | Anatomically interpretable replication family. |
| Replication mechanistic | ICVF and ISOVF in weighted tracts. | Mechanistic replication context. |
| Controls | WMH volume/count/mean volume. | Tests whether findings are explained by gross white-matter lesion burden. |
| Exploratory context | QSM, T2*, grey-white contrast, and structural fields. | Secondary biological context, not headline evidence. |

### Composite Scores

Composite scores were built by
[`../src/04b_derive_composites_and_pca.py`](../src/04b_derive_composites_and_pca.py).
Within each tract, component metrics were standardized and combined using
signed weights. A score required at least two available components.

| composite | signed components | intended meaning |
|---|---|---|
| demyelination-like | `+RD +MD -FA -ICVF +ISOVF` | A diffusion pattern often discussed as demyelination-relevant, while still indirect. |
| free-water-like | `+MD +ISOVF +L1 +L2 +L3` | Free-water or extracellular-water shift. |
| axonal-loss-like | `-L1 -FA -ICVF` | Lower axial/directional/neurite-density context. |
| tract-geometry-like | `-OD +FA +MO` | Geometry/coherence/orientation context. |

These composites are not external biomarkers. They are analytic summaries that
make related diffusion metrics easier to interpret together.

### ROI PCA

ROI PCA scores were also built by
[`../src/04b_derive_composites_and_pca.py`](../src/04b_derive_composites_and_pca.py).
For configured tracts, available diffusion metrics were z-scored, residual
missingness was mean-imputed after filtering, and PCA was fit with two possible
components. PC1 was retained; PC2 was retained only when it explained at least
10% of variance. Loadings are stored in `roi_pca_loadings.csv`.

### MS Case Phenotype

MS case status was derived by
[`../src/07_run_ms_risk.py`](../src/07_run_ms_risk.py) from UKB disease fields:
ICD-10 G35 hospital/all-source diagnosis, self-report code 1261, and
first-occurrence MS fields.

| source | count |
|---|---:|
| ICD-10 G35 | 2,190 |
| Self-report code 1261 | 1,856 |
| First occurrence | 5,217 |
| Union MS cases in basket | 5,217 |

## Statistical Models and Rationale

### Additive OLS

The broad imaging scan used additive OLS in the unrelated white-British subset:

```text
INRT(IDP) ~ dosage_A + age + age^2 + sex + age:sex + site
          + head_size + dMRI_motion + array + PC1..PC20 + error
```

Outliers were trimmed at `|z| > 6` before inverse-rank normalization. WMH
outcomes were log1p-transformed before inverse-rank normalization. The additive
model is the primary high-throughput scan because it is simple, interpretable,
and efficient across hundreds of IDPs.

### Robust and Genotype-Form Sensitivities

[`../src/05a_run_ols.py`](../src/05a_run_ols.py) also fit model variants:

| model | rationale |
|---|---|
| HC3 additive OLS | Checks sensitivity to heteroskedasticity without changing the point estimate. |
| 2-df genotypic model | Tests whether genotype categories depart from a strict additive trend. |
| AA vs GG, AG vs GG, AA vs AG pairwise contrasts | Makes genotype-group differences explicit and supports interpretation of AA/AG/GG trends. |
| Dominant A | Tests whether AA and AG jointly differ from GG. |
| Recessive A | Tests whether AA differs from AG/GG. |

### REGENIE Quantitative-Trait LMM

REGENIE was used for kinship-aware validation on a focused phenotype set:

```text
IDP ~ dosage_A + covariates + genetic relationship structure
```

Step 1 fit LOCO predictions from model SNPs, and step 2 tested rs2546890 only.
The final successful imaging REGENIE run focused on selected composite, ROI-PCA,
and WMH control phenotypes. It should be interpreted as focused validation, not
as a full replacement for the complete OLS scan.

### MS Disease-Risk Models

MS risk was assessed two ways:

| model | rationale |
|---|---|
| Logistic regression in unrelated subjects | Direct disease-risk sensitivity model matching the OLS unrelated sample philosophy. |
| REGENIE binary-trait/Firth | Larger kinship-tolerant disease-risk model with Firth correction for a low-prevalence binary phenotype. |

## Main Imaging Results

The most coherent single-IDP finding is in the **left cerebral peduncle**:
higher FA, lower RD, and lower MD per A allele.

| tract / IDP | model | n | beta per A | raw p | interpretation |
|---|---:|---:|---:|---:|---|
| Left cerebral peduncle FA | OLS | 32,560 | +0.0272 | 9.12e-05 | Strongest targeted single-IDP signal; higher directional diffusion. |
| Left cerebral peduncle RD | OLS | 32,558 | -0.0259 | 2.14e-04 | Lower radial diffusivity; paired with higher FA. |
| Left cerebral peduncle MD | OLS | 32,555 | -0.0163 | 0.0220 | Lower mean diffusivity; supports a restricted-diffusion pattern. |
| Left corticospinal tract FA | OLS | 32,565 | +0.0186 | 0.0110 | Motor/projection-tract extension. |
| Right uncinate fasciculus FA | OLS | 32,561 | +0.0168 | 0.0187 | Association-tract context. |
| Right posterior thalamic radiation ISOVF | OLS | 32,557 | +0.0158 | 0.0218 | Free-water context in a projection radiation. |
| Left posterior thalamic radiation L1 | OLS | 32,555 | -0.0134 | 0.0274 | Directional diffusivity context. |

![Targeted diffusion forest](../figures/targeted_diffusion_forest.png)

The heatmap shows that the association is not a uniform whole-brain shift. It
is concentrated in selected projection, cerebellar, callosal, and association
tract contexts.

![Effect heatmap by tract and metric](../figures/effect_heatmap_by_tract_metric.png)

## Composite and LMM Results

The focused REGENIE/LMM run gave its strongest evidence in the **superior
cerebellar peduncles**.

| phenotype | n | beta per A | raw p | note |
|---|---:|---:|---:|---|
| Demyelination-like, superior cerebellar peduncle left | 39,218 | -0.0220 | 0.000939 | Top LMM composite result. |
| Free-water-like, superior cerebellar peduncle right | 39,218 | -0.0220 | 0.00121 | Bilateral cerebellar-peduncle support. |
| Free-water-like, superior cerebellar peduncle left | 39,218 | -0.0216 | 0.00138 | Same tract family and similar direction. |
| Demyelination-like, superior cerebellar peduncle right | 39,218 | -0.0200 | 0.00272 | Bilateral composite support. |
| Superior cerebellar peduncle left PC1 | 39,218 | -0.0215 | 0.00120 | ROI-PCA validation. |
| Superior cerebellar peduncle right PC1 | 39,218 | -0.0188 | 0.00490 | Bilateral ROI-PCA validation. |

Because the demyelination-like composite is signed as `+RD +MD -FA -ICVF +ISOVF`,
a negative beta means the A allele is associated with a lower score on
that composite in these focused LMM phenotypes. This again argues against a
simple "the allele causes active demyelination visible as high RD" story. The
stronger story is that the allele marks tract-specific microstructure and MS
susceptibility, potentially reflecting baseline architecture, repair biology,
or vulnerability to later inflammatory injury.

![Composite effects forest](../figures/composite_effects_forest.png)

![Top ROI PCA loadings](../figures/top_roi_pca_loadings.png)

## Genotype-Group Results

The AA/AG/GG contrasts make the genotype pattern visible without assuming a
single additive line. The top AA-vs-GG diffusion contrasts include callosal and
forceps phenotypes.

| phenotype | contrast | n | beta | raw p | interpretation |
|---|---|---:|---:|---:|---|
| Splenium RD | AA vs GG | 30,126 | -0.0314 | 0.0346 | Lower RD in AA than GG. |
| Splenium FA | AA vs GG | 30,134 | +0.0326 | 0.0350 | Higher FA in AA than GG. |
| Forceps major FA | AA vs GG | 30,148 | +0.0318 | 0.0356 | Tractography-weighted callosal replication context. |
| Forceps major RD | AA vs GG | 30,145 | -0.0319 | 0.0365 | Directionally paired with forceps FA. |

This callosal/forceps pattern matters because interhemispheric white matter is
often implicated in demyelinating disease. It also supports the broader trend:
rs2546890-A is associated with selected long-tract diffusion differences rather
than with a generic MRI artifact.

![Genotype violin plots for top hits](../figures/genotype_violin_top3.png)

## MS Disease-Risk Results

The MS-risk evidence is modest but directionally consistent.

| model | n | cases / controls | OR per A | 95% CI | raw p |
|---|---:|---:|---:|---:|---:|
| REGENIE-BT Firth additive | 430,354 | 4,431 / 425,923 | 1.043 | 1.000-1.088 | 0.0492 |
| Logistic additive sensitivity | 399,161 | 4,101 / 395,060 | 1.048 | 1.003-1.095 | 0.0347 |
| Logistic dominant A | 399,161 | 4,101 / 395,060 | 1.089 | 1.011-1.174 | 0.0246 |
| Logistic AA vs GG | 200,287 | 2,031 / 198,256 | 1.104 | 1.011-1.206 | 0.0274 |
| Logistic AG vs GG | 291,325 | 2,960 / 288,365 | 1.081 | 0.999-1.170 | 0.0528 |
| Logistic AA vs AG | 306,710 | 3,211 / 303,499 | 1.019 | 0.947-1.096 | 0.6128 |

The disease-risk result is important because it anchors the imaging phenotype
to a clinical endpoint. A common allele with OR near 1.04 is not expected to
produce large single-tract MRI effects. The useful signal is the convergence:
small disease risk plus tract-specific white-matter microstructure differences.

## Control and Alternative Explanations

The focused LMM control result for WMH volume was null:

| control phenotype | n | beta per A | raw p | interpretation |
|---|---:|---:|---:|---|
| WMH volume | 38,995 | +0.000351 | 0.953 | No evidence that the variant's MRI story is explained by gross WMH burden. |

This does not exclude all vascular or inflammatory explanations. It does argue
against a simple model where the variant merely increases visible WMH lesion
burden in UKB and all diffusion findings follow from that.

## Biological Interpretation

The most defensible interpretation is a **tract-specific vulnerability model**:

1. rs2546890-A is associated with a small increase in MS risk.
2. The same allele is associated with diffusion MRI differences in long tracts
   that are anatomically relevant to MS.
3. The diffusion signature is not a direct lesion or direct myelin measure.
4. The pattern may reflect baseline tract architecture, axon/myelin packing,
   extracellular-water context, fiber coherence, or repair/injury-response
   biology.
5. This makes rs2546890-A a good candidate for follow-up in datasets with
   direct myelin-sensitive imaging, MS-enriched cohorts, lesion maps, or
   longitudinal conversion/progression outcomes.

![Top effect-size forest](../figures/effect_size_forest_top_hits.png)

## Limitations

- Raw p values are reported descriptively. No multiple-testing correction is
  used in this report.
- UKB diffusion IDPs are indirect microstructure phenotypes. They do not
  directly quantify myelin content.
- The strongest OLS single-IDP results and focused REGENIE/LMM composite
  results are convergent but not identical phenotype sets.
- The final successful imaging LMM was focused on selected composite/ROI-PCA
  and WMH-control phenotypes, not the full OLS panel.
- MS-risk effects are small and should be treated as candidate-variant evidence,
  not as proof of mechanism.
- Follow-up with MTsat, MTR, MWF, quantitative T1, lesion-level MRI, or
  MS-enriched cohorts would be needed to test direct demyelination mechanisms.

## Source Tables and Figures

Primary source tables:

- `association_results_primary.csv`
- `targeted_diffusion_check.csv`
- `association_results_composites.csv`
- `association_results_pairwise_genotype.csv`
- `association_results_model_matrix.csv`
- `association_results_composites_lmm.csv`
- `association_results_roi_pca_lmm.csv`
- `association_results_controls_lmm.csv`
- `association_ms_risk.csv`
- `association_ms_risk_lmm.csv`
- `phenotype_manifest.csv`
- `phenotype_tiers.csv`
- `roi_pca_loadings.csv`

Embedded figures:

- `../figures/targeted_diffusion_forest.png`
- `../figures/effect_heatmap_by_tract_metric.png`
- `../figures/composite_effects_forest.png`
- `../figures/top_roi_pca_loadings.png`
- `../figures/genotype_violin_top3.png`
- `../figures/effect_size_forest_top_hits.png`
