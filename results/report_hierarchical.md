# rs2546890 → UKBB MRI: hierarchical analysis report

_generated: 2026-05-19T20:18:17_

- rsID: **rs2546890**
- position (GRCh37): chr5:158759900
- effect allele: **A** / other: **G**
- INFO: 1.0, EAF: 0.5160327363671988, HWE p: 0.0003608555846144

## Required scientific framing

> Ordinary diffusion MRI signals (FA, MD, L1-L3, ICVF, OD, ISOVF) are *indirect* and should not be labeled as myelination effects unless supported by myelin-sensitive MRI (MTR, MTsat, MWF, qT1) or orthogonal validation. The current expected interpretation is that rs2546890-A appears more consistent with a free-water or tract-geometry signal than a canonical demyelination signal unless the composite analysis below shows otherwise.

## Sample funnel
| step | n | n dropped |
|---|---:|---:|
| genotype_input | 487,409 | 0 |
| phenotype_input | 502,357 | 0 |
| covariate_input | 488,221 | 0 |
| merged_inner_join | 487,150 | 0 |
| ancestry=White British_via_sqc | 430,518 | 56,632 |
| required_covars_nonmissing | 430,518 | 0 |
| ols_unrelated_king≥0.0884 | 399,161 | 31,357 |
| lmm_keep_kinship<0.354 | 430,354 | 164 |

## Per-family multiple-testing summary (LMM)
| family    |   n_tests |   alpha_bonf |   n_sig_bonf |   n_sig_fdr_05 |   n_sig_fdr_10 |   min_p |   q_at_top_hit |
|:----------|----------:|-------------:|-------------:|---------------:|---------------:|--------:|---------------:|
| controls  |         1 |       0.05   |            0 |              0 |              0 |   0.953 |          0.953 |
| secondary |         4 |       0.0125 |            0 |              0 |              0 |   0.696 |          0.903 |

## Family 1 — primary myelin-sensitive MRI
_(LMM results not present)_

## Family 2a — biological composites (LMM headline)
| column_name                                             | family    | panel     | modality            | metric             | region                          | source    |     n |   beta_per_A |      se |   chisq |     p |   log10p | engine   |   info |   eaf | status   | effect_allele   | other_allele   |   family_bonferroni |   family_fdr_bh |
|:--------------------------------------------------------|:----------|:----------|:--------------------|:-------------------|:--------------------------------|:----------|------:|-------------:|--------:|--------:|------:|---------:|:---------|-------:|------:|:---------|:----------------|:---------------|--------------------:|----------------:|
| demyelination_like__wm__acoustic_radiation_left         | secondary | secondary | dMRI_composite_wm   | demyelination_like | acoustic radiation (left)       | composite | 39216 |     -0.0037  | 0.0095  |  0.152  | 0.696 |   0.157  | regenie  |    nan | 0.481 | ok       | A               | G              |                   1 |           0.903 |
| demyelination_like__wm__acoustic_radiation_right        | secondary | secondary | dMRI_composite_wm   | demyelination_like | acoustic radiation (right)      | composite | 39216 |      0.00268 | 0.00972 |  0.0758 | 0.783 |   0.106  | regenie  |    nan | 0.481 | ok       | A               | G              |                   1 |           0.903 |
| demyelination_like__tbss__anterior_corona_radiata_left  | secondary | secondary | dMRI_composite_tbss | demyelination_like | anterior corona radiata (left)  | composite | 39218 |      0.00225 | 0.00968 |  0.0542 | 0.816 |   0.0883 | regenie  |    nan | 0.481 | ok       | A               | G              |                   1 |           0.903 |
| demyelination_like__tbss__anterior_corona_radiata_right | secondary | secondary | dMRI_composite_tbss | demyelination_like | anterior corona radiata (right) | composite | 39218 |      0.00118 | 0.00965 |  0.0148 | 0.903 |   0.0443 | regenie  |    nan | 0.481 | ok       | A               | G              |                   1 |           0.903 |

### OLS sensitivity for composites
| column_name                                                   |   field_id | family    | panel     | modality            | metric              | region                               | source    | effect_allele   | other_allele   |     n |   beta_per_A |      se |     t |       p | vcov_type   |   aa_vs_gg_additive_2beta | status   |   bonferroni |   fdr_bh |   family_bonferroni |   family_fdr_bh |
|:--------------------------------------------------------------|-----------:|:----------|:----------|:--------------------|:--------------------|:-------------------------------------|:----------|:----------------|:---------------|------:|-------------:|--------:|------:|--------:|:------------|--------------------------:|:---------|-------------:|---------:|--------------------:|----------------:|
| free_water_like__tbss__superior_cerebellar_peduncle_left      |        nan | secondary | secondary | dMRI_composite_tbss | free_water_like     | superior cerebellar peduncle (left)  | composite | A               | G              | 36352 |       0.0472 | 0.014   |  3.38 | 0.00072 | nonrobust   |                    0.0944 | ok       |        0.247 |    0.247 |               0.246 |           0.246 |
| demyelination_like__tbss__superior_cerebellar_peduncle_left   |        nan | secondary | secondary | dMRI_composite_tbss | demyelination_like  | superior cerebellar peduncle (left)  | composite | A               | G              | 36352 |       0.0329 | 0.0106  |  3.1  | 0.00192 | nonrobust   |                    0.0658 | ok       |        0.659 |    0.295 |               0.657 |           0.294 |
| free_water_like__tbss__superior_cerebellar_peduncle_right     |        nan | secondary | secondary | dMRI_composite_tbss | free_water_like     | superior cerebellar peduncle (right) | composite | A               | G              | 36352 |       0.0413 | 0.0141  |  2.93 | 0.00344 | nonrobust   |                    0.0826 | ok       |        1     |    0.295 |               1     |           0.294 |
| demyelination_like__tbss__superior_cerebellar_peduncle_right  |        nan | secondary | secondary | dMRI_composite_tbss | demyelination_like  | superior cerebellar peduncle (right) | composite | A               | G              | 36352 |       0.0282 | 0.0105  |  2.68 | 0.00739 | nonrobust   |                    0.0563 | ok       |        1     |    0.507 |               1     |           0.505 |
| axonal_loss_like__tbss__tapetum_right                         |        nan | secondary | secondary | dMRI_composite_tbss | axonal_loss_like    | tapetum (right)                      | composite | A               | G              | 36352 |      -0.0161 | 0.00685 | -2.35 | 0.0188  | nonrobust   |                   -0.0322 | ok       |        1     |    0.772 |               1     |           0.77  |
| tract_geometry_like__wm__posterior_thalamic_radiation_left    |        nan | secondary | secondary | dMRI_composite_wm   | tract_geometry_like | posterior thalamic radiation (left)  | composite | A               | G              | 36351 |      -0.0212 | 0.00916 | -2.32 | 0.0206  | nonrobust   |                   -0.0424 | ok       |        1     |    0.772 |               1     |           0.77  |
| tract_geometry_like__tbss__cerebral_peduncle_left             |        nan | secondary | secondary | dMRI_composite_tbss | tract_geometry_like | cerebral peduncle (left)             | composite | A               | G              | 36352 |       0.0234 | 0.0105  |  2.23 | 0.0254  | nonrobust   |                    0.0468 | ok       |        1     |    0.772 |               1     |           0.77  |
| axonal_loss_like__tbss__cerebral_peduncle_left                |        nan | secondary | secondary | dMRI_composite_tbss | axonal_loss_like    | cerebral peduncle (left)             | composite | A               | G              | 36352 |      -0.0172 | 0.00779 | -2.2  | 0.0277  | nonrobust   |                   -0.0343 | ok       |        1     |    0.772 |               1     |           0.77  |
| axonal_loss_like__tbss__corticospinal_tract_left              |        nan | secondary | secondary | dMRI_composite_tbss | axonal_loss_like    | corticospinal tract (left)           | composite | A               | G              | 36352 |      -0.0161 | 0.0074  | -2.18 | 0.0294  | nonrobust   |                   -0.0322 | ok       |        1     |    0.772 |               1     |           0.77  |
| tract_geometry_like__tbss__posterior_thalamic_radiation_right |        nan | secondary | secondary | dMRI_composite_tbss | tract_geometry_like | posterior thalamic radiation (right) | composite | A               | G              | 36352 |      -0.0216 | 0.0101  | -2.13 | 0.0332  | nonrobust   |                   -0.0431 | ok       |        1     |    0.772 |               1     |           0.77  |

## Family 2b — ROI PCA (LMM headline)
_(no ROI PCA LMM results)_

### OLS sensitivity for ROI PCA
| column_name                                        |   field_id | family    | panel     | modality     | metric   | region                                     | source   | effect_allele   | other_allele   |     n |   beta_per_A |     se |     t |       p | vcov_type   |   aa_vs_gg_additive_2beta | status   |   bonferroni |   fdr_bh |   family_bonferroni |   family_fdr_bh |
|:---------------------------------------------------|-----------:|:----------|:----------|:-------------|:---------|:-------------------------------------------|:---------|:----------------|:---------------|------:|-------------:|-------:|------:|--------:|:------------|--------------------------:|:---------|-------------:|---------:|--------------------:|----------------:|
| roi__superior_cerebellar_peduncle_left__PC1        |        nan | secondary | secondary | dMRI_roi_pca | PC1      | superior cerebellar peduncle (left)        | roi_pca  | A               | G              | 36349 |       0.0446 | 0.0151 |  2.96 | 0.00307 | nonrobust   |                    0.0893 | ok       |            1 |    0.295 |                   1 |           0.294 |
| roi__superior_cerebellar_peduncle_right__PC1       |        nan | secondary | secondary | dMRI_roi_pca | PC1      | superior cerebellar peduncle (right)       | roi_pca  | A               | G              | 36349 |       0.0363 | 0.0151 |  2.41 | 0.0159  | nonrobust   |                    0.0727 | ok       |            1 |    0.772 |                   1 |           0.77  |
| roi__superior_cerebellar_peduncle_right__PC2       |        nan | secondary | secondary | dMRI_roi_pca | PC2      | superior cerebellar peduncle (right)       | roi_pca  | A               | G              | 36349 |       0.0258 | 0.012  |  2.15 | 0.0317  | nonrobust   |                    0.0516 | ok       |            1 |    0.772 |                   1 |           0.77  |
| roi__anterior_limb_of_internal_capsule_left__PC2   |        nan | secondary | secondary | dMRI_roi_pca | PC2      | anterior limb of internal capsule (left)   | roi_pca  | A               | G              | 36349 |      -0.0217 | 0.0102 | -2.12 | 0.0338  | nonrobust   |                   -0.0435 | ok       |            1 |    0.772 |                   1 |           0.77  |
| roi__posterior_thalamic_radiation_left__PC2        |        nan | secondary | secondary | dMRI_roi_pca | PC2      | posterior thalamic radiation (left)        | roi_pca  | A               | G              | 36349 |       0.027  | 0.0136 |  1.99 | 0.0465  | nonrobust   |                    0.0541 | ok       |            1 |    0.888 |                   1 |           0.885 |
| roi__superior_cerebellar_peduncle_left__PC2        |        nan | secondary | secondary | dMRI_roi_pca | PC2      | superior cerebellar peduncle (left)        | roi_pca  | A               | G              | 36349 |       0.023  | 0.0121 |  1.9  | 0.0572  | nonrobust   |                    0.046  | ok       |            1 |    0.888 |                   1 |           0.885 |
| roi__posterior_limb_of_internal_capsule_left__PC2  |        nan | secondary | secondary | dMRI_roi_pca | PC2      | posterior limb of internal capsule (left)  | roi_pca  | A               | G              | 36349 |       0.0249 | 0.0132 |  1.89 | 0.0593  | nonrobust   |                    0.0499 | ok       |            1 |    0.888 |                   1 |           0.885 |
| roi__posterior_limb_of_internal_capsule_right__PC2 |        nan | secondary | secondary | dMRI_roi_pca | PC2      | posterior limb of internal capsule (right) | roi_pca  | A               | G              | 36349 |       0.0242 | 0.0131 |  1.85 | 0.0649  | nonrobust   |                    0.0485 | ok       |            1 |    0.888 |                   1 |           0.885 |
| roi__corticospinal_tract_left__PC2                 |        nan | secondary | secondary | dMRI_roi_pca | PC2      | corticospinal tract (left)                 | roi_pca  | A               | G              | 36349 |       0.0247 | 0.0144 |  1.72 | 0.086   | nonrobust   |                    0.0494 | ok       |            1 |    0.888 |                   1 |           0.885 |
| roi__posterior_thalamic_radiation_right__PC2       |        nan | secondary | secondary | dMRI_roi_pca | PC2      | posterior thalamic radiation (right)       | roi_pca  | A               | G              | 36349 |       0.023  | 0.0136 |  1.69 | 0.0908  | nonrobust   |                    0.0459 | ok       |            1 |    0.888 |                   1 |           0.885 |

## Family 3 — exploratory single-IDP PheWAS (LMM headline)

_Interpret single-IDP hits below only in light of the composite + ROI-PC results above._

_(no LMM single-IDP results)_

### OLS sensitivity for single IDPs (top 10 by OLS p)
| column_name                                                     |   field_id | family      | panel   | modality           | metric   | region                                   | source   | effect_allele   | other_allele   |     n |   beta_per_A |      se |    t |       p | vcov_type   |   aa_vs_gg_additive_2beta | status   |   bonferroni |   fdr_bh |   family_bonferroni |   family_fdr_bh |
|:----------------------------------------------------------------|-----------:|:------------|:--------|:-------------------|:---------|:-----------------------------------------|:---------|:----------------|:---------------|------:|-------------:|--------:|-----:|--------:|:------------|--------------------------:|:---------|-------------:|---------:|--------------------:|----------------:|
| dMRI_TBSS_L3_superior_cerebellar_peduncle_left                  |   2.53e+04 | exploratory | primary | dMRI_TBSS          | L3       | superior cerebellar peduncle (left)      | basket   | A               | G              | 36329 |       0.0272 | 0.0071  | 3.83 | 0.00013 | nonrobust   |                    0.0544 | ok       |       0.0799 |   0.0799 |              0.0799 |          0.0799 |
| dMRI_TBSS_MD_superior_cerebellar_peduncle_left                  |   2.51e+04 | exploratory | primary | dMRI_TBSS          | MD       | superior cerebellar peduncle (left)      | basket   | A               | G              | 36338 |       0.0226 | 0.00711 | 3.18 | 0.00145 | nonrobust   |                    0.0453 | ok       |       0.893  |   0.305  |              0.893  |          0.305  |
| dMRI_TBSS_ISOVF_superior_cerebellar_peduncle_left               |   2.55e+04 | exploratory | primary | dMRI_TBSS          | ISOVF    | superior cerebellar peduncle (left)      | basket   | A               | G              | 36338 |       0.0228 | 0.00718 | 3.18 | 0.00149 | nonrobust   |                    0.0456 | ok       |       0.914  |   0.305  |              0.914  |          0.305  |
| dMRI_TBSS_L2_superior_cerebellar_peduncle_left                  |   2.53e+04 | exploratory | primary | dMRI_TBSS          | L2       | superior cerebellar peduncle (left)      | basket   | A               | G              | 36325 |       0.0215 | 0.00709 | 3.03 | 0.00245 | nonrobust   |                    0.0429 | ok       |       1      |   0.322  |              1      |          0.322  |
| dMRI_TBSS_MD_superior_cerebellar_peduncle_right                 |   2.51e+04 | exploratory | primary | dMRI_TBSS          | MD       | superior cerebellar peduncle (right)     | basket   | A               | G              | 36335 |       0.0214 | 0.00719 | 2.98 | 0.0029  | nonrobust   |                    0.0428 | ok       |       1      |   0.322  |              1      |          0.322  |
| dMRI_TBSS_L3_superior_cerebellar_peduncle_right                 |   2.53e+04 | exploratory | primary | dMRI_TBSS          | L3       | superior cerebellar peduncle (right)     | basket   | A               | G              | 36334 |       0.0209 | 0.00718 | 2.91 | 0.00365 | nonrobust   |                    0.0418 | ok       |       1      |   0.322  |              1      |          0.322  |
| dMRI_TBSS_L1_tapetum_right                                      |   2.52e+04 | exploratory | primary | dMRI_TBSS          | L1       | tapetum (right)                          | basket   | A               | G              | 36345 |       0.0202 | 0.00694 | 2.91 | 0.00366 | nonrobust   |                    0.0403 | ok       |       1      |   0.322  |              1      |          0.322  |
| dMRI_weighted_mean_ISOVF_inferior_longitudinal_fasciculus_right |   2.57e+04 | exploratory | primary | dMRI_weighted_mean | ISOVF    | inferior longitudinal fasciculus (right) | basket   | A               | G              | 36334 |       0.0188 | 0.0069  | 2.73 | 0.00631 | nonrobust   |                    0.0377 | ok       |       1      |   0.434  |              1      |          0.434  |
| dMRI_TBSS_ISOVF_superior_cerebellar_peduncle_right              |   2.55e+04 | exploratory | primary | dMRI_TBSS          | ISOVF    | superior cerebellar peduncle (right)     | basket   | A               | G              | 36342 |       0.0196 | 0.00719 | 2.72 | 0.00646 | nonrobust   |                    0.0392 | ok       |       1      |   0.434  |              1      |          0.434  |
| dMRI_weighted_mean_OD_posterior_thalamic_radiation_left         |   2.57e+04 | exploratory | primary | dMRI_weighted_mean | OD       | posterior thalamic radiation (left)      | basket   | A               | G              | 36326 |       0.0189 | 0.00703 | 2.69 | 0.00706 | nonrobust   |                    0.0379 | ok       |       1      |   0.434  |              1      |          0.434  |

## Family 4 — controls (WMH / QC outcomes)
| column_name   | family   | panel     | modality   | metric     |   region | source   |     n |   beta_per_A |      se |   chisq |     p |   log10p | engine   |   info |   eaf | status   | effect_allele   | other_allele   |   family_bonferroni |   family_fdr_bh |
|:--------------|:---------|:----------|:-----------|:-----------|---------:|:---------|------:|-------------:|--------:|--------:|------:|---------:|:---------|-------:|------:|:---------|:----------------|:---------------|--------------------:|----------------:|
| WMH_volume    | controls | secondary | T2_FLAIR   | WMH_volume |      nan | basket   | 38995 |     0.000351 | 0.00599 | 0.00344 | 0.953 |   0.0208 | regenie  |    nan | 0.481 | ok       | A               | G              |               0.953 |           0.953 |

## LMM (REGENIE) vs OLS comparison

## Composite-driven interpretation

**Auto-classified signature: `inconclusive`**

Interpretation key:
- `canonical_demyelination` — `demyelination_like` composite dominates the FDR-significant set
- `free_water_dominant` — `free_water_like` composite dominates
- `tract_geometry_dominant` — `tract_geometry_like` composite dominates
- `axonal_loss_dominant` — `axonal_loss_like` composite dominates
- `inconclusive` — no composite has any tract hit at family FDR ≤ 0.10

## Limitations
- No myelin-sensitive MRI present in this UKB extract (family 1 empty).
- ROI PCA uses mean imputation for residual missing data → biased if missingness is informative.
- Family-FDR applied independently per family (not sequential gatekeeping).
- REGENIE step-1 LOCO model fit on the imaging subsample (~40K), not the full WB-MM cohort.
- BOLT-LMM (when used) relies on rs2546890 being in HapMap3; for imputed-dosage stats see lmm.bolt.regen_bgen documentation.

## Figures
- `figures/qqplot_pvalues.png`
- `figures/manhattan_like_idp_results.png`
- `figures/effect_size_forest_top_hits.png`
- `figures/effect_heatmap_by_tract_metric.png` (NEW)
- `figures/composite_effects_forest.png` (NEW)
- `figures/top_roi_pca_loadings.png` (NEW)

