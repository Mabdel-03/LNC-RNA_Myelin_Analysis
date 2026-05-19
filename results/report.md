# rs2546890 → UKBB MRI IDP association — report

_generated: 2026-05-19T18:28:40_

## Variant
- rsID: **rs2546890**
- position: chr5:159332892 (GRCh38)
- effect allele: **A** (modeled in additive dosage)
- other allele:  **G**
- INFO (MFI):    1.0
- est. EAF:      0.5160327363671988
- HWE p:         0.0003608555846144
- hardcall AA: 130414
- hardcall AG: 242210
- hardcall GG: 114785

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

## Phenotypes tested
- primary: 615
- secondary: 1

## Top 20 associations (primary additive model)
| column_name                                                     | panel   | modality           | metric   | region                                       |     n |   beta_per_A |      se |       p |   bonferroni |   fdr_bh |   aa_vs_gg_additive_2beta |
|:----------------------------------------------------------------|:--------|:-------------------|:---------|:---------------------------------------------|------:|-------------:|--------:|--------:|-------------:|---------:|--------------------------:|
| dMRI_TBSS_L3_superior_cerebellar_peduncle_left                  | primary | dMRI_TBSS          | L3       | superior cerebellar peduncle (left)          | 36329 |       0.0272 | 0.0071  | 0.00013 |       0.0799 |   0.0799 |                    0.0544 |
| dMRI_TBSS_MD_superior_cerebellar_peduncle_left                  | primary | dMRI_TBSS          | MD       | superior cerebellar peduncle (left)          | 36338 |       0.0226 | 0.00711 | 0.00145 |       0.893  |   0.305  |                    0.0453 |
| dMRI_TBSS_ISOVF_superior_cerebellar_peduncle_left               | primary | dMRI_TBSS          | ISOVF    | superior cerebellar peduncle (left)          | 36338 |       0.0228 | 0.00718 | 0.00149 |       0.914  |   0.305  |                    0.0456 |
| dMRI_TBSS_L2_superior_cerebellar_peduncle_left                  | primary | dMRI_TBSS          | L2       | superior cerebellar peduncle (left)          | 36325 |       0.0215 | 0.00709 | 0.00245 |       1      |   0.322  |                    0.0429 |
| dMRI_TBSS_MD_superior_cerebellar_peduncle_right                 | primary | dMRI_TBSS          | MD       | superior cerebellar peduncle (right)         | 36335 |       0.0214 | 0.00719 | 0.0029  |       1      |   0.322  |                    0.0428 |
| dMRI_TBSS_L3_superior_cerebellar_peduncle_right                 | primary | dMRI_TBSS          | L3       | superior cerebellar peduncle (right)         | 36334 |       0.0209 | 0.00718 | 0.00365 |       1      |   0.322  |                    0.0418 |
| dMRI_TBSS_L1_tapetum_right                                      | primary | dMRI_TBSS          | L1       | tapetum (right)                              | 36345 |       0.0202 | 0.00694 | 0.00366 |       1      |   0.322  |                    0.0403 |
| dMRI_weighted_mean_ISOVF_inferior_longitudinal_fasciculus_right | primary | dMRI_weighted_mean | ISOVF    | inferior longitudinal fasciculus (right)     | 36334 |       0.0188 | 0.0069  | 0.00631 |       1      |   0.434  |                    0.0377 |
| dMRI_TBSS_ISOVF_superior_cerebellar_peduncle_right              | primary | dMRI_TBSS          | ISOVF    | superior cerebellar peduncle (right)         | 36342 |       0.0196 | 0.00719 | 0.00646 |       1      |   0.434  |                    0.0392 |
| dMRI_weighted_mean_OD_posterior_thalamic_radiation_left         | primary | dMRI_weighted_mean | OD       | posterior thalamic radiation (left)          | 36326 |       0.0189 | 0.00703 | 0.00706 |       1      |   0.434  |                    0.0379 |
| dMRI_TBSS_L1_superior_fronto_occipital_fasciculus_lef           | primary | dMRI_TBSS          | L1       | superior fronto-occipital fasciculus (left)  | 36305 |       0.0164 | 0.00665 | 0.0136  |       1      |   0.656  |                    0.0328 |
| dMRI_TBSS_L1_posterior_limb_of_internal_capsule_left            | primary | dMRI_TBSS          | L1       | posterior limb of internal capsule (left)    | 36339 |       0.0175 | 0.00711 | 0.0142  |       1      |   0.656  |                    0.0349 |
| dMRI_TBSS_OD_superior_fronto_occipital_fasciculus_lef           | primary | dMRI_TBSS          | OD       | superior fronto-occipital fasciculus (left)  | 36327 |      -0.0179 | 0.00734 | 0.0146  |       1      |   0.656  |                   -0.0359 |
| dMRI_TBSS_L1_superior_fronto_occipital_fasciculus_rig           | primary | dMRI_TBSS          | L1       | superior fronto-occipital fasciculus (right) | 36305 |       0.0163 | 0.00671 | 0.015   |       1      |   0.656  |                    0.0326 |
| dMRI_TBSS_FA_cerebral_peduncle_left                             | primary | dMRI_TBSS          | FA       | cerebral peduncle (left)                     | 36337 |       0.0166 | 0.00699 | 0.0179  |       1      |   0.656  |                    0.0331 |
| dMRI_TBSS_L2_cerebral_peduncle_left                             | primary | dMRI_TBSS          | L2       | cerebral peduncle (left)                     | 36334 |      -0.0169 | 0.00717 | 0.0185  |       1      |   0.656  |                   -0.0337 |
| dMRI_TBSS_FA_corticospinal_tract_left                           | primary | dMRI_TBSS          | FA       | corticospinal tract (left)                   | 36346 |       0.0168 | 0.00715 | 0.0191  |       1      |   0.656  |                    0.0335 |
| dMRI_TBSS_FA_superior_cerebellar_peduncle_left                  | primary | dMRI_TBSS          | FA       | superior cerebellar peduncle (left)          | 36313 |      -0.0164 | 0.00701 | 0.0196  |       1      |   0.656  |                   -0.0327 |
| dMRI_weighted_mean_OD_middle_cerebellar_peduncle                | primary | dMRI_weighted_mean | OD       | middle cerebellar peduncle                   | 36305 |      -0.0166 | 0.00716 | 0.0205  |       1      |   0.656  |                   -0.0332 |
| dMRI_TBSS_L1_anterior_limb_of_internal_capsule_left             | primary | dMRI_TBSS          | L1       | anterior limb of internal capsule (left)     | 36330 |       0.0153 | 0.00668 | 0.0216  |       1      |   0.656  |                    0.0307 |

## FDR-q ≤ 0.1: 1 IDP(s)

### Sensitivity (additive 2β vs genotypic AA-vs-GG) for FDR hits
| column_name                                    |   beta_per_A |   aa_vs_gg_additive_2beta |   beta_AA_vs_GG |       p |   wald_p |
|:-----------------------------------------------|-------------:|--------------------------:|----------------:|--------:|---------:|
| dMRI_TBSS_L3_superior_cerebellar_peduncle_left |       0.0272 |                    0.0544 |          0.0521 | 0.00013 | 1.08e-05 |

## Caveats
- UKB imaging confounds file (Smith et al. 2020 Resource 1977) not supplied; only core imaging covariates (age, sex, site, head size, dMRI outlier slices) included.
- OLS additive model is primary; AA-vs-GG read off as 2 × β. A 2-df genotypic model provides the sensitivity check.
- Full panel run (test_mode disabled).
