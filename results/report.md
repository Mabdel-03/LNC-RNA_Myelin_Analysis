# rs2546890 → UKBB MRI IDP association — report

_generated: 2026-05-23T07:15:02_

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
- primary: 198
- exploratory: 104
- secondary: 3

## Top 20 associations (primary additive model)
| column_name                                                     | panel       | modality           | metric   | region                                                     |     n |   beta_per_A |      se |        p |    raw_p |   minus_log10_raw_p |   raw_p_rank_within_family |   aa_vs_gg_additive_2beta |
|:----------------------------------------------------------------|:------------|:-------------------|:---------|:-----------------------------------------------------------|------:|-------------:|--------:|---------:|---------:|--------------------:|---------------------------:|--------------------------:|
| dMRI_TBSS_FA_cerebral_peduncle_left                             | primary     | dMRI_TBSS          | FA       | cerebral peduncle (left)                                   | 32560 |       0.0272 | 0.00694 | 9.12e-05 | 9.12e-05 |                4.04 |                          1 |                    0.0543 |
| dMRI_TBSS_RD_cerebral_peduncle_left                             | primary     | dMRI_TBSS          | RD       | cerebral peduncle (left)                                   | 32558 |      -0.0259 | 0.007   | 0.000214 | 0.000214 |                3.67 |                          2 |                   -0.0519 |
| SWI_deep_gm_QSM_f24483                                          | exploratory | SWI_deep_gm        | QSM      | nan                                                        | 32535 |      -0.0209 | 0.00779 | 0.0074   | 0.0074   |                2.13 |                          1 |                   -0.0417 |
| SWI_T2star_T2star_median_t2star_in_pallidum_left                | exploratory | SWI_T2star         | T2star   | Median T2star in pallidum (left)                           | 32542 |      -0.0198 | 0.00776 | 0.0109   | 0.0109   |                1.96 |                          2 |                   -0.0395 |
| dMRI_TBSS_FA_corticospinal_tract_left                           | primary     | dMRI_TBSS          | FA       | corticospinal tract (left)                                 | 32565 |       0.0186 | 0.00731 | 0.011    | 0.011    |                1.96 |                          3 |                    0.0372 |
| SWI_T2star_T2star_median_t2star_in_putamen_right                | exploratory | SWI_T2star         | T2star   | Median T2star in putamen (right)                           | 32547 |      -0.0187 | 0.00762 | 0.0142   | 0.0142   |                1.85 |                          3 |                   -0.0374 |
| dMRI_weighted_mean_FA_uncinate_fasciculus_right                 | primary     | dMRI_weighted_mean | FA       | uncinate fasciculus (right)                                | 32561 |       0.0168 | 0.00714 | 0.0187   | 0.0187   |                1.73 |                          1 |                    0.0335 |
| dMRI_TBSS_ISOVF_posterior_thalamic_radiation_right              | primary     | dMRI_TBSS          | ISOVF    | posterior thalamic radiation (right)                       | 32557 |       0.0158 | 0.0069  | 0.0218   | 0.0218   |                1.66 |                          1 |                    0.0317 |
| dMRI_TBSS_MD_cerebral_peduncle_left                             | primary     | dMRI_TBSS          | MD       | cerebral peduncle (left)                                   | 32555 |      -0.0163 | 0.00713 | 0.022    | 0.022    |                1.66 |                          4 |                   -0.0326 |
| dMRI_weighted_mean_L1_posterior_thalamic_radiation_left         | primary     | dMRI_weighted_mean | L1       | posterior thalamic radiation (left)                        | 32555 |      -0.0134 | 0.00608 | 0.0274   | 0.0274   |                1.56 |                          1 |                   -0.0268 |
| T1_GWC_GWC_grey_white_contrast_in_cuneus_right_hemisphe         | exploratory | T1_GWC             | GWC      | Grey-white contrast in cuneus (right hemisphere)           | 32144 |      -0.0151 | 0.00697 | 0.0303   | 0.0303   |                1.52 |                          4 |                   -0.0302 |
| dMRI_weighted_mean_ICVF_parahippocampal_part_of_cingulum_left   | primary     | dMRI_weighted_mean | ICVF     | parahippocampal part of cingulum (left)                    | 32551 |       0.015  | 0.0071  | 0.0346   | 0.0346   |                1.46 |                          1 |                    0.03   |
| dMRI_weighted_mean_ISOVF_superior_longitudinal_fasciculus_right | primary     | dMRI_weighted_mean | ISOVF    | superior longitudinal fasciculus (right)                   | 32559 |       0.0138 | 0.00662 | 0.0366   | 0.0366   |                1.44 |                          2 |                    0.0277 |
| SWI_deep_gm_QSM_f24467                                          | exploratory | SWI_deep_gm        | QSM      | nan                                                        | 32538 |       0.015  | 0.00733 | 0.0403   | 0.0403   |                1.39 |                          5 |                    0.0301 |
| SWI_deep_gm_QSM_f24482                                          | exploratory | SWI_deep_gm        | QSM      | nan                                                        | 32537 |       0.0154 | 0.00772 | 0.0465   | 0.0465   |                1.33 |                          6 |                    0.0307 |
| T1_GWC_GWC_grey_white_contrast_in_isthmuscingulate_righ         | exploratory | T1_GWC             | GWC      | Grey-white contrast in isthmuscingulate (right hemisphere) | 32143 |      -0.014  | 0.00706 | 0.047    | 0.047    |                1.33 |                          7 |                   -0.028  |
| SWI_T2star_T2star_median_t2star_in_pallidum_right               | exploratory | SWI_T2star         | T2star   | Median T2star in pallidum (right)                          | 32543 |      -0.0152 | 0.00777 | 0.0506   | 0.0506   |                1.3  |                          8 |                   -0.0304 |
| dMRI_weighted_mean_ISOVF_posterior_thalamic_radiation_right     | primary     | dMRI_weighted_mean | ISOVF    | posterior thalamic radiation (right)                       | 32544 |       0.0129 | 0.00662 | 0.0514   | 0.0514   |                1.29 |                          3 |                    0.0258 |
| dMRI_weighted_mean_ICVF_uncinate_fasciculus_right               | primary     | dMRI_weighted_mean | ICVF     | uncinate fasciculus (right)                                | 32553 |       0.0128 | 0.00669 | 0.0552   | 0.0552   |                1.26 |                          4 |                    0.0257 |
| dMRI_weighted_mean_MD_posterior_thalamic_radiation_left         | primary     | dMRI_weighted_mean | MD       | posterior thalamic radiation (left)                        | 32547 |      -0.0114 | 0.00592 | 0.0552   | 0.0552   |                1.26 |                          2 |                   -0.0227 |

## Raw p < 0.05: 16 IDP(s)

### Sensitivity (additive 2β vs genotypic AA-vs-GG) for raw-p hits
| column_name                                                     |   beta_per_A |   aa_vs_gg_additive_2beta |   beta_AA_vs_GG |        p |   wald_p |
|:----------------------------------------------------------------|-------------:|--------------------------:|----------------:|---------:|---------:|
| SWI_deep_gm_QSM_f24483                                          |      -0.0209 |                   -0.0417 |         -0.0402 | 0.0074   | 0.00704  |
| SWI_T2star_T2star_median_t2star_in_pallidum_left                |      -0.0198 |                   -0.0395 |         -0.0371 | 0.0109   | 0.00126  |
| SWI_T2star_T2star_median_t2star_in_putamen_right                |      -0.0187 |                   -0.0374 |         -0.0363 | 0.0142   | 0.025    |
| T1_GWC_GWC_grey_white_contrast_in_cuneus_right_hemisphe         |      -0.0151 |                   -0.0302 |         -0.0312 | 0.0303   | 0.0484   |
| SWI_deep_gm_QSM_f24467                                          |       0.015  |                    0.0301 |          0.0291 | 0.0403   | 0.0651   |
| SWI_deep_gm_QSM_f24482                                          |       0.0154 |                    0.0307 |          0.0294 | 0.0465   | 0.0503   |
| T1_GWC_GWC_grey_white_contrast_in_isthmuscingulate_righ         |      -0.014  |                   -0.028  |         -0.0291 | 0.047    | 0.0622   |
| dMRI_TBSS_FA_cerebral_peduncle_left                             |       0.0272 |                    0.0543 |          0.0549 | 9.12e-05 | 0.000371 |
| dMRI_TBSS_RD_cerebral_peduncle_left                             |      -0.0259 |                   -0.0519 |         -0.0524 | 0.000214 | 0.000848 |
| dMRI_TBSS_FA_corticospinal_tract_left                           |       0.0186 |                    0.0372 |          0.0375 | 0.011    | 0.0372   |
| dMRI_TBSS_MD_cerebral_peduncle_left                             |      -0.0163 |                   -0.0326 |         -0.0328 | 0.022    | 0.0713   |
| dMRI_weighted_mean_L1_posterior_thalamic_radiation_left         |      -0.0134 |                   -0.0268 |         -0.0266 | 0.0274   | 0.0823   |
| dMRI_weighted_mean_ICVF_parahippocampal_part_of_cingulum_left   |       0.015  |                    0.03   |          0.0298 | 0.0346   | 0.105    |
| dMRI_weighted_mean_ISOVF_superior_longitudinal_fasciculus_right |       0.0138 |                    0.0277 |          0.0273 | 0.0366   | 0.101    |
| dMRI_weighted_mean_FA_uncinate_fasciculus_right                 |       0.0168 |                    0.0335 |          0.0319 | 0.0187   | 0.00934  |
| dMRI_TBSS_ISOVF_posterior_thalamic_radiation_right              |       0.0158 |                    0.0317 |          0.0315 | 0.0218   | 0.0713   |

## Caveats
- UKB imaging confounds file (Smith et al. 2020 Resource 1977) not supplied; only core imaging covariates (age, sex, site, head size, dMRI outlier slices) included.
- OLS additive model is primary; AA-vs-GG read off as 2 × β. A 2-df genotypic model provides the sensitivity check.
- Full panel run (test_mode disabled).
