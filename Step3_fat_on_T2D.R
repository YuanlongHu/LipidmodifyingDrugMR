library(tidyverse)
library(openxlsx)
library(TwoSampleMR)
#####################################################
exp_fat <- extract_instruments(outcomes = c("ebi-a-GCST90016675", "ebi-a-GCST90016673", "ebi-a-GCST90016672"), p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 10000)
out_fat_T2D <- extract_outcome_data(snps = exp_drugLDL$SNP, outcomes = "ebi-a-GCST006867", proxies = T, maf_threshold = 0.1, access_token = NULL)
data_fat_T2D <- harmonise_data(exposure_dat = exp_fat, outcome_dat = out_fat_T2D, action = 2)
res_fat_T2D <- mr(data_fat_T2D, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))