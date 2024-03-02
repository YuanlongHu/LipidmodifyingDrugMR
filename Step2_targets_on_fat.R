library(tidyverse)
library(openxlsx)
library(TwoSampleMR)
#####################################################
res_drugLDL_fat <- lapply(exp, function(x){
    
    exp_drugLDL <- extract_instruments(outcomes = "ieu-a-300", p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                    filter(chr.exposure == x$CHR) %>%
                    filter(pos.exposure > (x$START - (100 * 1000))) %>%
                    filter(pos.exposure < (x$STOP + (100 * 1000))) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 250)

    out_drugLDL_fat <- extract_outcome_data(snps = exp_drugLDL$SNP, outcomes = c("ebi-a-GCST90016675", "ebi-a-GCST90016673", "ebi-a-GCST90016672"), proxies = T, maf_threshold = 0.1, access_token = NULL)
    data_drugLDL_fat <- harmonise_data(exposure_dat = exp_drugLDL, outcome_dat = out_drugLDL_fat, action = 2)
    res <- mr(data_drugLDL_fat, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    return(res)
})

res_drugTG_fat <- lapply(exp2, function(x){
    
    exp_drugTG <- extract_instruments(outcomes = "ieu-a-302", p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                    filter(chr.exposure == x$CHR) %>%
                    filter(pos.exposure > (x$START - (100 * 1000))) %>%
                    filter(pos.exposure < (x$STOP + (100 * 1000))) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 250)

    out_drugTG_fat <- extract_outcome_data(snps = exp_drugTG$SNP, outcomes = c("ebi-a-GCST90016675", "ebi-a-GCST90016673", "ebi-a-GCST90016672"), proxies = T, maf_threshold = 0.1, access_token = NULL)
    data_drugTG_fat <- harmonise_data(exposure_dat = exp_drugTG, outcome_dat = out_drugTG_fat, action = 2)
    res <- mr(data_drugTG_fat, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    return(res)
})

res_drugHDL_fat <- lapply(exp3, function(x){
    
    exp_drugTG <- extract_instruments(outcomes = "ieu-a-299", p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                    filter(chr.exposure == x$CHR) %>%
                    filter(pos.exposure > (x$START - (100 * 1000))) %>%
                    filter(pos.exposure < (x$STOP + (100 * 1000))) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 250)

    out_drugTG_fat <- extract_outcome_data(snps = exp_drugTG$SNP, outcomes = c("ebi-a-GCST90016675", "ebi-a-GCST90016673", "ebi-a-GCST90016672"), proxies = T, maf_threshold = 0.1, access_token = NULL)
    data_drugTG_fat <- harmonise_data(exposure_dat = exp_drugTG, outcome_dat = out_drugTG_fat, action = 2)
    res <- mr(data_drugTG_fat, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    return(res)
})