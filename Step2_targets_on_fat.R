library(tidyverse)
library(openxlsx)
library(TwoSampleMR)
####################################################

out_drug_fat <- extract_outcome_data(snps = exp_drug$SNP, outcomes = c("ebi-a-GCST90016675", "ebi-a-GCST90016673", "ebi-a-GCST90016672"), proxies = T, maf_threshold = 0.1, access_token = NULL)
data_drug_fat <- harmonise_data(exposure_dat = exp_drug, outcome_dat = out_drug_fat, action = 2)
res_drug_fat <- mr(data_drug_fat, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

res_drug_fat <- res_drug_fat %>%
  mutate(b = ifelse(exposure == "CETP", b, -1*b)) %>%
  mutate(CImin = b-1.96*se,
         CImax = b+1.96*se)

###################################################
file_fat <- list.files("./GWAS_data/fat/all", pattern = "_bgen_stats.gz$", full.names = T)
res_drug_fat2 <- pbapply::pblapply(file_fat, function(x){
  
  pheno <- sub(".*/([^/]+)_bgen_stats\\.gz", '\\1', x)
  pheno <- sub("^\\d+_", "", pheno)
  
  out <- read_outcome_data(filename = x,
                      type = "outcome",
                      sep = "\t",
                      snps = NULL,
                      snp_col = "SNP",
                      beta_col = "BETA",
                      se_col = "SE",
                      eaf_col = "A1FREQ",
                      effect_allele_col = "ALLELE1",
                      other_allele_col = "ALLELE0",
                      pval_col = "P_LINREG",
                      chr_col = "CHR",
                      pos_col = "BP",
                      log_pval = FALSE)
  out$outcome <- pheno
  out$id.outcome <- pheno
  res <- harmonise_data(exposure_dat = exp_drug, 
                            outcome_dat = out, 
                            action = 2) %>%
    mr(method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", 
                       "mr_weighted_mode"))
  
  return(res)
  
})

res_drug_fat2 <- Reduce(rbind, res_drug_fat2)
res_drug_fat2 <- res_drug_fat2 %>%
  mutate(b = ifelse(exposure == "CETP", b, -1*b)) %>%
  mutate(CImin = b-1.96*se,
         CImax = b+1.96*se)