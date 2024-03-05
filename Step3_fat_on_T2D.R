library(tidyverse)
library(openxlsx)
library(TwoSampleMR)
#####################################################
exp_fat <- extract_instruments(outcomes = c("ebi-a-GCST90016675", "ebi-a-GCST90016673", "ebi-a-GCST90016672"), p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 10000)
out_fat_T2D <- extract_outcome_data(snps = exp_drugLDL$SNP, outcomes = "ebi-a-GCST006867", proxies = T, maf_threshold = 0.1, access_token = NULL)
data_fat_T2D <- harmonise_data(exposure_dat = exp_fat, outcome_dat = out_fat_T2D, action = 2)
res_fat_T2D <- mr(data_fat_T2D, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

#######################################################
# local clump
mr_clump_local2 <- function(x,clump_r2 = 0.01){
  
  d <- data.frame(rsid = x$SNP, pval = x$pval.exposure, 
                  id = x$id.exposure)
  out <- ld_clump_local(dat=d, clump_kb=10000, clump_r2=clump_r2, clump_p=1, bfile="./clump/EUR/EUR", plink_bin="./clump/plink_Windows.exe")
  
  keep <- paste(x$SNP, x$id.exposure) %in% paste(out$rsid, out$id)
  return(x[keep, ])
}

###################################################
file_fat <- list.files("./GWAS_data/fat/all", pattern = "_bgen_stats.gz$", full.names = T)

exp_fat_2 <- pbapply::pblapply(file_fat, function(x){
  
  pheno <- sub(".*/([^/]+)_bgen_stats\\.gz", '\\1', x)
  pheno <- sub("^\\d+_", "", pheno)
  
  exp <- read_exposure_data(filename = x,
                            type = "exposure",
                            sep = "\t",
                            snp_col = "SNP",
                            beta_col = "BETA",
                            se_col = "SE",
                            eaf_col = "A1FREQ",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            pval_col = "P_LINREG",
                            chr_col = "CHR",
                            pos_col = "BP") %>%
    filter(pval.exposure < 5e-08)
    exp$exposure <- pheno
    exp$id.exposure <- pheno
    exp_0.2 <- mr_clump_local2(exp,clump_r2 = 0.2)
    return(exp_0.2)
  
})

exp_fat_2 <- Reduce(rbind, exp_fat_2)
######################################################
res_fat_T2D_2 <- harmonise_data(exposure_dat = exp_fat_2, 
                                     outcome_dat = extract_outcome_data(snps = exp_fat_2$SNP, outcomes = "ebi-a-GCST006867", proxies = T, maf_threshold = 0.1, access_token = NULL), 
                                     action = 2) %>%
  mr(method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", 
                     "mr_weighted_mode"))