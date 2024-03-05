library(tidyverse)
library(openxlsx)
library(TwoSampleMR)
#####################################################
# Reduced LDL-C
exp <- list(
  LDLR = list(GENE = "LDLR", CHR= 19, START = 11200139, STOP = 11244496),
  HMGCR = list(GENE = "HMGCR", CHR=5, START = 74632993, STOP = 74657941),
  NPC1L1 = list(GENE = "NPC1L1",CHR=7, START =44552134, STOP = 44580929),
  PCSK9 = list(GENE = "PCSK9", CHR=1, START =55505221, STOP = 55530525),
  APOB = list(GENE = "APOB", CHR=2, START = 21224301, STOP = 21266945),
  ABCG5_ABCG8 = list(GENE="ABCG5/ABCG8", CHR = 2, START = 44039611, STOP = 44110127)
)
# Reduced TG
exp2 <- list(
  LPL = list(GENE = "LPL", CHR = 8, START = 19796764, STOP = 19824770),
  #PPARA=list(GENE="PPARA", CHR=22,START=46546429, STOP=46639652),
  ANGPTL3 = list(GENE = "ANGPTL3", CHR = 1, START = 63063191, STOP = 63071984),
  APOC3 = list(GENE = "APOC3", CHR = 11, START = 116700623, STOP = 116703788)
)
# Elevated HDL-C
exp3 <- list(
  CETP = list(GENE = "CETP", CHR = 16, START = 56995762, STOP = 57017757)
)
############################################################
exp_drugLDL <- lapply(exp, function(x){
    
    exp_drugLDL <- extract_instruments(outcomes = "ieu-a-300", p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                    filter(chr.exposure == x$CHR) %>%
                    filter(pos.exposure > (x$START - (100 * 1000))) %>%
                    filter(pos.exposure < (x$STOP + (100 * 1000))) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 250)
    exp_drugLDL$exposure <- x$GENE
    return(exp_drugLDL)
})

exp_drugTG <- lapply(exp2, function(x){
    
    exp_drugTG <- extract_instruments(outcomes = "ieu-a-302", p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                    filter(chr.exposure == x$CHR) %>%
                    filter(pos.exposure > (x$START - (100 * 1000))) %>%
                    filter(pos.exposure < (x$STOP + (100 * 1000))) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 250)
    exp_drugTG$exposure <- x$GENE
    return(exp_drugTG)
})

exp_drugHDL <- lapply(exp3, function(x){
    
    exp_drugHDL <- extract_instruments(outcomes = "ieu-a-299", p1 = 5e-08, clump = FALSE, access_token = NULL) %>%
                    filter(chr.exposure == x$CHR) %>%
                    filter(pos.exposure > (x$START - (100 * 1000))) %>%
                    filter(pos.exposure < (x$STOP + (100 * 1000))) %>%
                      clump_data(clump_r2 = 0.2, clump_kb = 250)
    exp_drugHDL$exposure <- x$GENE
    return(exp_drugHDL)
})


exp_drug <- Reduce(rbind, c(exp_drugLDL, exp_drugTG, exp_drugHDL))

############################################################

out_drug_T2D <- extract_outcome_data(snps = exp_drug$SNP, outcomes = "ebi-a-GCST006867", proxies = T, maf_threshold = 0.1, access_token = NULL)
data_drug_T2D <- harmonise_data(exposure_dat = exp_drug, outcome_dat = out_drug_T2D, action = 2)
res_drug_T2D <- mr(data_drugLDL_T2D, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

res_drug_T2D <- res_drug_T2D %>%
  mutate(b = ifelse(exposure == "CETP", b, -1*b)) %>%
  mutate(OR = exp(b),
         CImin = exp(b-1.96*se),
         CImax = exp(b+1.96*se))


