library(tidyverse)
#######################################
# Delta Method

Delta_method <- function(b_exp_M, se_exp_M,
                         b_M_out, se_M_out,
                         b_exp_out, se_exp_out){
    # indect effect
    EO <- b_exp_M * b_M_out
    prop_mediated <- EO/b_exp_out
    S <- sqrt(b_exp_M^2 * se_M_out^2 + b_M_out^2 * se_exp_M^2)
    indirect_se <- S
    Z <- (EO)/S
    P <- pnorm(q = abs(Z), lower.tail = F) * 2
 
    D_EO <- eval(stats::D(expression(EO/b_exp_out), "EO"))
    D_b_exp_out <- eval(stats::D(expression(EO/b_exp_out), "b_exp_out"))
    prop_mediated_se <- sqrt((D_EO^2) * indirect_se^2 + (D_b_exp_out^2) * se_exp_out^2)

res <- list(
    b_indect_effect = EO,                         # Indect effect
    se_indect_effect =  indirect_se,              # SE of indect effect
    P_indect_effect = P,                          # P-value of indect effect
    proportion_mediated = prop_mediated,          # Proportion of the mediated effect
    se_proportion_mediated = prop_mediated_se     # SE of mediated proportion
)
return(res)

                         }