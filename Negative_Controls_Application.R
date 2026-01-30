
##########   NEGATIVE CONTROL REAL-DATA ANALYSIS   ##########

## This contains R code to run the real-data analysis 
## for the paper "Negative control outcomes and 
## selection bias in Mendelian randomization".

## We use UK Biobank data from in Schoeler et al. (2023).
## Details on how to access the data can be found in:
## https://github.com/TabeaSchoeler/TS2021_UKBBweighting

## Set working directory.
setwd("...")
options(digits = 4)

## Load R packages.
library(ieugwasr)
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(usethis)
#remotes::install_github("n-mounier/MRlap")
library(MRlap)
library(openxlsx)
library(gt)

## The TwoSampleMR package needs authentication - see:
## https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication

##################################################

##########   ORIGINAL ANALYSIS   ##########

## Perform the paper's main analyses.

## ---------- SET UP THE ANALYSIS ---------- ##

## Names and files of the various traits considered.
trait_files <- data.frame(matrix(0, 20, 2))
colnames(trait_files) <- c("Trait", "Filename")
trait_files$Trait <- c("BMI", "Height", "Alcohol", "Smoking", "Coffee", "Fruit", 
                       "Vegetable", "Education", "Physical-activity", "LDL", "SBP",
                       "Cancer", "Non-cancer", "Diabetes", "Depression", 
                       "Loneliness", "Reaction-time", "Risk-taking", "Insomnia", "Sex")
trait_files$Filename <- c("bmi", "height", "alcfrequencyweekly", "smoking_status", 
                          "coffee_intake", "fruit_intake_con", "vegetable_intake_con", 
                          "education_age", "physical_activity", "cholesterol_ldl", 
                          "systolic_blood_pressure", "cancer", "noncancer_illness_number",
                          "diabetes", "depression_anxiety", "loneliness", "reaction_time", 
                          "risk_taking", "insomnia_con", "sex")

## Create tables where results will be stored.
MR_results_standard <- data.frame(rep(trait_files$Trait[1:19], each = 20),
                                  rep(trait_files$Trait, times = 19),
                                  rep(NA, 19*20), matrix(NA, 19*20, 15))
colnames(MR_results_standard) <- c("Exposure", "Outcome", "nSNP", "IVW_beta", "IVW_se", "IVW_p",
                                   "Egger_beta", "Egger_se", "Egger_p", "wMed_beta", "wMed_se", "wMed_p",
                                   "Mode_beta", "Mode_se", "Mode_p", "wMode_beta", "wMode_se", "wMode_p")
MR_results_weighted <- MR_results_standard

## Store F statistics here.
MR_results_fstat <- data.frame("Exposure" = rep(trait_files$Trait[1:19], each = 20),
                               "Outcome" = rep(trait_files$Trait, times = 19),
                               "F_stat_standard" = rep(NA, 19*20), 
                               "F_stat_weighted" = rep(NA, 19*20))

## ---------- RUN THE ANALYSIS ---------- ##

## Loop through exposures and run MR.
for (I in 1:19) {
  
  ## Report progress.
  print(paste("Analyzing data on exposure ", trait_files$Trait[I], ".", sep = ""))
  
  ## Load standard GWAS data for the exposure.
  exp_file <- paste(trait_files$Filename[I], "_standard.tsv", sep = "")
  gwas <- read.table(file = exp_file, sep = '\t', header = TRUE)
  
  ## Get GWAS hits, convert to TwoSampleMR format and clump.
  exp_data <- gwas[which(gwas$p_value < 5e-8), ]
  exp_data_tsmr <- format_data(exp_data, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
  exp_data_tsmr <- clump_data(exp_data_tsmr)
  gwas_hits <- as.character(exp_data_tsmr$SNP)
  exp_data <- exp_data[which(exp_data$variant_id %in% as.character(exp_data_tsmr$SNP)), ]
  
  ## Load weighted GWAS data for the exposure.
  exp_file_w <- paste(trait_files$Filename[I], "_weighted.tsv", sep = "")
  gwas <- read.table(file = exp_file_w, sep = '\t', header = TRUE)
  
  ## Get summary stats for the standard GWAS hits from the weighted GWAS.
  exp_data_w <- gwas[which(gwas$variant_id %in% gwas_hits), ]
  exp_data_tsmr_w <- format_data(exp_data_w, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
  
  ## Loop through outcomes for the given exposure.
  for (J in 1:20) {
    
    ## No point running MR of a trait on itself.
    if (I == J) next
    
    ## Report progress.
    print(paste("Running MR of ", trait_files$Trait[I], " on ", trait_files$Trait[J], ".", sep = ""))
    
    ## ---------- STANDARD GWAS ---------- ##
    
    ## Load standard GWAS data for the outcome.
    out_file <- paste(trait_files$Filename[J], "_standard.tsv", sep = "")
    gwas <- read.table(file = out_file, sep = '\t', header = TRUE)
    
    ## Extract SNP-outcome estimates for the exposure SNPs.
    out_data <- gwas[which(gwas$variant_id %in% gwas_hits), ]
    out_data_tsmr <- format_data(out_data, type = "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
    
    ## Harmonize data.
    both_data <- harmonise_data(exp_data_tsmr, out_data_tsmr)
    
    ## Run the MR analysis.
    mr_results <- mr(both_data)
    Fstat_standard <- mean( both_data$beta.exposure^2 / both_data$se.exposure^2 )
    
    ## Store results
    MR_results_standard[20 * (I-1) + J, 3] <- mr_results[3, 6]
    MR_results_standard[20 * (I-1) + J, 4:6] <- mr_results[3, 7:9]
    MR_results_standard[20 * (I-1) + J, 7:9] <- mr_results[1, 7:9]
    MR_results_standard[20 * (I-1) + J, 10:12] <- mr_results[2, 7:9]
    MR_results_standard[20 * (I-1) + J, 13:15] <- mr_results[4, 7:9]
    MR_results_standard[20 * (I-1) + J, 16:18] <- mr_results[5, 7:9]
    MR_results_fstat[20 * (I-1) + J, 3] <- Fstat_standard
    
    ## ---------- WEIGHTED GWAS ---------- ##
    
    ## Load standard GWAS data for the outcome.
    out_file_w <- paste(trait_files$Filename[J], "_weighted.tsv", sep = "")
    gwas <- read.table(file = out_file_w, sep = '\t', header = TRUE)
    
    ## Extract SNP-outcome estimates for the exposure SNPs.
    out_data_w <- gwas[which(gwas$variant_id %in% gwas_hits), ]
    out_data_tsmr_w <- format_data(out_data_w, type = "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
    
    ## Harmonize data.
    both_data_w <- harmonise_data(exp_data_tsmr_w, out_data_tsmr_w)
    
    ## Run the MR analysis.
    mr_results_w <- mr(both_data_w)
    Fstat_weighted <- mean( both_data_w$beta.exposure^2 / both_data_w$se.exposure^2 )
    
    ## Store results
    MR_results_weighted[20 * (I-1) + J, 3] <- mr_results_w[3, 6]
    MR_results_weighted[20 * (I-1) + J, 4:6] <- mr_results_w[3, 7:9]
    MR_results_weighted[20 * (I-1) + J, 7:9] <- mr_results_w[1, 7:9]
    MR_results_weighted[20 * (I-1) + J, 10:12] <- mr_results_w[2, 7:9]
    MR_results_weighted[20 * (I-1) + J, 13:15] <- mr_results_w[4, 7:9]
    MR_results_weighted[20 * (I-1) + J, 16:18] <- mr_results_w[5, 7:9]
    MR_results_fstat[20 * (I-1) + J, 4] <- Fstat_weighted
    
  }
}

## ---------- ASSESS IVW RESULTS ---------- ##

## Load the data, if needed.
#load("Application_Schoeler.RData")

## Restrict to exposures with at least 10 SNPs.
MR_results_fstat10 <- MR_results_fstat[which(MR_results_standard$nSNP > 10), ]
MR_results_standard10 <- MR_results_standard[which(MR_results_standard$nSNP > 10), ]
MR_results_weighted10 <- MR_results_weighted[which(MR_results_weighted$nSNP > 10), ]

## Check instrument strength.
MR_results_fstat10
## All above 10, though some just above. Using the 
## weighted GWAS reduces the F stats, which is not 
## unexpected since instruments are selected in the 
## unweighted one. Still, no cause for concern.

## Effects of the various traits on sex (standard GWAS).
MR_results_standard10[1:13 * 19, 1:6]
## BMI, alcohol consumption, education and reaction time 
## are the variables with IVW estimates at p < 0.05.
## The Bonferroni p-value is 0.05/13 = 0.003846, which
## is surpassed only by alcohol consumption.

## Effects of the various traits on sex (weighted GWAS).
MR_results_weighted10[1:13 * 19, 1:6]
## Smoking is borderline significant but doesn't pass 
## the Bonferroni correction threshold.

## Which MR analyses have the biggest difference between 
## estimates in the standard and weighted GWAS?
MR_diff <- cbind(MR_results_standard10[, 1:6], MR_results_weighted10[, 4:6])
colnames(MR_diff)[4:9] <- c("Standard_IVW_beta", "Standard_IVW_se", "Standard_IVW_p", 
                            "Weighted_IVW_beta", "Weighted_IVW_se", "Weighted_IVW_p")
MR_diff$beta_diff <- MR_diff[, 4] - MR_diff[, 7]
MR_diff$ci_overlap <- ifelse(MR_diff$beta_diff > 0, 
                             MR_diff[, 4] + qnorm(0.025) * MR_diff[, 5] < MR_diff[, 7] + qnorm(0.975) * MR_diff[, 8],
                             MR_diff[, 4] + qnorm(0.975) * MR_diff[, 5] > MR_diff[, 7] + qnorm(0.025) * MR_diff[, 8])

## Which MR analyses have non-overlapping CIs?
MR_diff[which(MR_diff$ci_overlap == FALSE), ]
## Only four, and all have education as the exposure.

## As a sensitivity analysis, use 75% confidence intervals.
MR_diff$ci_75 <- ifelse(MR_diff$beta_diff > 0, 
                        MR_diff[, 4] + qnorm(0.125) * MR_diff[, 5] < MR_diff[, 7] + qnorm(0.875) * MR_diff[, 8],
                        MR_diff[, 4] + qnorm(0.875) * MR_diff[, 5] > MR_diff[, 7] + qnorm(0.125) * MR_diff[, 8])
MR_diff[which(MR_diff$ci_75 == FALSE), ]
## 11 exposure-outcome pairs, 10 of which have an exposure
## associated with sex in the negative control analysis.

## ---------- ASSESS MEDIAN RESULTS ---------- ##

## The TwoSampleMR package we have used also implements
## standard pleiotropy-robust methods. Here, we present
## results from the weighted median method, as an added
## diagnostic. 

## Variables affecting sex.
MR_results_standard10[1:13 * 19, c(1:3, 10:12)]
MR_results_weighted10[1:13 * 19, c(1:3, 10:12)]
## BMI, Alcohol and reaction time in the unweighted
## analysis, all the effects vanish after weighting.

## In practice, we won't want to use pleiotropy robust
## methods in the negative control analyses, as the
## whole point is to detect biases and pleiotropy
## adjustments will work against that point.

## Get significant differences.
MR_med <- cbind(MR_results_standard10[, c(1:3, 10:12)], MR_results_weighted10[, 10:12])
MR_med$beta_diff <- MR_med[, 4] - MR_med[, 7]
MR_med$ci_overlap <- ifelse(MR_med$beta_diff > 0, 
                            MR_med[, 4] + qnorm(0.025) * MR_med[, 5] < MR_med[, 7] + qnorm(0.975) * MR_med[, 8],
                            MR_med[, 4] + qnorm(0.975) * MR_med[, 5] > MR_med[, 7] + qnorm(0.025) * MR_med[, 8])
MR_med$ci_75 <- ifelse(MR_med$beta_diff > 0, 
                       MR_med[, 4] + qnorm(0.125) * MR_med[, 5] < MR_med[, 7] + qnorm(0.875) * MR_med[, 8],
                       MR_med[, 4] + qnorm(0.875) * MR_med[, 5] > MR_med[, 7] + qnorm(0.125) * MR_med[, 8])

## Significant differences (95%).
MR_med[which(MR_med$ci_overlap == FALSE), ]
## Two pairs only, BMI-Education and Diabetes-SBP.

## Significant differences (75%).
MR_med[which(MR_med$ci_75 == FALSE), ]
## 15 pairs with plausible bias (plus Insomnia-sex). 

## This is hard to interpret. Depending on the DAG, e.g.
## whether X -> S vs G -> S, pleiotropy adjustments may 
## or may not help with selection bias. Our main purpose 
## here was to explore the impact of selection bias, so 
## we did not include these results in the manuscript - 
## they are reported in the supplementary data file.

##################################################

##########   SENSITIVITY ANALYSIS 1   ##########

## In theory our negative control analysis can go wrong if
## SNPs affect selection directly (via pathways not associated
## with the MR exposure or outcome). Here, we identify SNPs 
## that hit 5e-8 in the participation GWAS and exclude them 
## from our analyses. 

## ---------- SET UP THE ANALYSIS ---------- ##

## This identifies the variants associated with participation.
gwas <- read.table("PS_standard.tsv", sep = "\t", header = TRUE)
ps_snps <- gwas[which(gwas$p_value < 5e-8), 1]
length(ps_snps)   ## 2474.
min(gwas$p_value)   ## 6.728e-15.

## Now repeat our analysis without the participation SNPs.

## Create tables where results will be stored.
MR_results_standard_noP <- data.frame(rep(trait_files$Trait[1:19], each = 20),
                                      rep(trait_files$Trait, times = 19),
                                      rep(NA, 19*20), matrix(NA, 19*20, 15))
colnames(MR_results_standard_noP) <- c("Exposure", "Outcome", "nSNP", "IVW_beta", "IVW_se", "IVW_p",
                                       "Egger_beta", "Egger_se", "Egger_p", "wMed_beta", "wMed_se", "wMed_p",
                                       "Mode_beta", "Mode_se", "Mode_p", "wMode_beta", "wMode_se", "wMode_p")
MR_results_weighted_noP <- MR_results_standard_noP

## Loop through exposures and run MR.
for (I in 1:19) {
  
  ## Report progress.
  print(paste("Analyzing data on exposure ", trait_files$Trait[I], ".", sep = ""))
  
  ## Load standard GWAS data for the exposure.
  exp_file <- paste(trait_files$Filename[I], "_standard.tsv", sep = "")
  gwas <- read.table(file = exp_file, sep = '\t', header = TRUE)
  
  ## Get GWAS hits, convert to TwoSampleMR format and clump.
  exp_data <- gwas[which(gwas$p_value < 5e-8), ]
  exp_data_tsmr <- format_data(exp_data, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
  exp_data_tsmr <- clump_data(exp_data_tsmr)
  gwas_hits <- as.character(exp_data_tsmr$SNP)
  exp_data <- exp_data[which(exp_data$variant_id %in% as.character(exp_data_tsmr$SNP)), ]
  
  ## Check if any of these SNPs affect participation.
  snps_to_exclude <- which(exp_data$variant_id %in% ps_snps)
  if (length(snps_to_exclude) == 0) {
    
    ## If nothing is excluded, there is no point rerunning.
    MR_results_standard_noP[(20 * (I-1) + 1):(20 * I), ] <- MR_results_standard[(20 * (I-1) + 1):(20 * I), ]
    MR_results_weighted_noP[(20 * (I-1) + 1):(20 * I), ] <- MR_results_weighted[(20 * (I-1) + 1):(20 * I), ]
    
  } else {
    
    ## Otherwise, exclude and repeat the analysis.
    exp_data <- exp_data[- snps_to_exclude, ]
    exp_data_tsmr <- exp_data_tsmr[- snps_to_exclude, ]
    
    ## Load weighted GWAS data for the exposure.
    exp_file_w <- paste(trait_files$Filename[I], "_weighted.tsv", sep = "")
    gwas <- read.table(file = exp_file_w, sep = '\t', header = TRUE)
    
    ## Get summary stats for the standard GWAS hits from the weighted GWAS.
    exp_data_w <- gwas[which(gwas$variant_id %in% gwas_hits), ]
    exp_data_tsmr_w <- format_data(exp_data_w, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
    
    ## Loop through outcomes for the given exposure.
    for (J in 1:20) {
      
      ## No point running MR of a trait on itself.
      if (I == J) next
      
      ## Report progress.
      print(paste("Running MR of ", trait_files$Trait[I], " on ", trait_files$Trait[J], ".", sep = ""))
      
      ## ---------- STANDARD GWAS ---------- ##
      
      ## Load standard GWAS data for the outcome.
      out_file <- paste(trait_files$Filename[J], "_standard.tsv", sep = "")
      gwas <- read.table(file = out_file, sep = '\t', header = TRUE)
      
      ## Extract SNP-outcome estimates for the exposure SNPs.
      out_data <- gwas[which(gwas$variant_id %in% gwas_hits), ]
      out_data_tsmr <- format_data(out_data, type = "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
      
      ## Harmonize data.
      both_data <- harmonise_data(exp_data_tsmr, out_data_tsmr)
      
      ## Run the MR analysis.
      mr_results <- mr(both_data)
      
      ## Store results
      MR_results_standard_noP[20 * (I-1) + J, 3] <- mr_results[3, 6]
      MR_results_standard_noP[20 * (I-1) + J, 4:6] <- mr_results[3, 7:9]
      MR_results_standard_noP[20 * (I-1) + J, 7:9] <- mr_results[1, 7:9]
      MR_results_standard_noP[20 * (I-1) + J, 10:12] <- mr_results[2, 7:9]
      MR_results_standard_noP[20 * (I-1) + J, 13:15] <- mr_results[4, 7:9]
      MR_results_standard_noP[20 * (I-1) + J, 16:18] <- mr_results[5, 7:9]
      
      ## ---------- WEIGHTED GWAS ---------- ##
      
      ## Load standard GWAS data for the outcome.
      out_file_w <- paste(trait_files$Filename[J], "_weighted.tsv", sep = "")
      gwas <- read.table(file = out_file_w, sep = '\t', header = TRUE)
      
      ## Extract SNP-outcome estimates for the exposure SNPs.
      out_data_w <- gwas[which(gwas$variant_id %in% gwas_hits), ]
      out_data_tsmr_w <- format_data(out_data_w, type = "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
      
      ## Harmonize data.
      both_data_w <- harmonise_data(exp_data_tsmr_w, out_data_tsmr_w)
      
      ## Run the MR analysis.
      mr_results_w <- mr(both_data_w)
      
      ## Store results
      MR_results_weighted_noP[20 * (I-1) + J, 3] <- mr_results_w[3, 6]
      MR_results_weighted_noP[20 * (I-1) + J, 4:6] <- mr_results_w[3, 7:9]
      MR_results_weighted_noP[20 * (I-1) + J, 7:9] <- mr_results_w[1, 7:9]
      MR_results_weighted_noP[20 * (I-1) + J, 10:12] <- mr_results_w[2, 7:9]
      MR_results_weighted_noP[20 * (I-1) + J, 13:15] <- mr_results_w[4, 7:9]
      MR_results_weighted_noP[20 * (I-1) + J, 16:18] <- mr_results_w[5, 7:9]
      
    }
  }
}

## ---------- ASSESS RESULTS ---------- ##

## Load the data, if needed.
#load("Application_Schoeler.RData")

## Restrict to exposures with at least 10 SNPs.
MR_results_standard_noP10 <- MR_results_standard_noP[which(MR_results_standard_noP$nSNP > 10), ]
MR_results_weighted_noP10 <- MR_results_weighted_noP[which(MR_results_weighted_noP$nSNP > 10), ]

## Effects of the various traits on sex.
MR_results_standard_noP10[1:13 * 19, 1:6]
MR_results_weighted_noP10[1:13 * 19, 1:6]
## Not much changes. The effects of the stronger exposures 
## (e.g. education) attenuate as some SNPs affecting 
## participation through them may have been discarded.

## Analyses with evidence of bias.
MR_noP <- cbind(MR_results_standard_noP10[, 1:6], MR_results_weighted_noP10[, 4:6])
MR_noP$beta_diff <- MR_noP[, 4] - MR_noP[, 7]
MR_noP$ci_overlap <- ifelse(MR_noP$beta_diff > 0, 
                            MR_noP[, 4] + qnorm(0.025) * MR_noP[, 5] < MR_noP[, 7] + qnorm(0.975) * MR_noP[, 8],
                            MR_noP[, 4] + qnorm(0.975) * MR_noP[, 5] > MR_noP[, 7] + qnorm(0.025) * MR_noP[, 8])
MR_noP$ci_75 <- ifelse(MR_noP$beta_diff > 0, 
                       MR_noP[, 4] + qnorm(0.125) * MR_noP[, 5] < MR_noP[, 7] + qnorm(0.875) * MR_noP[, 8],
                       MR_noP[, 4] + qnorm(0.875) * MR_noP[, 5] > MR_noP[, 7] + qnorm(0.125) * MR_noP[, 8])
MR_noP[which(MR_noP$ci_overlap == FALSE), ]
MR_noP[which(MR_noP$ci_75 == FALSE), ]
## Again, results are similar to the original analysis.

##################################################

##########   SENSITIVITY ANALYSIS 2   ##########

## The above analysis can exclude SNPs associated 
## with selection via X or Y, and can include SNPs
## weakly associated with selection. Here we use 
## the mtCoJo method to compute conditional SNP-sex
## associations given each pair of traits and then
## exclude SNPs associated with sex from MR analyses.

## ---------- MTCOJO FUNCTIONS ---------- ##

## These are obtained from the GitHub page:
## https://github.com/yangq001/conditional
JointSum=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0,adj_Y=1,lam=0){
  if(adj_Y==1){
    return(YYX4(B1,S1,B2,S2,N,XX,YY0,lam))
  }
  else{
    return(YYX5(B1,S1,B2,S2,N,XX,YY0))    
  }
}
YYX4=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0,lam=0){ #x first
  B1=as.matrix(B1)
  B2=as.matrix(B2)
  S1=as.matrix(S1)
  S2=as.matrix(S2)
  n=max(N)
  
  nrx=1
  for(rx in 1:nrx){
    ny=ncol(B2)
    xx=XX[rx,rx]
    s1=S1[rx]
    b1=B1[rx]
    S21=S2[rx,]
    B21=B2[rx,]
    if(rx==1){
      yy1=N[1,rx]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=N[-1,rx]/n*(n-1)*xx*S21^2+xx*B21^2
    }
    else{
      yy1=yy1+N[1,rx]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=YY2+N[-1,rx]/n*(n-1)*xx*S21^2+xx*B21^2   
    }
  }
  yy1=yy1/nrx
  YY2=YY2/nrx
  
  YYself=c(yy1,YY2)
  xxd=diag(XX)
  XY1=xxd*B1       #vector
  XY2=t(t(B2)*xxd)   #matrix
  
  YYcross=matrix(nrow=ny+1,ncol=ny+1)
  for(i1 in 1:(ny+1)){
    for(j1 in 1:(ny+1)){
      if(i1==j1){
        YYcross[i1,j1]=YYself[i1]
      }
      else{
        YYcross[i1,j1]=YY0[i1,j1]*sqrt(YYself[i1]*YYself[j1])
      }
    }
  }
  
  nx=ncol(XX)
  A=matrix(nrow=ny+nx,ncol=ny+nx)
  A[1:nx,1:nx]=XX
  A[(nx+1):(ny+nx),1:nx]=t(XY2)
  A[1:nx,(nx+1):(ny+nx)]=XY2
  
  A[(nx+1):(ny+nx),(nx+1):(ny+nx)]=YYcross[2:(ny+1),2:(ny+1)]
  
  B=c(XY1,YYcross[1,2:(ny+1)])
  
  O=matrix(0,ncol=(ny+nx),nrow=(ny+nx))
  for(i in 1:(ny+nx)){
    O[i,i]=1
  }
  A=(A+lam*O)/(1+lam)
  
  beta=solve(A)%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx+ny))
  se=sigma2[1,1]*solve(A)   #cov
  
  if(sum(diag(se)<=0)==0){
    pvalue=c()
    for(i in 1:length(beta)){
      pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
    }
  }else{
    diag(se)=abs(diag(se))
    pvalue=beta*0-1 
  }
  
  
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}
YYX5=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0){ #x first
  B1=as.matrix(B1)
  B2=as.matrix(B2)
  S1=as.matrix(S1)
  S2=as.matrix(S2)
  n=max(N)
  ny=ncol(B2)
  xx=XX[1,1]
  s1=S1[1]
  b1=B1[1]
  S21=S2[1,]
  B21=B2[1,]
  yy1=N[1,1]/n*(n-1)*xx*s1^2+xx*b1^2
  YY2=N[-1,1]/n*(n-1)*xx*S21^2+xx*B21^2
  YYself=c(yy1,YY2)
  xxd=diag(XX)
  XY1=xxd*B1       #vector
  XY2=t(t(B2)*xxd)   #matrix
  nx=ncol(XX)
  A=matrix(nrow=nx,ncol=nx)
  A[1:nx,1:nx]=XX
  
  B=c(XY1)
  
  beta=solve(A)%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx))
  se=sigma2[1,1]*solve(A)   #cov
  
  pvalue=c()
  for(i in 1:length(beta)){
    pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
  }
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}

## ---------- COMPUTE YY0 ---------- ##

## The method requires one to approximate the 
## covariance matrix of SNPs by computing Z-scores
## from a list of SNPs not associated with any trait
## and computing the covariances of those Z-scores.

## To get SNPs for the Z scores, compute minimum p-values
## across the 21 traits and select those below 0.05.
gwas <- read.table("PS_standard.tsv", sep = "\t", header = TRUE)
min_p <- data.frame(rsid = gwas$variant_id, p_value = gwas$p_value)

## Loop through exposures, get minimum p-values.
for (I in 1:20) {
  exp_file <- paste(trait_files$Filename[I], "_standard.tsv", sep = "")
  gwas <- read.table(file = exp_file, sep = '\t', header = TRUE)
  new_p <- gwas$p_value[match(min_p$rsid, gwas$variant_id)]
  to_replace <- which(min_p$p_value > new_p)
  min_p$p_value[to_replace] <- new_p[to_replace]
}

## Now clump the SNPs with p-values higher than 0.05.
## Use the participation GWAS as a starting point.
gwas <- read.table("PS_standard.tsv", sep = "\t", header = TRUE)
z_data <- gwas[which(min_p$p_value > 0.05), ]
z_data_tsmr <- format_data(z_data, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
z_data_tsmr <- clump_data(z_data_tsmr)
z_snps <- z_data_tsmr$SNP
rm(z_data)

## Use summary statistics for z_snps to compute Z scores.
Z_scores <- matrix(0, length(z_snps), 21)
colnames(Z_scores) <- c("Participation", trait_files[, 1])
ind <- match(z_snps, gwas$variant_id)
Z_scores[, 1] <- gwas[ind, "beta"] / gwas[ind, "standard_error"]
for (I in 1:20) {
  print(paste("Analyzing data on exposure ", trait_files$Trait[I], ".", sep = ""))
  exp_file <- paste(trait_files$Filename[I], "_standard.tsv", sep = "")
  gwas <- read.table(file = exp_file, sep = '\t', header = TRUE)
  ind <- match(z_snps, gwas$variant_id)
  Z_scores[, I + 1] <- gwas[ind, "beta"] / gwas[ind, "standard_error"]
}
## This gives a NA result for entry 1561, so remove it.
Z_scores <- Z_scores[-1561, ]

## Get the Z-score genetic correlations.
YY0_est <- cor(Z_scores)

## ---------- RE-RUN THE ANALYSIS ---------- ##

## Now repeat our analysis without the participation SNPs.

## Create tables where results will be stored.
MR_results_standard_mtc <- data.frame(rep(trait_files$Trait[1:19], each = 20),
                                      rep(trait_files$Trait, times = 19),
                                      rep(NA, 19*20), matrix(NA, 19*20, 15))
colnames(MR_results_standard_mtc) <- c("Exposure", "Outcome", "nSNP", "IVW_beta", "IVW_se", "IVW_p",
                                       "Egger_beta", "Egger_se", "Egger_p", "wMed_beta", "wMed_se", "wMed_p",
                                       "Mode_beta", "Mode_se", "Mode_p", "wMode_beta", "wMode_se", "wMode_p")
MR_results_weighted_mtc <- MR_results_standard_mtc

## Load the participation GWAS once - this will take up
## a lot of memory but is faster than loading it per iteration.
prt_gwas_standard <- read.table(file = "PS_standard.tsv", header = TRUE)

## Loop through exposures and run MR.
for (I in 1:19) {
  
  ## Report progress.
  print(paste("Analyzing data on exposure ", trait_files$Trait[I], ".", sep = ""))
  
  ## Vegetable intake had no SNPs previously, so skip it.
  if (I == 7) next
  
  ## Load standard GWAS data for the exposure.
  exp_file <- paste(trait_files$Filename[I], "_standard.tsv", sep = "")
  gwas <- read.table(file = exp_file, sep = '\t', header = TRUE)
  
  ## Get GWAS hits, convert to TwoSampleMR format and clump.
  exp_data <- gwas[which(gwas$p_value < 5e-8), ]
  exp_data_tsmr <- format_data(exp_data, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
  exp_data_tsmr <- clump_data(exp_data_tsmr)
  gwas_hits <- as.character(exp_data_tsmr$SNP)
  exp_data <- exp_data[which(exp_data$variant_id %in% as.character(exp_data_tsmr$SNP)), ]
  
  ## Load weighted GWAS data for the exposure.
  exp_file_w <- paste(trait_files$Filename[I], "_weighted.tsv", sep = "")
  gwas <- read.table(file = exp_file_w, sep = '\t', header = TRUE)
  
  ## Get summary stats for the standard GWAS hits from the weighted GWAS.
  exp_data_w <- gwas[which(gwas$variant_id %in% gwas_hits), ]
  exp_data_tsmr_w <- format_data(exp_data_w, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
  
  ## Loop through outcomes for the given exposure.
  for (J in 1:20) {
    
    ## No point running MR of a trait on itself.
    if (I == J) next
    
    ## Report progress.
    print(paste("Running MR of ", trait_files$Trait[I], " on ", trait_files$Trait[J], ".", sep = ""))
    
    ## ---------- STANDARD GWAS ---------- ##
    
    ## Load standard GWAS data for the outcome.
    out_file <- paste(trait_files$Filename[J], "_standard.tsv", sep = "")
    gwas <- read.table(file = out_file, sep = '\t', header = TRUE)
    
    ## Extract SNP-outcome estimates for the exposure SNPs.
    out_data <- gwas[which(gwas$variant_id %in% gwas_hits), ]
    out_data_tsmr <- format_data(out_data, type = "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
    
    ## Harmonize data.
    both_data <- harmonise_data(exp_data_tsmr, out_data_tsmr)
    
    ## Get participation GWAS data for the MR SNPs.
    prt_sumstats <- prt_gwas_standard[which(prt_gwas_standard$variant_id %in% both_data$SNP), ]
    prt_sumstats <- prt_sumstats[order(prt_sumstats$variant_id), ]
    
    ## Run mtCoJo, identify SNPs affecting participation
    ## independently of the exposure and outcome.
    mtcojo <- JointSum(B1 = matrix(prt_sumstats$beta, nrow(both_data), 1), S1 = matrix(prt_sumstats$standard_error, nrow(both_data), 1),
                       B2 = both_data[, c("beta.exposure", "beta.outcome")], S2 = both_data[, c("se.exposure", "se.outcome")], 
                       N = cbind(prt_sumstats$N, both_data[, c("samplesize.exposure", "samplesize.outcome")]), 
                       XX = diag(nrow(both_data)), YY0 = YY0_est[c(1, I+1, J+1), c(1, I+1, J+1)])
    
    ## Identify SNPs that have a corrected p-value for 
    ## participation below the multiple testing threshold
    ## and exclude them from the MR analysis.
    both_data <- both_data[- which(mtcojo$pvalue < 0.05 / nrow(both_data)), ]
    
    ## Run the MR analysis.
    mr_results <- mr(both_data)
    
    ## Store results
    MR_results_standard_mtc[20 * (I-1) + J, 3] <- mr_results[3, 6]
    MR_results_standard_mtc[20 * (I-1) + J, 4:6] <- mr_results[3, 7:9]
    MR_results_standard_mtc[20 * (I-1) + J, 7:9] <- mr_results[1, 7:9]
    MR_results_standard_mtc[20 * (I-1) + J, 10:12] <- mr_results[2, 7:9]
    MR_results_standard_mtc[20 * (I-1) + J, 13:15] <- mr_results[4, 7:9]
    MR_results_standard_mtc[20 * (I-1) + J, 16:18] <- mr_results[5, 7:9]
    
    ## ---------- WEIGHTED GWAS ---------- ##
    
    ## Load standard GWAS data for the outcome.
    out_file_w <- paste(trait_files$Filename[J], "_weighted.tsv", sep = "")
    gwas <- read.table(file = out_file_w, sep = '\t', header = TRUE)
    
    ## Extract SNP-outcome estimates for the exposure SNPs.
    out_data_w <- gwas[which(gwas$variant_id %in% gwas_hits), ]
    out_data_tsmr_w <- format_data(out_data_w, type = "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error", eaf_col = "EAF", pval_col = "p_value", samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
    
    ## We don't apply mtCoJo here, as the SNPs selected 
    ## for the weighted analysis are the same as those for
    ## the unweighted analysis (after mtCoJo adjustment).
    
    ## Harmonize data.
    both_data_w <- harmonise_data(exp_data_tsmr_w, out_data_tsmr_w)
    
    ## Run the MR analysis.
    mr_results_w <- mr(both_data_w)
    
    ## Store results
    MR_results_weighted_mtc[20 * (I-1) + J, 3] <- mr_results_w[3, 6]
    MR_results_weighted_mtc[20 * (I-1) + J, 4:6] <- mr_results_w[3, 7:9]
    MR_results_weighted_mtc[20 * (I-1) + J, 7:9] <- mr_results_w[1, 7:9]
    MR_results_weighted_mtc[20 * (I-1) + J, 10:12] <- mr_results_w[2, 7:9]
    MR_results_weighted_mtc[20 * (I-1) + J, 13:15] <- mr_results_w[4, 7:9]
    MR_results_weighted_mtc[20 * (I-1) + J, 16:18] <- mr_results_w[5, 7:9]
    
  }
}

## Discard the participation GWAS.
rm(prt_gwas_standard)

## ---------- ASSESS RESULTS ---------- ##

## Load the data, if needed.
#load("Application_Schoeler.RData")

## Restrict to exposures with at least 10 SNPs.
MR_results_standard_mtc10 <- MR_results_standard_mtc[which(MR_results_standard_noP$nSNP > 10), ]
MR_results_weighted_mtc10 <- MR_results_weighted_mtc[which(MR_results_weighted_noP$nSNP > 10), ]
## Fruit intake had 10+ SNPs in previous analyses but just
## 10 SNPs here - we keep it in the analysis for comparison.

## Effects of the various traits on sex.
MR_results_standard_mtc10[1:13 * 19, 1:6]
MR_results_weighted_mtc10[1:13 * 19, 1:6]
## Again, not much changes. Here education and reaction time
## remain significant, while alcohol and BMI see their effects
## attenuate and their p-values are just above 0.05.

## Analyses with evidence of bias.
MR_mtc <- cbind(MR_results_standard_mtc10[, 1:6], MR_results_weighted_mtc10[, 4:6])
MR_mtc$beta_diff <- MR_mtc[, 4] - MR_mtc[, 7]
MR_mtc$ci_overlap <- ifelse(MR_mtc$beta_diff > 0, 
                            MR_mtc[, 4] + qnorm(0.025) * MR_mtc[, 5] < MR_mtc[, 7] + qnorm(0.975) * MR_mtc[, 8],
                            MR_mtc[, 4] + qnorm(0.975) * MR_mtc[, 5] > MR_mtc[, 7] + qnorm(0.025) * MR_mtc[, 8])
MR_mtc$ci_75 <- ifelse(MR_mtc$beta_diff > 0, 
                       MR_mtc[, 4] + qnorm(0.125) * MR_mtc[, 5] < MR_mtc[, 7] + qnorm(0.875) * MR_mtc[, 8],
                       MR_mtc[, 4] + qnorm(0.875) * MR_mtc[, 5] > MR_mtc[, 7] + qnorm(0.125) * MR_mtc[, 8])
MR_mtc[which(MR_mtc$ci_overlap == FALSE), ]
MR_mtc[which(MR_mtc$ci_75 == FALSE), ]
## Again, results are similar to the original analysis - 
## mainly analyses with education and BMI as exposures 
## are identified as biased.

##################################################

##########   MR-LAP IMPLEMENTATION   ##########

## We also implement MRlap to adjust for weak instrument bias,
## winner's curse bias and sample overlap in the analysis.

## As with the Median, this will be less clear to interpret,
## but we do it as these biases are likely in our application.

## ---------- PREPARE THE ANALYSIS ---------- ##

## Create tables where results will be stored.
MR_results_standard_mrlap <- data.frame(rep(trait_files$Trait[1:19], each = 20),
                                        rep(trait_files$Trait, times = 19),
                                        rep(NA, 19*20), matrix(NA, 19*20, 6))
colnames(MR_results_standard_mrlap) <- c("Exposure", "Outcome", "nSNP", "MRlap_obs_beta", "MRlap_obs_se", "MRlap_obs_p",
                                         "MRlap_corr_beta", "MRlap_corr_se", "MRlap_corr_p")
MR_results_weighted_mrlap <- MR_results_standard_mrlap

## Loop through exposures and run MR.
for (I in 1:19) {
  
  ## Report progress.
  print(paste("Analyzing data on exposure ", trait_files$Trait[I], ".", sep = ""))
  
  ## Load standard GWAS data for the exposure.
  exp_file <- paste(trait_files$Filename[I], "_standard.tsv", sep = "")
  exp_gwas <- read.table(file = exp_file, sep = '\t', header = TRUE)
  colnames(exp_gwas) <- c("snp", "chr", "pos", "a1", "a0", "beta", "se", "pval", "N", "eaf")
  
  ## Load weighted GWAS data for the exposure.
  exp_file_w <- paste(trait_files$Filename[I], "_weighted.tsv", sep = "")
  exp_gwas_w <- read.table(file = exp_file_w, sep = '\t', header = TRUE)
  colnames(exp_gwas_w) <- c("snp", "chr", "pos", "a1", "a0", "beta", "se", "pval", "N", "eaf")
  
  ## Loop through outcomes for the given exposure.
  for (J in 1:20) {
    
    ## No point running MR of a trait on itself.
    if (I == J) next
    
    ## Report progress.
    print(paste("Running MR of ", trait_files$Trait[I], " on ", trait_files$Trait[J], ".", sep = ""))
    
    ## ---------- STANDARD GWAS ---------- ##
    
    ## Load standard GWAS data for the outcome.
    out_file <- paste(trait_files$Filename[J], "_standard.tsv", sep = "")
    out_gwas <- read.table(file = out_file, sep = '\t', header = TRUE)
    colnames(out_gwas) <- c("snp", "chr", "pos", "a1", "a0", "beta", "se", "pval", "N", "eaf")
    
    ## Run MRlap.
    mrlap_results <- try(MRlap(exposure =  exp_gwas, exposure_name = "Exp", outcome = out_gwas, outcome_name = "Out", MR_threshold = 5e-8, MR_pruning_LD = 0.001,
                               ld = ".../eur_w_ld_chr",
                               hm3 = ".../w_hm3.snplist"),
                         silent = TRUE)
    
    ## Store results
    if(!(inherits(mrlap_results, "try-error"))) {
      MR_results_standard_mrlap[20 * (I-1) + J, 3] <- mrlap_results$MRcorrection$m_IVs
      MR_results_standard_mrlap[20 * (I-1) + J, 4] <- mrlap_results$MRcorrection$observed_effect
      MR_results_standard_mrlap[20 * (I-1) + J, 5] <- mrlap_results$MRcorrection$observed_effect_se
      MR_results_standard_mrlap[20 * (I-1) + J, 6] <- mrlap_results$MRcorrection$observed_effect_p
      MR_results_standard_mrlap[20 * (I-1) + J, 7] <- mrlap_results$MRcorrection$corrected_effect
      MR_results_standard_mrlap[20 * (I-1) + J, 8] <- mrlap_results$MRcorrection$corrected_effect_se
      MR_results_standard_mrlap[20 * (I-1) + J, 9] <- mrlap_results$MRcorrection$corrected_effect_p
    }
    
    ## ---------- WEIGHTED GWAS ---------- ##
    
    ## Load standard GWAS data for the outcome.
    out_file_w <- paste(trait_files$Filename[J], "_weighted.tsv", sep = "")
    out_gwas_w <- read.table(file = out_file_w, sep = '\t', header = TRUE)
    colnames(out_gwas_w) <- c("snp", "chr", "pos", "a1", "a0", "beta", "se", "pval", "N", "eaf")
    
    ## Run MRlap.
    mrlap_results_w <- try(MRlap(exposure =  exp_gwas_w, exposure_name = "Exp", outcome = out_gwas_w, outcome_name = "Out", MR_threshold = 5e-8, MR_pruning_LD = 0.001,
                                 ld = ".../eur_w_ld_chr",
                                 hm3 = ".../w_hm3.snplist"),
                           silent = TRUE)
    
    ## Store results
    if(!(inherits(mrlap_results_w, "try-error"))) {
      MR_results_weighted_mrlap[20 * (I-1) + J, 3] <- mrlap_results_w$MRcorrection$m_IVs
      MR_results_weighted_mrlap[20 * (I-1) + J, 4] <- mrlap_results_w$MRcorrection$observed_effect
      MR_results_weighted_mrlap[20 * (I-1) + J, 5] <- mrlap_results_w$MRcorrection$observed_effect_se
      MR_results_weighted_mrlap[20 * (I-1) + J, 6] <- mrlap_results_w$MRcorrection$observed_effect_p
      MR_results_weighted_mrlap[20 * (I-1) + J, 7] <- mrlap_results_w$MRcorrection$corrected_effect
      MR_results_weighted_mrlap[20 * (I-1) + J, 8] <- mrlap_results_w$MRcorrection$corrected_effect_se
      MR_results_weighted_mrlap[20 * (I-1) + J, 9] <- mrlap_results_w$MRcorrection$corrected_effect_p
    }
    
    ## Clear the results.
    rm(mrlap_results, mrlap_results_w)
    
  }
}
## Note: some of these runs return an error, especially for 
## exposures without many SNPs to instrument them.

## Save results.
save(MR_results_standard, MR_results_weighted, MR_results_fstat, trait_files,
     ps_snps, MR_results_standard_noP, MR_results_weighted_noP, 
     MR_results_standard_mtc, MR_results_weighted_mtc, z_data_tsmr, YY0_est,
     MR_results_standard_mrlap, MR_results_weighted_mrlap, file = "Negative_Controls_Application.RData")

## ---------- DIAGNOSTICS ---------- ##

## Load the data, if needed.
#load("Negative_Controls_Application.RData")

## Informal diagnostics: differences between unadjusted and
## adjusted MRlap estimates for each exposure-oucome pair.
max(MR_results_standard_mrlap$MRlap_obs_beta - MR_results_standard_mrlap$MRlap_corr_beta, na.rm = TRUE)
min(MR_results_standard_mrlap$MRlap_obs_beta - MR_results_standard_mrlap$MRlap_corr_beta, na.rm = TRUE)
mean(MR_results_standard_mrlap$MRlap_obs_beta - MR_results_standard_mrlap$MRlap_corr_beta, na.rm = TRUE)
mean(abs(MR_results_standard_mrlap$MRlap_obs_beta - MR_results_standard_mrlap$MRlap_corr_beta), na.rm = TRUE)
hist(MR_results_standard_mrlap$MRlap_obs_beta - MR_results_standard_mrlap$MRlap_corr_beta)


## Note that because MRlap does SNP selection internally,
## the standard and weighted GWAS are conducted using
## different numbers of SNPs (fewer SNPs for the weighted one).

## Restrict to exposures with at least 10 SNPs.
MR_results_standard_mrlap10 <- MR_results_standard_mrlap[which(MR_results_standard_mrlap$nSNP > 10), ]
MR_results_weighted_mrlap10 <- MR_results_weighted_mrlap[which(MR_results_standard_mrlap$nSNP > 10), ]
## Fruit intake and Insomnia have >10 SNPs in the original analysis
## but fewer than 10 SNPs when using MRlap. We leave these two
## exposures out of the MRlap results, especially since the method
## can fail with few SNPs. Some traits also have >10 SNPs in the
## unweighted and <10 in the weighted Mlap runs, but we keep those.

## Effects of the various traits on sex.
MR_results_standard_mrlap[which(MR_results_standard_mrlap$Outcome == "Sex"), ]
MR_results_weighted_mrlap[which(MR_results_weighted_mrlap$Outcome == "Sex"), ]
## Both the original analysis and the corrected one identify 
## BMI, alcohol and education as causes of participation.

## Analyses with evidence of bias - unadjusted analyses first.
MR_mrlap <- cbind(MR_results_standard_mrlap[, 1:6], MR_results_weighted_mrlap[, 3:6])
MR_mrlap <- MR_mrlap[which(MR_mrlap$nSNP >= 10), ]
MR_mrlap$beta_diff <- MR_mrlap[, 4] - MR_mrlap[, 8]
MR_mrlap$ci_overlap <- ifelse(MR_mrlap$beta_diff > 0, 
                              MR_mrlap[, 4] + qnorm(0.025) * MR_mrlap[, 5] < MR_mrlap[, 8] + qnorm(0.975) * MR_mrlap[, 9],
                              MR_mrlap[, 4] + qnorm(0.975) * MR_mrlap[, 5] > MR_mrlap[, 8] + qnorm(0.025) * MR_mrlap[, 9])
MR_mrlap$ci_75 <- ifelse(MR_mrlap$beta_diff > 0, 
                         MR_mrlap[, 4] + qnorm(0.125) * MR_mrlap[, 5] < MR_mrlap[, 8] + qnorm(0.875) * MR_mrlap[, 9],
                         MR_mrlap[, 4] + qnorm(0.875) * MR_mrlap[, 5] > MR_mrlap[, 8] + qnorm(0.125) * MR_mrlap[, 9])
MR_mrlap[which(MR_mrlap$ci_overlap == FALSE), ]
MR_mrlap[which(MR_mrlap$ci_75 == FALSE), ]
## Again, results are similar to the original analysis,
## with most biased pairs having BMI, alcohol or education
## as the exposure. Height also has quite a few analyses 
## suggestive of bias - this may be because it has by far 
## the most SNPs in the weighted GWAS, hence higher precision.

## Now repeat this for the adjusted analyses.
MR_mrlap2 <- cbind(MR_results_standard_mrlap[, c(1:3, 7:9)], MR_results_weighted_mrlap[, c(3, 7:9)])
MR_mrlap2 <- MR_mrlap2[which(MR_mrlap2$nSNP >= 10), ]
MR_mrlap2$beta_diff <- MR_mrlap2[, 4] - MR_mrlap2[, 8]
MR_mrlap2$ci_overlap <- ifelse(MR_mrlap2$beta_diff > 0, 
                               MR_mrlap2[, 4] + qnorm(0.025) * MR_mrlap2[, 5] < MR_mrlap2[, 8] + qnorm(0.975) * MR_mrlap2[, 9],
                               MR_mrlap2[, 4] + qnorm(0.975) * MR_mrlap2[, 5] > MR_mrlap2[, 8] + qnorm(0.025) * MR_mrlap2[, 9])
MR_mrlap2$ci_75 <- ifelse(MR_mrlap2$beta_diff > 0, 
                          MR_mrlap2[, 4] + qnorm(0.125) * MR_mrlap2[, 5] < MR_mrlap2[, 8] + qnorm(0.875) * MR_mrlap2[, 9],
                          MR_mrlap2[, 4] + qnorm(0.875) * MR_mrlap2[, 5] > MR_mrlap2[, 8] + qnorm(0.125) * MR_mrlap2[, 9])
MR_mrlap2[which(MR_mrlap2$ci_overlap == FALSE), ]
MR_mrlap2[which(MR_mrlap2$ci_75 == FALSE), ]
## Very little changes.

##################################################

##########   SUMMARIZE RESULTS   ##########

## Create Tables for the main part of the paper. This 
## includes two tables, variables' effects on sex and 
## pairwise MR results for pairs with evidence of bias.

## ---------- MAIN TABLES ---------- ##

## This computes confidence interval limits from a Table.
compute_ci <- function (Table, mean_col, se_col, alpha = 0.05) {
  ci_lower <- Table[, mean_col] + qnorm(alpha/2) * Table[, se_col]
  ci_upper <- Table[, mean_col] + qnorm(1 - alpha/2) * Table[, se_col]
  cbind(ci_lower, ci_upper)
}

## Effects of traits on sex in the negative control analysis.
Main_Table1 <- cbind(MR_results_standard10[which(MR_results_standard10$Outcome == "Sex"), 1:6],
                     MR_results_weighted10[which(MR_results_weighted10$Outcome == "Sex"), 4:6])
Main_Table1 <- cbind(Main_Table1, compute_ci(Main_Table1, 4, 5), compute_ci(Main_Table1, 7, 8))
colnames(Main_Table1) <- paste(c(rep("", 3), rep("Unw_", 3), rep("Wei_", 3), rep("Unw_", 2), rep("Wei_", 2)), colnames(Main_Table1), sep = "")
Main_Table1$Unw_CI <- paste("(", round(Main_Table1$Unw_ci_lower, 3), ", ", round(Main_Table1$Unw_ci_upper, 3), ")", sep = "")
Main_Table1$Wei_CI <- paste("(", round(Main_Table1$Wei_ci_lower, 3), ", ", round(Main_Table1$Wei_ci_upper, 3), ")", sep = "")
Main_Table1 <- Main_Table1[, c(1:5, 14, 6:8, 15, 9)]

## Pairs of traits with significant effects.
Main_Table2 <- MR_diff[which(MR_diff$ci_75 == FALSE), ]
colnames(Main_Table2) <- paste(c(rep("", 3), rep("Unw_", 3), rep("Wei_", 3), rep("", 3)), colnames(Main_Table2), sep = "")
Main_Table2 <- Main_Table2[order(Main_Table2$ci_overlap), ]
Main_Table2 <- cbind(Main_Table2, compute_ci(Main_Table2, 4, 5), compute_ci(Main_Table2, 4, 5, alpha = 0.25), compute_ci(Main_Table2, 7, 8), compute_ci(Main_Table2, 7, 8, alpha = 0.25))
Main_Table2$Sex_pval <- rep(NA, nrow(Main_Table2))
for (i in 1:nrow(Main_Table2)) Main_Table2$Sex_pval[i] <- Main_Table1$Unw_IVW_p[which(Main_Table1$Exposure == Main_Table2$Exposure[i])]
colnames(Main_Table2)[11:20] <- c("ci95_overlap", "ci75_overlap", "Unw_ci95_lower", "Unw_ci95_upper", "Unw_ci75_lower", "Unw_ci75_upper", 
                                  "Wei_ci95_lower", "Wei_ci95_upper", "Wei_ci75_lower", "Wei_ci75_upper")
Main_Table2$Unw_CI95 <- paste("(", round(Main_Table2$Unw_ci95_lower, 3), ", ", round(Main_Table2$Unw_ci95_upper, 3), ")", sep = "")
Main_Table2$Unw_CI75 <- paste("(", round(Main_Table2$Unw_ci75_lower, 3), ", ", round(Main_Table2$Unw_ci75_upper, 3), ")", sep = "")
Main_Table2$Wei_CI95 <- paste("(", round(Main_Table2$Wei_ci95_lower, 3), ", ", round(Main_Table2$Wei_ci95_upper, 3), ")", sep = "")
Main_Table2$Wei_CI75 <- paste("(", round(Main_Table2$Wei_ci75_lower, 3), ", ", round(Main_Table2$Wei_ci75_upper, 3), ")", sep = "")
Main_Table2 <- Main_Table2[, c(1:6, 22:23, 7:9, 24:25, 10:12, 21)]

## Save Tables for export in MS Word.
gt::gtsave(data = fmt_number(gt(Main_Table1), decimals = 3), file = "Main_Table1.docx")
gt::gtsave(data = fmt_number(gt(Main_Table2), decimals = 3), file = "Main_Table2.docx")

## ---------- SUPPLEMENTARY TABLES ---------- ##

## Create tables for the paper's supplement. This
## is the same two Tables as before but for the 
## analysis without the 2474 participation SNPs.

## Trait-sex effects in the analysis with 2474 SNPs excluded.
Supp_Table1 <- cbind(MR_results_standard_noP10[which(MR_results_standard_noP10$Outcome == "Sex"), 1:6],
                     MR_results_weighted_noP10[which(MR_results_weighted_noP10$Outcome == "Sex"), 4:6])
Supp_Table1 <- cbind(Supp_Table1, compute_ci(Supp_Table1, 4, 5), compute_ci(Supp_Table1, 7, 8))
colnames(Supp_Table1) <- paste(c(rep("", 3), rep("Unw_", 3), rep("Wei_", 3), rep("Unw_", 2), rep("Wei_", 2)), colnames(Supp_Table1), sep = "")
Supp_Table1$Unw_CI <- paste("(", round(Supp_Table1$Unw_ci_lower, 3), ", ", round(Supp_Table1$Unw_ci_upper, 3), ")", sep = "")
Supp_Table1$Wei_CI <- paste("(", round(Supp_Table1$Wei_ci_lower, 3), ", ", round(Supp_Table1$Wei_ci_upper, 3), ")", sep = "")
Supp_Table1 <- Supp_Table1[, c(1:5, 14, 6:8, 15, 9)]

## Trait-sex effects in the mtCOJO analysis.
Supp_Table2 <- cbind(MR_results_standard_mtc10[which(MR_results_standard_mtc10$Outcome == "Sex"), 1:6],
                     MR_results_weighted_mtc10[which(MR_results_weighted_mtc10$Outcome == "Sex"), 4:6])
Supp_Table2 <- cbind(Supp_Table2, compute_ci(Supp_Table2, 4, 5), compute_ci(Supp_Table2, 7, 8))
colnames(Supp_Table2) <- paste(c(rep("", 3), rep("Unw_", 3), rep("Wei_", 3), rep("Unw_", 2), rep("Wei_", 2)), colnames(Supp_Table2), sep = "")
Supp_Table2$Unw_CI <- paste("(", round(Supp_Table2$Unw_ci_lower, 3), ", ", round(Supp_Table2$Unw_ci_upper, 3), ")", sep = "")
Supp_Table2$Wei_CI <- paste("(", round(Supp_Table2$Wei_ci_lower, 3), ", ", round(Supp_Table2$Wei_ci_upper, 3), ")", sep = "")
Supp_Table2 <- Supp_Table2[, c(1:5, 14, 6:8, 15, 9)]

## Save Tables for export in MS Word.
gt::gtsave(data = fmt_number(gt(Supp_Table1), decimals = 3), file = "Supp_Table1.docx")
gt::gtsave(data = fmt_number(gt(Supp_Table2), decimals = 3), file = "Supp_Table2.docx")

## ---------- MRLAP TABLES ---------- ##

## Now create two more tables with the MR-lap results.

## Effects of traits on sex in the negative control analysis.
Supp_Table3 <- cbind(MR_results_standard_mrlap10[which(MR_results_standard_mrlap10$Outcome == "Sex"), c(1:3, 7:9)],
                     MR_results_weighted_mrlap10[which(MR_results_weighted_mrlap10$Outcome == "Sex"), c(3, 7:9)])
Supp_Table3 <- cbind(Supp_Table3, compute_ci(Supp_Table3, 4, 5), compute_ci(Supp_Table3, 8, 9))
colnames(Supp_Table3) <- paste(c(rep("", 2), rep("Unw_", 4), rep("Wei_", 4), rep("Unw_", 2), rep("Wei_", 2)), colnames(Supp_Table3), sep = "")
Supp_Table3$Unw_CI <- paste("(", round(Supp_Table3$Unw_ci_lower, 3), ", ", round(Supp_Table3$Unw_ci_upper, 3), ")", sep = "")
Supp_Table3$Wei_CI <- paste("(", round(Supp_Table3$Wei_ci_lower, 3), ", ", round(Supp_Table3$Wei_ci_upper, 3), ")", sep = "")
Supp_Table3 <- Supp_Table3[, c(1:5, 15, 6:9, 16, 10)]

## Pairs of traits with significant effects.
Supp_Table4 <- MR_mrlap2[which(MR_mrlap2$ci_75 == FALSE), ]
colnames(Supp_Table4) <- paste(c(rep("", 2), rep("Unw_", 4), rep("Wei_", 4), rep("", 3)), colnames(Supp_Table4), sep = "")
Supp_Table4 <- Supp_Table4[order(Supp_Table4$ci_overlap), ]
Supp_Table4 <- cbind(Supp_Table4, compute_ci(Supp_Table4, 4, 5), compute_ci(Supp_Table4, 4, 5, alpha = 0.25), compute_ci(Supp_Table4, 8, 9), compute_ci(Supp_Table4, 8, 9, alpha = 0.25))
Supp_Table4$Sex_pval <- rep(NA, nrow(Supp_Table4))
for (i in 1:nrow(Supp_Table4)) Supp_Table4$Sex_pval[i] <- Supp_Table3$Unw_MRlap_corr_p[which(Supp_Table3$Exposure == Supp_Table4$Exposure[i])]
colnames(Supp_Table4)[12:21] <- c("ci95_overlap", "ci75_overlap", "Unw_ci95_lower", "Unw_ci95_upper", "Unw_ci75_lower", "Unw_ci75_upper", 
                                  "Wei_ci95_lower", "Wei_ci95_upper", "Wei_ci75_lower", "Wei_ci75_upper")
Supp_Table4$Unw_CI95 <- paste("(", round(Supp_Table4$Unw_ci95_lower, 3), ", ", round(Supp_Table4$Unw_ci95_upper, 3), ")", sep = "")
Supp_Table4$Unw_CI75 <- paste("(", round(Supp_Table4$Unw_ci75_lower, 3), ", ", round(Supp_Table4$Unw_ci75_upper, 3), ")", sep = "")
Supp_Table4$Wei_CI95 <- paste("(", round(Supp_Table4$Wei_ci95_lower, 3), ", ", round(Supp_Table4$Wei_ci95_upper, 3), ")", sep = "")
Supp_Table4$Wei_CI75 <- paste("(", round(Supp_Table4$Wei_ci75_lower, 3), ", ", round(Supp_Table4$Wei_ci75_upper, 3), ")", sep = "")
Supp_Table4 <- Supp_Table4[, c(1:6, 23:24, 7:10, 25:26, 11:13, 22)]

## Save Tables for export in MS Word.
gt::gtsave(data = fmt_number(gt(Supp_Table3), decimals = 3), file = "Supp_Table3.docx")
gt::gtsave(data = fmt_number(gt(Supp_Table4), decimals = 3), file = "Supp_Table4.docx")

## ---------- EXCEL FILES ---------- ##

## The full list of results will be submitted as a
## supplementary data file. Here, create the file.

## Create an Excel file replica to store the data.
IAS <- createWorkbook()
addWorksheet(IAS, "MR Results (Unw)")
addWorksheet(IAS, "MR Results (Wei)")
addWorksheet(IAS, "Excl Sex SNPs (Unw)")
addWorksheet(IAS, "Excl Sex SNPs (Wei)")
addWorksheet(IAS, "mtCOJO (Unw)")
addWorksheet(IAS, "mtCOJO (Wei)")
addWorksheet(IAS, "MR-lap (Unw)")
addWorksheet(IAS, "MR-lap (Wei)")

## Add the Tables.
writeData(IAS, 1, MR_results_standard10)
writeData(IAS, 2, MR_results_weighted10)
writeData(IAS, 3, MR_results_standard_noP10)
writeData(IAS, 4, MR_results_weighted_noP10)
writeData(IAS, 5, MR_results_standard_mtc10)
writeData(IAS, 6, MR_results_weighted_mtc10)
writeData(IAS, 7, MR_results_standard_mrlap10)
writeData(IAS, 8, MR_results_weighted_mrlap10)

## Create the Excel file.
saveWorkbook(IAS, "Supplementary_Data.xlsx")

##################################################
