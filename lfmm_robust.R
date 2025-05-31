###############################################
### Running a robust LFMM Credit:           ###
### This is based on the lfmm2 tutorial     ###
### Modified for robust analysis by         ###
### Nikos Tourvas                           ###
### nikostourvas@gmail.com                  ###
###############################################

#### Load packages ####
#install.packages("BiocManager")
#BiocManager::install("LEA")
library(LEA) # https://bcm-uga.github.io/lea/index.html
#install.packages("viridis")
library(viridis)
#install.packages("scales")
library(scales)
#BiocManager::install("qvalue")
library(qvalue)
#install.packages("MASS")
library(MASS) # for robust regression
rm(list = ls())

# load simulated data
data("offset_example")
# 200 diploid individuals genotyped at 510 SNP
Y <- offset_example$geno
# 4 environmental variables
X <- offset_example$env
mod.lfmm2 <- lfmm2(input = Y, env = X[,2], K = 2)

stat.lfmm2 <- lfmm2.test(mod.lfmm2, input=Y, env = X[,2], full=TRUE, 
                         genomic.control = TRUE)

# I will recreate lfmm2.test function here

  # Suppose:
  # Y      : n x M genotype matrix (rows: individuals, columns: SNPs)
  # X[,2]  : the environmental variable of interest (a vector of length n)
  # mod.lfmm2@U : latent factors from the LFMM2 model (n x K matrix)
  
  # Combine the environmental variable and latent factors into one data frame:
  covariates <- data.frame(env = X[,2], mod.lfmm2@U)
  # Give names to the latent factors:
  colnames(covariates)[-1] <- paste0("LF", 1:ncol(mod.lfmm2@U))
  
  # Initialize matrices to store the statistics.
  M <- ncol(Y)
  pvals <- numeric(M)
  zvals <- numeric(M)
  
  # Loop over each SNP (each column of Y)
  for (j in 1:M) {
    # Fit the regression model for SNP j
    fit <- lm(Y[, j] ~ ., data = covariates)
    
    # Get the summary and extract the coefficient for the 'env' predictor.
    sfit <- summary(fit)
    
    # In the design matrix, the first coefficient is the intercept.
    # The second coefficient corresponds to env.
    # Extract the t value and p value for the environmental variable.
    zvals[j] <- sfit$coefficients["env", "t value"]
    pvals[j] <- sfit$coefficients["env", "Pr(>|t|)"]
  }
  
  # Calculate GIF using the squared z-values
  gif <- median(zvals^2) / qchisq(0.5, df = 1, lower.tail = FALSE)
  
  # Adjusted p-values:
  adjusted_pvals <- pchisq(zvals^2 / gif, df = 1, lower.tail = FALSE)
  
  # Non GIF-adjusted p-values:
  raw_pvals <- pchisq(zvals^2, df = 1, lower.tail = FALSE)
  
  # Reproduced results
  rep_results <- list(pvalues = adjusted_pvals, 
                      zscores = zvals, 
                      gif = gif)
  
  # Does the re-creation of lfmm2.test work properly?
  # Compare p-values. They should be almost identical.
  cor.test(stat.lfmm2$pvalues, rep_results$pvalues, method = "pearson")
  plot(stat.lfmm2$pvalues, rep_results$pvalues, 
       xlab = "p-values from LFMM2", ylab = "p-values from lm",
       main = "Comparison of p-values")

# Implement robust regression with the rlm function from MASS package
  # Initialize matrices to store the statistics.
  rm(M, pvals, zvals, gif, adjusted_pvals)
  M <- ncol(Y)
  pvals <- numeric(M)
  zvals <- numeric(M)
  
  # Loop over each SNP (each column of Y)
  for (j in 1:M) {
    # Fit the regression model for SNP j using robust regression
    fit <- rlm(Y[, j] ~ ., data = covariates,
               maxit=1000)
    
    # Get the summary and extract the coefficient for the 'env' predictor.
    sfit <- summary(fit)
    
    # In the design matrix, the first coefficient is the intercept.
    # The second coefficient corresponds to env.
    # Extract the t value for the environmental variable.
    zvals[j] <- sfit$coefficients["env", "t value"]
  }
  
  # Calculate GIF using the squared z-values
  gif <- median(zvals^2) / qchisq(0.5, df = 1, lower.tail = FALSE)
  # Adjusted p-values:
  adjusted_pvals_robust <- pchisq(zvals^2 / gif, df = 1, lower.tail = FALSE)
  # Non GIF-adjusted p-values:
  raw_pvals <- pchisq(zvals^2, df = 1, lower.tail = FALSE)
  # Reproduced results
  rep_results_robust <- list(pvalues = adjusted_pvals_robust, 
                             zscores = zvals, 
                             gif = gif)
  
# Compare p-values from two methods
  cor.test(rep_results$pvalues, rep_results_robust$pvalues, method = "pearson")
  plot(rep_results$pvalues, rep_results_robust$pvalues, 
       xlab = "p-values from LFMM2", ylab = "p-values from robust regression",
       main = "Comparison of p-values")
  
# Implement robust regression with robustbase::lmrob
# This didn't work for me, but I will leave it here for reference.
  # library(robustbase)
  # # Initialize matrices to store the statistics.
  # rm(M, pvals, zvals, gif, adjusted_pvals)
  # M <- ncol(Y)
  # pvals <- numeric(M)
  # zvals <- numeric(M)
  # 
  # # Loop over each SNP (each column of Y)
  # for (j in 1:M) {
  #   # Fit the regression model for SNP j using robust regression
  #   fit <- lmrob(Y[, j] ~ ., data = covariates, setting="KS2014")
  #   
  #   # Get the summary and extract the coefficient for the 'env' predictor.
  #   sfit <- summary(fit)
  #   
  #   # In the design matrix, the first coefficient is the intercept.
  #   # The second coefficient corresponds to env.
  #   # Extract the t value and p value for the environmental variable.
  #   zvals[j] <- sfit$coefficients["env", "t value"]
  #   pvals[j] <- sfit$coefficients["env", "Pr(>|t|)"]
  # }
  # 
  # # Calculate GIF using the squared z-values
  # gif <- median(zvals^2) / qchisq(0.5, df = 1, lower.tail = FALSE)
  # # Adjusted p-values:
  # adjusted_pvals <- pchisq(zvals^2 / gif, df = 1, lower.tail = FALSE)
  # # Reproduced results
  # rep_results_robustbase <- list(pvalues = adjusted_pvals, 
  #                 zscores = zvals, 
  #                 gif = gif)
  