##################################################################
##               Obtain ATE estimates for Nethery               ##
##                    Use handful of methods                    ##
##################################################################

# Requires functions defined in 
source("data_gen_processes_Zhu2023.R")
library(bcf)
library(CausalGAM)
library(grf)

# Five methods
# Linear model, IPW-Linear model, BART, BART-PS, GRF
# Store the ATE.k, Psi.k, Bias.k, SD.k, COVER.k
L=250
raw.results <- matrix(NA, nrow = L*5, ncol=7) # include c, method, metrics,
colnames(raw.results) <- c("c", "method", "ATE.k",
                           "Psi.k","Bias.k","SD.k","COVER.k")

cset=0
  for(ell in 1:L){
  window=(ell-1)*5
  #################################################################
  ##                        Generate data                        ##
  #################################################################
  dat1 <- gen_dat_Nethery(c=cset)
  # (sample) true ATE, averaged over the 500 individual causal effects 
  ATE.true.k <- mean(dat1$tau.true)
  
  ##################################################################
  ##                         Linear model                         ##
  ##################################################################
  
  ols.mod <- lm(Y ~ x1 + x2 + trt , data=dat1)
  Psi.k <- coef(ols.mod)["trt"] 
  
  # Bias = Psi_k - ATE_true_k
  BIAS.k <- Psi.k -ATE.true.k
  # Standard deviation of Psi
  SD.k <- sqrt(vcov(ols.mod)["trt","trt"])
  # Coverage indicator
  COVER.k <- ATE.true.k >= confint(ols.mod)["trt",1] && ATE.true.k <= confint(ols.mod)["trt",2] 
  
  raw.results[window+1, ] <- c(cset, "OLS", ATE.true.k,
                               Psi.k, BIAS.k, SD.k, COVER.k)
  
  
  
  #################################################################
  ##          Inverse probability weighted linear model          ##
  #################################################################
  
  ipw.mod <- lm(Y ~ x1 + x2 + trt , data=dat1, weights = 0.5/dat1$ps.true) # use stabilized weights
  Psi.k <- coef(ipw.mod)["trt"] 
  
  # Bias = Psi_k - ATE_true_k
  BIAS.k <- Psi.k - ATE.true.k
  # Standard deviation of Psi
  SD.k <- sqrt(vcov(ipw.mod)["trt","trt"])
  # Coverage indicator
  COVER.k <- ATE.true.k >= confint(ipw.mod)["trt",1] && ATE.true.k <= confint(ipw.mod )["trt",2] 
  
  raw.results[window+2, ] <- c(cset, "IPW-OLS", ATE.true.k,
                               Psi.k, BIAS.k, SD.k, COVER.k)
  
  ##################################################################
  ##                         vanilla BART                         ##
  ##################################################################
  
  bart.mod <- bartc(data = dat1, response = Y, treatment = trt,
        confounders = x1 + x2, estimand = "ate",
        p.scoreAsCovariate = F)
  bart.res <- summary(bart.mod)$estimates
  
  Psi.k <- bart.res$estimate
  SD.k <- bart.res$sd
  COVER.k <- ATE.true.k >= bart.res[,3] && ATE.true.k <= bart.res[,4]
  
  raw.results[window+3, ] <- c(cset, "BART", ATE.true.k,
                               Psi.k, BIAS.k, SD.k, COVER.k)
  
  #################################################################
  ##                           BART-ps                           ##
  #################################################################
  bart.ps.mod <- bartc(data = dat1, response = Y, treatment = trt,
                    confounders = x1 + x2, estimand = "ate",
                    p.scoreAsCovariate = T, p.score=dat1$ps.true)
  bart.ps.res <- summary(bart.ps.mod)$estimates
  
  Psi.k <- bart.ps.res$estimate
  SD.k <- bart.ps.res$sd
  COVER.k <- ATE.true.k >= bart.ps.res[,3] && ATE.true.k <= bart.ps.res[,4]
    
  raw.results[window+4, ] <- c(cset, "BART-ps", ATE.true.k,
                               Psi.k, BIAS.k, SD.k, COVER.k)
  
  #################################################################
  ##                  Generalized Random Forest                  ##
  #################################################################
  grf.mod <- causal_forest(X=as.matrix(dat1[, c("x1", "x2")]),
                           Y=dat1$Y, W=dat1$trt, num.trees = 4000)
  grf.res <- average_treatment_effect(grf.mod,target.sample = "all", method = "AIPW")
  Psi.k <- grf.res[1]
  SD.k <- grf.res[2] # standard error...
  COVER.k <- ATE.true.k >= Psi.k-SD.k*1.96 && ATE.true.k <= Psi.k+SD.k*1.96

  raw.results[window+5, ] <- c(cset, "GRF", ATE.true.k,
                               Psi.k, BIAS.k, SD.k, COVER.k)
  }


#save(raw.results, file="raw_results_0.Rdata")


#################################################################
##                     Analysis of Results                     ##
#################################################################

#L=250

#load("raw_results_0.Rdata")
raw.results <- as.data.frame(raw.results)
raw.results[,-c(2, 7)] <- apply(raw.results[,-c(2, 7)], 2, as.numeric)
raw.results <- raw.results %>% mutate(COVER.k = if_else(COVER.k=="FALSE", 0, 1))
raw.results %>% select(-c) %>% group_by(method) %>% 
  summarise(ATE=mean(ATE.k), 
            Bias=mean(Psi.k-ATE.k), 
            PercentBias=mean((Psi.k-ATE.k)/abs(ATE.k)*100),
            SD.bar = mean(SD.k),
            SE = sqrt(sum((Psi.k-ATE)^2)/(L-1)),
            MSE = sqrt(sum((Psi.k-ATE.k)^2)/(L-1)),
            Coverage= mean(COVER.k)
            )
  
