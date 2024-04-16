#################################################################
##    Obtain ATE estimates from the linear response surface    ##
##            DGP from Zhu et al (2023), denoted Y1            ##
#################################################################

# Requires the two functions defined in 
source("data_gen_processes_Zhu2023.R")

library(bcf)
library(CausalGAM)
library(grf)


#################################################################
##             ATE for linear response surface, Y1             ##
##                            TRUE                             ##
#################################################################
# Generate data
# Perfect Overlap (not in paper)
data0 <- gen_dat_Zhu(mu1 = 0, mu2 = 2, P=0.4)
# Estimate propensity score using bart
ps.bart0 <- dbarts::bart2(formula = trt ~ x1 + x2 + x3, data = data0, n.chains = 10L,
                          verbose = FALSE)
data0 <- data0 %>% group_by(id) %>%  bind_cols(ps.estimate=fitted(ps.bart0))


##################################################################
##                         vanilla BART                         ##
##################################################################
bart.mod <- bart.ps.mod <- bartc(data = data0, response = Y1, treatment = trt,
                                 confounders = x1 + x2 + x3, estimand = "ate",
                                 p.scoreAsCovariate = F)
# Some recommended diagnostics
summary(bart.mod)
# plot_sigma(bart.mod)


#################################################################
##                           BART-ps                           ##
#################################################################
bart.ps.mod <- bartc(data = data0, response = Y1, treatment = trt,
                     confounders = x1 + x2 + x3, estimand = "ate",
                     p.scoreAsCovariate = T)
# Some recommended diagnostics
summary(bart.ps.mod)
# plot_sigma(bart.ps.mod)

##################################################################
##                    Adjusted  Linear model                    ##
##################################################################

# unweighted
fit0.unw  <- lm(Y1 ~ x1 + x2 + x3 + trt, data = data0)
summary(fit0.unw)
confint(fit0.unw)["trt",]

# Inverse probability weighted observations
fit0 <- lm(Y1 ~ x1 + x2 + x3 + trt, data = data0, weights = 1/data0$ps.estimate)
summary(fit0)
confint(fit0)["trt",]


#################################################################
##                  Generalized Random Forest                  ##
#################################################################
# https://grf-labs.github.io/grf/articles/grf_guide.html#efficiently-estimating-summaries-of-the-cates-1

cf0 <- causal_forest(X=as.matrix(data0[, c("x1", "x2", "x3")]), Y=data0$Y1, W=data0$trt)
average_treatment_effect(cf0, target.sample="treated")




##################################################################
##                    Bayesian Causal Forest                    ##
##################################################################

# Note: This may take some time to run. 
attach(data0)
colnames(data0)
bcf_fit <- bcf(y = Y1, z = trt, x_control = matrix(c(x1=x1, x2=x2, x3=x3), ncol = 3),
               x_moderate = matrix(c(x1=x1, x2=x2, x3=x3), ncol = 3),
               pihat = ps.estimate, nburn = 2000, nsim = 1000, nthin = 5)

detach(data0)

# ?summary.bcf
summary(bcf_fit)

tau_post = bcf_fit$tau
hist(tau_post)
# sd_post = bcf_fit$sigma
# tauhat = colMeans(tau_post)
# sdhat = mean(sd_post)

