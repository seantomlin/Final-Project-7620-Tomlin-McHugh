############################################################################
############################################################################
###                                                                      ###
###            CODE TO REPLICATE PORTIONS OF ZHU ET AL (2023)            ###
###                                                                      ###
############################################################################
############################################################################

# https://doi-org.proxy.lib.ohio-state.edu/10.1002/sim.9600

# load packages
library(dplyr)
library(ggplot2)
theme_set(theme_bw()) # set theme to BW for ggplots
library(tibble)
library(knitr)
library(gridExtra)
library(grid)
library(ggpubr)

# Some packs for the causal inference methods
library(bartCause)
library(plotBart)
library(CausalGAM)
library(grf)

##################################################################
##        Three covariates, Y1 := Linear response surface       ##
##      Three covariates, Y2 := Nonlinear response surface      ##
##################################################################
#' @description 
#' This replicates the simulation scenario from page 40 in Zhu et al (2023)
#'  @params
#' mu1 is mean of x1
#' mu2 is mean of X2
#' P is mean of X3, which is Bernoulli probability
#' N.study is number of samples 
#' @returns 
#' N=500 x 9 tibble with treatment, two Normal covariates, one binary covariate
#' observed outcome for linear response surface Y1
#' observed outcome for nonlinear response surface Y2
gen_dat_Zhu <- function(mu1, mu2, P, N.study=500, prob.Z=0.5){
  dat.Zhu <-  expand.grid(id = 1:N.study) %>% mutate(trt = rbinom(N.study, 1, prob.Z)) %>% 
    group_by(id, trt) %>% 
    mutate(x1 = if_else(trt==1, rnorm(1, mu1, 1), rnorm(1, 0, 1)), 
           x2=if_else(trt==1, rnorm(1, mu2, 1), rnorm(1, 2, 1))) %>% 
    mutate(prob.X3 =if_else(trt==1, P, 0.4)) %>% 
    mutate(x3 = rbinom(n=1, size=1, prob = prob.X3)) %>% 
    mutate(Y1 = 1 - 2*x1 + x2 - 1.2*x3 + 2*trt + rnorm(1), 
           Y2 = -3 - 2.5*x1 + 2*(x1^2)*trt + 
             exp(1.4 - x2*trt) + x2*x3 - 1.2*x3 
           - 2*x3*trt + 2*trt + rnorm(1)) %>% 
    mutate(Z=factor(trt))
  dat.Zhu
}

##########################################################################
##  DGP adapted from Nethery et al (2019), is used in Zhu et al (2023)  ##
##                Two covariates, True PS from Bayes rule               ##
##                    Potential outcomes Y(0), Y(1)                     ##
##  ATE true is really a sample ATE over the individual causal effects  ##
##########################################################################

#' @description 
#' This replicates the simulation scenario from page 41 in Zhu et al (2023)
#' and originally appears in Nethery ey al (2019)
#'  @params
#' Overlap parameter c = {0, 0.35, 0.75} 
#' @returns 
#' N=500 x 11 tibble with treatment, three covariates, true propensity score, 
#' potential outcomes, observed outcome, true individual treatment effect.
gen_dat_Nethery <- function(c, N=500, prob.Z=0.5){
  # apply this rowwise to x1, x2, prob.x1  
  true_ps <- function(x1, x2){
    A <- dnorm(x2, mean=2+c, sd=1.25+0.1*c)*dbinom(x1, 1, prob = 0.5)
    B <- dnorm(x2, mean=2+c, sd=1.25+0.1*c)*dbinom(x1, 1, prob = 0.5) + 
      dnorm(x2, mean=1, sd=1)*dbinom(x1, 1, prob = 0.4)
    A/B
  }
  # Treatment, covariates, potential outcomes
  dat <-  expand.grid(id = 1:N) %>% mutate(trt = rbinom(N, 1, prob.Z)) %>% 
    group_by(id, trt) %>% 
    mutate(prob.X1 =if_else(trt==1, 0.5, 0.4)) %>% 
    mutate(x1 = rbinom(n=1, size=1, prob = prob.X1)) %>% 
    mutate(x2 = if_else(trt==1, rnorm(1, mean=2+c, sd=1.25+0.1*c), rnorm(1, mean=1, sd=1))) %>% # change with trueps 
    mutate(ps.true = true_ps(x1, x2)) %>% 
    mutate(Y.0 = -1.5*x2, Y.1 = (-3/(1+exp(-10*(x2-1))) + 0.25*x1 - x1*x2)) %>% 
    mutate(Y = if_else(trt==1, Y.1, Y.0)) %>% 
    mutate(tau.true = Y.1- Y.0)  %>% 
    mutate(Z=factor(trt))
  dat
}



