##################################################################
##             Demo of each data generating process             ##
##################################################################

# Requires the two functions defined in 
# source("data_gen_processes_Zhu2023.R")


#################################################################
##       Generate data from Nethery et al (2019) process       ##
#################################################################
# Reproducibility
set.seed(1000) 

# Overlap settings
c1 <- 0 # good overlap
c2 <- 0.35 # some nonoverlap
c3 <- 0.7 # Substantial nonoverlap

dat1 <- gen_dat_Nethery(c=c1)
# (sample) true ATE, averaged over the 500 individual causal effects 
mean(dat1$tau.true)
p1 <- dat1 %>% ggplot(aes(ps.true, fill=Z, group=Z)) +
  geom_histogram(alpha=0.5) + ggtitle("c=0") + xlab("True Propensity Score") + xlim(c(0,1))

# REPEAT WITH C2
dat2 <- gen_dat_Nethery(c=c2)
p2 <- dat2 %>% ggplot(aes(ps.true, fill=Z, group=Z)) + 
  geom_histogram(alpha=0.5) + ggtitle("c=0.35") + xlab("True Propensity Score") + xlim(c(0,1))

# repeat with C3
dat3 <- gen_dat_Nethery(c=c3)
p3 <- dat3 %>% ggplot(aes(ps.true, fill=Z, group=Z)) +
  geom_histogram(alpha=0.5) + ggtitle("c=0.7") + xlab("True Propensity Score") + xlim(c(0,1))

# plot the three levels of overlap
ggarrange(p1, p2, p3, nrow=3, common.legend = TRUE, legend="bottom")

# True potential outcome models for each overlap setting c=0, 0.35, 0.7
tau.p1 <- dat1 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") +
  ggtitle("c=0") + ylab("Individual causal effect") + geom_rug(sides="b")

tau.p2 <- dat2 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") + 
  ggtitle("c=0.35") + ylab("Individual causal effect") + geom_rug(sides="b")

tau.p3 <- dat3 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") + 
  ggtitle("c=0.7") + ylab("Individual causal effect") + geom_rug(sides="b")

ggarrange(tau.p1, tau.p2, tau.p3, nrow=3, common.legend = TRUE, legend="bottom") 



#################################################################
##         Generate data from Zhu et al (2023) process         ##
#################################################################
# Reproducibility
set.seed(2000)
# Perfect Overlap (not in paper, but corresponds to setting where trt, ctrl have same distributions)
data0 <- gen_dat_Zhu(mu1 = 0, mu2 = 2, P=0.4)
# some nonoverlap
data1 <- gen_dat_Zhu(mu1 = 1, mu2 = 2, P=0.5)
# substantial nonoverlap
data2 <- gen_dat_Zhu(mu1 = 1, mu2 = 3, P=0.6)

#' Plot the overlap for Vars
#' The estimated propensity scores are obtained using "plotBart" package. 
#' They use BART to estimate the propensity scores.
#' Sourced from line 6.
g0 <- plot_overlap_pScores(.data = data0, treatment = "trt", response = "Y1", 
                           confounders = c("x1", "x2", "x3"))  + ggtitle("Good overlap") +labs(subtitle = "mu1 = 0, mu2 = 2, P=0.4")
v0 <- plot_overlap_vars(.data = data0, treatment = "trt", 
                        confounders = c("x1", "x2", "x3"))  + ggtitle("Good overlap") +labs(subtitle = "mu1 = 0, mu2 = 2, P=0.4")

g1 <- plot_overlap_pScores(.data = data1, treatment = "trt", response = "Y1", 
                           confounders = c("x1", "x2", "x3"))  + ggtitle("Some nonoverlap") +labs(subtitle = "mu1 = 1, mu2 = 2, P=0.5")
v1 <-plot_overlap_vars(.data = data1, treatment = "trt",  
                       confounders = c("x1", "x2", "x3"))  + ggtitle("Some nonoverlap")+labs(subtitle = "mu1 = 1, mu2 = 2, P=0.5")

g2 <- plot_overlap_pScores(.data = data2, treatment = "trt", response = "Y1", 
                           confounders = c("x1", "x2", "x3")) + ggtitle("Substantial nonoverlap") +labs(subtitle = "mu1 = 1, mu2 = 3, P=0.6")
v2 <-plot_overlap_vars(.data = data2, treatment = "trt",  
                       confounders = c("x1", "x2", "x3")) + ggtitle("Substantial nonoverlap") +labs(subtitle = "mu1 = 1, mu2 = 3, P=0.6")

# overlap in PS, vars by treatment status, with increasing degrees of nonoverlap. 
grid.arrange(g0,g1,g2, nrow=3)
grid.arrange(v0,v1,v2, nrow=3)





