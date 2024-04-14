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

sdy <- sd(dat1[["Y"]])

burn_in <- 500
sims <- 1000
chains <- 1
bcf_3 <- bcf(y                = dat1[["Y"]],
               z                = dat1[["trt"]],
               x_control        = as.matrix(data.frame(dat1)[,c("x1", "x2")]),
               x_moderate       = as.matrix(data.frame(dat1)[,c("x2")]),
               pihat            = dat1[["ps.true"]],
               nburn            = burn_in,
               nsim             = sims,
               n_chains         = chains,
               random_seed      = 1,
               update_interval  = 1, 
               no_output        = TRUE,
               ntree_control = 200,
               ntree_moderate = 100,
               sigq = 0.75,
               use_tauscale = F,
               use_muscale = F,
               sd_control = 3*sdy,
               sd_moderate = 3*sdy,
               power_moderate = 1,
               base_moderate = 0.95,
               power_control = 1,
               base_control = 0.95,
               include_pi = "both")

plot(1:(sims), bcf_3$tau[1:sims,2], type="l")
plot((sims+1):(2*sims), bcf_3$tau[(sims+1):(2*sims),2], type="l")

hist(bcf_3$tau[1:sims,1])
hist(bcf_3$yhat[1:sims,1])


true_ate <- mean(dat1$tau.true)

tau_ests <- data.frame(Mean  = colMeans(bcf_3$tau),
                       Low95 = apply(bcf_3$tau, 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(bcf_3$tau, 2, function(x) quantile(x, 0.975)))

ggplot(NULL, aes(x = dat1$x2)) +
  geom_pointrange(aes(y = tau_ests$Mean, ymin = tau_ests$Low95, ymax = tau_ests$Up95), color = "forestgreen", alpha = 0.5) +
  geom_smooth(aes(y = dat1$tau.true), se = FALSE) +
  xlab("x2") +
  ylab("Estimated TE") +
  xlim(-2, 6) +
  geom_segment(aes(x = 3, xend = 4, y = 0.2, yend = 0.2), color = "blue", alpha = 0.9) +
  geom_text(aes(x = 4.5, y = 0.2, label = "Truth"), color = "black") +
  geom_segment(aes(x = 3, xend = 4, y = 0.1, yend = 0.1), color = "forestgreen", alpha = 0.7) +
  geom_text(aes(x = 5.2, y = 0.1, label = "Estimate (95% CI)"), color = "black")


x11_tau_ests <- data.frame(Mean  = colMeans(bcf_3$tau[,which(dat1$x1 == 1)]),
                       Low95 = apply(bcf_3$tau[,which(dat1$x1 == 1)], 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(bcf_3$tau[,which(dat1$x1 == 1)], 2, function(x) quantile(x, 0.975)))

ggplot(NULL, aes(x = dat1[which(dat1$x1 == 1),]$x2)) +
  geom_pointrange(aes(y = x11_tau_ests$Mean, ymin = x11_tau_ests$Low95, ymax = x11_tau_ests$Up95), color = "forestgreen", alpha = 0.5) +
  geom_smooth(aes(y = dat1[which(dat1$x1 == 1),]$tau.true), se = FALSE) +
  xlab("x2") +
  ylab("Estimated TE") +
  xlim(-2, 6) +
  geom_segment(aes(x = 3, xend = 4, y = 0.2, yend = 0.2), color = "blue", alpha = 0.9) +
  geom_text(aes(x = 4.5, y = 0.2, label = "Truth"), color = "black") +
  geom_segment(aes(x = 3, xend = 4, y = 0.1, yend = 0.1), color = "forestgreen", alpha = 0.7) +
  geom_text(aes(x = 5.2, y = 0.1, label = "Estimate (95% CI)"), color = "black")

z1_tau_ests <- data.frame(Mean  = colMeans(bcf_3$tau[,which(dat1$trt == 1)]),
                           Low95 = apply(bcf_3$tau[,which(dat1$trt == 1)], 2, function(x) quantile(x, 0.025)),
                           Up95  = apply(bcf_3$tau[,which(dat1$trt == 1)], 2, function(x) quantile(x, 0.975)))

ggplot(NULL, aes(x = dat1[which(dat1$trt == 1),]$x2)) +
  geom_pointrange(aes(y = z1_tau_ests$Mean, ymin = z1_tau_ests$Low95, ymax = z1_tau_ests$Up95), color = "forestgreen", alpha = 0.5) +
  geom_smooth(aes(y = dat1[which(dat1$trt == 1),]$tau.true), se = FALSE) +
  xlab("x2") +
  ylab("Estimated TE") +
  xlim(-2, 6) +
  geom_segment(aes(x = 3, xend = 4, y = 0.2, yend = 0.2), color = "blue", alpha = 0.9) +
  geom_text(aes(x = 4.5, y = 0.2, label = "Truth"), color = "black") +
  geom_segment(aes(x = 3, xend = 4, y = 0.1, yend = 0.1), color = "forestgreen", alpha = 0.7) +
  geom_text(aes(x = 5.2, y = 0.1, label = "Estimate (95% CI)"), color = "black")


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





