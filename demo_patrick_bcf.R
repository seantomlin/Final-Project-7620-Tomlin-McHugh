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
N = 500

dat1 <- gen_dat_Nethery_witherror(c=c1)
newdat1 <- gen_dat_Nethery_witherror(c=c1)

# (sample) true ATE, averaged over the 500 individual causal effects 
mean(dat1$tau.true)
p1 <- dat1 %>% ggplot(aes(ps.true, fill=Z, group=Z)) +
  geom_histogram(alpha=0.5) + ggtitle("c=0") + xlab("True Propensity Score") + xlim(c(0,1))

# REPEAT WITH C2
dat2 <- gen_dat_Nethery_witherror(c=c2)
newdat2 <- gen_dat_Nethery_witherror(c=c2)
p2 <- dat2 %>% ggplot(aes(ps.true, fill=Z, group=Z)) + 
  geom_histogram(alpha=0.5) + ggtitle("c=0.35") + xlab("True Propensity Score") + xlim(c(0,1))

# repeat with C3
dat3 <- gen_dat_Nethery_witherror(c=c3)
newdat3 <- gen_dat_Nethery_witherror(c=c3)
p3 <- dat3 %>% ggplot(aes(ps.true, fill=Z, group=Z)) +
  geom_histogram(alpha=0.5) + ggtitle("c=0.7") + xlab("True Propensity Score") + xlim(c(0,1))

# plot the three levels of overlap
# ggarrange(p1, p2, p3, nrow=3, common.legend = TRUE, legend="bottom")
# 
# # True potential outcome models for each overlap setting c=0, 0.35, 0.7
# tau.p1 <- dat1 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") +
#   ggtitle("c=0") + ylab("Individual causal effect") + geom_rug(sides="b")
# 
# tau.p2 <- dat2 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") + 
#   ggtitle("c=0.35") + ylab("Individual causal effect") + geom_rug(sides="b")
# 
# tau.p3 <- dat3 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") + 
#   ggtitle("c=0.7") + ylab("Individual causal effect") + geom_rug(sides="b")
# 
# ggarrange(tau.p1, tau.p2, tau.p3, nrow=3, common.legend = TRUE, legend="bottom") 

tau.x2.p1 <- dat1 %>% ggplot(aes(x2, tau.true, color=Z)) + geom_point() + xlab("x2") +
  ggtitle("c=0") + ylab("Individual causal effect") + geom_rug(sides="b")

sdy <- sd(dat1[["Y"]])

burn_in <- 1200
sims <- 1000
chains <- 1
bcf_3 <- bcf(y                = dat1[["Y"]],
             z                = dat1[["trt"]],
             x_control        = as.matrix(data.frame(dat1)[,c("x1", "x2")]),
             x_moderate       = as.matrix(data.frame(dat1)[,c("x1", "x2")]),
             pihat            = dat1[["ps.true"]],
             nburn            = burn_in,
             nsim             = sims,
             n_chains         = chains,
             random_seed      = 1,
             update_interval  = 1, 
             no_output        = F,
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
plot(1:(sims), bcf_3$tau[1:sims,497], type="l")
#plot((sims+1):(2*sims), bcf_3$tau[(sims+1):(2*sims),2], type="l")

#hist(bcf_3$tau[1:sims,1])
#hist(bcf_3$yhat[1:sims,1])


true_ate <- mean(dat1$tau.true)

tau_ests <- data.frame(Mean  = colMeans(bcf_3$tau),
                       Low95 = apply(bcf_3$tau, 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(bcf_3$tau, 2, function(x) quantile(x, 0.975)))

c0plot <- ggplot(NULL, aes(x = dat1$x2, color = as.factor(dat1$trt))) +
  geom_pointrange(aes(y = tau_ests$Mean, ymin = tau_ests$Low95, ymax = tau_ests$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
c0plot

# x11_tau_ests <- data.frame(Mean  = colMeans(bcf_3$tau[,which(dat1$x1 == 1)]),
#                            Low95 = apply(bcf_3$tau[,which(dat1$x1 == 1)], 2, function(x) quantile(x, 0.025)),
#                            Up95  = apply(bcf_3$tau[,which(dat1$x1 == 1)], 2, function(x) quantile(x, 0.975)))
# 
# ggplot(NULL, aes(x = dat1[which(dat1$x1 == 1),]$x2)) +
#   geom_pointrange(aes(y = x11_tau_ests$Mean, ymin = x11_tau_ests$Low95, ymax = x11_tau_ests$Up95), color = "forestgreen", alpha = 0.5) +
#   geom_point(aes(y = dat1[which(dat1$x1 == 1),]$tau.true), se = FALSE) +
#   xlab("x2") +
#   ylab("Estimated TE") +
#   xlim(-2, 6) +
#   geom_segment(aes(x = 3, xend = 4, y = 0.2, yend = 0.2), color = "blue", alpha = 0.9) +
#   geom_text(aes(x = 4.5, y = 0.2, label = "Truth"), color = "black") +
#   geom_segment(aes(x = 3, xend = 4, y = 0.1, yend = 0.1), color = "forestgreen", alpha = 0.7) +
#   geom_text(aes(x = 5.2, y = 0.1, label = "Estimate (95% CI)"), color = "black")
# 
# z1_tau_ests <- data.frame(Mean  = colMeans(bcf_3$tau[,which(dat1$trt == 1)]),
#                           Low95 = apply(bcf_3$tau[,which(dat1$trt == 1)], 2, function(x) quantile(x, 0.025)),
#                           Up95  = apply(bcf_3$tau[,which(dat1$trt == 1)], 2, function(x) quantile(x, 0.975)))
# 
# ggplot(NULL, aes(x = dat1[which(dat1$trt == 1),]$x2)) +
#   geom_pointrange(aes(y = z1_tau_ests$Mean, ymin = z1_tau_ests$Low95, ymax = z1_tau_ests$Up95), color = "forestgreen", alpha = 0.5) +
#   geom_smooth(aes(y = dat1[which(dat1$trt == 1),]$tau.true), se = FALSE) +
#   xlab("x2") +
#   ylab("Estimated TE") +
#   xlim(-2, 6) +
#   geom_segment(aes(x = 3, xend = 4, y = 0.2, yend = 0.2), color = "blue", alpha = 0.9) +
#   geom_text(aes(x = 4.5, y = 0.2, label = "Truth"), color = "black") +
#   geom_segment(aes(x = 3, xend = 4, y = 0.1, yend = 0.1), color = "forestgreen", alpha = 0.7) +
#   geom_text(aes(x = 5.2, y = 0.1, label = "Estimate (95% CI)"), color = "black")


#################################################################
##         Generate data from Zhu et al (2023) process         ##
#################################################################
# Reproducibility
#' set.seed(2000)
#' # Perfect Overlap (not in paper, but corresponds to setting where trt, ctrl have same distributions)
#' data0 <- gen_dat_Zhu(mu1 = 0, mu2 = 2, P=0.4)
#' # some nonoverlap
#' data1 <- gen_dat_Zhu(mu1 = 1, mu2 = 2, P=0.5)
#' # substantial nonoverlap
#' data2 <- gen_dat_Zhu(mu1 = 1, mu2 = 3, P=0.6)
#' 
#' #' Plot the overlap for Vars
#' #' The estimated propensity scores are obtained using "plotBart" package. 
#' #' They use BART to estimate the propensity scores.
#' #' Sourced from line 6.
#' g0 <- plot_overlap_pScores(.data = data0, treatment = "trt", response = "Y1", 
#'                            confounders = c("x1", "x2", "x3"))  + ggtitle("Good overlap") +labs(subtitle = "mu1 = 0, mu2 = 2, P=0.4")
#' v0 <- plot_overlap_vars(.data = data0, treatment = "trt", 
#'                         confounders = c("x1", "x2", "x3"))  + ggtitle("Good overlap") +labs(subtitle = "mu1 = 0, mu2 = 2, P=0.4")
#' 
#' g1 <- plot_overlap_pScores(.data = data1, treatment = "trt", response = "Y1", 
#'                            confounders = c("x1", "x2", "x3"))  + ggtitle("Some nonoverlap") +labs(subtitle = "mu1 = 1, mu2 = 2, P=0.5")
#' v1 <-plot_overlap_vars(.data = data1, treatment = "trt",  
#'                        confounders = c("x1", "x2", "x3"))  + ggtitle("Some nonoverlap")+labs(subtitle = "mu1 = 1, mu2 = 2, P=0.5")
#' 
#' g2 <- plot_overlap_pScores(.data = data2, treatment = "trt", response = "Y1", 
#'                            confounders = c("x1", "x2", "x3")) + ggtitle("Substantial nonoverlap") +labs(subtitle = "mu1 = 1, mu2 = 3, P=0.6")
#' v2 <-plot_overlap_vars(.data = data2, treatment = "trt",  
#'                        confounders = c("x1", "x2", "x3")) + ggtitle("Substantial nonoverlap") +labs(subtitle = "mu1 = 1, mu2 = 3, P=0.6")
#' 
#' # overlap in PS, vars by treatment status, with increasing degrees of nonoverlap. 
#' grid.arrange(g0,g1,g2, nrow=3)
#' grid.arrange(v0,v1,v2, nrow=3)


sd2 <- sd(dat2[["Y"]])

burn_in <- 1000
sims <- 1000
chains <- 2
bcf_dat2 <- bcf(y                = dat2[["Y"]],
                z                = dat2[["trt"]],
                x_control        = as.matrix(data.frame(dat2)[,c("x1", "x2")]),
                x_moderate       = as.matrix(data.frame(dat2)[,c("x1", "x2")]),
                pihat            = dat2[["ps.true"]],
                nburn            = burn_in,
                nsim             = sims,
                n_chains         = chains,
                random_seed      = 1,
                update_interval  = 1, 
                no_output        = F,
                ntree_control = 200,
                ntree_moderate = 100,
                sigq = 0.75,
                use_tauscale = F,
                use_muscale = F,
                sd_control = 4*sd2,
                sd_moderate = 4*sd2,
                power_moderate = 1,
                base_moderate = 0.95,
                power_control = 1,
                base_control = 0.95,
                include_pi = "both")

plot(1:(sims), bcf_dat2$tau[1:sims,2], type="l")
plot(1:(sims), bcf_dat2$tau[1:sims,453], type="l")
#plot((sims+1):(2*sims), bcf_3$tau[(sims+1):(2*sims),2], type="l")

#hist(bcf_3$tau[1:sims,1])
#hist(bcf_3$yhat[1:sims,1])



tau_ests2 <- data.frame(Mean  = colMeans(bcf_dat2$tau),
                        Low95 = apply(bcf_dat2$tau, 2, function(x) quantile(x, 0.025)),
                        Up95  = apply(bcf_dat2$tau, 2, function(x) quantile(x, 0.975)))

c035plot <- ggplot(NULL, aes(x = dat2$x2, color = as.factor(dat2$trt))) +
  geom_pointrange(aes(y = tau_ests2$Mean, ymin = tau_ests2$Low95, ymax = tau_ests2$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat2$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0.35)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
c035plot




sd3 <- sd(dat3[["Y"]])

burn_in <- 1000
sims <- 1000
chains <- 1
bcf_dat3 <- bcf(y                = dat3[["Y"]],
                z                = dat3[["trt"]],
                x_control        = as.matrix(data.frame(dat3)[,c("x1", "x2")]),
                x_moderate       = as.matrix(data.frame(dat3)[,c("x1", "x2")]),
                pihat            = dat3[["ps.true"]],
                nburn            = burn_in,
                nsim             = sims,
                n_chains         = chains,
                random_seed      = 1,
                update_interval  = 1, 
                no_output        = F,
                ntree_control = 200,
                ntree_moderate = 100,
                sigq = 0.75,
                use_tauscale = F,
                use_muscale = F,
                sd_control = 4*sd2,
                sd_moderate = 4*sd2,
                power_moderate = 1,
                base_moderate = 0.95,
                power_control = 1,
                base_control = 0.95,
                include_pi = "both")

plot(1:(sims), bcf_dat3$tau[1:sims,2], type="l")
plot(1:(sims), bcf_dat3$tau[1:sims,380], type="l")
#plot((sims+1):(2*sims), bcf_3$tau[(sims+1):(2*sims),380], type="l")

#hist(bcf_3$tau[1:sims,1])
#hist(bcf_3$yhat[1:sims,1])



tau_ests3 <- data.frame(Mean  = colMeans(bcf_dat3$tau),
                        Low95 = apply(bcf_dat3$tau, 2, function(x) quantile(x, 0.025)),
                        Up95  = apply(bcf_dat3$tau, 2, function(x) quantile(x, 0.975)))

c07plot <- ggplot(NULL, aes(x = dat3$x2, color = as.factor(dat3$trt))) +
  geom_pointrange(aes(y = tau_ests3$Mean, ymin = tau_ests3$Low95, ymax = tau_ests3$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat3$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0.7)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
c07plot

library(patchwork)
patch <- c0plot / c035plot / c07plot
print(patch)



pred_c0 = predict(object=bcf_3,
                   x_predict_control=as.matrix(data.frame(newdat1)[,c("x1", "x2")]),
                   x_predict_moderate=as.matrix(data.frame(newdat1)[,c("x1", "x2")]),
                   pi_pred=newdat1[["ps.true"]],
                   z_pred=newdat1[["trt"]],
                   n_cores = 1,
                   save_tree_directory = '.')


tau_ests_pred_c0 <- data.frame(Mean  = colMeans(pred_c0$tau),
                       Low95 = apply(pred_c0$tau, 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(pred_c0$tau, 2, function(x) quantile(x, 0.975)))

c0_pred_plot <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_pointrange(aes(y = tau_ests_pred_c0$Mean, ymin = tau_ests_pred_c0$Low95, ymax = tau_ests_pred_c0$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
c0_pred_plot

c0_patch <- c0plot / c0_pred_plot
c0_patch

c0_mse_train <- sum((colMeans(bcf_3$tau) - dat1[["tau.true"]])^2)/N
c0_mse_test <- sum((colMeans(pred_c0$tau) - newdat1[["tau.true"]])^2)/N



pred_c035 = predict(object=bcf_dat2,
                  x_predict_control=as.matrix(data.frame(newdat1)[,c("x1", "x2")]),
                  x_predict_moderate=as.matrix(data.frame(newdat1)[,c("x1", "x2")]),
                  pi_pred=newdat1[["ps.true"]],
                  z_pred=newdat1[["trt"]],
                  n_cores = 1,
                  save_tree_directory = '.')


tau_ests_pred_c035 <- data.frame(Mean  = colMeans(pred_c035$tau),
                               Low95 = apply(pred_c035$tau, 2, function(x) quantile(x, 0.025)),
                               Up95  = apply(pred_c035$tau, 2, function(x) quantile(x, 0.975)))

c035_pred_plot <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_pointrange(aes(y = tau_ests_pred_c035$Mean, ymin = tau_ests_pred_c035$Low95, ymax = tau_ests_pred_c035$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0.35)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
c035_pred_plot


c035_mse_train <- sum((colMeans(bcf_dat2$tau) - dat2[["tau.true"]])^2)/N
c035_mse_test <- sum((colMeans(pred_c035$tau) - newdat1[["tau.true"]])^2)/N

pred_c07 = predict(object=bcf_dat3,
                    x_predict_control=as.matrix(data.frame(newdat1)[,c("x1", "x2")]),
                    x_predict_moderate=as.matrix(data.frame(newdat1)[,c("x1", "x2")]),
                    pi_pred=newdat1[["ps.true"]],
                    z_pred=newdat1[["trt"]],
                    n_cores = 1,
                    save_tree_directory = '.')


tau_ests_pred_c07 <- data.frame(Mean  = colMeans(pred_c07$tau),
                                 Low95 = apply(pred_c07$tau, 2, function(x) quantile(x, 0.025)),
                                 Up95  = apply(pred_c07$tau, 2, function(x) quantile(x, 0.975)))

c07_pred_plot <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_point(aes(y = tau_ests_pred_c07$Mean), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0.7)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
c07_pred_plot


c07_mse_train <- sum((colMeans(bcf_dat3$tau) - dat3[["tau.true"]])^2)/N
c07_mse_test <- sum((colMeans(pred_c07$tau) - newdat1[["tau.true"]])^2)/N

mse <- data.frame(c = c(0, 0.35, 0.7), train = c(c0_mse_train, c035_mse_train, c07_mse_train), test =  c(c0_mse_test, c035_mse_test, c07_mse_test))






lm_c0 <- lm(Y ~ x1*x2*trt , data=dat1)
dat1z1 <- mutate(dat1, trt=c(1))
dat1z0 <- mutate(dat1, trt=c(0))
lm_pred_tau_train <- predict(lm_c0, newdata = dat1z1) - predict(lm_c0, newdata = dat1z0)
lm_mse_train_c0 <- mean((dat1$tau.true - lm_pred_tau_train)^2)

c07_pred_plot_lm <- ggplot(NULL, aes(x = dat1$x2, color = as.factor(dat1$trt))) +
  geom_point(aes(y = lm_pred_tau_train), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0)", subtitle="Black: Linear Estimate; Red/Blue: True Treatment Effect", col="z")
c07_pred_plot_lm

newdat1z1 <- mutate(newdat1, trt=c(1))
newdat1z0 <- mutate(newdat1, trt=c(0))
lm_pred_tau_test <- predict(lm_c0, newdata = newdat1z1) - predict(lm_c0, newdata = newdat1z0)
lm_mse_test_c0 <- mean((newdat1$tau.true - lm_pred_tau_test)^2)

c07_pred_plot_lm_test <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_point(aes(y = lm_pred_tau_test), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0)", subtitle="Black: Linear Estimate; Red/Blue: True Treatment Effect", col="z")
c07_pred_plot_lm_test


lm_c035 <- lm(Y ~ (x1+x2+trt)^3 , data=dat2)
dat2z1 <- mutate(dat2, trt=c(1))
dat2z0 <- mutate(dat2, trt=c(0))
lm_pred_tau_train <- predict(lm_c035, newdata = dat2z1) - predict(lm_c035, newdata = dat2z0)
lm_mse_train_c035 <- mean((dat2$tau.true - lm_pred_tau_train)^2)

c035_pred_plot_lm <- ggplot(NULL, aes(x = dat2$x2, color = as.factor(dat2$trt))) +
  geom_point(aes(y = lm_pred_tau_train), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat2$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0.35)", subtitle="Black: Linear Estimate; Red/Blue: True Treatment Effect", col="z")
c035_pred_plot_lm

lm_pred_tau_test <- predict(lm_c035, newdata = newdat1z1) - predict(lm_c035, newdata = newdat1z0)
lm_mse_test_c035 <- mean((newdat1$tau.true - lm_pred_tau_test)^2)

c035_pred_plot_lm_test <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_point(aes(y = lm_pred_tau_test), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0.35)", subtitle="Black: Linear Estimate; Red/Blue: True Treatment Effect", col="z")
c035_pred_plot_lm_test

lm_c07 <- lm(Y ~ (x1+x2+trt)^3 , data=dat3)
dat3z1 <- mutate(dat3, trt=c(1))
dat3z0 <- mutate(dat3, trt=c(0))
lm_pred_tau_train <- predict(lm_c07, newdata = dat3z1) - predict(lm_c07, newdata = dat3z0)
lm_mse_train_c07 <- mean((dat3$tau.true - lm_pred_tau_train)^2)

c07_pred_plot_lm <- ggplot(NULL, aes(x = dat3$x2, color = as.factor(dat3$trt))) +
  geom_point(aes(y = lm_pred_tau_train), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat3$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0.7)", subtitle="Black: Linear Estimate; Red/Blue: True Treatment Effect", col="z")
c07_pred_plot_lm

lm_pred_tau_test <- predict(lm_c07, newdata = newdat1z1) - predict(lm_c07, newdata = newdat1z0)
lm_mse_test_c07 <- mean((newdat1$tau.true - lm_pred_tau_test)^2)

c07_pred_plot_lm_test <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_point(aes(y = lm_pred_tau_test), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0.7)", subtitle="Black: Linear Estimate; Red/Blue: True Treatment Effect", col="z")
c07_pred_plot_lm_test


d <- dat3
dz1 <- mutate(d, trt=c(1))
dz0 <- mutate(d, trt=c(0))
bart.mod <- bartc(data = d, response = Y, treatment = trt,
                  confounders = x1 + x2, estimand = "ate",
                  p.scoreAsCovariate = F, keepTrees=T)
bart_pred_tau_train <- colMeans(predict(bart.mod, newdata = dz1)) - colMeans(predict(bart.mod, newdata = dz0))
bart_mse_train <- mean((d$tau.true - bart_pred_tau_train)^2)
pred_plot_bart_train <- ggplot(NULL, aes(x = d$x2, color = as.factor(d$trt))) +
  geom_point(aes(y = bart_pred_tau_train), color = "black", alpha = 0.3) +
  geom_point(aes(y = d$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Train Set (c = 0.7)", subtitle="Black: BART Estimate; Red/Blue: True Treatment Effect", col="z")
pred_plot_bart_train

bart_pred_tau_test <- colMeans(predict(bart.mod, newdata = newdat1z1)) - colMeans(predict(bart.mod, newdata = newdat1z0))
bart_mse_test <- mean((newdat1$tau.true - bart_pred_tau_test)^2)

pred_plot_bart_test <- ggplot(NULL, aes(x = newdat1$x2, color = as.factor(newdat1$trt))) +
  geom_point(aes(y = bart_pred_tau_test), color = "black", alpha = 0.3) +
  geom_point(aes(y = newdat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation, Test Set (c = 0.7)", subtitle="Black: BART Estimate; Red/Blue: True Treatment Effect", col="z")
pred_plot_bart_test





burn_in <- 1000
sims <- 2000
chains <- 1
bcf_good <- bcf(y                = dat1[["Y"]],
                z                = dat1[["trt"]],
                x_control        = as.matrix(data.frame(dat1)[,c("x1", "x2")]),
                x_moderate       = as.matrix(data.frame(dat1)[,c("x1", "x2")]),
                pihat            = dat1[["ps.true"]],
                nburn            = burn_in,
                nsim             = sims,
                n_chains         = chains,
                random_seed      = 1,
                update_interval  = 1, 
                no_output        = TRUE,
                ntree_control = 200,
                ntree_moderate = 150,
                sigq = 0.75,
                use_tauscale = F,
                use_muscale = F,
                sd_control = 3*sdy,
                sd_moderate = 3*sdy,
                power_moderate = 1,
                base_moderate = 0.95,
                power_control = 1,
                base_control = 0.95,
                include_pi = "control")

tau_ests_good <- data.frame(Mean  = colMeans(bcf_good$tau),
                       Low95 = apply(bcf_good$tau, 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(bcf_good$tau, 2, function(x) quantile(x, 0.975)))

ggplot(NULL, aes(x = dat1$x2, color = as.factor(dat1$trt))) +
  geom_pointrange(aes(y = tau_ests_good$Mean, ymin = tau_ests_good$Low95, ymax = tau_ests_good$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation (c = 0)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")





tau_ests <- data.frame(Mean  = colMeans(bcf_2$tau),
                            Low95 = apply(bcf_2$tau, 2, function(x) quantile(x, 0.025)),
                            Up95  = apply(bcf_2$tau, 2, function(x) quantile(x, 0.975)))

ggplot(NULL, aes(x = dat1$x2, color = as.factor(dat1$trt))) +
  geom_pointrange(aes(y = tau_ests$Mean, ymin = tau_ests$Low95, ymax = tau_ests$Up95), color = "black", alpha = 0.3) +
  geom_point(aes(y = dat1$tau.true)) +
  xlab("x2") +
  ylab("Treatment Effect") +
  xlim(-2, 6) +
  labs(title = "Treatment Effect Estimation (c = 0)", subtitle="Black: BCF Estimate/95% CI", col="z")
