#################################################################
##            Obtain ATE estimates from Nethery DGP            ##
#################################################################

# Requires the two functions defined in 
source("data_gen_processes_Zhu2023.R")

library(bcf)
# library(CausalGAM)
# library(grf)


dat1 <- gen_dat_Nethery(c=0)
# (sample) true ATE, averaged over the 500 individual causal effects 
mean(dat1$tau.true)
p1 <- dat1 %>% ggplot(aes(ps.true, fill=Z, group=Z)) +
  geom_histogram(alpha=0.5) + ggtitle("c=0") + xlab("True Propensity Score") + xlim(c(0,1))

tau.p1 <- dat1 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") +
  ggtitle("c=0") + ylab("Individual causal effect") + geom_rug(sides="b")

grid.arrange(p1, tau.p1)


##################################################################
##                    Bayesian Causal Forest                    ##
##################################################################

attach(dat1)
colnames(dat1)
start_time <- Sys.time()
bcf_fit <- bcf(y = Y, z = trt, x_control = matrix(c(x1=x1, x2=x2), ncol = 2),
               x_moderate = matrix(c(x1=x1, x2=x2), ncol = 2),
               pihat = ps.true, nburn = 10000, nsim = 5000, nthin = 5)
end_time <- Sys.time()
end_time - start_time
detach(dat1)

# ?summary.bcf
summary(bcf_fit)

tau_post = bcf_fit$tau
hist(tau_post)

sd_post = apply(bcf_fit$sigma, 1, mean)

# posterior distribution of y for each individual (i in columns)
str(bcf_fit$yhat)



# Make predicted values plot by propensity score. 
#load("C:/Users/User/OneDrive - The Ohio State University/Spring 2024/STAT 7620/Final Project/Final-Project-7620-Tomlin-Mchugh/bcf_nethery_fit_onehour.RData")

tauhat.i <-apply(bcf_fit$tau, 2,mean) 
tauhat.lower <-apply(bcf_fit$tau, 2,quantile, probs=c(0.025))  # lower quantile
tauhat.upper <-apply(bcf_fit$tau, 2,quantile, probs=c(0.975)) # upper quantile


tau.true.i <- dat1$tau.true 
data.frame(tauhat.i, tau.true.i) %>% ggplot(aes(tau.true.i, tauhat.i)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1)

dat1 <- dat1 %>% group_by(id) %>% bind_cols(tau.hat=tauhat.i, tauhat.lower=tauhat.lower, tauhat.upper=tauhat.upper)

TE.plot1 <- dat1 %>%  ggplot(aes(ps.true, tau.hat)) + geom_point(col="red") +
  geom_point(aes(y=tau.true), color="black") + 
  geom_errorbar(aes(ymin=tauhat.lower, ymax=tauhat.upper), color="red") +
  geom_rug(aes(ps.true, color=Z), sides="b", linewidth=2, alpha=0.5) + 
  ylab("Treatment effect") + xlab("True propensity score") + 
  ggtitle("True and estimated individual treatment effects by True PS")


TE.plot2 <- dat1 %>%  ggplot(aes(x2, tau.hat)) + geom_point(col="red") + 
  geom_point(aes(y=tau.true), color="black") + 
  geom_errorbar(aes(ymin=tauhat.lower, ymax=tauhat.upper), color="red") + 
  geom_rug(aes(x2, color=Z), sides="b", linewidth=2, alpha=0.5) + 
  ylab("Treatment effect") + xlab("X2") + 
  ggtitle("True and estimated individual treatment effects by X2")


grid.arrange(TE.plot1, TE.plot2)


#dat1 %>% ggplot(aes(x2, fill=Z)) + geom_histogram(alpha=3)




#################################################################
##                       Good overlap ###                      ##
#################################################################
set.seed(999)
dat1 <- gen_dat_Nethery(c=0)
# (sample) true ATE, averaged over the 500 individual causal effects 
mean(dat1$tau.true)
p1 <- dat1 %>% ggplot(aes(ps.true, fill=Z, group=Z)) +
  geom_histogram(alpha=0.5) + ggtitle("c=0") + xlab("True Propensity Score") + xlim(c(0,1))

tau.p1 <- dat1 %>% ggplot(aes(ps.true, tau.true, color=Z)) + geom_point() + xlab("True propensity score") +
  ggtitle("c=0") + ylab("Individual causal effect") + geom_rug(sides="b")

grid.arrange(p1, tau.p1)

attach(dat1)
colnames(dat1)
start_time <- Sys.time()
bcf_fit <- bcf(y = Y, z = trt, x_control = matrix(c(x1=x1, x2=x2), ncol = 2),
               x_moderate = matrix(c(x1=x1, x2=x2), ncol = 2),
               pihat = ps.true, nburn = 10000, nsim = 5000, nthin = 5)
end_time <- Sys.time()
end_time - start_time
detach(dat1)


summary(bcf_fit)



tauhat.i <-apply(bcf_fit$tau, 2,mean) 
tauhat.lower <-apply(bcf_fit$tau, 2,quantile, probs=c(0.025))  # lower quantile
tauhat.upper <-apply(bcf_fit$tau, 2,quantile, probs=c(0.975)) # upper quantile
muhat.i <- apply(bcf_fit$mu, 2,mean) # prognostic function mean


tau.true.i <- dat1$tau.true 
data.frame(tauhat.i, tau.true.i) %>% ggplot(aes(tau.true.i, tauhat.i)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1)

dat1 <- dat1 %>% group_by(id) %>% bind_cols(tau.hat=tauhat.i, tauhat.lower=tauhat.lower, tauhat.upper=tauhat.upper)
dat1 <- dat1 %>% bind_cols(muhat=muhat.i)


TE.plot1 <- dat1 %>%  ggplot(aes(ps.true, tau.hat)) + geom_point(col="red") +
  geom_point(aes(y=tau.true), color="black") + 
  geom_errorbar(aes(ymin=tauhat.lower, ymax=tauhat.upper), color="red") +
  geom_rug(aes(ps.true, color=Z), sides="b", linewidth=2, alpha=0.5) + 
  ylab("Treatment effect") + xlab("True propensity score") + 
  ggtitle("True and estimated individual treatment effects by True PS")


TE.plot2 <- dat1 %>%  ggplot(aes(x2, tau.hat)) + geom_point(col="red") + 
  geom_point(aes(y=tau.true), color="black") + 
  geom_errorbar(aes(ymin=tauhat.lower, ymax=tauhat.upper), color="red") + 
  geom_rug(aes(x2, color=Z), sides="b", linewidth=2, alpha=0.5) + 
  ylab("Treatment effect") + xlab("X2") + 
  ggtitle("True and estimated individual treatment effects by X2")


grid.arrange(TE.plot1, TE.plot2)
#save.image("C:/Users/User/OneDrive - The Ohio State University/Spring 2024/STAT 7620/Final Project/Final-Project-7620-Tomlin-Mchugh/bcf_fit_nethery_c0_onehour.RData")


plot(1:5000, bcf_fit$sigma[1:5000], type="l")


prog.plot <- dat1 %>%  ggplot(aes(ps.true, muhat)) + geom_point(col="red") + 
  ylab("Estimated prognostic function") + xlab("True PS") + 
  ggtitle("Estimated prognostic function by True PS", subtitle = "Evidence of targeted selection") 
prog.plot




