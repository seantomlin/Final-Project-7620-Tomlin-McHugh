# bcf ATE estimation

set.seed(20240419) 
source("data_gen_processes_Zhu2023.R")


# Overlap settings
c1 <- 0 # good overlap
c2 <- 0.35 # some nonoverlap
c3 <- 0.7 # Substantial nonoverlap
overlaps <- c(c1, c2, c3)
N = 500
n_reps <- 90

results <- matrix(NA, nrow = n_reps*3, ncol=9)


start_time <- Sys.time()

for (i in 1:n_reps){
  
  dat1 <- gen_dat_Nethery(c=c1)
  dat2 <- gen_dat_Nethery(c=c2)
  dat3 <- gen_dat_Nethery(c=c3)
  
  
  for (j in 1:3){
    
    if (j == 1){
      d <- dat1
    } else if (j == 2){
      d <- dat2
    } else if (j == 3){
      d <- dat3
    }
    
    sdy <- sd(d[["Y"]])
    
    burn_in <- 1000
    sims <- 1000
    chains <- 1
    bcf_model <- bcf(y                = d[["Y"]],
                 z                = d[["trt"]],
                 x_control        = as.matrix(data.frame(d)[,c("x1", "x2")]),
                 x_moderate       = as.matrix(data.frame(d)[,c("x1", "x2")]),
                 pihat            = d[["ps.true"]],
                 nburn            = burn_in,
                 nsim             = sims,
                 n_chains         = chains,
                 random_seed      = 1,
                 update_interval  = 1, 
                 no_output        = T,
                 ntree_control = 200,
                 ntree_moderate = 100,
                 sigq = 0.75,
                 use_tauscale = F,
                 use_muscale = F,
                 sd_control = 4*sdy,
                 sd_moderate = 4*sdy,
                 power_moderate = 1,
                 base_moderate = 0.95,
                 power_control = 1,
                 base_control = 0.95,
                 include_pi = "both")
    
    # tau_ests <- data.frame(Mean  = colMeans(bcf_model$tau),
    #                        Low95 = apply(bcf_model$tau, 2, function(x) quantile(x, 0.025)),
    #                        Up95  = apply(bcf_model$tau, 2, function(x) quantile(x, 0.975)))
    # 
    # c0plot <- ggplot(NULL, aes(x = d$x2, color = as.factor(d$trt))) +
    #   geom_pointrange(aes(y = tau_ests$Mean, ymin = tau_ests$Low95, ymax = tau_ests$Up95), color = "black", alpha = 0.3) +
    #   geom_point(aes(y = d$tau.true)) +
    #   xlab("x2") +
    #   ylab("Treatment Effect") +
    #   xlim(-2, 6) +
    #   labs(title = "Treatment Effect Estimation, Train Set (c = 0)", subtitle="Black: BCF Estimate/95% CI; Red/Blue: True Treatment Effect", col="z")
    # print(c0plot)
    
    ate_posterior <- rowMeans(bcf_model$tau)
    
    est_ate <- mean(ate_posterior)
    true_ate <- mean(d$tau.true)
    ci_95 <- quantile(ate_posterior, c(.025, .975))
    ci_width <- ci_95[2] - ci_95[1]
    covered <- ifelse(true_ate < ci_95[2] & true_ate > ci_95[1], 1, 0)
    tau_hat_se <- sd(ate_posterior)
    
    
    new_results <- c(i, overlaps[j], true_ate, est_ate, tau_hat_se, ci_95[1], ci_95[2], ci_width, covered)
    results[(3*(i-1)+j),] <- new_results
    
    cat( "\n \n \n \n Completed rep", i,"-", j, " \n \n \n \n")
    
  }
  
  
}

colnames(results) <- c("rep", "c", "true ATE", "estimated ATE", "se(tau-hat)", "95%CI_lower", "95%CI_upper", "CI_width", "covered")


end_time <- Sys.time()
total_time <- end_time - start_time
time_per_rep <- total_time / n_reps
cat("Total Time", total_time)
cat("Time per rep", time_per_rep)



write.csv(results, paste("C:/Users/pmchu/OneDrive - The Ohio State University/7620/project/bcf_results2.csv", sep=""),  row.names = F)



