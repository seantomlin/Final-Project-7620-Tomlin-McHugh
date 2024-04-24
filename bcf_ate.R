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

results <- data.frame(results)
results$bias <- results$estimated.ATE - results$true.ATE
results$pct.bias <- (results$bias / abs(results$true.ATE))*100
results$sqderr <- (results$estimated.ATE - results$true.ATE)^2

ate_c0 <- mean(results[which(results[,"c"] == 0), "estimated.ATE"])
ate_c035 <- mean(results[which(results[,"c"] == 0.35), "estimated.ATE"])
ate_c07 <- mean(results[which(results[,"c"] == 0.7), "estimated.ATE"])

# Coverage
mean(results[which(results[,"c"] == 0), "covered"])
mean(results[which(results[,"c"] == 0.35), "covered"])
mean(results[which(results[,"c"] == 0.7), "covered"])

# ATE
mean(results[which(results[,"c"] == 0), "true.ATE"])
mean(results[which(results[,"c"] == 0.35), "true.ATE"])
mean(results[which(results[,"c"] == 0.7), "true.ATE"])

# Bias
mean(results[which(results[,"c"] == 0), "bias"])
mean(results[which(results[,"c"] == 0.35), "bias"])
mean(results[which(results[,"c"] == 0.7), "bias"])

# Percent bias
# mean(results[which(results[,"c"] == 0), "bias"]) / mean(results[which(results[,"c"] == 0), "true.ATE"])
# mean(results[which(results[,"c"] == 0.35), "bias"]) / mean(results[which(results[,"c"] == 0.35), "true.ATE"])
# mean(results[which(results[,"c"] == 0.7), "bias"]) / mean(results[which(results[,"c"] == 0.7), "true.ATE"])
mean(results[which(results[,"c"] == 0), "pct.bias"])
mean(results[which(results[,"c"] == 0.35), "pct.bias"])
mean(results[which(results[,"c"] == 0.7), "pct.bias"])

# SD-bar (std dev of posterior distribution of tau)
mean(results[which(results[,"c"] == 0), "se.tau.hat."])
mean(results[which(results[,"c"] == 0.35), "se.tau.hat."])
mean(results[which(results[,"c"] == 0.7), "se.tau.hat."])

# MSE 
mean(results[which(results[,"c"] == 0), "sqderr"])
mean(results[which(results[,"c"] == 0.35), "sqderr"])
mean(results[which(results[,"c"] == 0.7), "sqderr"])

# SE
sqrt(sum((results[which(results[,"c"] == 0), "estimated.ATE"] - ate_c0)^2)/89)
sqrt(sum((results[which(results[,"c"] == 0.35), "estimated.ATE"] - ate_c035)^2)/89)
sqrt(sum((results[which(results[,"c"] == 0.7), "estimated.ATE"] - ate_c07)^2)/89)
