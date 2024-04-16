# Sim1 results

source("data_gen_processes_Zhu2023.R")

L=250

##################################################################
##                   Load results and combine                   ##
##################################################################

# load("raw_results_c000.Rdata")
# raw.results.c000 <- as.data.frame(raw.results.c000)
# raw.results.c000[,-c(2, 7)] <- apply(raw.results.c000[,-c(2, 7)], 2, as.numeric)
# raw.results.c000 <- raw.results.c000 %>% mutate(COVER.k = if_else(COVER.k=="FALSE", 0, 1))
# 
# 
# load("raw_results_c035.Rdata")
# raw.results.c035 <- as.data.frame(raw.results.c035)
# raw.results.c035[,-c(2, 7)] <- apply(raw.results.c035[,-c(2, 7)], 2, as.numeric)
# raw.results.c035 <- raw.results.c035 %>% mutate(COVER.k = if_else(COVER.k=="FALSE", 0, 1))
# 
# 
# load("raw_results_c070.Rdata")
# raw.results.c070 <- as.data.frame(raw.results.c070)
# raw.results.c070[,-c(2, 7)] <- apply(raw.results.c070[,-c(2, 7)], 2, as.numeric)
# raw.results.c070 <- raw.results.c070 %>% mutate(COVER.k = if_else(COVER.k=="FALSE", 0, 1))
# 
# all_results <- bind_rows(raw.results.c000, raw.results.c035, raw.results.c070)

load("all_results.Rdata")

table.1 <- all_results %>% group_by(c, method) %>%
  summarise(ATE=mean(ATE.k),
            Bias=mean(Psi.k-ATE.k),
            PercentBias=mean((Psi.k-ATE.k)/abs(ATE.k)*100),
            SD.bar = mean(SD.k),
            SE = sqrt(sum((Psi.k-ATE)^2)/(L-1)),
            MSE = sqrt(sum((Psi.k-ATE.k)^2)/(L-1)),
            Coverage= mean(COVER.k), 
            ci.cover = 1.96*sd(COVER.k)/sqrt(L)
  )


table.1


results.plot <- table.1 %>% 
  ggplot(aes(x = method, y = Coverage, group=method)) + 
  geom_hline(yintercept = .95, linetype = "dashed") + 
  geom_point(aes(color=method)) + 
  geom_errorbar(aes(ymin=Coverage-ci.cover, ymax=Coverage+ci.cover, color=method), width=0.5) + 
  facet_grid(. ~ c, scales = "free", labeller=labeller(label_both)) + 
  labs(x = "Method", y = "ATE Coverage (nominal 0.95)", 
       caption = "FIGURE 1. Coverage rate of five methods for ATE estimation. Dashed lines indicate the nominal confidence level.") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90), legend.position = "bottom",
        plot.caption = element_text(hjust = 0), legend.title = element_blank()) 
grid.arrange(results.plot, top="Overlap setting c")








