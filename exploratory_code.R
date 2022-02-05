library(ggplot2)
library(stringr)
lambda <- 0.2
theor_mean <- 1/lambda
theor_stdev <- 1/lambda
n <- 40
nsims <- 1000

samples <- matrix(rexp(n*nsims, rate = lambda), nrow = nsims, ncol = n)
means <- apply(samples, MARGIN = 1, FUN = mean)

avg_of_means <- mean(means) #as expected
stdev_of_means <- sd(means)  #as expected - se
var_of_means <- var(means) #this variance multiplied by sample size (40) approaches the population variance (5^2 = 25)

random_sample <-matrix(rexp(nsims, rate = lambda), nrow = nsims, ncol = 1)
avg_of_sample <- mean(random_sample)

png(filename="sample_avg_sim_dist.png")
ggplot(data = as.data.frame(means), mapping = aes(means)) + 
    geom_histogram(binwidth = 0.25, fill = "steelblue2", color = "black") +
    geom_vline(aes(xintercept = theor_mean, color = "Theoretical Mean"), size = 0.5, linetype = 2) +
    geom_vline(aes(xintercept = avg_of_means, color = "Mean of Sample Averages"), size = 0.5, linetype = 4) +
    xlab("Sample Average") +
    ylab("Count") +
    labs(title = str_wrap("Distribution of Sample Averages (n=40) from 1,000 Simulations of Samples from an Exponential Distribution (lambda = 0.2)",70))+
    scale_color_manual(name = "",
                       breaks = c('Theoretical Mean', 'Mean of Sample Averages'),
                       values = c('red', 'navyblue'))+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(2,8))
dev.off()

png(filename="sample_value_sim_dist.png")
ggplot(data = as.data.frame(random_sample), mapping = aes(random_sample)) + 
    geom_histogram(binwidth = 0.25, fill = "steelblue2", color = "black") +
    geom_vline(aes(xintercept = theor_mean, color = "Theoretical Mean"), size = 0.5, linetype = 2) +
    geom_vline(aes(xintercept = avg_of_sample, color = "Mean of Sample Values"), size = 0.5, linetype = 4) +
    xlab("Sampled Value") +
    ylab("Count") +
    labs(title = str_wrap("Distribution of Sample Values from 1,000 Simulations of Draws from an Exponential Distribution (lambda = 0.2)",70))+
    scale_color_manual(name = "",
                       breaks = c('Theoretical Mean', 'Mean of Sample Values'),
                       values = c('red', 'navyblue'))+
    theme(legend.position = "bottom")
dev.off()

install.packages("crosstable")
library(crosstable)
library(tidyverse)
library(knitr)
tooth_data <- ToothGrowth 
summary(tooth_data) #no missing values
str(tooth_data)
crosstable(ToothGrowth, supp, by = dose) %>% as_flextable(keep_id=FALSE)
# 60 obs total, 10 per dose/supp combo

#Check normality of length var
avg_len <- mean(ToothGrowth$len)
ggplot(data = ToothGrowth, mapping = aes(len)) + geom_histogram(binwidth = 5, fill = "steelblue", color = "black") + geom_vline(xintercept = avg_len)

#Check dist of length by dose/supp
ToothGrowth %>% group_by(dose, supp) %>% summarize(median(len), mean(len), sd(len))
#Maybe some differences between OJ/VC - particularly for lower dose, but maybe not all doses
#Higher doses do seem to have larger teeth length, particuarly from 0.5 to 1, but may differ by supp
#Conduct 9 pairwise t-tests
 # 3 comparing OJ and VC - 1 within each dose level
 # 3 comparing each dose pair within OJ
 # 3 comparing each dose pair within VC

#Two independent sample t-tests, assume unequal variance
doses <- unique(ToothGrowth$dose)
supps <- unique(ToothGrowth$supp)

#test_results <- data.frame(matrix(ncol = 5, nrow = 9))
#colnames(test_results) <- c("reference_group", "comparator_group", "lower_ci", "upper_ci", "raw_pvalue")
reference_group <- vector(length = 9)
comparator_group <- vector(length = 9)
lower_ci <- vector(length = 9)
upper_ci <- vector(length = 9)
raw_pvalue <- vector(length = 9)

n_tests <- 1
#loop over doses
for(d in doses){
    vc_dose <- ToothGrowth %>% subset(supp == 'VC' & dose == d, select = len)
    oj_dose <- ToothGrowth %>% subset(supp == 'OJ' & dose == d, select = len)
    
    reference_group[[n_tests]] <- paste0("vc_", d)
    comparator_group[[n_tests]] <- paste0("oj_", d)
    raw_pvalue[[n_tests]] <- t.test(oj_dose, vc_dose, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value
    lower_ci[[n_tests]] <- t.test(oj_dose, vc_dose, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[1]
    upper_ci[[n_tests]] <- t.test(oj_dose, vc_dose, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[2]
    n_tests <- n_tests+1
}

#loop over supps
for(s in supps){
    supp_0.5 <- ToothGrowth %>% subset(supp == s & dose == 0.5, select = len)
    supp_1 <- ToothGrowth %>% subset(supp == s & dose == 1, select = len)
    supp_2 <- ToothGrowth %>% subset(supp == s & dose == 2, select = len)

    reference_group[[n_tests]] <- paste0(s, "_0.5")
    comparator_group[[n_tests]] <- paste0(s, "_1")
    raw_pvalue[[n_tests]] <- t.test(supp_1, supp_0.5, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value
    lower_ci[[n_tests]] <- t.test(supp_1, supp_0.5, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[1]
    upper_ci[[n_tests]] <- t.test(supp_1, supp_0.5, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[2]
    n_tests <- n_tests+1
    
    reference_group[[n_tests]] <- paste0(s, "_0.5")
    comparator_group[[n_tests]] <- paste0(s, "_2")
    raw_pvalue[[n_tests]] <- t.test(supp_2, supp_0.5, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value
    lower_ci[[n_tests]] <- t.test(supp_2, supp_0.5, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[1]
    upper_ci[[n_tests]] <- t.test(supp_2, supp_0.5, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[2]
    n_tests <- n_tests+1
    
    reference_group[[n_tests]] <- paste0(s, "_1")
    comparator_group[[n_tests]] <- paste0(s, "_2")
    raw_pvalue[[n_tests]] <- t.test(supp_2, supp_1, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value
    lower_ci[[n_tests]] <- t.test(supp_2, supp_1, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[1]
    upper_ci[[n_tests]] <- t.test(supp_2, supp_1, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)$conf.int[2]
    n_tests <- n_tests+1
}

test_results <- data.frame(reference_group, comparator_group, lower_ci, upper_ci, raw_pvalue)
test_results$BH_pvalue <- p.adjust(test_results$raw_pvalue, method = "BH")
kable(test_results)
