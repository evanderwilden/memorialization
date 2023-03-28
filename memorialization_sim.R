#########################
# Ethan vanderWilden
# Simulating Data
# Memorialization and the Far Right (with Laia Balcells)

# March 2023
#########################

############## 1. Load in packages ####
# List of packages
pkg = c("tidyverse", "gridExtra", "haven", "modelsummary")

# Checks if they are installed, install if not
if (length(setdiff(pkg, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkg, rownames(installed.packages())))
}
# Load
lapply(pkg, library, character.only = TRUE)

############ 2. Load and clean the data ####
setwd('~/Projects/memorialization')
toy <- read.csv("memorialization_toy.csv")
attributes(toy)$subheaders <- toy[1,]
toy <- toy[-c(1:5), -c(1:17)]


########## 3. Validating baseline assumptions ####
cis <-read_sav('dec19.sav')

# from CIS, what is approval for Abascal, what is SD?
#around 2 and 3 seem appropriate
cis$B29_3 <- ifelse(cis$B29_3>10, NA, cis$B29_3)
mean(cis$B29_3, na.rm = T)
sd(cis$B29_3, na.rm = T)


#cis$C8 <- ifelse(cis$C8>95, NA, cis$C8)
#cis$C8 <- ifelse(cis$C8 == 18, 1, 0)
#mean(cis$C8, na.rm = T)
#sd(cis$C8, na.rm = T)


############# 4. Simulate Data ####

#write function to simulate data
simulate_data <- function(N, effects){
  #initialize data frame
  df <- data.frame(matrix(ncol = 35, nrow = N))
  names(df) <- names(toy)
  
  #covariates
  df$sex <- sample(c(1,2), N, replace = T)
  df$age <- sample(c(1:6), N, replace = T)
  df$education <- sample(c(1:6), N, replace = T)
  df$region_1 <- sample(c(1:19), N, replace = T)
  df$ideology_1 <- sample(c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7), N, replace = T)
  
  #assign treatment
  df$treat <- sample(c(0:6), N, replace = T)
  
  for (i in 1:nrow(df)){
    if (df$treat[i] == 0){
      df$inf[i] <- 0
      df$edu[i] <- 0
      df$pd[i] <- 0
      df$additive[i] <- 0
      df$removal[i] <- 0
    } else if (df$treat[i] == 1){
      df$inf[i] <- 1
      df$edu[i] <- 0
      df$pd[i] <- 0
      df$additive[i] <- 1
      df$removal[i] <- 0
    } else if (df$treat[i] == 2){
      df$inf[i] <- 1
      df$edu[i] <- 0
      df$pd[i] <- 0
      df$additive[i] <- 0
      df$removal[i] <- 1
    } else if (df$treat[i] == 3){
      df$inf[i] <- 0
      df$edu[i] <- 1
      df$pd[i] <- 0
      df$additive[i] <- 1
      df$removal[i] <- 0
    } else if (df$treat[i] == 4){
      df$inf[i] <- 0
      df$edu[i] <- 1
      df$pd[i] <- 0
      df$additive[i] <- 0
      df$removal[i] <- 1
    } else if (df$treat[i] == 5){
      df$inf[i] <- 0
      df$edu[i] <- 0
      df$pd[i] <- 1
      df$additive[i] <- 1
      df$removal[i] <- 0
    } else {
      df$inf[i] <- 0
      df$edu[i] <- 0
      df$pd[i] <- 1
      df$additive[i] <- 0
      df$removal[i] <- 1
    }
  }
  
  #Simulate treatment effects
  ### hypothesize average sympathy for vox is 3
  base <- effects[1]
  # hypothesize average treatment effect of additive policy
  additive_tau <- effects[2]
  # hypothesize average treatment effect of removal policy
  removal_tau <- effects[3]
  # hypothesize standard deviation of the outcome
  sd_outcome <- effects[4]
  
  #simulate main outcome
  for (i in 1:nrow(df)){
    df$sympathy_4[i] <- rnorm(1, mean = (base + additive_tau*df$additive[i] + removal_tau*df$removal[i]), sd = sd_outcome)
  }
  
  #keep all outcomes within bounds of possibility (0:10)
  for (i in 1:nrow(df)){
    if(df$sympathy_4[i] < 0){
      df$sympathy_4[i] <- 0
    } else if(df$sympathy_4[i] > 10){
      df$sympathy_4[i] <- 10
    } else{
      df$sympathy_4[i]
    }
  }
  
  return(df)
}

set.seed(1234)


#Simulate one effect size
effects <- c(3, -0.4, 0.4, 2)

#test simulation
sim <- simulate_data(2800, effects)

#Extract data and plot (looks like things are working)
sympathy_coef <- as.data.frame(coef(summary(lm(sympathy_4 ~ additive + removal, 
                                 data = subset(sim, sim$treat <=2)))))

names(sympathy_coef) <- c("estimate", "sd", "tval", "pval")
sympathy_coef$variable <- row.names(sympathy_coef)

ggplot(data = sympathy_coef[-c(1),]) +
  geom_pointrange(aes(x = estimate, xmin = estimate-1.96*sd, xmax = estimate + 1.96*sd, y = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(title = "Regression coefficients (single group compared to control)",
       y = "variable", 
       x = "treatment effect")


################## 5. Power analysis ####

#create df to store results
simulations <- data.frame(
  n = seq(100, 5000, 500),
  sig_results_4  = rep(NA, length(n)),
  pooled_AMCE_4 = rep(NA, length(n)),
  sig_results_2 = rep(NA, length(n)),
  pooled_AMCE_2 = rep(NA, length(n))
)

#hypothesize [LARGER] effect sizes
effects <- c(3, -0.8, 0.8, 2)

#simulate results, calculate power
for(i in 1:nrow(simulations)){
  sig_ind <- 0
  sig_pooled <- 0
  for (j in 1:100){ #note: just running 100 simulations for now
    sim <- simulate_data(simulations[i,1], effects)
    sig_ind <- sig_ind + ifelse(coef(summary(lm(sympathy_4 ~ additive, data = subset(sim, sim$treat <=1))))[2,4] < 0.10, 1, 0)
    sig_pooled <- sig_pooled + ifelse(coef(summary(lm(sympathy_4 ~ additive, data = sim)))[2,4] < 0.10, 1, 0)
    }
  simulations[i,2] <- sig_ind/100
  simulations[i,3] <- sig_pooled/100
}


#hypothesize [smaller] effect sizes
effects <- c(3, -0.4, 0.4, 2)

#simulate results, calculate power
for(i in 1:nrow(simulations)){
  sig_ind <- 0
  sig_pooled <- 0
  for (j in 1:100){ #note: just running 100 simulations for now
    sim <- simulate_data(simulations[i,1], effects)
    sig_ind <- sig_ind + ifelse(coef(summary(lm(sympathy_4 ~ additive, data = subset(sim, sim$treat <=1))))[2,4] < 0.10, 1, 0)
    sig_pooled <- sig_pooled + ifelse(coef(summary(lm(sympathy_4 ~ additive, data = sim)))[2,4] < 0.10, 1, 0)
  }
  simulations[i,4] <- sig_ind/100
  simulations[i,5] <- sig_pooled/100
}


############# 6. Plot power calculations ####

cols <- c("Individual Comparison [H1] \n (d = 0.4)"= "red",
         "Pooled Comparison [H1] \n (d = 0.4)"= "black",
         "Individual Comparison [H2/3] \n (d = 0.2)"= "blue",
         "Pooled Comparison [H2/3] \n (d = 0.2)"= "orange")

ggplot(data = simulations) +
  geom_point(aes(x = n, y = sig_results_4, color = "Individual Comparison [H1] \n (d = 0.4)"), shape = 16, size = 3) + 
  geom_smooth(aes(x = n, y = sig_results_4, color ="Individual Comparison [H1] \n (d = 0.4)"), se = F, method = "loess") +
  
  geom_point(aes(x = n, y = pooled_AMCE_4, color = "Pooled Comparison [H1] \n (d = 0.4)"), shape = 16, size = 3) + 
  geom_smooth(aes(x = n, y = pooled_AMCE_4, se = F, color ="Pooled Comparison [H1] \n (d = 0.4)"), se = F, method = "loess") +
  
  geom_point(aes(x = n, y = sig_results_2, color = "Individual Comparison [H2/3] \n (d = 0.2)"), shape = 16, size = 3) + 
  geom_smooth(aes(x = n, y = sig_results_2, se = F, color ="Individual Comparison [H2/3] \n (d = 0.2)"), se = F,  method = "loess") +
  
  geom_point(aes(x = n, y = pooled_AMCE_2, color = "Pooled Comparison [H2/3] \n (d = 0.2)"), shape = 16, size = 3) + 
  geom_smooth(aes(x = n, y = pooled_AMCE_2, se = F, color ="Pooled Comparison [H2/3] \n (d = 0.2)"), se = F,  method = "loess") +
  
  geom_hline(yintercept = 0.8, linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(0, 4800), breaks = seq(0, 4500, 500))+
  scale_color_manual(values = cols)+
  theme_bw()+
  labs(title = "Simulating significance for H1-3",
       x = "Sample Size Required (total)", y = "Power") +
  theme(legend.title = element_blank())




