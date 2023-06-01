#########################
# Ethan vanderWilden
# Simulating Data
# Memorialization and the Far Right (with Laia Balcells)

# May 2023
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
setwd('~/memorialization')
toy <- read.csv("power_sim.csv")
attributes(toy)$subheaders <- toy[1,]
toy <- toy[-c(1:2), -c(1:17)]


########## 3. Validating baseline assumptions ####
cis <-read_sav('replications/dec19.sav')


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

#write function to keep data within bounds
inbounds <- function(issue, min, max) {
  for (i in 1:length(issue)){
    if (issue[i] < min){
      issue[i] <- min + sample(seq(0, 1, 0.01), 1)
    } else if (issue[i] > max){
      issue[i] <- max - sample(seq(0, 1, 0.01), 1)
    }
  }
  return(issue)
}

#write function to simulate data
simulate_data <- function(N, r_effects, l_effects){
  #initialize data frame
  df <- data.frame(matrix(ncol = 34, nrow = N))
  names(df) <- names(toy)
  
  #covariates
  df$sex <- sample(c(1,2), N, replace = T)
  df$age <- sample(c(1:6), N, replace = T)
  df$education <- sample(c(1:6), N, replace = T)
  df$region_1 <- sample(c(1:19), N, replace = T)
  df$ideology_1 <- inbounds(rnorm(N, 4, 1.5), 1, 7) #normal distribution corrected for bounds
  
  #assign treatment
  df$treat <- sample(c(0:4), N, replace = T)
  
  for (i in 1:nrow(df)){
    if (df$treat[i] == 0){
      df$inf[i] <- 0
      df$edu[i] <- 0
      df$victim[i] <- 0
      df$perpetrator[i] <- 0
      df$control[i] <- 1
      df$inf_victim[i] <- 0
      df$inf_perp[i] <- 0
      df$edu_victim[i] <- 0
      df$edu_perp[i] <- 0
    } else if (df$treat[i] == 1){
      df$inf[i] <- 1
      df$edu[i] <- 0
      df$victim[i] <- 1
      df$perpetrator[i] <- 0
      df$control[i] <- 0
      df$inf_victim[i] <- 1
      df$inf_perp[i] <- 0
      df$edu_victim[i] <- 0
      df$edu_perp[i] <- 0
    } else if (df$treat[i] == 2){
      df$inf[i] <- 1
      df$edu[i] <- 0
      df$victim[i] <- 0
      df$perpetrator[i] <- 1
      df$control[i] <- 0
      df$inf_victim[i] <- 0
      df$inf_perp[i] <- 1
      df$edu_victim[i] <- 0
      df$edu_perp[i] <- 0
    } else if (df$treat[i] == 3){
      df$inf[i] <- 0
      df$edu[i] <- 1
      df$victim[i] <- 1
      df$perpetrator[i] <- 0
      df$control[i] <- 0
      df$inf_victim[i] <- 0
      df$inf_perp[i] <- 0
      df$edu_victim[i] <- 1
      df$edu_perp[i] <- 0
    } else {
      df$inf[i] <- 0
      df$edu[i] <- 1
      df$victim[i] <- 0
      df$perpetrator[i] <- 1
      df$control[i] <- 0
      df$inf_victim[i] <- 0
      df$inf_perp[i] <- 0
      df$edu_victim[i] <- 0
      df$edu_perp[i] <- 1
    }
  }
  
  #Simulate treatment effects (right leaning)
  ### hypothesize average sympathy for vox is 3
  base_r <- r_effects[1]
  
  # hypothesize average treatment effect of infrastructure victim
  inf_victim_r <- r_effects[2]
  # hypothesize average treatment effect of infrastructure perpetrator
  inf_perp_r <- r_effects[3]
  
  # hypothesize average treatment effect of infrastructure victim
  edu_victim_r <- r_effects[4]
  # hypothesize average treatment effect of infrastructure perpetrator
  edu_perp_r <- r_effects[5]
  
  # hypothesize standard deviation of the outcome
  sd_outcome_r <- r_effects[6]
  
  #Simulate treatment effects (left leaning)
  ### hypothesize average sympathy for vox is 3
  base_l <- l_effects[1]
  
  # hypothesize average treatment effect of infrastructure victim
  inf_victim_l <- l_effects[2]
  # hypothesize average treatment effect of infrastructure perpetrator
  inf_perp_l <- l_effects[3]
  
  # hypothesize average treatment effect of infrastructure victim
  edu_victim_l <- l_effects[4]
  # hypothesize average treatment effect of infrastructure perpetrator
  edu_perp_l <- l_effects[5]
  
  # hypothesize standard deviation of the outcome
  sd_outcome_l <- l_effects[6]
  
  
  ### Subset data by left and right
  right <- subset(df, df$ideology>4)
  left <- subset(df, df$ideology <4)
  
  
  #simulate main outcome - RIGHTISTS
  for (i in 1:nrow(right)){
    right$sympathy_3[i] <- rnorm(1, mean = (base_r +
                                              inf_victim_r*right$inf_victim[i] +
                                              inf_perp_r*right$inf_perp[i] +
                                              edu_victim_r*right$edu_victim[i] +
                                              edu_perp_r*right$edu_perp[i]), sd = sd_outcome_r)
  }
  
  #simulate main outcome - LEFTISTS
  for (i in 1:nrow(left)){
    left$sympathy_3[i] <- rnorm(1, mean = (base_l +
                                              inf_victim_l*left$inf_victim[i] +
                                              inf_perp_l*left$inf_perp[i] +
                                              edu_victim_l*left$edu_victim[i] +
                                              edu_perp_l*left$edu_perp[i]), sd = sd_outcome_r)
  }
  
  
  
  #keep all outcomes within bounds of possibility (0:100)
  right$sympathy_3 <- inbounds(right$sympathy_3, 0, 100)
  left$sympathy_3 <- inbounds(left$sympathy_3, 0, 100)
  
  return(list(right = right, left = left))
}

set.seed(1234)

# c(baseline, inf-victim, inf-perp, edu-victim, edu-perp, sd)
r_effects <- c(65, 5, 9, 8, 12, 15)
l_effects <- c(20, -5, -9, -8, -12, 15)
simulate_data(1000, r_effects, l_effects)

#### SO: We want to ...
# (1) simulate these
# (2) correct for multiple hypotheses
# (3) re-enter when we do the pilot


################## 5. Power analysis ####
#create df to store results
simulations <- data.frame(
  n = seq(100, 6000, 500),
  all_sig_2 = rep(NA, length(n))
)

#hypothesize d = 0.25 (for smallest effect, and for both ACMEs)
r_effects <- c(65, 4, 8, 8, 12, 16)
l_effects <- c(20, -5, -9, -8, -12, 15) #doesn't matter here


#see if all six estimates are significant (return 1 if all are sig, 0 if else)
calc_full <- function(data){
  #store p-values
  results <- 
    c((coef(summary(lm(sympathy_3 ~ inf_victim, 
                       data = subset(data$right, data$right$treat %in% c(0,1)))))[2,4]),
      (coef(summary(lm(sympathy_3 ~ inf_perp, 
                       data = subset(data$right, data$right$treat %in% c(0,2)))))[2,4]),
      (coef(summary(lm(sympathy_3 ~ edu_victim, 
                       data = subset(data$right, data$right$treat %in% c(0,3)))))[2,4]),
      (coef(summary(lm(sympathy_3 ~ edu_perp, 
                       data = subset(data$right, data$right$treat %in% c(0,4)))))[2,4]),
      (coef(summary(lm(sympathy_3 ~ perpetrator, 
                       data = subset(data$right, data$right$treat %in% c(1:4)))))[2,4]),
      (coef(summary(lm(sympathy_3 ~ edu, 
                       data = subset(data$right, data$right$treat %in% c(1:4)))))[2,4]))
  #sort from low to high
  results <- sort(results)
  
  #add Bonferroni correction
  ifelse(length(which(results < c(.1/6, .1/5, .1/4, .1/3, .1/2, .1))) == 6, 1, 0)
}

#simulate 100 data draws
corrected_full <- function(simulations, first, last){
  #simulate results, calculate power
  for(i in first:last){
    sig_number <- 0
    for (j in 1:100){ #note: just running 100 simulations for now
      sim <- simulate_data(simulations[i,1], r_effects, l_effects)
      sig_number <- sig_number + calc_full(sim)
    }
    simulations[i,2] <- sig_number/100 #calculate proportion that was significant
  }
  
  return(simulations)
}

#run simulations
simulations <- corrected_full(simulations, 1, 4) #first 4 simulations
simulations <- corrected_full(simulations, 5, 8) #next 4 simulations
simulations <- corrected_full(simulations, 9, 12) #final 4 simulations


#Plot simulations
ggplot(data = simulations) +
  geom_point(aes(x = n, y = all_sig_2), shape = 16, size = 3) +
  geom_line(aes(x = n, y = all_sig_2)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(0, 6100), breaks = seq(0, 6000, 1000))+
  theme_bw()+
  labs(title = "Simulating significance for H1, 3, and 4", 
       subtitle = "(BH correction for multiple hypotheses)",
       x = "Sample Size Required (total)", y = "Power") +
  theme(legend.title = element_blank())



