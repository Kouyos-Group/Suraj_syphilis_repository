########################################################################################################
##################################### TO OPTIMIZE THE VALUE OF tau ####################################
########################################################################################################

#AUTHOR: Suraj Balakrishna

#RATIONALE: MAIN MODEL

##################################################  INPUT  ##############################################
rm(list = ls())

## Setting working directory
setwd("~/Downloads/STM_and_SRM/output")

## Setting up required libraries
library(deSolve)
library(plotly)
library(parallel)
library(Rcpp)
library(FME)
library(tidyr)

##  Load the Rcpp program i.e., the mathematical model written in C++
sourceCpp("~/Downloads/STM_and_SRM/programs/STM9_b_model_birth.cpp")

## Load required data
# Data with the number of HIV diagnosed MSM enrolled in the SHCS and diagnosed with syphilis
inc <- read.csv('STM2_c_rs_inc_data_new.csv')

# Syphilis notification data from the FOPH
inc_neg <- read.csv('STM2_f_inc_data_neg.csv')
inc <- inc[inc$year > 2005,]

# Population data from the Federal statistical office and inferred MSM population
pop_data <- read.csv("~/Downloads/STM_and_SRM/Literature/BAG_pop/pop_data.csv")
pop_data <- pop_data[pop_data$Year %in% c(2006:2018),]

# Yearly rate of switching from MSM without nsP to MSM with nsP for HIV diagnosed (sigma)
risk <- read.csv('STM2_e_survival_p_occas_sigma.csv')

# Yearly rate of switching from MSM with nsP to MSM with nsP for HIV diagnosed (kappa)
risk2 <- read.csv('STM2_e_survival_p_occas_kappa.csv')

# Proportion of MSM with nsP among MSM without HIV diagnosis
risk_neg <- read.csv('STM2_f_rs_french_doc.csv')

risk <- merge(risk, risk2, by='year')

risk <- risk[names(risk) %in% c('year', 'sigma', 'kappa')]
risk <- risk[risk$year %in% c(2006:2018),]

#Plot population dynamics
if(F){
  png("STM9_fpop_plot_MSM.png", width = 14, height = 8, units = 'in', res = 100)
  layout(matrix(nrow=2,ncol=2, 1:4),widths=rep(6,6),heights=rep(6,6))
  par(mfrow=c(2,2))
  par(mar=c(6, 18, 1, 1))
  plot(pop_data$Year, pop_data$MSM_pop, xlab='Year', ylab='', ylim = c(70000, 90000), type='l',cex.lab =2, lwd=3, col='red', cex.axis=1.5, cex=1.5, las=1)
  mtext('Estimated MSM\nin Switzerland',side=2,las=1, line=11, cex = 1.5, adj = 0.5)
  plot(pop_data$Year, pop_data$Birth, xlab='Year', ylab='', type='l', ylim=c(0, 1000), cex.lab =2, lwd=3, col='red', cex.axis=1.5, cex=1.5, las=1)
  mtext('Estimated inflow \nof MSM in \nSwitzerland (E)',side=2,las=1, line=10, cex = 1.5, adj = 0.5)
  plot(pop_data$Year, pop_data$Net_migration, xlab='Year', ylab='', type='l', ylim=c(0, 1000), cex.lab =2, lwd=3, col='red', cex.axis=1.5, cex=1.5, las=1)
  mtext('Estimated \nnet-migration \nof MSM in \nSwitzerland (M)',side=2,las=1, line=10, cex = 1.5, adj = 0.5)
  plot(pop_data$Year, pop_data$Death, xlab='Year', ylab='', type='l', ylim=c(0, 1000), cex.lab =2, lwd=3, col='red', cex.axis=1.5, cex=1.5, las=1)
  mtext('Estimated outflow \nof MSM in \nSwitzerland (O)',side=2,las=1, line=10, cex = 1.5, adj = 0.5)
  dev.off()
}

#Additional parameters
num_latin <- 20
startyear <- 2006
stopyear <- 2018
timestep <- 52
cores_used <- 1
overest_factor_value <- 1
new_corr_factor <- 1/0.8
inc <- inc[inc$year  %in% c((startyear):(stopyear)),]
inc_neg <- inc_neg[inc_neg$year %in% c((startyear):(stopyear)),]
numyears <- stopyear-startyear+1
inc$diag_all_pos <- inc$diag_first + inc$diag_second
inc_neg$msm_neg <- inc_neg$syph_msm - inc$diag_all_pos
inc$diag_all <- inc_neg$syph_msm
inc$diag_neg <- inc$diag_all - inc$diag_all_pos
risk$kappa_neg <- mean(risk$kappa)

##############################################  FUNCTIONS  ##############################################

#Function containing the parameters used in the model
parameter_values_generator <- function(beta_hiv, beta, risk_sort, overest_factor, theta){
  parameters <- c(beta_hiv=beta_hiv, beta=beta, risk_sort=risk_sort, 
                  overest_factor=overest_factor, theta=theta, kappa_neg=mean(risk$kappa),
                  deltaID_low_pos  = 16,  deltaLD_low_pos  = 1,    deltaIL_low_pos  = 12/3,  lambda_low_pos  = 12/1,
                  deltaID_high_pos = 16,  deltaLD_high_pos = 1,    deltaIL_high_pos = 12/3,  lambda_high_pos = 12/1,
                  deltaID_low_neg  = 16,  deltaLD_low_neg  = 0.5,  deltaIL_low_neg  = 12/3,  lambda_low_neg  = 12/1, 
                  deltaID_high_neg = 16,  deltaLD_high_neg = 0.5,  deltaIL_high_neg = 12/3,  lambda_high_neg = 12/1,
                  HR=3.167, sero_sort_pos=0.34, sero_sort_neg=0.91)
  return(parameters)
} 

#Function with starting state of the compartments of the model
starting_values_generator <- function(pos_inc, neg_inc, sero_naive, overest_factor){
  
  starting_values <- c(S1_low_pos  = new_corr_factor*(inc$pop_first_low[inc$year==2006]   - inc$diag_first_low[inc$year==2006]),                                                        I1_low_pos  = pos_inc,                                                                           L1_low_pos  = 0,  D1_low_pos  = 0,  
                       S2_low_pos  = new_corr_factor*(inc$pop_second_low[inc$year==2006]  - inc$diag_second_low[inc$year==2006]),                                                       I2_low_pos  = pos_inc*inc$diag_second_low[inc$year==2006] /inc$diag_first_low[inc$year==2006],   L2_low_pos  = 0,  D2_low_pos  = 0,   
                       S1_high_pos = new_corr_factor*(inc$pop_first_high[inc$year==2006]  - inc$diag_first_high[inc$year==2006]),                                                       I1_high_pos = pos_inc*inc$diag_first_high[inc$year==2006] /inc$diag_first_low[inc$year==2006],   L1_high_pos = 0,  D1_high_pos = 0,  
                       S2_high_pos = new_corr_factor*(inc$pop_second_high[inc$year==2006] - inc$diag_second_high[inc$year==2006]),                                                      I2_high_pos = pos_inc*inc$diag_second_high[inc$year==2006]/inc$diag_first_low[inc$year==2006],   L2_high_pos = 0,  D2_high_pos = 0, 
                       
                       S1_low_neg  = (pop_data$MSM_pop[pop_data$Year == 2006] - new_corr_factor*inc$all[inc$year==2006])*   sero_naive *(1-risk_neg$neg_high_prop[1]*overest_factor), I1_low_neg  = neg_inc,                                                                           L1_low_neg  = 0,  D1_low_neg  = 0,  
                       S2_low_neg  = (pop_data$MSM_pop[pop_data$Year == 2006] - new_corr_factor*inc$all[inc$year==2006])*(1-sero_naive)*(1-risk_neg$neg_high_prop[1]*overest_factor), I2_low_neg  = neg_inc,                                                                           L2_low_neg  = 0,  D2_low_neg  = 0,  
                       S1_high_neg = (pop_data$MSM_pop[pop_data$Year == 2006] - new_corr_factor*inc$all[inc$year==2006])*   sero_naive *   risk_neg$neg_high_prop[1]*overest_factor,  I1_high_neg = neg_inc,                                                                           L1_high_neg = 0,  D1_high_neg = 0, 
                       S2_high_neg = (pop_data$MSM_pop[pop_data$Year == 2006] - new_corr_factor*inc$all[inc$year==2006])*(1-sero_naive)*   risk_neg$neg_high_prop[1]*overest_factor,  I2_high_neg = neg_inc,                                                                           L2_high_neg = 0,  D2_high_neg = 0) 
  return(starting_values)
}

#Function which generates the initial parameter sets using a latin hypercube sampling algorithm.
random_initial_set_generator <- function(){
  
  parRange <- data.frame(min = c(1E-3, 1E-3, 1E-10, 1E-3, 1E-3,  0.5), max = c(30, 5, 0.1, 10, 10, 1))
  rownames(parRange) <- c("beta_hiv", "beta", "theta", "pos_inc", "neg_inc", "risk_sort")
  parameters_set <- Latinhyper(parRange, num_latin)
  return(parameters_set)
}

#Function to compute sigma for MSM without HIV diagnosis using spline function
sigma_neg_value_generator <- function(overest_factor){
  risk_neg$neg_high_prop_temp <<- risk_neg$neg_high_prop*overest_factor
  risk_neg$neg_low_prop_temp  <<- 1 - risk_neg$neg_high_prop_temp
  riskSpline <- splinefun(risk_neg$year, risk_neg$neg_low_prop_temp)
  risk$sigma_neg <<- risk$kappa_neg*(1- riskSpline(risk$year))/riskSpline(risk$year) -
    riskSpline(risk$year, deriv = T)/riskSpline(risk$year)
  return(risk$sigma_neg)
}

#Plotting the sigma_neg
if(F){
  sigma_neg_for_plot <- sigma_neg_value_generator(overest_factor = 1)
  sigma_neg_for_plot <- as.data.frame(sigma_neg_for_plot)
  sigma_neg_for_plot$year <- 2006:2018
  sigma_neg_for_plot$risk_neg <- risk$kappa_neg
  sigma_neg_for_plot <- sigma_neg_for_plot[sigma_neg_for_plot$year %in% c(2006:2017),]
  
  png('STM9_f_sigma_kappa.png', width = 14, height = 8, units = 'in', res = 100)
  par(mfrow=c(1,2))
  par(mar=c(6, 8, 1, 1))
  plot(sigma_neg_for_plot$year, sigma_neg_for_plot$sigma_neg_for_plot, xlab='Year', ylab='', ylim = c(0, 2), type='l',cex.lab =2, lwd=3, col='red', cex.axis=1.5, cex=1.5, las=1)
  mtext('sigma',side=2,las=1, line=5, cex = 1.5, adj = 0.5)
  plot(sigma_neg_for_plot$year, sigma_neg_for_plot$risk_neg, xlab='Year', ylab='', ylim = c(0, 2), type='l',cex.lab =2, lwd=3, col='red', cex.axis=1.5, cex=1.5, las=1)
  mtext('kappa',side=2,las=1, line=5, cex = 1.5, adj = 0.5)
  dev.off()
}

#Function returning the optimised output/ model prediction  
predictor_output_generator <- function(parameters, dataset){
  
  parms <- parameter_values_generator(beta_hiv = parameters['beta_hiv'], beta = parameters['beta'], risk_sort = parameters['risk_sort'], overest_factor = overest_factor_value, theta = parameters['theta'])
  
  inits <- starting_values_generator(pos_inc = parameters['pos_inc'], neg_inc = parameters['neg_inc'], sero_naive = 0.93, overest_factor = overest_factor_value)           
  
  sigma_value <- risk$sigma
  kappa_value <- risk$kappa
  sigma_neg_value <- sigma_neg_value_generator(overest_factor = overest_factor_value)
  
  birth <- pop_data$Birth
  death <- pop_data$Death
  migration <- pop_data$Net_migration
  
  times <- seq(startyear, stopyear, length= timestep*(numyears-1)+1)
  
  { sim1 <- bmv_modelQSS(timestep,
                         inits,
                         
                         parms[1],parms[2],parms[3],
                         parms[4],parms[5],parms[6],
                         parms[7],parms[8],parms[9],
                         parms[10],parms[11],parms[12],
                         parms[13],parms[14],parms[15],
                         parms[16],parms[17],parms[18],
                         parms[19],parms[20],parms[21],
                         parms[22],parms[23],parms[24],
                         parms[25],
                         
                         sigma_value,
                         kappa_value,
                         sigma_neg_value,
                         
                         birth,
                         death,
                         migration
                         
  )}
  
  add_on_names <- c('Diag_first_low', 'Diag_second_low', 'Diag_first_high', 'Diag_second_high',
                    'Diag_first', 'Diag_second', 'Diag_pos', 'Diag_neg', 'Diag_all',
                    'Diag_low_pos', 'Diag_high_pos', 'Diag_low_neg', 'Diag_high_neg')
  sim1 <- as.data.frame(t(sim1))
  #sim1 <- exp(sim1)
  colnames(sim1) <- c(names(inits), add_on_names)
  sim1$time <- times
  sim1$SHCS <- rowSums(sim1[names(sim1) %in%  c('S1_low_pos', 'I1_low_pos', 'L1_low_pos', 'D1_low_pos', 'S2_low_pos', 'I2_low_pos', 'L2_low_pos', 'D2_low_pos',
                                                'S1_high_pos', 'I1_high_pos', 'L1_high_pos', 'D1_high_pos', 'S2_high_pos', 'I2_high_pos', 'L2_high_pos', 'D2_high_pos')])
  
  sim1$denom_pos_low  <- rowSums(sim1[names(sim1) %in%  c('S1_low_pos',  'I1_low_pos',  'L1_low_pos',  'D1_low_pos',  'S2_low_pos',  'I2_low_pos',  'L2_low_pos',  'D2_low_pos')])
  sim1$denom_pos_high <- rowSums(sim1[names(sim1) %in%  c('S1_high_pos', 'I1_high_pos', 'L1_high_pos', 'D1_high_pos', 'S2_high_pos', 'I2_high_pos', 'L2_high_pos', 'D2_high_pos')])
  sim1$denom_neg_low  <- rowSums(sim1[names(sim1) %in%  c('S1_low_neg',  'I1_low_neg',  'L1_low_neg',  'D1_low_neg',  'S2_low_neg',  'I2_low_neg',  'L2_low_neg',  'D2_low_neg')])
  sim1$denom_neg_high <- rowSums(sim1[names(sim1) %in%  c('S1_high_neg', 'I1_high_neg', 'L1_high_neg', 'D1_high_neg', 'S2_high_neg', 'I2_high_neg', 'L2_high_neg', 'D2_high_neg')])
  
  sim1$total_pop <- rowSums(sim1[names(sim1) %in%  c('denom_pos_low', 'denom_pos_high', 'denom_neg_low', 'denom_neg_high')])
  
  yearIndice <- seq(1, (numyears-1)*timestep+1, timestep)
  sim1 <- sim1[yearIndice,]
  predicted_dataset_temp <- cbind(sim1[,grepl('Diag', names(sim1))])
  predicted_dataset <- predicted_dataset_temp[0,]
  for(i in 2:length(predicted_dataset_temp$Diag_first_low)){
    predicted_dataset[i-1,] <- predicted_dataset_temp[i,] - predicted_dataset_temp[i-1,]
  }
  predicted_dataset[length(predicted_dataset$Diag_first_low)+1,] <-   predicted_dataset[length(predicted_dataset$Diag_first_low)+1,]
  predicted_dataset$year <- startyear:(stopyear)
  predicted_dataset$SHCS <- rowSums(sim1[names(sim1) 
                                         %in%  c('S1_low_pos', 'I1_low_pos', 'L1_low_pos', 'D1_low_pos', 'S2_low_pos', 'I2_low_pos', 'L2_low_pos', 'D2_low_pos',
                                                 'S1_high_pos', 'I1_high_pos', 'L1_high_pos', 'D1_high_pos', 'S2_high_pos', 'I2_high_pos', 'L2_high_pos', 'D2_high_pos')])
  predicted_dataset$denom_pos_low  <- rowSums(sim1[names(sim1) %in%  c('S1_low_pos',  'I1_low_pos',  'L1_low_pos',  'D1_low_pos',  'S2_low_pos',  'I2_low_pos',  'L2_low_pos',  'D2_low_pos')])
  predicted_dataset$denom_pos_high <- rowSums(sim1[names(sim1) %in%  c('S1_high_pos', 'I1_high_pos', 'L1_high_pos', 'D1_high_pos', 'S2_high_pos', 'I2_high_pos', 'L2_high_pos', 'D2_high_pos')])
  predicted_dataset$denom_neg_low  <- rowSums(sim1[names(sim1) %in%  c('S1_low_neg',  'I1_low_neg',  'L1_low_neg',  'D1_low_neg',  'S2_low_neg',  'I2_low_neg',  'L2_low_neg',  'D2_low_neg')])
  predicted_dataset$denom_neg_high <- rowSums(sim1[names(sim1) %in%  c('S1_high_neg', 'I1_high_neg', 'L1_high_neg', 'D1_high_neg', 'S2_high_neg', 'I2_high_neg', 'L2_high_neg', 'D2_high_neg')])
  
  predicted_dataset$total_denom <- rowSums(sim1[names(sim1) %in%  c('S1_low_pos', 'I1_low_pos', 'L1_low_pos', 'D1_low_pos', 'S2_low_pos', 'I2_low_pos', 'L2_low_pos', 'D2_low_pos',
                                                                    'S1_high_pos', 'I1_high_pos', 'L1_high_pos', 'D1_high_pos', 'S2_high_pos', 'I2_high_pos', 'L2_high_pos', 'D2_high_pos',
                                                                    'S1_low_neg', 'I1_low_neg', 'L1_low_neg', 'D1_low_neg', 'S2_low_neg', 'I2_low_neg', 'L2_low_neg', 'D2_low_neg',
                                                                    'S1_high_neg', 'I1_high_neg', 'L1_high_neg', 'D1_high_neg', 'S2_high_neg', 'I2_high_neg', 'L2_high_neg', 'D2_high_neg')])
  
  predicted_dataset$inc_low_pos_per_1000 <- 1000*predicted_dataset$Diag_low_pos/predicted_dataset$denom_pos_low
  predicted_dataset$inc_high_pos_per_1000 <- 1000*predicted_dataset$Diag_high_pos/predicted_dataset$denom_pos_high
  
  predicted_dataset$inc_low_neg_per_1000 <- 1000*predicted_dataset$Diag_low_neg/predicted_dataset$denom_neg_low
  predicted_dataset$inc_high_neg_per_1000 <- 1000*predicted_dataset$Diag_high_neg/predicted_dataset$denom_neg_high
  
  predicted_dataset$inc_pos_per_1000 <- 1000*predicted_dataset$Diag_pos/(predicted_dataset$denom_pos_high + predicted_dataset$denom_pos_low)
  predicted_dataset$inc_neg_per_1000 <- 1000*predicted_dataset$Diag_neg/(predicted_dataset$denom_neg_high + predicted_dataset$denom_neg_low)
  
  predicted_dataset <- predicted_dataset[predicted_dataset$year %in% c(startyear:(stopyear-1)),]
  
  new_name_order <- c("year",   "Diag_first_low",  "Diag_second_low", "Diag_first_high",  "Diag_second_high",
                      "Diag_first",      "Diag_second",     "Diag_pos",        "Diag_neg",       
                      "Diag_all",        "Diag_low_pos",    "Diag_high_pos",   "Diag_low_neg",   
                      "Diag_high_neg",            "SHCS",            "denom_pos_low",  
                      "denom_pos_high",  "denom_neg_low",   "denom_neg_high",  "total_denom",
                      'inc_low_pos_per_1000', 'inc_high_pos_per_1000', 'inc_low_neg_per_1000', 'inc_high_neg_per_1000',
                      'inc_pos_per_1000', 'inc_neg_per_1000')
  predicted_dataset <- predicted_dataset[,(new_name_order)]
  
  #80% of all HIV-positive MSM are in the SHCS.
  predicted_dataset[ , c("Diag_first_low",  "Diag_second_low", "Diag_first_high",  "Diag_second_high",
                         "Diag_first",      "Diag_second",     "Diag_pos",      
                         "Diag_low_pos",    "Diag_high_pos",  "SHCS", 
                         "denom_pos_low",  "denom_pos_high")] <-   
    predicted_dataset[ , c("Diag_first_low",  "Diag_second_low", "Diag_first_high",  "Diag_second_high",
                           "Diag_first",      "Diag_second",     "Diag_pos",         
                           "Diag_low_pos",    "Diag_high_pos",  "SHCS", 
                           "denom_pos_low",  "denom_pos_high")]/new_corr_factor
  
  
  return(predicted_dataset)
}

#Function returns -2*loglikelihood based on Poisson distribution
loglik_mini <- function(x,y){
  sum(-2*(y*log(x) - x -lfactorial(y)))
}

#Function computes total of -2*loglikelihood between observed and predicted datapoints 
loglik <- function(parameters){
  predicted <- predictor_output_generator(parameters)
  out <- loglik_mini(predicted$Diag_first_low[1:(2017-startyear+1)],  inc$diag_first_low[1:(2017-startyear+1)]) + 
    loglik_mini(predicted$Diag_first_high[1:(2017-startyear+1)], inc$diag_first_high[1:(2017-startyear+1)]) + 
    loglik_mini(predicted$Diag_second[1:(2017-startyear+1)],     inc$diag_second[1:(2017-startyear+1)]) +
    loglik_mini(predicted$Diag_all[1:(2017-startyear+1)],        inc$diag_all[1:(2017-startyear+1)]) +
    loglik_mini(predicted$SHCS[1:(2017-startyear+1)],            inc$all[1:(2017-startyear+1)])
  #print(c(parameters,out))
  if(is.nan(out)){out <- 100000}
  if(out==Inf){out <- 100000}
  return(out) #loglik value is returned
}

#Computes the sum of squared weighted residuals between observed and predicted datapoints 
mod_cost <- function(parameters){
  predicted <- predictor_output_generator(parameters)
  predicted <- predicted[predicted$year %in% c(2006:2017),]
  observed <- inc[inc$year %in% c(2006:2017),]
  observed <- observed[names(observed) %in% c("diag_first_low",  "diag_first_high",  "diag_second",
                                              "diag_all", "all", "year")]
  #observed$all <- observed$all/20
  
  names(predicted) <- c("year", "diag_first_low",  "diag_second_low",  "diag_first_high",  "diag_second_high",
                        "diag_first",       "diag_second",      "diag_pos",         "diag_neg",
                        "diag_all",        'Diag_low_pos', 'Diag_high_pos', 'Diag_low_neg', 'Diag_high_neg',
                        "all")
  #predicted$all <- predicted$all/20
  predicted <- predicted[names(predicted) %in% c("diag_first_low",  "diag_first_high",  "diag_second",
                                                 "diag_all", "all", "year")]
  predicted <- predicted[, names(observed)]
  
  predicted_long <- gather(predicted, name, val, diag_second:diag_all)
  predicted_long <- predicted_long[, c('name', 'year', 'val')]
  #names(predicted_long) <- c('name', 'time', 'val')
  #predicted_long <- predicted_long[names(predicted_long) %in% c('year', 'val')]
  
  
  observed_long <- gather(observed, name, val, diag_second:diag_all)
  observed_long <- observed_long[, c('name', 'year', 'val')]
  #names(observed_long) <- c('name', 'time', 'val')
  #observed_long <- observed_long[names(observed_long) %in% c('year', 'val')]
  observed_long$err <- NA
  observed_long$err[observed_long$name == 'diag_second'] <- sqrt(predicted$diag_second)
  observed_long$err[observed_long$name == 'diag_first_high'] <- sqrt(predicted$diag_first_high)
  observed_long$err[observed_long$name == 'diag_first_low'] <- sqrt(predicted$diag_first_low)
  observed_long$err[observed_long$name == 'diag_all'] <- sqrt(predicted$diag_all)
  observed_long$err[observed_long$name == 'all'] <- sqrt(predicted$all)
  
  cost <- modCost(model=predicted, obs=observed_long, x = 'year', y='val', err='err')
  #print(c(parameters, cost[["model"]]),scientific = F)
  return(cost)
}

#Function to fit the model
mod_function <- function(parameters){
  modFit(f=mod_cost, parameters, method = 'L-BFGS-B',
         lower = c(1E-10, 1E-10, 1E-10, 1E-3, 1E-3, 0.5), upper = c(30, 5, 0.1, 10, 10, 1),
         control = list(pgtol=1E-5, maxit=30000))
}

####################### Creating MCMC CODA CHAIN  ##############################################

load('STM9_f_optim_par_mc.Rdata')
load('STM9_f_optim_par_set.Rdata')

lower_bound <- c(1E-3, 1E-3, 1E-10, 1E-3, 1E-3, 0.5)
upper_bound <- c(30, 5, 0.1, 10, 10, 1)
jump_value <- (upper_bound - lower_bound)/1000

mc_list <- list()
mc_list[[1]] <- (mc)

for(i in 1:2){
  mc <- modMCMC(f=mod_cost, p = data_16$par,  
                updatecov = 10 , 
                ntrydr = 2,
                jump = jump_value,
                var0 = data_16$var_ms_unweighted,
                lower = lower_bound, upper = upper_bound,
                niter = 10000)
  mc_list[[i+1]] <- mc
  cat(i, '\n')
}

mc_coda1 <- as.mcmc(mc_list[[1]]$pars)
mc_coda2 <- as.mcmc(mc_list[[2]]$pars)
mc_coda3 <- as.mcmc(mc_list[[3]]$pars)
mc_coda_list <- list(mc_coda1, mc_coda2, mc_coda3)
mc_coda_list <- as.mcmc.list(mc_coda_list)

save(mc_coda_list,    file='STM9_gelman_mc_coda_list.Rdata')


#####################      KEY OUTPUT OF THE CODA MODEL   #######################################
load('STM9_gelman_mc_coda_list.Rdata')

#MAIN RESULT
gelman.diag(mc_coda_list)

png("~/Downloads/STM_and_SRM/output/STM9_gelman_mc_coda_trace.png", width = 12, height = 8, units = 'in', res = 100)
plot(mc_coda_list, density = F, col=c('black', 'red', 'yellow'))
dev.off()

plot(mc_coda_list, trace = F, col=c('black', 'red', 'yellow'))

