###############################################################################
#                                                                    April 2017
#           Stan Model Fitting for Lees Ferry State Space Model 
#
#  Notes:
#  * 
#
###############################################################################
setwd('C:/Users/mdodrill/Desktop/FB_Git/LF_Drift')
rm(list = ls(all = TRUE))
library(dplyr)
# library(ggplot2)

# load functions stored in R.script
source("C:/Users/mdodrill/Desktop/FB_Git/Drift/FB_Drift_Functions_Git_V1.R", chdir = F)
#-----------------------------------------------------------------------------#
# see look -- these are the species to do summaries for....
# Right now the functions are only 'wired' for LUM, NZM, GAM, SIM, CHI
sp_num = unique(look[,1])    

# Data to import (these are global and not passed in)
drift_specimen = 'DRIFT_SPECIMEN_01_18_17.csv'
drift_sample   = 'DRIFT_SAMPLE_01_18_17.csv'

data = tbl_df(drift_summary(c(2012,2013,2014,2015), sp_num))

# need the section for getting/prepping the Lees Ferry data....


#-----------------------------------------------------------------------------#
params <- c("beta_o") 


data.in = list( int Nmd;                 # // Number of drift months from start   
                int Nmb;                 # // Number of benthic months from start
                int Nmonth;              # // Number of months (12)
                int Nsamps;              # // Number of drift samples 
                vector[Nmd] d_trip;       #// Index for month from start (drift)
                int month[Nmd];            #// Index for month in year (1 - 12)
                int DC[Nsamps, Nmd];        #// Drift conc. 
                matrix[Nsamps, Nmd] log_Q;   #// log Q)

## MCMC settings
ni <- 2000
nt <- 1
nb <- 500
nc <- 3

fit <- stan("LF_State_Space_V1.stan", 
            data = data.in, 
            pars = params,
            # init = inits,
            # control = list(adapt_delta = .9),
            # control = list(max_treedepth = 12),
            control = list(max_treedepth = 14, adapt_delta = .95),
            chains = nc, thin = nt, iter = ni, warmup = nb)

#-----------------------------------------------------------------------------#
# print(get_elapsed_time(fit))/60/60
# Rhat = rstan::summary(fit)$summary[,"Rhat"]
# 
# bad = Rhat[which(Rhat > 1.2)]
# length(bad)
# length(Rhat)

# windows(record = T)
# source("C:/Users/mdodrill/Desktop/Multinomial/stan_plotting.R", chdir = F)

# my_plots(fit, "beta_sz", 1:7)
# my_plots(fit, "fix_lprod",  list(rep(1,7),c(1:7)))
# 
# library("shinystan")
# launch_shinystan(fit)