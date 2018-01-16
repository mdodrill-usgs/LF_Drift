################################################################################
#                                                                     
#
#         Fit count models (neg.binom) explore flow effects
#     
#
#                                                                     
################################################################################
setwd('C:/Users/mdodrill/Desktop/FOODBASE/LF_Drift/Data/')
rm(list = ls(all = TRUE))
# library(MASS)   # for glm.nb
# library(bbmle)  # AICtab - delta AIC table 
library(foodbase)
library(dplyr)
library(ggplot2)
library(xts)     # align.time (round to 15min) library(zoo)
library(lme4)
library(R2admb)
# install.packages("glmmADMB",repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")),type="source")
library(glmmADMB)


#-----------------------------------------------------------------------------#

all.dat = readDB()

all.dat$year = substr(all.dat$Date, 1, 4)

# don't need this because those are the only years in the db
# years = c('2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017')
# all.dat.2 = all.dat[all.dat$year %in% years,]

all.dat.2 = all.dat[which(all.dat$Reach == "CRLeesFerry"),]

dat = sampspec(samp = all.dat.2)




#-----------------------------------------------------------------------------#

# make jan 2008 --> 2007 
idx = which(substr(data4$DriftDate,1,7) == "2008-01") 

idx2 = which(substr(data4$DriftDate,1,10) == "2008-03-03")

data4$year = as.numeric(substr(data4$DriftDate,1,4))

data4[idx,which(names(data4) == "year")] = "2007"
data4[idx2,which(names(data4) == "year")] = "2007"

data3 = data4






