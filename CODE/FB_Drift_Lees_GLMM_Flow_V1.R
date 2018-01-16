################################################################################
#                                                                       Jan 2018
#
#         Fit count models (neg.binom) explore flow effects
#     
#
#                                                                     
################################################################################
# setwd('C:/Users/mdodrill/Desktop/FOODBASE/LF_Drift/Data/')
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

#--------------------------------------
keep.samp = c("BarcodeID", "TripID", "Date", "TimeBegin", "RiverMile", "DepthTotal", "DepthSample",
              "DepthIntegrated", "TimeDay", "Volume")

samp = dat$Samples

# cut down the sample table
ltl.samp = samp[,which(names(samp) %in% keep.samp)] 

#--------------------------------------
# only the taxa that we want
spec = dat$Specimens

taxa = c("LUMB", "NZMS", "GAM", "SIMA", "SIMP", "SIML", "CHIA", "CHIP", "CHIL")

ltl.spec = spec[which(spec$SpeciesID %in% taxa),]

#--------------------------------------
# merge the sample table with the specimen counts

tmp.dat = left_join(ltl.spec, ltl.samp, by = "BarcodeID")

lf.dat = tmp.dat[,c(1,24:32,2:23)]

lf.dat$SpeciesID = as.character(lf.dat$SpeciesID)

rm(list = setdiff(ls(), "lf.dat"))

#-----------------------------------------------------------------------------#
# make jan 2008 --> 2007 
idx = which(substr(lf.dat$Date, 1, 7) == "2008-01") 

idx2 = which(substr(lf.dat$Date, 1, 10) == "2008-03-03")

lf.dat$year = as.numeric(substr(lf.dat$Date, 1, 4))

lf.dat[idx, which(names(lf.dat) == "year")] = "2007"
lf.dat[idx2, which(names(lf.dat) == "year")] = "2007"

#------------------------------------------------------------------------------#
# Group the life stages of midges and blackflies together

lf.dat$SpeciesID = ifelse(lf.dat$SpeciesID == "SIMA", "SIML", lf.dat$SpeciesID)
lf.dat$SpeciesID = ifelse(lf.dat$SpeciesID == "SIMP", "SIML", lf.dat$SpeciesID)

lf.dat$SpeciesID = ifelse(lf.dat$SpeciesID == "CHIA", "CHIL", lf.dat$SpeciesID)
lf.dat$SpeciesID = ifelse(lf.dat$SpeciesID == "CHIP", "CHIL", lf.dat$SpeciesID)

unique(lf.dat$SpeciesID)

#-----------------------------------------------------------------------------#
# adding in the discharge data 

lf.q = read.csv('LFQ_2007_2016.csv', header = T)  

# round the drift sampling time to 15 mins 
time = as.POSIXct(strptime(as.character(lf.dat$TimeBegin), "%H:%M"))  
tmp.time = substr(align.time(time, n = 60 * 15), 11, 19) # n is in seconds (rounds up)
lf.dat$date.time = as.POSIXct(paste(lf.dat$Date, tmp.time), format = "%Y-%m-%d %H:%M:%S")

# format the discharge time column to match the bug data
lf.q$time = as.POSIXct(strptime(as.character(lf.q$time),"%m/%d/%Y  %H:%M"))
lf.q$year = substr(lf.q$time, 1, 4)

# add a new col. with the coresponding discharges.
lf.dat$q = lf.q[match(lf.dat$date.time, lf.q$time),2]
lf.dat$q = lf.dat$q * .028317

# !!! some of the discharges are NA !!! ---> Need to fix this....

#-----------------------------------------------------------------------------#
# add in some other stuff... organize for the model 

# trying to avoid referencing the cols by position, instead get the names
cols = which(!is.na(as.numeric(gsub("[^0-9]", "", names(lf.dat)))))

# total counts 
lf.dat$count.sum = apply(lf.dat[,cols], 1, sum)

lf.dat$month = substr(lf.dat$Date, 6, 7)

lf.dat$RM = as.factor(lf.dat$RiverMile)

lf.dat$year = as.factor(lf.dat$year)

#-----------------------------------------------------------------------------#
# fit some models to look at the effect of steady vs. fluc. flows
st.yr = c(2008, 2009, 2010, 2011, 2012)

dat.1 = lf.dat[which(lf.dat$SpeciesID == "CHIL"),]  

dat.2 = dat.1[which(dat.1$month %in% c("08", "09", "10")),]

# dat.3 = dat.2[which(dat.2$year %in% st.yr),]
# dat.4 = dat.2[which(!dat.2$year %in% st.yr),]
# dat.3$mo.num = as.numeric(dat.3$month) - 9
# dat.4$mo.num = as.numeric(dat.4$month) - 9


dat.2$steady = as.factor(ifelse(dat.2$year %in% st.yr, 1, 0))
dat.2$mo.num = as.numeric(dat.2$month) - 9



fit.1 = glmmadmb(count.sum ~ 0 + (1|year) + steady:mo.num + (1|RM) + offset(log(Volume)), 
                 data = dat.2, family = 'nbinom')












