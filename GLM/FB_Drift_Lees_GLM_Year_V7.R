################################################################################
#                                                                       11/12/14
#
#   Fit count model (neg.binom) to the drift data to look at yearly effects 
#     
#
#                                                                     
################################################################################
#setwd('C:/Users/mdodrill/Desktop/FOODBASE/DRIFT_LEE_GLM')
setwd('P:/BIOLOGICAL/Flyco/DRIFT_LEE_GLM/DATA')
rm(list=ls(all=TRUE))
# library(MASS)   # for glm.nb
# library(bbmle)  # AICtab - delta AIC table 
library(dplyr)
library(ggplot2)
library(xts)     # align.time (round to 15min) library(zoo)
library(lme4)
library(R2admb)
library(glmmADMB)


# load functions stored in R.script
source("P:/BIOLOGICAL/Flyco/DRIFT_COMMON/FB_Drift_Functions_Common_V1.R", chdir = F)

lf.q = read.csv('LFQ_2007_2016.csv', header = T)  

#-----------------------------------------------------------------------------#
# see look -- these are the species to do summaries for....
# Right now the functions are only 'wired' for LUM, NZM, GAM, SIM, CHI
sp_num = unique(look[,1])      

# Data to import (these are global and not passed in)
drift_specimen = 'DRIFT_SPECIMEN_01_18_17.csv'
drift_sample   = 'DRIFT_SAMPLE_01_18_17.csv'


# This is the main function that takes a year (or years) and species of interest as 
# arguments. Creates a summary (acconting for the structure of the data, what you get
# from the functions 'summarize' & 'derive') and associates the sampling info (stuff
# from the DRIFT_SAMPLE table in Access).

data = tbl_df(drift_summary(c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016), sp_num)   )

# 2007 has a different set of gear !
data[which(substr(data$DriftDate,1,4) == "2007"),]$GearID
#-----------------------------------------------------------------------------#
# filter out any data that you don't want 
# (make sure this matches "FB_Drift_to_BioE_Control_1mm")
#data2 = filter(data, purpose =="monthly"| purpose =="flood")
data3 = filter(data,DriftRM<0.1)
#data2 = filter(data, GearID == 4 | GearID == 3 |GearID ==2 |GearID ==1, DriftRM == 0.0 | DriftRM == -0.2)
#data3 = filter(data2, DepthIntegrated == TRUE)


tmp2 = data3[which(substr(data3$DriftDate,1,10) == "2008-03-03"),]

# for 2008 -- not March or April -- flood
#tmp = data3[which(substr(data3$DriftDate,1,7) != "2008-03" & substr(data3$DriftDate,1,7) != "2012-11" 
#            & substr(data3$DriftDate,1,7) != "2013-11"),]


data4 = data3
#data4 = rbind(data3, tmp2)
#unique(data4$DriftDate)

# make jan 2008 --> 2007 
idx = which(substr(data4$DriftDate,1,7) == "2008-01") 

idx2 = which(substr(data4$DriftDate,1,10) == "2008-03-03")

data4$year = as.numeric(substr(data4$DriftDate,1,4))

data4[idx,which(names(data4) == "year")] = "2007"
data4[idx2,which(names(data4) == "year")] = "2007"

data3 = data4

#------------------------------------------------------------------------------#
# Group the life stages of midges and blackflies together

tmp.ID = ifelse(data3$SpeciesID == 19, 21, data3$SpeciesID)
tmp.ID = ifelse(tmp.ID == 20, 21, tmp.ID)

tmp.ID = ifelse(tmp.ID == 5, 7, tmp.ID)
tmp.ID = ifelse(tmp.ID == 6, 7, tmp.ID)
unique(tmp.ID)

data3$SpeciesID = tmp.ID

unique(substr(data3$DriftDate,1,4))

#------------------------------------------------------------------------------#
#Bug meat. Skip below to do GLM models on counts
data3 = filter(data4,Dtot>=0) #gets rid of barcode D02986 with negative vol filtered
data3 = filter(data3,SpeciesID!=2&SpeciesID!=1) #removes NZMS and worms

out = dplyr::select(data3, DriftSampleID, Dtot) %>%
  group_by(DriftSampleID) %>%
  summarise(tot.mass = sum(Dtot)) #sums all the Dtots for a Sample 


tmp = data3[match(out$DriftSampleID, data3$DriftSampleID),26:42] #grabs the other info for the sample

out2 = cbind(out,tmp) 

out2$flood = ifelse(substr(out2$DriftDate,1,7) == "2008-03" | substr(out2$DriftDate,1,7) == "2012-11" 
                            | substr(out2$DriftDate,1,7) == "2013-11" ,1,0)

drift<-filter(out2,flood==0)
#drift<-filter(counts.data4,purpose=="monthly" | purpose=="flood")
drift<-filter(drift,DriftRM==-11.5 |DriftRM==-11|DriftRM==-8|DriftRM==-4.9|DriftRM==-3.5|DriftRM==-3.4|DriftRM==-2.1|DriftRM==-0.2
              |DriftRM==0)

time = as.POSIXct(strptime(as.character(drift$TimeBegin),"%H:%M"))  
tmp.time = substr(align.time(time, n = 60 * 15), 11, 19) # n is in seconds (rounds up)
drift$date.time = as.POSIXct(paste(drift$DriftDate, tmp.time), format = "%Y-%m-%d %H:%M:%S")


# format the discharge time column to match the bug data
lf.q$time = as.POSIXct(strptime(as.character(lf.q$time),"%m/%d/%Y  %H:%M"))
lf.q$year<-substr(lf.q$time,1,4)
# add a new col. with the coresponding discharges.
drift$q = lf.q[match(drift$date.time, lf.q$time),2]
#drift$q = drift$q * .028317 converts to m3/s

counts.data3 = counts.data2[which(!is.na(counts.data2$q)),]
counts.data4 = filter(counts.data3, Volume_Sampled_M.3>10 & Velocity_m.s>0.1)



drift = filter(drift, Volume_Sampled_M.3>10 & Velocity_m.s>0.1)
drift$month = as.factor(substr(drift$date.time, 6, 7))
drift$yr.mo = as.factor(paste(drift$year, drift$month))
drift$year = as.factor(drift$year) 
drift$jul = (as.numeric(drift$DriftDate - drift$DriftDate[1])/86400)
drift$tlq<- log(drift$q)-mean(log(drift$q))
drift$RM = as.factor(drift$DriftRM)
drift$mo.num <-as.numeric(substr(drift$date.time, 6, 7))
drift2 <- filter(drift, mo.num>3 & mo.num<10)#|month==5|month==6|month==7|month==8|month==9)

################################################
##########models of bug meat

#models = NULL  # variable to save model fits 
#sp. = unique(drift$SpeciesID)
#big.out = NULL

#for(i in 1:length(sp.)){
  
#  sp.sub = filter(drift, SpeciesID == sp.[i])
  
   bugmeat = lmer(log(tot.mass) ~ 0 + year + (1|RM) + (1|month), data = drift)
   bugmeat2 = lmer(log(tot.mass) ~ 0 + year + (1|RM) + (1|month), data = drift2)
  tmp.len = length(coef(fm_1))
  
  tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
  tmp.out[,1] = rep(sp.[i], tmp.len)
  tmp.out[,2] = as.numeric(paste(unique(drift$year)))  # this should match the order of exp.count (coefs in model)
  tmp.out[,3] = exp(coef(fm_1)[1:10])*mean(exp(ranef(fm_1)$month)) 
  tmp.out[,4] = exp(confint(fm_1)[1:10,1])*mean(exp(ranef(fm_1)$month))
  tmp.out[,5] = exp(confint(fm_1)[1:10,2])*mean(exp(ranef(fm_1)$month)) 
  
  # coef(fm_1)
  # confint(fm_1)
  # ranef(fm_1)
  
  # Model predictions 
  # newdata = expand.grid(year = unique(sp.sub$year), month = unique(sp.sub$month), vol = 1)
  # newdata$yhatmean <- predict(fm_1, newdata = newdata, type = "response", se = T)   # predicted mean rate (#/m^3) 
  
  
  models[[i]] = fm_1
  
  big.out[[i]] = tmp.out
}

model.counts = do.call(rbind, big.out)
colnames(model.counts) = list("sp", "year", "est", "ll", "ul")

###########################
# make some YEAR plots 
windows()
model.counts$sp2 = look[match(model.counts[,1], look[,1]),2]


p = ggplot(model.counts, aes(x = factor(year), y = est)) +
  geom_point() +  # aes(color = year)
  geom_errorbar(aes(ymax = ul, ymin = ll)) + #, color = year  
  facet_wrap(~ sp, scales = "free_y") + 
  labs(title = "Long-term drift monitoring",
       y = "Count / m3",x ="Year")  # (1|yr.mo) or (1|month)
p

p2 = p + theme_bw()
p3 = p2 +theme(panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"))
p3











#------------------------------------------------------------------------------#
# to counts
counts.data = data_to_counts(data3)
# add zeros
counts.data2 = tbl_df(add_zeros(counts.data))



#-----------------------------------------------------------------------------#
# adding in the discharge data 


# round the drift sampling time to 15 mins 
time = as.POSIXct(strptime(as.character(counts.data2$TimeBegin),"%H:%M"))  
tmp.time = substr(align.time(time, n = 60 * 15), 11, 19) # n is in seconds (rounds up)
counts.data2$date.time = as.POSIXct(paste(counts.data2$DriftDate, tmp.time), format = "%Y-%m-%d %H:%M:%S")


# format the discharge time column to match the bug data
lf.q$time = as.POSIXct(strptime(as.character(lf.q$time),"%m/%d/%Y  %H:%M"))
lf.q$year<-substr(lf.q$time,1,4)
# add a new col. with the coresponding discharges.
counts.data2$q = lf.q[match(counts.data2$date.time, lf.q$time),2]
counts.data2$q = counts.data2$q * .028317

counts.data3 = counts.data2[which(!is.na(counts.data2$q)),]
counts.data4 = filter(counts.data3, Volume_Sampled_M.3>10 & Velocity_m.s>0.1)

# above excludes DriftID 1303,1816,2363,2374,2563, 2942, 4552, 1303, 1816, 2363, 2374, 2563, 2942,4552
####### ALSO NEED TO CHECK DriftID 1938 and 4693, both of which have velocity <0.1 m/s.  
########  Did not exclude these bc Volume Sampled is reasonable.  Suggest problem with duration of sampling. 
#-----------------------------------------------------------------------------#
# set up models for each species 


counts.data4$month = as.factor(substr(counts.data4$date.time, 6, 7))
counts.data4$yr.mo = as.factor(paste(counts.data4$year, counts.data4$month))
counts.data4$year = as.factor(counts.data4$year) 
counts.data4$jul = (as.numeric(counts.data4$DriftDate - counts.data4$DriftDate[1])/86400)
counts.data4$tlq<- log(counts.data4$q)-mean(log(counts.data4$q))
#counts.data4$bin3_6<-rowSums(counts.data4[,7:10])
#counts.data4$Station = as.factor(counts.data4$Station)
counts.data4$RM = as.factor(counts.data4$DriftRM)
################################################

counts.data4$flood = ifelse(substr(counts.data4$DriftDate,1,7) == "2008-03" | substr(counts.data4$DriftDate,1,7) == "2012-11" 
                  | substr(counts.data4$DriftDate,1,7) == "2013-11" ,1,0)

drift<-filter(counts.data4,flood==0)
#drift<-filter(counts.data4,purpose=="monthly" | purpose=="flood")
drift<-filter(drift,DriftRM==-11.5 |DriftRM==-11|DriftRM==-8|DriftRM==-4.9|DriftRM==-3.5|DriftRM==-3.4|DriftRM==-2.1|DriftRM==-0.2
              |DriftRM==0)

#drift1<-filter(drift,purpose!="flood")
#drift2<-filter(drift1,Station==3 & (DriftRM==0| DriftRM==-0.2))

##############
models = NULL  # variable to save model fits 
sp. = unique(drift$SpeciesID)
big.out = NULL
drift1$DriftRM<-as.factor(drift$DriftRM)

for(i in 1:length(sp.)){

  sp.sub = filter(drift, SpeciesID == sp.[i])

  fm_1 = glmmadmb(count.sum ~ 0 + year + (1|RM) + (1|month) + offset(log(Volume_Sampled_M.3)), 
                  data = sp.sub, family = 'nbinom')
  
   tmp.len = length(coef(fm_1))
  
   tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
   tmp.out[,1] = rep(sp.[i], tmp.len)
   tmp.out[,2] = as.numeric(paste(unique(drift$year)))  # this should match the order of exp.count (coefs in model)
   tmp.out[,3] = exp(coef(fm_1)[1:10])*mean(exp(ranef(fm_1)$month)) 
   tmp.out[,4] = exp(confint(fm_1)[1:10,1])*mean(exp(ranef(fm_1)$month))
   tmp.out[,5] = exp(confint(fm_1)[1:10,2])*mean(exp(ranef(fm_1)$month)) 
  
  # coef(fm_1)
  # confint(fm_1)
  # ranef(fm_1)
  
  # Model predictions 
  # newdata = expand.grid(year = unique(sp.sub$year), month = unique(sp.sub$month), vol = 1)
  # newdata$yhatmean <- predict(fm_1, newdata = newdata, type = "response", se = T)   # predicted mean rate (#/m^3) 
  
  
  models[[i]] = fm_1
  
  big.out[[i]] = tmp.out
}

model.counts = do.call(rbind, big.out)
colnames(model.counts) = list("sp", "year", "est", "ll", "ul")

###########################
# make some YEAR plots 
windows()
model.counts$sp2 = look[match(model.counts[,1], look[,1]),2]


p = ggplot(model.counts, aes(x = factor(year), y = est)) +
  geom_point() +  # aes(color = year)
  geom_errorbar(aes(ymax = ul, ymin = ll)) + #, color = year  
  facet_wrap(~ sp, scales = "free_y") + 
  labs(title = "Long-term drift monitoring",
       y = "Count / m3",x ="Year")  # (1|yr.mo) or (1|month)
p

p2 = p + theme_bw()
p3 = p2 +theme(panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"))
p3



# need to source Facet_Wrap_Labels_Function 
# source ("C:/Users/mdodrill/Desktop/FOODBASE/DRIFT_LEE_GLM/Facet_Wrap_Labels_Function.R", chdir = F)
facetAdjust(p3)
#########################################################################
gpp<-read.csv("glen_canyon_metab_export_may14.csv",header=T)
gpp$Time<-as.POSIXct(gpp$Time,format="%m/%d/%Y")
gpp$Year<-substr(gpp$Time,1,4)
gpp$Month<-substr(gpp$Time,6,7)
gppsum<-aggregate(gpp$gpp,by=list(gpp$Month,gpp$Year),FUN=mean)




#############################################################Month models 
models.mo = NULL  # variable to save model fits 
sp. = unique(drift$SpeciesID)
big.out.mo = NULL



for(i in 1:length(sp.)){
  
  sp.sub = filter(drift, SpeciesID == sp.[i])
  fm_1 = glmmadmb(count.sum ~ 0 + month + (1|year) + offset(log(Volume_Sampled_M.3)),  
                  data = sp.sub, family = 'nbinom')
  
  
  tmp.len = length(coef(fm_1))
  
  tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
  tmp.out[,1] = rep(sp.[i], tmp.len)
  tmp.out[,2] = sort(as.numeric(paste(unique(drift$month))))  # this should match the order of exp.count (coefs in model)
  tmp.out[,3] = exp(coef(fm_1))*mean(exp(ranef(fm_1)$year)) 
  tmp.out[,4] = exp(confint(fm_1)[,1])*mean(exp(ranef(fm_1)$year))
  tmp.out[,5] = exp(confint(fm_1)[,2])*mean(exp(ranef(fm_1)$year)) 
  
  
  models.mo[[i]] = fm_1
  
  big.out.mo[[i]] = tmp.out
}

model.counts.mo = do.call(rbind, big.out.mo)
colnames(model.counts.mo) = list("sp", "month", "est", "ll", "ul")
#model.counts2 = data.frame(matrix(NA,dim(model.counts)[1],dim(model.counts)[2]))
###########################
# make some month plots 

windows()
model.counts.mo$sp2 = look[match(model.counts.mo[,1], look[,1]),2]




p = ggplot(model.counts.mo, aes(x = factor(month), y = est)) +
  geom_point() +  # aes(color = year)
  geom_errorbar(aes(ymax = ul, ymin = ll)) + #, color = year  
  facet_wrap(~ sp, scales = "free_y") + 
  labs(title = "Long-term drift monitoring",
       y = "Count / m3",x ="Month" )  # (1|yr.mo) or (1|month)
p

p2 = p + theme_bw()
p3 = p2 +theme(panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"))
p3


#################
models.q = NULL  # variable to save model fits 
sp. = unique(counts.data4$SpeciesID)
big.out.q = NULL



for(i in 1:length(sp.)){
  
  sp.sub = filter(counts.data4, SpeciesID == sp.[i])
  fm_1 = glmmadmb(count.sum ~ 0 + tlq + month + year + (1|yr.mo) + offset(log(Volume_Sampled_M.3)),  
                  data = sp.sub, family = 'nbinom')
  
  tmp.len = length(coef(fm_1))+1
  
  tmp.out = data.frame(matrix(NA, tmp.len, length(unique(counts.data4$SpeciesID))))
  tmp.out[,1] = rep(sp.[i], tmp.len)
  tmp.out[,2] = c('tlq',seq(1,12,1),seq(2007,2014,1))  # this should match the order of exp.count (coefs in model)
  tmp.out[1,3] = coef(fm_1)[1]
  tmp.out[2:13,3] = exp(coef(fm_1))[2:13]*mean(exp(ranef(fm_1)$yr.mo))
  tmp.out[14,3] = 0
  tmp.out[15:21,3] = exp(coef(fm_1))[14:20]*mean(exp(ranef(fm_1)$yr.mo))
  tmp.out[,4] = exp(confint(fm_1)[,1])*mean(exp(ranef(fm_1)$year))
  tmp.out[,5] = exp(confint(fm_1)[,2])*mean(exp(ranef(fm_1)$year)) 
  
  
  models.mo[[i]] = fm_1
  
  big.out.mo[[i]] = tmp.out
}

model.counts.mo = do.call(rbind, big.out.mo)
colnames(model.counts.mo) = list("sp", "month", "est", "ll", "ul")
#model.counts2 = data.frame(matrix(NA,dim(model.counts)[1],dim(model.counts)[2]))







p = ggplot(drift, aes(x = sp., y = count.sum)) +
  geom_boxplot() +  # aes(color = year)
  geom_errorbar(aes(ymax = ul, ymin = ll)) + #, color = year  
  facet_wrap(~ sp, scales = "free_y") + 
  labs(title = "Long-term drift monitoring",
       y = "Count / m3",x ="Month" )  # (1|yr.mo) or (1|month)
p





###########nzms by size bin
i<-2
sp. = unique(counts.data4$SpeciesID)
sp.sub = filter(counts.data4, SpeciesID == sp.[i])

fm_1 = glmmadmb(sp.sub[,4] ~ 0 + sp.sub$year + sp.sub$tlq + (1|sp.sub$month) + (1|sp.sub$yr.mo) + offset(log(sp.sub$Volume_Sampled_M.3)),  #  (1|yr.mo)  or (1|month)
                family = 'nbinom')

fm_2 = glmmadmb(sp.sub[,5] ~ 0 + sp.sub$year + sp.sub$tlq + (1|sp.sub$month) + (1|sp.sub$yr.mo) + offset(log(sp.sub$Volume_Sampled_M.3)),  #  (1|yr.mo)  or (1|month)
                family = 'nbinom')

fm_3 = glmmadmb(sp.sub[,6] ~ 0 + sp.sub$year + sp.sub$tlq + (1|sp.sub$month) + (1|sp.sub$yr.mo) + offset(log(sp.sub$Volume_Sampled_M.3)),  #  (1|yr.mo)  or (1|month)
                family = 'nbinom')

fm_4 = glmmadmb(sp.sub$bin3_6 ~ 0 + sp.sub$year + sp.sub$tlq + (1|sp.sub$month) + (1|sp.sub$yr.mo) + offset(log(sp.sub$Volume_Sampled_M.3)),  #  (1|yr.mo)  or (1|month)
                family = 'nbinom')




