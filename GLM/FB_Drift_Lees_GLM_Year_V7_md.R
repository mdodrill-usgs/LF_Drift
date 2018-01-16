################################################################################
#                                                                     
#
#   Fit count model (neg.binom) to the drift data to look at yearly effects 
#     
#
#                                                                     
################################################################################
setwd('C:/Users/mdodrill/Desktop/FOODBASE/LF_Drift/Data/')
# setwd('C:/Users/tkennedy/Documents/Repositories/LF_Drift/Data')
rm(list=ls(all=TRUE))
# library(MASS)   # for glm.nb
# library(bbmle)  # AICtab - delta AIC table 
library(dplyr)
library(ggplot2)
library(xts)     # align.time (round to 15min) library(zoo)
library(lme4)
library(R2admb)
# install.packages("glmmADMB",repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")),type="source")
library(glmmADMB)




# load functions stored in R.script
# source("P:/BIOLOGICAL/Flyco/DRIFT_COMMON/FB_Drift_Functions_Common_V1.R", chdir = F)
# source("C:/Users/tkennedy/Documents/Repositories/Drift/FB_Drift_Functions_Git_V1.R", chdir = F)
source("C:/Users/mdodrill/Desktop/FB_Git/Drift/FB_Drift_Functions_Git_V2.R", chdir = F)

lf.q = read.csv('LFQ_2007_2016.csv', header = T)  

#-----------------------------------------------------------------------------#
# see look -- these are the species to do summaries for....
# Right now the functions are only 'wired' for LUM, NZM, GAM, SIM, CHI
sp_num = unique(look[,1])      

# Data to import (these are global and not passed in)
# drift_specimen = 'DRIFT_SPECIMEN_05_30_2017.csv'
# drift_sample   = 'DRIFT_SAMPLE_05_30_2017.csv'
drift_specimen = 'DriftSpecimen.csv'
drift_sample   = 'DriftSample.csv'



# This is the main function that takes a year (or years) and species of interest as 
# arguments. Creates a summary (acconting for the structure of the data, what you get
# from the functions 'summarize' & 'derive') and associates the sampling info (stuff
# from the DRIFT_SAMPLE table in Access).

# data = drift_summary(c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016), sp_num, update.data = TRUE)
data = drift_summary(c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017), sp_num)

# 2007 has a different set of gear !
data[which(substr(data$DriftDate,1,4) == "2007"),]$GearID
#-----------------------------------------------------------------------------#
# filter out any data that you don't want 
# (make sure this matches "FB_Drift_to_BioE_Control_1mm")
#data2 = filter(data, purpose =="monthly"| purpose =="flood")
data3 = filter(data,DriftRM<0.1)
#data2 = filter(data, GearID == 4 | GearID == 3 |GearID ==2 |GearID ==1, DriftRM == 0.0 | DriftRM == -0.2)
#data3 = filter(data2, DepthIntegrated == TRUE)


# tmp2 = data3[which(substr(data3$DriftDate,1,10) == "2008-03-03"),]

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
#make some graphs
# lf.q$m3<-lf.q[,2]*0.028317
# windows(width = 10, height = 10, record = TRUE)
# plot(lf.q[,1],lf.q[,4], ylim = c(0,1300),type='l',xlab="DATE",ylab="Discharge (m3/s)",cex.lab=1.3,cex.axis=1.3)
# drift.dates<-unique(counts.data4$DriftDate)
# points(drift.dates,rep(100,178),cex=1.5)




#-----------------------------------------------------------------------------#
# set up models for each species 


counts.data4$month = as.factor(substr(counts.data4$date.time, 6, 7))
counts.data4$yr.mo = as.factor(paste(counts.data4$year, counts.data4$month))
counts.data4$year = as.factor(counts.data4$year) 
counts.data4$jul = (as.numeric(counts.data4$DriftDate - counts.data4$DriftDate[1])/86400)
#counts.data4$tlq<- log(counts.data4$q)-mean(log(counts.data4$q))
#counts.data4$bin3_6<-rowSums(counts.data4[,7:10])
#counts.data4$Station = as.factor(counts.data4$Station)
counts.data4$DriftRM = as.factor(counts.data4$DriftRM)
################################################

counts.data4$flood = ifelse(substr(counts.data4$DriftDate,1,7) == "2008-03" | substr(counts.data4$DriftDate,1,7) == "2012-11" 
                            | substr(counts.data4$DriftDate,1,7) == "2013-11" ,1,0)

#drift<-filter(counts.data4,purpose=="monthly" | purpose=="flood")
counts.data4<-filter(counts.data4, DriftRM==-11.5 |DriftRM==-11|DriftRM==-8|DriftRM==-4.9|DriftRM==-3.5|DriftRM==-3.4|DriftRM==-2.1|DriftRM==-0.2
                     |DriftRM==0)
flood.drift<-filter(counts.data4,SpeciesID!=1)
drift<-filter(counts.data4,flood==0)
drift<-filter(drift,SpeciesID!=1)
##############
models = NULL  # variable to save model fits 
sp. = unique(drift$SpeciesID)
big.out = NULL


for(i in 1:length(sp.)){
  
  sp.sub = filter(drift, SpeciesID == sp.[i])
  
  fm_1 = glmmadmb(count.sum ~ 0 + year + (1|month)  + (1|DriftRM) + offset(log(Volume_Sampled_M.3)), 
                  data = sp.sub, family = 'nbinom')
  
  tmp.len = length(coef(fm_1))
  
  tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
  tmp.out[,1] = rep(sp.[i], tmp.len)
  
  tmp.out[,2] = sort(as.numeric(paste(unique(drift$year))))
  
  # tmp.out[,3] = exp(coef(fm_1)[1:10])*mean(exp(ranef(fm_1)$month)) 
  #   tmp.out[,3] = exp(coef(fm_1)[1:10]) * exp(ranef(fm_1)$month[5]) 
  #   tmp.out[,4] = exp(confint(fm_1)[1:10,1])*exp(ranef(fm_1)$month[5])
  #   tmp.out[,5] = exp(confint(fm_1)[1:10,2])*exp(ranef(fm_1)$month[5]) 
  tmp.out[,3] = exp(coef(fm_1)[1:tmp.len]) #* exp(ranef(fm_1)$month[5]) 
  tmp.out[,4] = exp(confint(fm_1)[1:tmp.len,1])#*exp(ranef(fm_1)$month[5])
  tmp.out[,5] = exp(confint(fm_1)[1:tmp.len,2])#*exp(ranef(fm_1)$month[5]) 
  
  # coef(fm_1)
  # confint(fm_1)
  # ranef(fm_1)
  
  # Model predictions 
  # newdata = expand.grid(year = unique(sp.sub$year), month = unique(sp.sub$month), Volume_Sampled_M.3 = 1)
  # newdata = expand.grid(year = unique(sp.sub$year), month = as.factor(5), Volume_Sampled_M.3 = 1)
  # newdata$yhatmean <- predict(fm_1, newdata = newdata, type = "response")#, se = T)   # predicted mean rate (#/m^3) 
  # test <- predict(fm_1, newdata = newdata, type = "response", se = T)   # predicted mean rate (#/m^3) 
  
  
  models[[i]] = fm_1
  
  big.out[[i]] = tmp.out
}

model.counts = do.call(rbind, big.out)
colnames(model.counts) = list("sp", "year", "est", "ll", "ul")

###########################
# make some YEAR plots 
windows(width = 10, height = 10, record = TRUE)
model.counts$sp2 = look[match(model.counts[,1], look[,1]),2]

model.counts$year2 = as.character(model.counts$year)
model.counts$year2 = as.numeric(model.counts$year)

model.counts$year2 = seq(2007, 2017,1)


# p = ggplot(model.counts, aes(x = factor(year), y = est)) +
p = ggplot(model.counts, aes(x = year2, y = est)) +
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



#-----------------------------------------------------------------------------#
st.yr = c(2008,2009,2010,2011,2012)

dat.1 = drift[which(drift$SpeciesID == 21),]

dat.2 = dat.1[which(dat.1$month %in% c("08", "09", "10")),]

dat.3 = dat.2[which(dat.2$year %in% st.yr),]

dat.4 = dat.2[which(!dat.2$year %in% st.yr),]


dat.3$mo.num = as.numeric(dat.3$month) - 9
dat.4$mo.num = as.numeric(dat.4$month) - 9


dat.2$steady = as.factor(ifelse(dat.2$year %in% st.yr, 1, 0))
dat.2$mo.num = as.numeric(dat.2$month) - 9





fit.1 = glmmadmb(count.sum ~ 0 + (1|year) + steady:mo.num + (1|DriftRM) + offset(log(Volume_Sampled_M.3)), 
                 data = dat.2, family = 'nbinom')

ranef(fit.1)




# fit.1 = glmmadmb(count.sum ~ 0 + (1|year) + month  + (1|DriftRM) + offset(log(Volume_Sampled_M.3)),
#                 data = dat.3, family = 'nbinom')
# 
# fit.2 = glmmadmb(count.sum ~ 0 + (1|year) + month  + (1|DriftRM) + offset(log(Volume_Sampled_M.3)),
#                 data = dat.4, family = 'nbinom')


tmp.len = length(coef(fit.1))

tmp.out = data.frame(matrix(NA, tmp.len, ncol = 4))
# tmp.out[,1] = rep(sp.[i], tmp.len)

tmp.out[,1] = rep("steady", tmp.len)

tmp.out[,2] = c("08", "09", "10")

tmp.out[,3] = exp(coef(fit.1)[1:tmp.len]) 
tmp.out[,4] = exp(confint(fit.1)[1:tmp.len,1])
tmp.out[,5] = exp(confint(fit.1)[1:tmp.len,2])

tmp.1 = tmp.out

tmp.out = data.frame(matrix(NA, tmp.len, ncol = 4))
# tmp.out[,1] = rep(sp.[i], tmp.len)

tmp.out[,1] = rep("fluc", tmp.len)

tmp.out[,2] = c("08", "09", "10")

tmp.out[,3] = exp(coef(fit.2)[1:tmp.len]) 
tmp.out[,4] = exp(confint(fit.2)[1:tmp.len,1])
tmp.out[,5] = exp(confint(fit.2)[1:tmp.len,2])


all = rbind(tmp.1, tmp.out)

colnames(all) = list("flow", "month", "est", "ll", "ul")



p = ggplot(all, aes(x = month, y = est)) +
  geom_point(aes(color = month)) +  # aes(color = year)
  geom_errorbar(aes(ymax = ul, ymin = ll, color = month)) + #, color = year  
  facet_wrap(~ flow) +
  labs(title = "Long-term drift monitoring",
       y = "Count / m3",x ="month")  # (1|yr.mo) or (1|month)
p

p2 = p + theme_bw()
p3 = p2 +theme(panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"))
p3

windows()

#-----------------------------------------------------------------------------#

# # need to source Facet_Wrap_Labels_Function 
# # source ("C:/Users/mdodrill/Desktop/FOODBASE/DRIFT_LEE_GLM/Facet_Wrap_Labels_Function.R", chdir = F)
# facetAdjust(p3)
# 
# #############################################################Month models 
# models.mo = NULL  # variable to save model fits 
# sp. = unique(drift$SpeciesID)
# big.out.mo = NULL
# 
# 
# 
# for(i in 1:length(sp.)){
#   
#   sp.sub = filter(drift, SpeciesID == sp.[i])
#   fm_1 = glmmadmb(count.sum ~ 0 + month + (1|year) + offset(log(Volume_Sampled_M.3)),  
#                   data = sp.sub, family = 'nbinom')
#   
#   
#   tmp.len = length(coef(fm_1))
#   
#   tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
#   tmp.out[,1] = rep(sp.[i], tmp.len)
#   tmp.out[,2] = sort(as.numeric(paste(unique(drift$month))))  # this should match the order of exp.count (coefs in model)
#   tmp.out[,3] = exp(coef(fm_1))#*mean(exp(ranef(fm_1)$year)) 
#   tmp.out[,4] = exp(confint(fm_1)[,1])#*mean(exp(ranef(fm_1)$year))
#   tmp.out[,5] = exp(confint(fm_1)[,2])#*mean(exp(ranef(fm_1)$year)) 
#   
#   
#   models.mo[[i]] = fm_1
#   
#   big.out.mo[[i]] = tmp.out
# }
# 
# model.counts.mo = do.call(rbind, big.out.mo)
# colnames(model.counts.mo) = list("sp", "month", "est", "ll", "ul")
# #model.counts2 = data.frame(matrix(NA,dim(model.counts)[1],dim(model.counts)[2]))
# ###########################
# # make some month plots 
# 
# windows(width = 10, height = 10, record = TRUE)
# model.counts.mo$sp2 = look[match(model.counts.mo[,1], look[,1]),2]
# 
# 
# 
# 
# p = ggplot(model.counts.mo, aes(x = factor(month), y = est)) +
#   geom_point() +  # aes(color = year)
#   geom_errorbar(aes(ymax = ul, ymin = ll)) + #, color = year  
#   facet_wrap(~ sp, scales = "free_y") + 
#   labs(title = "Long-term drift monitoring",
#        y = "Count / m3",x ="Month" )  # (1|yr.mo) or (1|month)
# p
# 
# p2 = p + theme_bw()
# p3 = p2 +theme(panel.grid.major = element_line(colour = "white"),
#                panel.grid.minor = element_line(colour = "white"))
# p3
# 
# 
# #################
# 
