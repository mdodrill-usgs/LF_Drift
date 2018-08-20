################################################################################
#                                                                       Aug 2018
#
#   Fit count model (neg.binom) to the drift data to look at yearly and monthly
#   effects
#     
#  Notes:
#  * 
#
#  To Do:
#  * 
#                                                                     
################################################################################
# setwd('C:/Users/mdodrill/Desktop/FOODBASE/LF_Drift/Data/')
setwd('C:\\Users\\mdodrill\\Desktop\\FB_Git\\LF_Drift\\DATA')
# setwd('C:/Users/tkennedy/Documents/Repositories/LF_Drift/Data')
rm(list=ls(all=TRUE))
# library(MASS)   # for glm.nb
# library(bbmle)  # AICtab - delta AIC table 
library(foodbase)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(xts)     # align.time (round to 15min) library(zoo)
# library(lme4)
library(glmmTMB)

#-----------------------------------------------------------------------------#

# all.dat = readDB(updater = TRUE)
all.dat = readDB()

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

taxa = c("LUMB", "NZMS", "GAMM", "SIMA", "SIMP", "SIML", "CHIA", "CHIP", "CHIL")

ltl.spec = spec[which(spec$SpeciesID %in% taxa),]

#--------------------------------------
# merge the sample table with the specimen counts

tmp.dat = left_join(ltl.spec, ltl.samp, by = "BarcodeID")

lf.dat = tmp.dat[,c(1,24:32,2:23)]

lf.dat$SpeciesID = as.character(lf.dat$SpeciesID)

# rm(list = setdiff(ls(), "lf.dat"))

#-----------------------------------------------------------------------------#
data4 = lf.dat

# make jan 2008 --> 2007 
idx = which(substr(data4$Date,1,7) == "2008-01") 

idx2 = which(substr(data4$Date,1,10) == "2008-03-03")

data4$year = as.numeric(substr(data4$Date,1,4))

data4[idx,which(names(data4) == "year")] = "2007"
data4[idx2,which(names(data4) == "year")] = "2007"

data3 = data4

#------------------------------------------------------------------------------#
# Group the life stages of midges and blackflies together

tmp.ID = ifelse(data3$SpeciesID == "CHIA", "CHIL", data3$SpeciesID)
tmp.ID = ifelse(tmp.ID == "CHIP", "CHIL", tmp.ID)

tmp.ID = ifelse(tmp.ID == "SIMA", "SIML", tmp.ID)
tmp.ID = ifelse(tmp.ID == "SIMP", "SIML", tmp.ID)
unique(tmp.ID)

data3$SpeciesID = tmp.ID

unique(substr(data3$Date,1,4))

#------------------------------------------------------------------------------#
# to counts

# the new data structure has the zeros...
counts.data2 = data3
#-----------------------------------------------------------------------------#
# adding in the discharge data 
lf.q = read.csv('LFQ_2007_2016.csv', header = T)  

# round the drift sampling time to 15 mins 
time = as.POSIXct(strptime(as.character(counts.data2$TimeBegin),"%H:%M"))  
tmp.time = substr(align.time(time, n = 60 * 15), 11, 19) # n is in seconds (rounds up)
counts.data2$date.time = as.POSIXct(paste(counts.data2$Date, tmp.time), format = "%Y-%m-%d %H:%M:%S")

# format the discharge time column to match the bug data
lf.q$time = as.POSIXct(strptime(as.character(lf.q$time),"%m/%d/%Y  %H:%M"))
lf.q$year<-substr(lf.q$time,1,4)

# add a new col. with the coresponding discharges.
counts.data2$q = lf.q[match(counts.data2$date.time, lf.q$time),2]
counts.data2$q = counts.data2$q * .028317

counts.data3 = counts.data2[which(!is.na(counts.data2$q)),]
# counts.data4 = filter(counts.data3, Volume_Sampled_M.3>10 & Velocity_m.s>0.1)

# maybe do with dat[[1]]$FlagStrange ? 
# don't have the updated Q, so dont drop these for now....!
counts.data4 = counts.data2


#-----------------------------------------------------------------------------#
# set up models for each species 

counts.data4$month = as.factor(substr(counts.data4$date.time, 6, 7))
counts.data4$yr.mo = as.factor(paste(counts.data4$year, counts.data4$month))
counts.data4$year = as.factor(counts.data4$year) 
counts.data4$jul = (as.numeric(counts.data4$Date - counts.data4$Date[1])/86400)
#counts.data4$tlq<- log(counts.data4$q)-mean(log(counts.data4$q))
#counts.data4$bin3_6<-rowSums(counts.data4[,7:10])
#counts.data4$Station = as.factor(counts.data4$Station)
counts.data4$RiverMile = as.factor(counts.data4$RiverMile)

#-----------------------------------------------------------------------------#
counts.data4$flood = ifelse(substr(counts.data4$Date,1,7) == "2008-03" | substr(counts.data4$Date,1,7) == "2012-11"
                            | substr(counts.data4$Date,1,7) == "2013-11" ,1,0)

#drift<-filter(counts.data4,purpose=="monthly" | purpose=="flood")
counts.data4<-filter(counts.data4, RiverMile==-11.5 |RiverMile==-11|RiverMile==-8|
                       RiverMile==-4.9|RiverMile==-3.5|RiverMile==-3.4|RiverMile==-2.1|
                       RiverMile==-0.2|RiverMile==0)

# flood.drift<-filter(counts.data4, SpeciesID!=1)
drift<-filter(counts.data4,flood==0)
# drift<-filter(drift,SpeciesID!=1)

drift$count.sum = rowSums(drift[,12:32])

#-----------------------------------------------------------------------------#
# run some models for each taxa looking at year effects
models = NULL  # variable to save model fits 
sp. = unique(drift$SpeciesID)
big.out = NULL


for(i in 1:length(sp.)){
  
  sp.sub = filter(drift, SpeciesID == sp.[i])
  
  start.time <- Sys.time()  # start timer
  fm_1 <- glmmTMB(count.sum ~ 0 + year + (1|month) + (1|RiverMile) + offset(log(Volume)),
  # fm_1 <- glmmTMB(count.sum ~ 0 + year + (1|month)  + offset(log(Volume)),
                 data = sp.sub,
                 family = nbinom2)
  
  end.time = Sys.time()
  time.taken = end.time - start.time
  print(round(time.taken,2))
  
  tmp.len = length(fixef(fm_1)[[1]])
  
  tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
  tmp.out[,1] = rep(sp.[i], tmp.len)
  
  tmp.out[,2] = sort(as.numeric(paste(unique(drift$year))))
  
  tmp.out[,3] = exp(fixef(fm_1)[[1]])
  tmp.out[,4] = exp(confint(fm_1)[1:tmp.len,1])
  tmp.out[,5] = exp(confint(fm_1)[1:tmp.len,2]) 
  
  models[[i]] = fm_1
  
  big.out[[i]] = tmp.out
}

model.counts = do.call(rbind, big.out)
colnames(model.counts) = list("sp", "year", "est", "ll", "ul")


#-----------------------------------------------------------------------------#
# make some plots of the year effect from the fitting above 
windows(width = 10, height = 10, record = TRUE)

name.key = data.frame(n1 = c("CHIL", "GAMM", "NZMS", "SIML"),
                      n2 = c("Midges", "Gammarus", "New Zealand Mud Snail", "Black Flies"))

model.counts$sp2 = name.key[match(model.counts[,1], name.key[,1]),2]

model.counts$year2 = as.character(model.counts$year)
model.counts$year2 = as.numeric(model.counts$year)

model.counts$year2 = seq(2007, 2018,1)


p = ggplot(model.counts, aes(x = year2, y = est)) +
  geom_point() +
  geom_errorbar(aes(ymax = ul, ymin = ll)) +
  facet_wrap(~ sp2, scales = "free_y") + 
  scale_x_continuous(labels = substr(as.character(seq(2007,2018,1)),3,4), 
                     breaks = c(2007:2018)) +
  labs(title = "Long-Term Drift Monitoring",
       y = expression(paste('Count / m'^' 3')), x = "Year")  
# p + theme_base()

G = p + theme(axis.title.x = element_text(size = 14, vjust = -.1),
              axis.title.y = element_text(size = 14, vjust = 1),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              title = element_text(size = 16),
              panel.background = element_rect(fill = "white"),
              panel.grid.minor = element_line(colour = "white"),
              panel.grid.major = element_line(colour = "white"),
              panel.border = element_rect(colour = "black", fill = NA),
              # panel.spacing = unit(c(1,1,1,1), "lines"),
              strip.background = element_blank(),
              strip.text = element_text(size = 14, vjust = 1),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.title.align = .5)
G


#-----------------------------------------------------------------------------#
# run some models for each taxa looking at month effects
models.2 = NULL  # variable to save model fits 
sp. = unique(drift$SpeciesID)
big.out.2 = NULL


for(i in 1:length(sp.)){
  
  sp.sub = filter(drift, SpeciesID == sp.[i])
  
  start.time <- Sys.time()  # start timer
  fm_1 <- glmmTMB(count.sum ~ 0 + (1|year) + month + (1|RiverMile) + offset(log(Volume)),
                  data = sp.sub,
                  family = nbinom2)
  
  end.time = Sys.time()
  time.taken = end.time - start.time
  print(round(time.taken,2))
  
  tmp.len = length(fixef(fm_1)[[1]])
  
  tmp.out = data.frame(matrix(NA, tmp.len, length(unique(drift$SpeciesID))))
  tmp.out[,1] = rep(sp.[i], tmp.len)
  
  tmp.out[,2] = sort(as.numeric(paste(unique(drift$month))))
  
  tmp.out[,3] = exp(fixef(fm_1)[[1]])
  tmp.out[,4] = exp(confint(fm_1)[1:tmp.len,1])
  tmp.out[,5] = exp(confint(fm_1)[1:tmp.len,2]) 
  
  models.2[[i]] = fm_1
  
  big.out.2[[i]] = tmp.out
}

model.counts.2 = do.call(rbind, big.out.2)
colnames(model.counts.2) = list("sp", "month", "est", "ll", "ul")


#-----------------------------------------------------------------------------#
# make some month plots 
windows(width = 10, height = 10, record = TRUE)

name.key = data.frame(n1 = c("CHIL", "GAMM", "NZMS", "SIML"),
                      n2 = c("Midges", "Gammarus", "New Zealand Mud Snail", "Black Flies"))

model.counts.2$sp2 = name.key[match(model.counts.2[,1], name.key[,1]),2]



p = ggplot(model.counts.2, aes(x = month, y = est)) +
  geom_point() +
  geom_errorbar(aes(ymax = ul, ymin = ll)) + 
  scale_x_continuous(breaks = c(1:12), labels = as.character(1:12)) +
  facet_wrap(~ sp2, scales = "free_y") + 
  labs(title = "Long-Term Drift Monitoring",
       y = expression(paste('Count / m'^' 3')), x = "Month")  #
# p + theme_base()

G = p + theme(axis.title.x = element_text(size = 14, vjust = -.1),
              axis.title.y = element_text(size = 14, vjust = 1),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              title = element_text(size = 16),
              panel.background = element_rect(fill = "white"),
              panel.grid.minor = element_line(colour = "white"),
              panel.grid.major = element_line(colour = "white"),
              panel.border = element_rect(colour = "black", fill = NA),
              # panel.spacing = unit(c(1,1,1,1), "lines"),
              strip.background = element_blank(),
              strip.text = element_text(size = 14, vjust = 1),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.title.align = .5)
G


#-----------------------------------------------------------------------------#
# End
#-----------------------------------------------------------------------------#