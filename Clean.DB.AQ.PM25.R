###############
#LIBS
###############
library(lme4)
library(reshape)
library(foreign) 
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(Hmisc)
library(mgcv)
library(gdata)
library(car)
library(dplyr)
library(ggmap)
library(broom)
library(splines)
library(DataCombine)

#load data

# AQUA data
mod1 <-readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM25_Daily/mod1.AQ.2003_2015.PM25_Daily_re.rds")
mod1=as.data.table(mod1)

## How many stations in each days? Before cleaning

stn_num <- mod1 %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(num=length(unique(stn)))
stn_num =as.data.table(stn_num)
summary(stn_num)
hist(stn_num$num)
nrow(stn_num[num==1,])
nrow(stn_num[num==1,])/nrow(stn_num)

# Raw cleaned correlation
mod1=as.data.table(mod1)
mod1<-filter(mod1,PM25 > 0) 
# mod1<-filter(mod1,RelAZ < 90) # Explore how using only forward scattering AOD obsrvations
mod1<-filter(mod1,UN < 0.04 & UN > 0) # Remove AOD with uncertentity values that are out of the suggested range
mod1<-filter(mod1,aod_047_mean < 3.5) # Remove AOD with very high values

# plot(mod1$PM25~ mod1$aod_047_mean, ylim=c(0,800), xlim=c(0,4))

# Use quantiles thresholds to remove PM~ aod removed records with conflicting AOD and PM10 values

x<-dplyr::select(mod1,aod_047,stn)
x$c<-1
x <- x %>%
  dplyr::group_by (stn) %>%
  dplyr::summarise(saod=sum(c))

x=as.data.table(x)
mod1=as.data.table(mod1)

#merge back count
setkey(x,stn)
setkey(mod1,stn)
mod1 <- merge(mod1,x, all.x = T)

mod1$exobs<-0
mod1<-mod1[aod_047_mean < quantile(aod_047, c(.50)) & PM25 >  quantile(PM25, c(.90)), exobs := 2]
mod1<-mod1[aod_047_mean > quantile(aod_047, c(.90)) & PM25 <  quantile(PM25, c(.50)), exobs := 3]

mod1<-mod1[saod < 10 , exobs := 5]

#take out bad exobs
mod1<-filter(mod1,exobs==0)

#based mixed model
m1.formula <- as.formula(PM25 ~ aod_047
                         #temporal
                         +(1+aod_047|day))  

#stage 1
mod1fit <- lmer(m1.formula,data=mod1)
summary(mod1fit)
mod1$pred.m1 <- predict(mod1fit)
print(summary(lm(PM25~pred.m1,data=mod1))$r.squared)
plot(mod1$PM25~mod1$aod_047)

## Filter outliers - this was not included in the final cleaned file, because it did not improve the model 
## Except for 2004

# Cleaning problematic values for 2004

m1.2004=dplyr::filter(mod1,year=="2004")
m1.2004=dplyr::filter(m1.2004,PM25<150)
m1.2004=dplyr::filter(m1.2004,aod_047<1)

mod1=filter(mod1,year!=2004)
mod1=rbind(mod1,m1.2004)

# Remove outliers in 2006
# mod1.06=filter(mod1,year == 2006)
# mod1.06=filter(mod1.06,aod_047_mean <0.9)
# mod1.06=filter(mod1.06,PM25 <200)
# 
# # Remove outliers in 2007
# mod1.07=filter(mod1,year == 2007)
# mod1.07=filter(mod1.07,aod_047_mean <1.2)
# 
# # Remove outliers in 2008
# mod1.08=filter(mod1,year == 2008)
# mod1.08=filter(mod1.08,aod_047_mean <1)
# mod1.08=filter(mod1.08,PM25 <200)
# 
# # Remove outliers in 2011
# mod1.11=filter(mod1,year == 2011)
# mod1.11=filter(mod1.11,aod_047_mean <1.2)
# mod1.11=filter(mod1.11,PM25 <250)
# 
# # Remove outliers in 2012
# mod1.12=filter(mod1,year == 2012)
# mod1.12=filter(mod1.12,PM25 <250)
# 
# # Bind all the data together
# mod1=filter(mod1,year != 2006)
# mod1=filter(mod1,year != 2007)
# mod1=filter(mod1,year != 2008)
# mod1=filter(mod1,year != 2011)
# mod1=filter(mod1,year != 2012)
# 
# mod1=rbind(mod1,mod1.06,mod1.07,mod1.08,mod1.11,mod1.12)
# mod1=as.data.table(mod1)
# mod1$day=as.Date(mod1$day)

plot(mod1$PM25~ mod1$aod_047_mean, ylim=c(0,800), xlim=c(0,4))

## Remove observations with missing data
mod1=mod1[!is.na(mod1$NO2_D.s),]

#take out station with wildly diff PM from surrounding stations (taken from Meytar's code)
mod1=as.data.table(mod1)
neveruse <- c("REM","HEF","AGR")
mod1 <-mod1[!stn %in% neveruse]

## cleaning mod1 - by meytar -These procesures were not implemented due to low observations left after applying them

#take out station with wildly diff PM from surrounding stations

# take out stn with co located PM10/25 with very high ratios
#calculate meanPM per grid per day to each station (excluding first station)
PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Pollution_stn_May16/PM25_D.csv")
PM25$date<-paste(PM25$Day,PM25$Month,PM25$Year,sep="/")
PM25[, day:=as.Date(date, "%d/%m/%Y")]
PM25[, c := as.numeric(format(day, "%Y")) ]
PM25[,c("Year","Month","Day","date"):=NULL]
PM25 <- PM25[X != 'NaN']
PM25<-PM25[!is.na(PM25)]
PM25<-PM25[PM25 > 0.000000000001 & PM25 < 900 ]
#clear non continous stations
setnames(PM25,"X","x_stn_ITM")
setnames(PM25,"Y","y_stn_ITM")
#calculate meanPM per grid per day to each station (excluding first station)
PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Pollution_stn_May16/PM10_D.csv")
PM10$date<-paste(PM10$Day,PM10$Month,PM10$Year,sep="/")
PM10[, day:=as.Date(date, "%d/%m/%Y")]
PM10[, c := as.numeric(format(day, "%Y")) ]
PM10[,c("Year","Month","Day","date"):=NULL]
PM10 <- PM10[X != 'NaN']
PM10<-PM10[!is.na(PM10)]
PM10<-PM10[PM10 > 0.000000000001 & PM10 <  2000 ]
#clear non continous stations
setnames(PM10,"X","x_stn_ITM")
setnames(PM10,"Y","y_stn_ITM")
setkey(PM10,stn,day)
setkey(PM25,stn,day)
PM.j=merge(PM10,PM25,by=c("stn","day"))
#leave only stations with both PM2.5 and PM 10 measurements
PM.j=na.omit(PM.j)
PM.j$ratio=PM.j[,PM25]/PM.j[,PM10]
PM.j[,badstn := paste(stn,day,sep="-")]
#################BAD STN
mod1[,badstn := paste(stn,day,sep="-")]
PM.j<- PM.j[ratio > 0.95]
####Take out bad stations
mod1 <- mod1[!(mod1$badstn %in% PM.j$badstn), ] 

################# clean BAD STN PM25 and check if improved model?
raWDaf <- ddply(mod1, c("stn","m"), 
                function(x) {
                  mod1 <- lm(PM25 ~ aod_047, data=x)
                  data.frame(R2 = round(summary(mod1)$r.squared, 5), 
                             nsamps = length(summary(mod1)$resid))
                })
raWDaf
raWDaf<-as.data.table(raWDaf)
bad<- raWDaf[R2 <= 0.1]
bad[,badid := paste(stn,m,sep="-")]

#################BAD STN

mod1[,badid := paste(stn,m,sep="-")]
####Take out bad stations
mod1 <- mod1[!(mod1$badid %in% bad$badid), ] 

saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM25_Daily/mod1.AQ.2003_2015.PM25_Daily_Re_Clean.rds")

## How many stations in each days? After cleaning

stn_num_ac <- mod1 %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(num=length(unique(stn)))

stn_num_ac =as.data.table(stn_num_ac)
summary(stn_num_ac)
hist(stn_num_ac$num)
nrow(stn_num_ac[num==1,])
nrow(stn_num_ac[num==1,])/nrow(stn_num_ac)
