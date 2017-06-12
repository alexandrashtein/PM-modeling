### This code builds the 3 input databases (mod1, mod2, mod3) for the model from the new MAIAC version (08.2016) 

# rm(list=ls()) - clean the environment 

#load libraries 
library(lme4)
library(reshape)
library(foreign) 
library(ggplot2)
library(plyr)
library(data.table)
library(Hmisc)
library(mgcv)
library(gdata)
library(car)
library(dplyr)
library(ggmap)
library(broom)
library(splines)
library(DataCombine)
library(readr)
library(bit64)
library(devtools)
install_github("allanjust/aodlur", dependencies = TRUE)
library("aodlur")
library(sp)
library(rgdal)
library(stringi)

# load aod data from new MAIAC data (08.2016)
maiac=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/MAIAC_data_082016/AQUA/MAIACAAOT_Israel_2015.csv")

# # cutting the data according to the project area
pol=readOGR(dsn="/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/General/Project_border/Project_aoi","Project_border_latlon")
# pol=readOGR(dsn="/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/General/Project_border/Project_aoi","Tel_aviv")
# Convert to data.frame
maiac = as.data.frame(maiac)
# Spatial subset
coordinates(maiac) = ~ lon + lat
proj4string(maiac) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
maiac = maiac[pol, ]
# Convert back to data.table
maiac = as.data.table(maiac)

# # creat id field
maiac$aodid=paste(formatC(round(maiac$lon,3),format='f',3),formatC(round(maiac$lat,3),format='f',3),sep="-")
#
# # set names abbriviations
setnames(maiac,"AOT_Uncertainty","UN")
setnames(maiac,"AOT_QA","QA")
setnames(maiac,"Optical_Depth_047","aod_047")
setnames(maiac,"Optical_Depth_055","aod_055")
setnames(maiac,"date","day")
setnames(maiac,"lon","long_aod")
setnames(maiac,"lat","lat_aod")
maiac$date=maiac$day

saveRDS(maiac,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/raw_data/AOD.AQ.2015.RDS")

# # Use the QA data to remove problematic observations
# # The QA have to be used bofore the next stage of creating a single aod point per aodid per day

system.time(maiac[, CloudMask := as.factor(sapply(QA, function(x){paste(rev(as.integer(intToBits(x))[1:3]), collapse = "")}))])
system.time(maiac[, MaskLandWaterSnow := as.factor(sapply(QA, function(x){paste(rev(as.integer(intToBits(x))[4:5]), collapse = "")}))])
system.time(maiac[, MaskAdjacency := as.factor(sapply(QA, function(x){paste(rev(as.integer(intToBits(x))[6:8]), collapse = "")}))])
system.time(maiac[, CloudDetection := as.factor(sapply(QA, function(x){paste(rev(as.integer(intToBits(x))[9:12]), collapse = "")}))])
system.time(maiac[, GlintMask := as.factor(sapply(QA, function(x){paste(rev(as.integer(intToBits(x))[13]), collapse = "")}))])
system.time(maiac[, AerosolModel := as.factor(sapply(QA, function(x){paste(rev(as.integer(intToBits(x))[14:15]), collapse = "")}))])

# Make sure that the values that were created are reasonable
# summary(maiac$CloudMask)
# summary(maiac$MaskLandWaterSnow)
# summary(maiac$CloudDetection)
# summary(maiac$AerosolModel)
# summary(maiac$GlintMask)
# summary(maiac$MaskAdjacency)

# remove cloudy QA
maiac=filter(maiac,CloudMask!="011")
maiac=filter(maiac,CloudMask!="010")
# remove observatiobs surrounded  by more than 8 cloudy pixels QA
maiac=filter(maiac,MaskAdjacency!="010")
# # remove water QA
maiac=filter(maiac,MaskLandWaterSnow!="01")
# # remove Adjacent to snow QA
maiac=filter(maiac,MaskAdjacency!="100")
#
# # create single aod point per aodid per day
maiac <-maiac %>%
  dplyr::group_by(aodid,day) %>%
  dplyr::summarise(long_aod=mean(long_aod,na.rm=TRUE),lat_aod=mean(lat_aod,na.rm=TRUE),aod_047=mean(aod_047,na.rm=TRUE),aod_055=mean(aod_055,na.rm=TRUE),UN=mean(UN,na.rm=TRUE),RelAZ=mean(RelAZ,na.rm=TRUE))

## Add General variables

# creating a filter field of the forward scattering (FS=1) and the backward scaterring (BS=0 or else)
maiac$FS_BS=1
# # First option for data devision be Azimuth angle:
maiac <- maiac[RelAZ> 90, FS_BS := 0]

#add season
maiac$month <- as.numeric(format(maiac$day, "%m"))
#1-winter, 2-spring,3-summer,4-autum
maiac$season<-car::recode(maiac$month,"1=1;2=1;3=2;4=2;5=2;6=3;7=3;8=3;9=4;10=4;11=4;12=1")
#1-winter, 2-summer
maiac$seasonSW<-car::recode(maiac$month,"1=1;2=1;3=1;4=2;5=2;6=2;7=2;8=2;9=2;10=1;11=1;12=1")

saveRDS(maiac,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/raw_data/AOD_after_QA.AQ.2015.RDS")

# maiac=readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/raw_data/AOD_after_QA.AQ.2015.RDS")

################ add Spatial Variables
#load clipped/LU grid 
lu<-fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Keytables/1km_grid/1km_MAIAC_grid.csv")
lu$V1=NULL

# replacing aodids not indetical to the lu aodids
replace=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Temp/correcting_grid/aod_to_replace.csv")
all(replace$aodid %in% lu$aodid)
setkey(lu,aodid)
setkey(replace,aodid)
lu=merge(lu,replace,all.x = T)
lu$correct[is.na(lu$correct)] =lu$aodid[is.na(lu$correct)]
lu=as.data.table(lu)
lu[,aodid:=NULL]
setnames(lu,"correct","aodid")

all(lu$aodid %in% maiac$aodid)
mean(lu$aodid %in% maiac$aodid)
maiac <- maiac[maiac$aodid %in% lu$aodid, ]

#create full LU-aod TS
days<-seq.Date(from = as.Date("2015-01-01"), to = as.Date("2015-12-31"), 1)
#create date range
days2015 <- data.table(expand.grid(aodid = lu[, unique(aodid)], day = days))
days2015 <- data.table(expand.grid(aodid = maiac[, unique(aodid)], day = days))
days2015$aodid<-as.character(days2015$aodid)

#merge maiac data
setkey(maiac,aodid,day)
setkey(days2015,aodid,day)
db2015 <- merge(days2015,maiac, all.x = T)

# precentage of NA in the data
sum(is.na(db2015$aod_055))*100/length(db2015$aod_055)

#add land use data
setkey(db2015,aodid)
setkey(lu,aodid)
db2015 <- merge(db2015, lu, all.x = T)
summary(db2015)
gc()

db2015[,.(aodid,long_aod,lon,lat_aod,lat)]
db2015[,c("long_aod","lat_aod"):=NULL]
setnames(db2015,"lon","lon_aod")
setnames(db2015,"lat","lat_aod")

summary(db2015)

#save maiac aod data for 2015
saveRDS(db2015,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/maiac_aod/AQ.AOD.2015.data.rds")
# db2015=readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/maiac_aod/AQ.AOD.2015.data.rds")

# saving as shapefile only unique id
# db2015=as.data.table(db2015)
# setkey(db2015,aodid)
# maiac_grid=db2015[!duplicated(aodid)]
# coordinates(maiac_grid) = ~ lon_lu + lat_lu
# proj4string(maiac_grid) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# writeOGR(maiac_grid,dsn="/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Temp/", layer="db2015_all", driver="ESRI Shapefile")

################ add TEMPORAL Variables

#### add ndvi

#import NDVI
ndvi<-readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/General/NDVI/MODIS/ndvi_2000_2015.rds")
ndvi=as.data.table(ndvi)
ndvi=filter(ndvi,c=="2015")
ndvi$m<-stri_pad_left(str=ndvi$m, 2, pad="0")

db2015$m=substr(db2015$day,6,7)

# add ndviid to db2015
#join actual NDVI to aod
ndvi=as.data.table(ndvi)
setkey(ndvi, ndviid, m)
setkey(db2015,ndviid, m)
db2015<- merge(db2015, ndvi,all.x = T)

#delete unnecessery columns
db2015[,c("lat_ndvi","long_ndvi","c.y"):=NULL]

###### Add Pbl

# add daily average PBL
hpbl=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/HPBL_Israel/newmodel.2003_2015_dailyavg.csv")
hpbl$year=substr(hpbl$date,1,4)
hpbl=filter(hpbl,hpbl$year=="2015")
hpbl=as.data.table(hpbl)
setnames(hpbl,"PBLid","pblid")
setnames(hpbl,"date","day")
hpbl$day=as.Date(hpbl$day)
hpbl[,c("V1","year"):=NULL]

setkey(db2015,pblid,day)
setkey(hpbl,pblid,day)
db2015 <- merge(db2015,hpbl,all.x = T)

#add overpass PBL 
pbl<-fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/HPBL_Israel/newmodel.2003_2015_11_12am.csv")
# pbl<-fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/HPBL_Israel/newmodel.2003_2015_8_9am.csv")
pbl$pblid=paste(formatC(round(pbl$lon,3),format='f',3),formatC(round(pbl$lat,3),format='f',3),sep="-")
setnames(pbl,"date","day")
pbl$day=as.Date(pbl$day)

pbl=as.data.table(pbl)
#join pbl to aod
setkey(pbl, pblid, day )
setkey(db2015,  pblid, day)
db2015 <- merge(db2015, pbl, all.x = T)

db2015[,c("V1.y","lon","lat","time"):=NULL]
db2015=as.data.table(db2015)
summary(db2015$pbl)

# setnames(db2015,"pbl","pbl_08")
setnames(db2015,"hpbl","pbl_11")

#add night PBL         
pbl<-fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/HPBL_Israel/newmodel.2003_2015_2_3pm.csv")
pbl$pblid=paste(formatC(round(pbl$lon,3),format='f',3),formatC(round(pbl$lat,3),format='f',3),sep="-")
setnames(pbl,"date","day")
pbl$day=as.Date(pbl$day)
pbl=as.data.table(pbl)

#join pbl to aod
setkey(pbl, pblid, day )
setkey(db2015,  pblid, day)
db2015 <- merge(db2015, pbl, all.x = T)

db2015=as.data.table(db2015)
db2015[,c("V1.y","lon","lat","time","PBLid","V1"):=NULL]

setnames(db2015,"hpbl","pbl_02")
summary(db2015$pbl_02)

###### Add Temperature

## Hourly Temperature (for AQUA- average of 11:30 to 14:30)
Temp <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/Hourly_Temp_aqua_IMS_Pollution_stn.csv")
## Hourly Temperature (for TERRA- average of 09:00 to 11:30)
# Temp <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/Hourly_Temp_terra_IMS_Pollution_stn.csv")
Temp$date<-paste(Temp$Day,Temp$Month,Temp$Year,sep="/")
Temp[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
Temp[, c := as.numeric(format(day, "%Y")) ]
Temp[,c("Year","Month","Day","date"):=NULL]
Temp <- Temp[X != 'NaN']
Temp <- Temp[Temp != 'NaN']
Temp <- Temp[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = Temp, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = Temp, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "Temp", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,Temp,aodid)], all.x = T)
head(db2015)
summary(db2015$Temp)
setnames(db2015,"Temp","Temp_H") # Hourly temperature

## Daily Temperature 
Temp <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Daily_Temp_data_IMS_Pollution.csv")
Temp$date<-paste(Temp$Day,Temp$Month,Temp$Year,sep="/")
Temp[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
Temp[, c := as.numeric(format(day, "%Y")) ]
Temp[,c("Year","Month","Day","date"):=NULL]
Temp <- Temp[X != 'NaN']
Temp <- Temp[Temp != 'NaN']
Temp <- Temp[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = Temp, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = Temp, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "Temp", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,Temp,aodid)], all.x = T)
head(db2015)
summary(db2015$Temp)
setnames(db2015,"Temp","Temp_D") # Daily temperature

###### Add Hourly WS
# for AQUA
WS <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_July16/WS_H.csv")
# for TERRA
# WS <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/Terra_Hourly_data_July16/WS_H.csv")

WS$date<-paste(WS$Day,WS$Month,WS$Year,sep="/")
WS[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
WS[, c := as.numeric(format(day, "%Y")) ]
WS[,c("Year","Month","Day","date"):=NULL]
WS <- WS[X != 'NaN']
WS <- WS[WS != 'NaN']
WS <- WS[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = WS, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = WS, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "WS", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,WS,aodid)], all.x = T)
head(db2015)
summary(db2015$WS)
setnames(db2015,"WS","WS_H")

###### Add daily WS
WS <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/IMS_stn_July16/WS_D.csv")

WS$date<-paste(WS$Day,WS$Month,WS$Year,sep="/")
WS[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
WS[, c := as.numeric(format(day, "%Y")) ]
WS[,c("Year","Month","Day","date"):=NULL]
WS <- WS[X != 'NaN']
WS <- WS[WS != 'NaN']
WS <- WS[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = WS, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = WS, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "WS", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,WS,aodid)], all.x = T)
head(db2015)
summary(db2015$WS)
setnames(db2015,"WS","WS_D")

# ## ADD ventilation coefficient
db2015$vc_D=c(db2015$WS_D/(db2015$daily_hpbl*1000))
db2015$vc_H=c(db2015$WS_H/(db2015$pbl_11*1000))

###### Add Hourly RH
# for AQUA
RH <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_July16/RH_H.csv")
# for TERRA
# RH <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/Terra_Hourly_data_July16/RH_H.csv")
RH$date<-paste(RH$Day,RH$Month,RH$Year,sep="/")
RH[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
RH[, c := as.numeric(format(day, "%Y")) ]
RH[,c("Year","Month","Day","date"):=NULL]
RH <- RH[X != 'NaN']
RH <- RH[RH != 'NaN']
RH <- RH[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = RH, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = RH, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "RH", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = FALSE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,RH,aodid)], all.x = T)
head(db2015)
summary(db2015$RH)
setnames(db2015,"RH","RH_H")

###### Add Daily RH
RH <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/IMS_stn_July16/RH_D.csv")
RH$date<-paste(RH$Day,RH$Month,RH$Year,sep="/")
RH[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
RH[, c := as.numeric(format(day, "%Y")) ]
RH[,c("Year","Month","Day","date"):=NULL]
RH <- RH[X != 'NaN']
RH <- RH[RH != 'NaN']
RH <- RH[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = RH, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = RH, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "RH", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = FALSE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,RH,aodid)], all.x = T)
head(db2015)
summary(db2015$RH)
setnames(db2015,"RH","RH_D")

###### Add hourly Rain
## for AQUA
Rain <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/Hourly_Rain_AQUA_IMS_Pollution_stn.csv")
## for TERRA
# Rain <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/Hourly_Rain_terra_IMS_Pollution_stn.csv")

Rain$date<-paste(Rain$Day,Rain$Month,Rain$Year,sep="/")
Rain[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
Rain[, c := as.numeric(format(day, "%Y")) ]
Rain[,c("Year","Month","Day","date"):=NULL]
Rain <- Rain[X != 'NaN']
Rain<- Rain[Rain != 'NaN']
Rain<- Rain[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = Rain, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = Rain, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "Rain", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,Rain,aodid)], all.x = T)
head(db2015)
summary(db2015$Rain)
setnames(db2015,"Rain","Rain_H")

## Daily Rain
Rain <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Daily_Rain_Sum_IMS_Pollutants.csv")
Rain$date<-paste(Rain$Day,Rain$Month,Rain$Year,sep="/")
Rain[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
Rain[, c := as.numeric(format(day, "%Y")) ]
Rain[,c("Year","Month","Day","date"):=NULL]
Rain <- Rain[X != 'NaN']
Rain<- Rain[Rain != 'NaN']
Rain<- Rain[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = Rain, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = Rain, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "Rain", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,Rain,aodid)], all.x = T)
head(db2015)
summary(db2015$Rain)
setnames(db2015,"Rain","Rain_D")

###### Add Hourly NO2
# for AQUA
NO2 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_May16/NO2_H.csv")
# for TERRA
# NO2 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/TERRA_Hourly_data_May16/NO2_H.csv")
NO2$date<-paste(NO2$Day,NO2$Month,NO2$Year,sep="/")
NO2[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
NO2[, c := as.numeric(format(day, "%Y")) ]
NO2[,c("Year","Month","Day","date"):=NULL]
NO2 <- NO2[X != 'NaN']
NO2<- NO2[NO2 != 'NaN']
NO2<- NO2[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = NO2, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = NO2, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "NO2", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = FALSE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,NO2,aodid)], all.x = T)
head(db2015)
summary(db2015$NO2)
setnames(db2015,"NO2","NO2_H")

###### Add Daily NO2

NO2 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Pollution_stn_May16/NO2_D.csv")
NO2$date<-paste(NO2$Day,NO2$Month,NO2$Year,sep="/")
NO2[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
NO2[, c := as.numeric(format(day, "%Y")) ]
NO2[,c("Year","Month","Day","date"):=NULL]
NO2 <- NO2[X != 'NaN']
NO2<- NO2[NO2 != 'NaN']
NO2<- NO2[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = NO2, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = NO2, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "NO2", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,NO2,aodid)], all.x = T)
head(db2015)
summary(db2015$NO2)
setnames(db2015,"NO2","NO2_D")

###### Add Hourly SO2
# for AQUA 
SO2 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_May16/SO2_H.csv")
# for TERRA
# SO2 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/TERRA_Hourly_data_May16/SO2_H.csv")
SO2$date<-paste(SO2$Day,SO2$Month,SO2$Year,sep="/")
SO2[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
SO2[, c := as.numeric(format(day, "%Y")) ]
SO2[,c("Year","Month","Day","date"):=NULL]
SO2 <- SO2[X != 'NaN']
SO2<- SO2[SO2 != 'NaN']
SO2<- SO2[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = SO2, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = SO2, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "SO2", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,SO2,aodid)], all.x = T)
head(db2015)
summary(db2015$SO2)
setnames(db2015,"SO2","SO2_H")

###### Add Daily SO2

SO2 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Pollution_stn_May16/SO2_D.csv")
SO2$date<-paste(SO2$Day,SO2$Month,SO2$Year,sep="/")
SO2[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
SO2[, c := as.numeric(format(day, "%Y")) ]
SO2[,c("Year","Month","Day","date"):=NULL]
SO2 <- SO2[X != 'NaN']
SO2<- SO2[SO2 != 'NaN']
SO2<- SO2[c == 2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = SO2, 
                                xvar = "X", yvar = "Y", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = SO2, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "SO2", 
                        knearest = 15, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,SO2,aodid)], all.x = T)
head(db2015)
summary(db2015$SO2)
setnames(db2015,"SO2","SO2_D")

#### Add Hourly mean PM2.5 
# for aqua
PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_May16/PM25_H.csv")

# for terra
# PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/TERRA_Hourly_data_May16/PM25_H.csv")
PM25$date<-paste(PM25$Day,PM25$Month,PM25$Year,sep="/")
PM25[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM25[, c := as.numeric(format(day, "%Y")) ]
PM25[,c("Year","Month","Day","date"):=NULL]
PM25 <- PM25[X != 'NaN']
PM25<-PM25[!is.na(PM25)]
#clear non continous stations
setnames(PM25,"X","x_stn_ITM")
setnames(PM25,"Y","y_stn_ITM")
pmall2015<- PM25[c==2015]

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = pmall2015, 
                                xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = pmall2015, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "PM25", 
                        knearest = 9, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,PM25,aodid)], all.x = T)
head(db2015)
summary(db2015$PM25)

setnames(db2015,"PM25","PM25_H_mean")

## Add Daily PM25

# Add Daily closest PM2.5

PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Pollution_stn_May16/PM25_D.csv")
PM25$date<-paste(PM25$Day,PM25$Month,PM25$Year,sep="/")
PM25[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM25[, c := as.numeric(format(day, "%Y")) ]
PM25[,c("Year","Month","Day","date"):=NULL]
PM25 <- PM25[X != 'NaN']
PM25<-PM25[!is.na(PM25)]
#clear non continous stations
setnames(PM25,"X","x_stn_ITM")
setnames(PM25,"Y","y_stn_ITM")
pmall2015<- PM25[c==2015]


# Join the closest PM2.5 value for each day

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = pmall2015, 
                                xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = pmall2015, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "PM25", 
                        knearest = 9, maxdistance = 100000, 
                        nearestmean = FALSE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,PM25,aodid)], all.x = T)
head(db2015)
summary(db2015$PM25)

setnames(db2015,"PM25","PM25_D_closest")

# Add daily mean PM2.5 

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = pmall2015, 
                                xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = pmall2015, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "PM25", 
                        knearest = 9, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,PM25,aodid)], all.x = T)
head(db2015)
summary(db2015$PM25)

setnames(db2015,"PM25","PM25_D_mean")

## Add IDW PM2.5

# calculate IDW for PM2.5

for(i in unique(db2015$day)) {
  
  x<-pmall2015[pmall2015$day==i, ]
  y= db2015[db2015$day==i, ]
  
  library(gstat)
  #defaults to idw (gstat)
  library(sp)
  coordinates(x) = ~ x_stn_ITM + y_stn_ITM
  coordinates(y) = ~ x_aod_ITM + y_aod_ITM
  #location statment uneeded since we defined coordinates
  inter = gstat(formula = PM25 ~ 1,  data =x)
  z<-predict(object = inter, newdata = y)
  # head(z)
  db2015$pred[db2015$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:100)
}

setnames(db2015,"pred","PM25_IDW")

#### ADD Hourly PM10 
# for AQUA
PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_May16/PM10_H.csv")
# for TERRA
# PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/TERRA_Hourly_data_May16/PM10_H.csv")
PM10$date<-paste(PM10$Day,PM10$Month,PM10$Year,sep="/")
PM10[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM10[, c := as.numeric(format(day, "%Y")) ]
PM10[,c("Year","Month","Day","date"):=NULL]
PM10 <- PM10[X != 'NaN']
PM10<-PM10[!is.na(PM10)]
#clear non continous stations
setnames(PM10,"X","x_stn_ITM")
setnames(PM10,"Y","y_stn_ITM")
pm_10all2015<- PM10[c==2015]

# Join the closest PM10 value for each day

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = pm_10all2015, 
                                xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = pm_10all2015, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "PM10", 
                        knearest = 9, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,PM10,aodid)], all.x = T)
head(db2015)
summary(db2015$PM10)

setnames(db2015,"PM10","PM10_H_mean")

#### ADD Daily PM10 

PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Pollution_stn_May16/PM10_D.csv")
PM10$date<-paste(PM10$Day,PM10$Month,PM10$Year,sep="/")
PM10[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM10[, c := as.numeric(format(day, "%Y")) ]
PM10[,c("Year","Month","Day","date"):=NULL]
PM10 <- PM10[X != 'NaN']
PM10<-PM10[!is.na(PM10)]
#clear non continous stations
setnames(PM10,"X","x_stn_ITM")
setnames(PM10,"Y","y_stn_ITM")
pm_10all2015<- PM10[c==2015]

# Join the closest PM10 value for each day

jointo.pt <- makepointsmatrix(datatable = db2015, 
                              xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinfrom.pt <- makepointsmatrix(datatable = pm_10all2015, 
                                xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = db2015, joinfrom = pm_10all2015, 
                        jointovarname = "aodid", joinfromvarname = "stn", 
                        joinprefix = "nearest", valuefield = "PM10", 
                        knearest = 9, maxdistance = 100000, 
                        nearestmean = TRUE, verbose = T)

setkey(db2015,aodid,day)
setkey(joinout,aodid,day)
db2015 <- merge(db2015, joinout[,list(day,PM10,aodid)], all.x = T)
head(db2015)
summary(db2015$PM10)

setnames(db2015,"PM10","PM10_D_closest")
setnames(db2015,"PM10","PM10_D_mean")

## Add daily IDW PM10

PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Daily_Data_Yuval/Pollution_stn_May16/PM10_D.csv")
PM10$date<-paste(PM10$Day,PM10$Month,PM10$Year,sep="/")
PM10[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM10[, c := as.numeric(format(day, "%Y")) ]
PM10[,c("Year","Month","Day","date"):=NULL]
PM10 <- PM10[X != 'NaN']
PM10<-PM10[!is.na(PM10)]
#clear non continous stations
setnames(PM10,"X","x_stn_ITM")
setnames(PM10,"Y","y_stn_ITM")
pmall2015<- PM10[c==2015]

# calculate IDW for PM10

for(i in unique(db2015$day)) {
  
  x<-pmall2015[pmall2015$day==i, ]
  y= db2015[db2015$day==i, ]
  
  library(gstat)
  #defaults to idw (gstat)
  library(sp)
  coordinates(x) = ~ x_stn_ITM + y_stn_ITM
  coordinates(y) = ~ x_aod_ITM + y_aod_ITM
  #location statment uneeded since we defined coordinates
  inter = gstat(formula = PM10 ~ 1,  data =x)
  z<-predict(object = inter, newdata = y)
  # head(z)
  db2015$pred[db2015$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:100)
}

setnames(db2015,"pred","PM10_IDW")

#take out uneeded

############################################### Save MOD3 ##########################################
gc()
# saveRDS(db2015,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/TR.MAIAC.2015.mod3.rds")
# x1db2015<- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/TR.PM25.2015.mod3.rds")

saveRDS(db2015,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/AQ.MAIAC.2015.mod3.rds")
x1db2015<- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/AQ.MAIAC.2015.mod3.rds")

#calculate weights
x1db2015[, m := as.numeric(format(day, "%m")) ]
x1db2015<-x1db2015[,obs:=1]
x1db2015[is.na(aod_055), obs:= 0]
ws.2015<-dplyr::select(x1db2015,obs,Elev,daily_hpbl,m,Temp_D,aodid,day)

#to save memory
gc()

w1 <- glm(obs ~ Elev+Temp_D+pbl_11+as.factor(m),family=binomial,data=ws.2015)
ws.2015$prob <- predict(w1 ,type = c("response"))  
ws.2015$wt <- 1/ws.2015$prob
ws.2015$normwt <- ws.2015$wt/mean(ws.2015$wt)
ws.2015[, c("prob", "wt","obs","Elev", "pbl_11" , "m","Temp_D"  ) := NULL]
gc()

setkey(x1db2015,aodid,day)
setkey(ws.2015,aodid,day)
x1db2015 <- merge(x1db2015,ws.2015,all.x = T)
x1db2015[,c("m","obs"):=NULL]
# saveRDS(x1db2015,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/TR.MAIAC.2015.mod3.rds")
saveRDS(x1db2015,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/AQ.MAIAC.2015.mod3.rds")

############################################### Create MOD2 ##########################################


#SPLIT the DATA
#create mod 2 file
db2015.m2 <- db2015[!is.na(aod_055)]
#rm db2015
rm(x1db2015)
gc()

# Scale mod2

db2015.m2=readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod2/AQUA/Data_buliding_Dec16/AQ.MAIAC.2015.mod2.rds")

names=c("daily_hpbl","ndvi","Elev","dis_inventory","Dis_Mroads","road_den","Pop_dens","Dis_Rd1_2012","Dis_Rd2_2012","Dist_Railw",
        "Dist_WB","Temp_D","P_In_Min_2014","P_OS_2014","P_Ur_2014","P_Ag_2014", "P_In_2004","P_OS_2004","P_Ur_2004","P_Ag_2004",
        "WS_D","RH_D","Rain_D","NO2_D" ,"SO2_D","pbl_02","pbl_11")

db2015.m2 = db2015.m2 %>% as.data.frame 
scaled = db2015.m2[,names] %>% dplyr::mutate_each(funs(scale))
colnames(scaled) = paste0(colnames(scaled), ".s")
db2015.m2= cbind(db2015.m2, scaled)
names(db2015.m2)
db2015.m2=as.data.table(db2015.m2)

#save mod2
# saveRDS(db2015.m2,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod2/TR.MAIAC.2015.mod2.rds")
saveRDS(db2015.m2,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod2/AQUA/AQ.MAIAC.2015.mod2.rds")
gc()

############################################# building model 1 -MOD1##################################

### Create daily PM2.5 mod1
# daily database
PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Pollution_stn_May16/PM25_D.csv")

PM25$date<-paste(PM25$Day,PM25$Month,PM25$Year,sep="/")
PM25[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM25[, c := as.numeric(format(day, "%Y")) ]
PM25[,c("Year","Month","Day","date"):=NULL]
PM25 <- PM25[X != 'NaN']
PM25<-PM25[!is.na(PM25)]
PM25<-PM25[PM25 > 0.000000000001 & PM25 < 1000 ]

# Add field classification (General or transportation)
PM_Type <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM_monitors_classification/PM_monitors.csv")
setnames(PM_Type,"Code","stn")
PM_Type$stn=substr(PM_Type$stn, 2, 4)
PM25=dplyr::left_join(PM25,PM_Type,by="stn")
PM25=as.data.table(PM25)
PM25[,c("Name","Region","X.y","Y.y","Long","Lat","HASL","HAGL","Parameters"):=NULL]
PM25$stn_type<-0
PM25[Type=="'Gener'",stn_type:=1]
PM25[Type=="'Trans'",stn_type:=0]
PM25[Type=="'NaN'",stn_type:=2]

#clear non continous stations
pmall2015<- PM25[c==2015]
setnames(pmall2015,"X.x","x_stn_ITM")
setnames(pmall2015,"Y.x","y_stn_ITM")

# ADD AOD 055 to PM25 mod 1
jointo.pt <- makepointsmatrix(datatable = pmall2015, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015.m2, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = pmall2015 , joinfrom =db2015.m2, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_055", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(pmall2015,stn,day)
setkey(joinout,stn,day)
PM25.m1 <- merge(pmall2015, joinout, all.x = T)

PM25.m1<-PM25.m1[!is.na(aod_055)]
setnames(PM25.m1,"nearestmean", "aod_055_mean")

PM25.m1[,nearestknn:=NULL]
PM25.m1[,nearestnobs:=NULL]
PM25.m1[,c.y:=NULL]
setnames(PM25.m1,"c.x", "year")

# ADD AOD 047
db2015.m2_s=db2015.m2[, c("aodid","x_aod_ITM","y_aod_ITM","aod_047","day"), with = FALSE]

jointo.pt <- makepointsmatrix(datatable = PM25.m1, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable =db2015.m2_s, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = PM25.m1 , joinfrom = db2015.m2_s, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_047", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(PM25.m1,stn,day)
setkey(joinout,stn,day)
PM25.m1<- merge(PM25.m1, joinout, all.x = T)

setnames(PM25.m1,"nearestmean", "aod_047_mean")
setnames(PM25.m1,"aod_047.x", "aod_047")
setnames(PM25.m1,"x_aod_ITM.x", "x_aod_ITM")
setnames(PM25.m1,"y_aod_ITM.x", "y_aod_ITM")
PM25.m1[,c("nearest.x","nearestknn","nearestnobs","x_aod_ITM.y", "y_aod_ITM.y", "aod_047.y","nearest.y"):=NULL]

# Join 200 m spatial variables
# add 200 m key field to the database
key_field=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM25_stn_200m_keytable_id/PM25_stn_200m_keytable_id.csv")
key_field=as.data.table(key_field)
key_field$Key200_id <- paste0(key_field$POINT_X,"-",key_field$POINT_Y)
setnames(key_field, "Code","stn")
key_field$stn=substr(key_field$stn,2,4)
key_field=key_field[,.(stn,Key200_id)]

setkey(PM25.m1,stn)
setkey(key_field,stn)
PM25.m1= merge(PM25.m1,key_field,all.x = T)

lu_200m=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Keytables/200m_grid/200m_grid_spatial_Data.csv")
setnames(lu_200m,"X_Y","Key200_id")
lu_200m$V1=NULL
colnames(lu_200m) <- paste(colnames(lu_200m),"200m", sep = "_")
setnames(lu_200m,"Key200_id_200m","Key200_id")

setkey(lu_200m,Key200_id)
setkey(PM25.m1,Key200_id)
PM25.m1 <- merge(PM25.m1, lu_200m, all.x = T)
PM25.m1[,c("X_ITM_200m","Y_ITM_200m","V1_200m"):=NULL]
setnames(PM25.m1,"aod_047.x" ,"aod_047")

PM25.m1_D=PM25.m1
# # delete hourly meteorological variables
PM25.m1_D[,c("Temp_H","WS_H" ,"RH_H","Rain_H","NO2_H" ,"SO2_H","PM25_D_closest","PM25_D_mean","PM25_IDW","PM10_D_closest","PM10_IDW","PM10_H_mean","vc_H","PM25_H_mean"):=NULL]

summary(PM25.m1_D)
# Save RDS files
saveRDS(PM25.m1_D,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM25_Daily/mod1.AQ.2015.PM25_Daily.rds")

##### Create Hourly PM2.5 mod1

# hourly terra database
# PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/TERRA_Hourly_data_May16/PM25_H.csv")
# hourly aqua database
PM25 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_May16/PM25_H.csv")

PM25$date<-paste(PM25$Day,PM25$Month,PM25$Year,sep="/")
PM25[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM25[, c := as.numeric(format(day, "%Y")) ]
PM25[,c("Year","Month","Day","date"):=NULL]
PM25 <- PM25[X != 'NaN']
PM25<-PM25[!is.na(PM25)]
PM25<-PM25[PM25 > 0.000000000001 & PM25 < 1000 ]

# Add field classification (General or transportation)
PM_Type <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM_monitors_classification/PM_monitors.csv")
setnames(PM_Type,"Code","stn")
PM_Type$stn=substr(PM_Type$stn, 2, 4)
PM25=dplyr::left_join(PM25,PM_Type,by="stn")
PM25=as.data.table(PM25)
PM25[,c("Name","Region","X.y","Y.y","Long","Lat","HASL","HAGL","Parameters"):=NULL]
PM25$stn_type<-0
PM25[Type=="'Gener'",stn_type:=1]
PM25[Type=="'Trans'",stn_type:=0]
PM25[Type=="'NaN'",stn_type:=2]

#clear non continous stations
pmall2015<- PM25[c==2015]
setnames(pmall2015,"X.x","x_stn_ITM")
setnames(pmall2015,"Y.x","y_stn_ITM")

# ADD AOD 055 to PM25 mod 1
jointo.pt <- makepointsmatrix(datatable = pmall2015, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015.m2, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = pmall2015 , joinfrom = db2015.m2, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_055", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(pmall2015,stn,day)
setkey(joinout,stn,day)
PM25.m1 <- merge(pmall2015, joinout, all.x = T)

PM25.m1<-PM25.m1[!is.na(aod_055)]
setnames(PM25.m1,"nearestmean", "aod_055_mean")

PM25.m1[,nearestknn:=NULL]
PM25.m1[,nearestnobs:=NULL]
PM25.m1[,c.y:=NULL]
setnames(PM25.m1,"c.x", "year")

# ADD AOD 047
db2015.m2_s=db2015.m2[, c("aodid","x_aod_ITM","y_aod_ITM","aod_047","day"), with = FALSE]

jointo.pt <- makepointsmatrix(datatable = PM25.m1, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015.m2_s, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = PM25.m1 , joinfrom = db2015.m2_s, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_047", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(PM25.m1,stn,day)
setkey(joinout,stn,day)
PM25.m1<- merge(PM25.m1, joinout, all.x = T)

setnames(PM25.m1,"nearestmean", "aod_047_mean")
setnames(PM25.m1,"aod_047.x", "aod_047")
setnames(PM25.m1,"x_aod_ITM.x", "x_aod_ITM")
setnames(PM25.m1,"y_aod_ITM.x", "y_aod_ITM")
PM25.m1[,c("nearest.x","nearestknn","nearestnobs","x_aod_ITM.y", "y_aod_ITM.y", "aod_047.y","nearest.y"):=NULL]

# Join 200 m spatial variables
# add 200 m key field to the database
key_field=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM25_stn_200m_keytable_id/PM25_stn_200m_keytable_id.csv")
key_field=as.data.table(key_field)
key_field$Key200_id <- paste0(key_field$POINT_X,"-",key_field$POINT_Y)
setnames(key_field, "Code","stn")
key_field$stn=substr(key_field$stn,2,4)
key_field=key_field[,.(stn,Key200_id)]

setkey(PM25.m1,stn)
setkey(key_field,stn)
PM25.m1= merge(PM25.m1,key_field,all.x = T)

lu_200m=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Keytables/200m_grid/200m_grid_spatial_Data.csv")
setnames(lu_200m,"X_Y","Key200_id")
lu_200m$V1=NULL
colnames(lu_200m) <- paste(colnames(lu_200m),"200m", sep = "_")
setnames(lu_200m,"Key200_id_200m","Key200_id")

setkey(lu_200m,Key200_id)
setkey(PM25.m1,Key200_id)
PM25.m1 <- merge(PM25.m1, lu_200m, all.x = T)
PM25.m1[,c("X_ITM_200m","Y_ITM_200m","V1_200m"):=NULL]
setnames(PM25.m1,"aod_047.x" ,"aod_047")

PM25.m1_H=PM25.m1
# delete daily meteorological variables
PM25.m1_H[,c("c.y","Temp_D","WS_D" ,"RH_D","Rain_D","NO2_D","SO2_D","PM25_D_closest","PM25_D_mean","PM25_IDW","PM10_D_closest","PM10_IDW","PM10_H_mean","vc_D"):=NULL ] 
names(PM25.m1_H)
summary(PM25.m1_H)

# Save RDS files
saveRDS(PM25.m1_H,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM25_Hourly/mod1.AQ.2015.PM25_Hourly.rds")
# saveRDS(PM25.m1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.TR.2015.PM25_Daily.rds")
# saveRDS(PM25.m1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.TR.2015.PM25_Hourly.rds")

# Create mod 1 for PM10

### PM10 mod1 Daily
# daily database
PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Daily_Data/Pollution_stn_May16/PM10_D.csv")

PM10$date<-paste(PM10$Day,PM10$Month,PM10$Year,sep="/")
PM10[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM10[, c := as.numeric(format(day, "%Y")) ]
PM10[,c("Year","Month","Day","date"):=NULL]
PM10 <- PM10[X != 'NaN']
PM10<-PM10[!is.na(PM10)]
PM10<-PM10[PM10 > 0.000000000001 & PM10 <  2000 ]

# Add field classification (General or transportation)
PM_Type <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM_monitors_classification/PM_monitors.csv")
setnames(PM_Type,"Code","stn")
PM_Type$stn=substr(PM_Type$stn, 2, 4)
PM10=left_join(PM10,PM_Type,by="stn")
PM10=as.data.table(PM10)
PM10[,c("Name","Region","X.y","Y.y","Long","Lat","HASL","HAGL","Parameters"):=NULL]
PM10$stn_type<-0
PM10[Type=="'Gener'",stn_type:=1]
PM10[Type=="'Trans'",stn_type:=0]
PM10[Type=="'NaN'",stn_type:=2]

#clear non continous stations
pmall2015<- PM10[c==2015]
setnames(pmall2015,"X.x","x_stn_ITM")
setnames(pmall2015,"Y.x","y_stn_ITM")

# ADD AOD 055 to MOD1
jointo.pt <- makepointsmatrix(datatable = pmall2015, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015.m2, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = pmall2015 , joinfrom = db2015.m2, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_055", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(pmall2015,stn,day)
setkey(joinout,stn,day)
PM10.m1 <- merge(pmall2015, joinout, all.x = T)

PM10.m1<-PM10.m1[!is.na(aod_055)]
setnames(PM10.m1,"nearestmean", "aod_055_mean")

PM10.m1[,nearestknn:=NULL]
PM10.m1[,nearestnobs:=NULL]
PM10.m1[,c.y:=NULL]
setnames(PM10.m1,"c.x", "year")

# ADD AOD 047
db2015.m2_s=db2015.m2[, c("aodid","x_aod_ITM","y_aod_ITM","aod_047","day"), with = FALSE]

jointo.pt <- makepointsmatrix(datatable = PM10.m1, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015.m2_s, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = PM10.m1 , joinfrom = db2015.m2_s, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_047", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(PM10.m1,stn,day)
setkey(joinout,stn,day)
PM10.m1<- merge(PM10.m1, joinout, all.x = T)

setnames(PM10.m1,"nearestmean", "aod_047_mean")
setnames(PM10.m1,"aod_047.x", "aod_047")
PM10.m1[,c("nearest.x","nearestknn","nearestnobs","x_aod_ITM.y", "y_aod_ITM.y", "aod_047.y","nearest.y"):=NULL]

# Join 200 m spatial variables
# add 200 m key field to the database
key_field=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM25_stn_200m_keytable_id/PM25_stn_200m_keytable_id.csv")
key_field=as.data.table(key_field)
key_field$Key200_id <- paste0(key_field$POINT_X,"-",key_field$POINT_Y)
setnames(key_field, "Code","stn")
key_field$stn=substr(key_field$stn,2,4)
key_field=key_field[,.(stn,Key200_id)]

setkey(PM10.m1,stn)
setkey(key_field,stn)
PM10.m1= merge(PM10.m1,key_field,all.x = T)

lu_200m=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Keytables/200m_grid/200m_grid_spatial_Data.csv")
setnames(lu_200m,"X_Y","Key200_id")
lu_200m$V1=NULL
colnames(lu_200m) <- paste(colnames(lu_200m),"200m", sep = "_")
setnames(lu_200m,"Key200_id_200m","Key200_id")

setkey(lu_200m,Key200_id)
setkey(PM10.m1,Key200_id)
PM10.m1 <- merge(PM10.m1, lu_200m, all.x = T)
PM10.m1[,c("X_ITM_200m","Y_ITM_200m"):=NULL]
setnames(PM10.m1,"aod_047.x" ,"aod_047")

PM10.m1_D=PM10.m1
# delete unneeded variables from daily database
PM10.m1_D[,c("aod_047.x","aod_055","nearest.x","Temp_H","WS_H" ,"RH_H","Rain_H","NO2_H" ,"SO2_H","PM25_D_closest","PM25_D_mean","PM25_IDW","PM25_H_mean","PM10_H_mean","PM10_D_closest","PM10_IDW","vc_H","V1_200m","c.y","lon_200m","lat_200m","m"):=NULL] 
setnames(PM10.m1_D,"c.x","c")
names(PM10.m1_D)
summary(PM10.m1_D)

# Save RDS files
saveRDS(PM10.m1_D,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM10_Daily/mod1.AQ.2015.PM10_Daily.rds")
# saveRDS(PM10.m1_D,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.TR.2015.PM10_Daily.rds")

### PM10 mod1 Hourly

# hourly terra database
# PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/TERRA_Hourly_data_May16/PM10_H.csv")
# hourly aqua database
PM10 <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Meteorological_Data/Hourly_data/AQUA_Hourly_data_May16/PM10_H.csv")

PM10$date<-paste(PM10$Day,PM10$Month,PM10$Year,sep="/")
PM10[, day:=as.Date(strptime(date, "%d/%m/%Y",tz="GMT"))]
PM10[, c := as.numeric(format(day, "%Y")) ]
PM10[,c("Year","Month","Day","date"):=NULL]
PM10 <- PM10[X != 'NaN']
PM10<-PM10[!is.na(PM10)]
PM10<-PM10[PM10 > 0.000000000001 & PM10 <  2000 ]


# Add field classification (General or transportation)
PM_Type <- fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM_monitors_classification/PM_monitors.csv")
setnames(PM_Type,"Code","stn")
PM_Type$stn=substr(PM_Type$stn, 2, 4)
PM10=left_join(PM10,PM_Type,by="stn")
PM10=as.data.table(PM10)
PM10[,c("Name","Region","X.y","Y.y","Long","Lat","HASL","HAGL","Parameters"):=NULL]
PM10$stn_type<-0
PM10[Type=="'Gener'",stn_type:=1]
PM10[Type=="'Trans'",stn_type:=0]
PM10[Type=="'NaN'",stn_type:=2]

#clear non continous stations
pmall2015<- PM10[c==2015]
setnames(pmall2015,"X.x","x_stn_ITM")
setnames(pmall2015,"Y.x","y_stn_ITM")

#--------->mod1
#PM10
# ADD AOD 055
jointo.pt <- makepointsmatrix(datatable = pmall2015, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015.m2, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = pmall2015 , joinfrom = db2015.m2, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_055", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(pmall2015,stn,day)
setkey(joinout,stn,day)
PM10.m1 <- merge(pmall2015, joinout, all.x = T)

PM10.m1<-PM10.m1[!is.na(aod_055)]
setnames(PM10.m1,"nearestmean", "aod_055_mean")

PM10.m1[,nearestknn:=NULL]
PM10.m1[,nearestnobs:=NULL]
PM10.m1[,c.y:=NULL]
setnames(PM10.m1,"c.x", "year")

# ADD AOD 047

jointo.pt <- makepointsmatrix(datatable = PM10.m1, 
                              xvar = "x_stn_ITM", yvar = "y_stn_ITM", idvar = "stn") 

joinfrom.pt <- makepointsmatrix(datatable = db2015_s.m2, 
                                xvar = "x_aod_ITM", yvar = "y_aod_ITM", idvar = "aodid") 

joinout <- nearestbyday(jointo.pts = jointo.pt, joinfrom.pts = joinfrom.pt, 
                        jointo = PM10.m1 , joinfrom = db2015_s.m2, 
                        jointovarname = "stn", joinfromvarname = "aodid", 
                        joinprefix = "nearest", valuefield = "aod_047", 
                        knearest = 9, maxdistance = 1500, 
                        nearestmean = TRUE, verbose = T)

setkey(PM10.m1,stn,day)
setkey(joinout,stn,day)
PM10.m1<- merge(PM10.m1, joinout, all.x = T)

setnames(PM10.m1,"nearestmean", "aod_047_mean")
setnames(PM10.m1,"aod_047.x", "aod_047")
PM10.m1[,c("nearest.x","nearestknn","nearestnobs","x_aod_ITM.y", "y_aod_ITM.y", "aod_047.y","nearest.y"):=NULL]

# Join 200 m spatial variables
# add 200 m key field to the database
key_field=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Qgis/Joins/PM25_stn_200m_keytable_id/PM25_stn_200m_keytable_id.csv")
key_field=as.data.table(key_field)
key_field$Key200_id <- paste0(key_field$POINT_X,"-",key_field$POINT_Y)
setnames(key_field, "Code","stn")
key_field$stn=substr(key_field$stn,2,4)
key_field=key_field[,.(stn,Key200_id)]

setkey(PM10.m1,stn)
setkey(key_field,stn)
PM10.m1= merge(PM10.m1,key_field,all.x = T)

lu_200m=fread("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Keytables/200m_grid/200m_grid_spatial_Data.csv")
setnames(lu_200m,"X_Y","Key200_id")
lu_200m$V1=NULL
colnames(lu_200m) <- paste(colnames(lu_200m),"200m", sep = "_")
setnames(lu_200m,"Key200_id_200m","Key200_id")

setkey(lu_200m,Key200_id)
setkey(PM10.m1,Key200_id)
PM10.m1 <- merge(PM10.m1, lu_200m, all.x = T)
PM10.m1[,c("X_ITM_200m","Y_ITM_200m"):=NULL]
setnames(PM10.m1,"aod_047.x" ,"aod_047")

PM10.m1_H=PM10.m1
# delete unneeded variables from hourly database
PM10.m1_H[,c("aod_047.x","aod_055","nearest.x","Temp_D","WS_D" ,"RH_D","Rain_D","NO2_D","SO2_D","PM25_D_closest","PM25_D_mean","PM25_IDW","PM25_H_mean","PM10_H_mean","PM10_D_closest","PM10_IDW","vc_D")] =NULL 
summary(PM10.m1_H)

saveRDS(PM10.m1_H,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM10_Hourly/mod1.AQ.201.PM10_Hourly.rds")
# saveRDS(PM10.m1_H,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.TR.2015.PM10_Hourly.rds")




