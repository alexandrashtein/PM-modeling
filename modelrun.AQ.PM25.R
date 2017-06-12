##### This code runs the Model for PM2.5 estimation

##### install packages

# rm(list = ls())
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
library(splines)
library(data.table)

#sourcing
source("/media/qnap/Data/code/R_functions/CV_splits.r") 
source("/media/qnap/Data/code/R_functions/rmspe.r")

#-------------------->> RES TABLE
res <- matrix(nrow=1, ncol=32)
res <- data.frame(res)
colnames(res) <- c(
  "m1.raw","m1.R2","m1.rmspe","m1.R2.space","m1.R2.time","m1.rmspe.space" #mod1 Full
  ,"m1cv.R2","m1cv.I","m1cv.Ise","m1cv.slope","m1cv.slopese","m1cv.rmspe","m1cv.R2.space","m1cv.R2.time","m1cv.rmspe.space" #mod1 CV
  ,"m1cvloc.R2","m1cvloc.Ise","m1cvloc.slope","m1cvloc.slopese","m1cvloc.rmspe","m1cvloc.R2.space","m1cvloc.R2.time" #loc m1
  #,"m2.R2" #mod2
  ,"m3.R2","m3.rmspe","m3.R2.space","m3.R2.time","m3.rmspe.space" #mod3
  ,"m3.I","m3.Ise","m3.slope","m3.slopese",
  "type")#Extra

res$type <- c("PM25")

#load data

# AQUA data
mod1 <-readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1/mod1.AQ.2003_2015.PM25_Daily/mod1.AQ.2003_2015.PM25_Daily_Re_Clean.rds")

mod1=as.data.table(mod1)

#Raw full correlation
m1.formula <- as.formula(PM25  ~ aod_047_mean)
mod1fit = lm(m1.formula,data=mod1)
summary(mod1fit)
mod1$pred.m1 <- predict(mod1fit)
res[res$type=="PM25", 'm1.raw'] <- print(summary(lm(PM25~pred.m1,data=mod1))$r.squared)
plot(mod1$PM25~ mod1$aod_047_mean, ylim=c(0,800), xlim=c(0,4))

#based mixed model
m1.formula <- as.formula(PM25 ~ aod_047_mean
                         #temporal
                         +(1+aod_047_mean|day))  

#stage 1
mod1fit <- lmer(m1.formula,data=mod1)
summary(mod1fit)
mod1$pred.m1 <- predict(mod1fit)
res[res$type=="PM25", 'm1.R2']=print(summary(lm(PM25~pred.m1,data=mod1))$r.squared)

#RMSPE
res[res$type=="PM25", 'm1.rmspe'] <- print(rmse(residuals(mod1fit)))

#spatial
spatialall<-mod1 %>%
  group_by(stn) %>%
  dplyr::summarise(barpm = mean(PM25, na.rm=TRUE), barpred = mean(pred.m1, na.rm=TRUE)) 
m1.fit.all.s <- lm(barpm ~ barpred, data=spatialall)
res[res$type=="PM25", 'm1.R2.space'] <-print(summary(lm(barpm ~ barpred, data=spatialall))$r.squared)
res[res$type=="PM25", 'm1.rmspe.space'] <- print(rmse(residuals(m1.fit.all.s)))

#temporal
#temporal (take out daily PM from yearly mean)
tempoall<-left_join(mod1,spatialall)
tempoall$delpm <-tempoall$PM25-tempoall$barpm
tempoall$delpred <-tempoall$pred.m1-tempoall$barpred
mod_temporal <- lm(delpm ~ delpred, data=tempoall)
res[res$type=="PM25", 'm1.R2.time']<- print(summary(lm(delpm ~ delpred, data=tempoall))$r.squared)

## add weights
#  mod3 <- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/mod3.AQ.2003_2015.PM25_Daily.rds")
#  
#  mod3 = mod3 %>% as.data.frame
#  scaled = mod3[,c("Elev","Temp_D","daily_hpbl")] %>% dplyr::mutate_each(funs(scale))
#  colnames(scaled) = paste0(colnames(scaled), ".s")
#  mod3 = cbind(mod3, scaled)
#  names(mod3)
#  
# mod3=as.data.table(mod3)
# # ##calculate weights
# mod3[, m := as.numeric(format(day, "%m")) ]
# mod3<-mod3[,obs:=1]
# mod3[is.na(aod_047), obs:= 0]
# mod3<-dplyr::select(mod3,aodid,obs,Elev.s,daily_hpbl.s,m,Temp_D.s,aodid,day)
#  
# ###to save memory
# gc()
# 
# w1 <- glm(obs ~ Elev.s+Temp_D.s+daily_hpbl.s+as.factor(m),family=binomial,data=mod3)
# mod3$prob <- predict(w1 ,type = c("response"))
# mod3$wt <- 1/mod3$prob
# mod3$normwt <- mod3$wt/mean(mod3$wt)
# mod3<-dplyr::select(mod3,aodid,day,normwt)
# names(mod3)
# saveRDS(mod3, "/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/mod3.AQ.2003_2015.PM25_Daily_weights.rds")
mod3=readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/mod3.AQ.2003_2015.PM25_Daily_weights.rds")
mod3<-dplyr::select(mod3,aodid,day,normwt)
mod1$aodid=paste(formatC(round(mod1$lon_aod,3),format='f',3),formatC(round(mod1$lat_aod,3),format='f',3),sep="-")
setkey(mod3,aodid,day)
setkey(mod1,aodid,day)
mod1=merge(mod1,mod3,all.x = T)
rm(mod3)
mod1=filter(mod1,!is.na(mod1$normwt))

# Mixed effects model

m1.formula <- as.formula(PM25 ~ aod_047
                         # +aod_047_mean*c +aod_047_mean*FS_BS+FS_BS+stn_type+aod_047_mean*stn_type
                         #spatial
                         +Elev.s +ndvi.s+Pop_dens.s+Dis_Rd1_2012.s
                         # +dis_inventory.s 
                         # +road_den.s
                         # +Dist_Railw.s
                         # +Dist_WB.s
                         +P_In.s 
                         +P_Ur.s
                         +P_Ag.s+P_OS.s
                         #temporal
                         +daily_hpbl.s+Temp_D.s+Rain_D.s+RH_D.s+NO2_D.s
                         # +pbl_02.s+vc_D.s+WS_D.s+O3_D+SO2_D.s
                         +(1+aod_047|day/metreg))


#stage 1
# mod1fit <- lmer(m1.formula,data=mod1,weights=normwt)
mod1fit <- lmer(m1.formula,data=mod1)
summary(mod1fit)
mod1$pred.m1 <- predict(mod1fit)
print(summary(lm(PM25~pred.m1,data=mod1))$r.squared)
res[res$type=="PM25", 'm1.R2']=print(summary(lm(PM25~pred.m1,data=mod1))$r.squared)

#RMSPE
res[res$type=="PM25", 'm1.rmspe'] <- print(rmse(residuals(mod1fit)))

#spatial
spatialall<-mod1 %>%
  group_by(stn) %>%
  dplyr::summarise(barpm = mean(PM25, na.rm=TRUE), barpred = mean(pred.m1, na.rm=TRUE)) 
m1.fit.all.s <- lm(barpm ~ barpred, data=spatialall)
res[res$type=="PM25", 'm1.R2.space'] <-print(summary(lm(barpm ~ barpred, data=spatialall))$r.squared)
res[res$type=="PM25", 'm1.rmspe.space'] <- print(rmse(residuals(m1.fit.all.s)))

#temporal
#temporal (take out daily PM from yearly mean)
tempoall<-left_join(mod1,spatialall)
tempoall$delpm <-tempoall$PM25-tempoall$barpm
tempoall$delpred <-tempoall$pred.m1-tempoall$barpred
mod_temporal <- lm(delpm ~ delpred, data=tempoall)
res[res$type=="PM25", 'm1.R2.time']<- print(summary(lm(delpm ~ delpred, data=tempoall))$r.squared)

# #check linearity for daily mod1
# 
# x<-gam(PM25~ s(Elev.s)+s(ndvi)+s(Dis_Rd1_2012.s)+ s(Pop_dens.s) +s(P_In.s)+s(P_Ur.s) +s(P_Ag.s)+s(P_OS.s)+s(daily_hpbl.s),data=mod1)
# plot(x)
# summary(x)

#save MOD1 
saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_Daily.AQ.2003_2015.PM25_clean.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_Daily.AQ.2012.PM25.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_Daily.AQ.2012.PM25_clean.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_hourly.AQ.2012.PM25_clean.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_hourly.AQ.2012.PM25.predm1.rds")

# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_Daily.TR.2012.PM25.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_Daily.TR.2012.PM25_clean.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_hourly.TR.2012.PM25_clean.predm1.rds")
# saveRDS(mod1,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_hourly.TR.2012.PM25.predm1.rds")

# Save AQUA results
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.AQ.2012.PM25_Daily.rds")
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.AQ.2012.PM25_Daily_clean.rds")
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.AQ.2012.PM25_Hourly_clean.rds")
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.AQ.2012.PM25_Hourly.rds")  
# Save TERRA results
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.TR.2012.PM25_Daily.rds")
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.TR.2012.PM25_Daily_clean.rds")
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.TR.2012.PM25_Hourly_clean.rds")
# saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.TR.2012.PM25_Hourly.rds")

## Building cross-validation (cv) dataset

library(foreach)
library(doMC)
registerDoMC(cores=10)

test_final = foreach(i = 1:10) %dopar% {
  
  splits_s <- splitdf(mod1)
  test_s <- splits_s$testset
  train_s <- splits_s$trainset
  # out_train_s <- lmer(m1.formula,data =  train_s,weights=normwt )
  out_train_s <- lmer(m1.formula,data =  train_s)
  test_s$pred.m1.cv <- predict(object=out_train_s ,newdata=test_s,allow.new.levels=TRUE,re.form=NULL )
  test_s$iter <- i
  test_s
}

#BIND 1 dataset
mod1.cv = do.call(rbind, test_final)

# Save AQUA
saveRDS(mod1.cv,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_CV/mod1.AQ.2003_2015.PM25_Daily_clean_CV.rds")

#table updates
m1.fit.all.cv<-lm(PM25~pred.m1.cv,data=mod1.cv)
res[res$type=="PM25", 'm1cv.R2'] <- print(summary(lm(PM25~pred.m1.cv,data=mod1.cv))$r.squared)
res[res$type=="PM25", 'm1cv.I'] <-print(summary(lm(PM25~pred.m1.cv,data=mod1.cv))$coef[1,1])
res[res$type=="PM25", 'm1cv.Ise'] <-print(summary(lm(PM25~pred.m1.cv,data=mod1.cv))$coef[1,2])
res[res$type=="PM25", 'm1cv.slope'] <-print(summary(lm(PM25~pred.m1.cv,data=mod1.cv))$coef[2,1])
res[res$type=="PM25", 'm1cv.slopese'] <-print(summary(lm(PM25~pred.m1.cv,data=mod1.cv))$coef[2,2])
#RMSPE
res[res$type=="PM25", 'm1cv.rmspe'] <- print(rmse(residuals(m1.fit.all.cv)))

# spatial cv.R2 & RMSPE
spatialall.cv<-mod1.cv %>%
  dplyr::group_by(stn) %>%
  dplyr::summarise(barpm = mean(PM25, na.rm=TRUE), barpred = mean(pred.m1.cv, na.rm=TRUE)) 
m1.fit.all.cv.s <- lm(barpm ~ barpred, data=spatialall.cv)
res[res$type=="PM25", 'm1cv.R2.space'] <-  print(summary(lm(barpm ~ barpred, data=spatialall.cv))$r.squared)
res[res$type=="PM25", 'm1cv.rmspe.space'] <- print(rmse(residuals(m1.fit.all.cv.s)))

# temporal cv.R2 & RMSPE
tempoall.cv<-left_join(mod1.cv,spatialall.cv)
tempoall.cv$delpm <-tempoall.cv$PM25-tempoall.cv$barpm
tempoall.cv$delpred <-tempoall.cv$pred.m1.cv-tempoall.cv$barpred
mod_temporal.cv <- lm(delpm ~ delpred, data=tempoall.cv)
res[res$type=="PM25", 'm1cv.R2.time'] <-  print(summary(lm(delpm ~ delpred, data=tempoall.cv))$r.squared)

# Save AQUA results
saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2003_2015/results.mod1.AQ.2003_2015.PM25_Daily_clean.rds")

#local stage
head(mod1)

# Explore svm
library(e1071)
mod1.cv$res.m1<- mod1.cv$PM25-mod1.cv$pred.m1.cv

# Scale local variables
names=c("Elev_200m" ,"Pop_dens_200m","Road_dis_200m","sp_ndvi_L_200m","Dis_Railways_200m","P_In_200m","P_Ur_200m")
mod1.cv = mod1.cv %>% as.data.frame 
scaled = mod1.cv[,names] %>% dplyr::mutate_each(funs(scale))
colnames(scaled) = paste0(colnames(scaled), ".s")
mod1.cv= cbind(mod1.cv, scaled)
names(mod1.cv)
mod1.cv=as.data.table(mod1.cv)
mod1.cv<-filter(mod1.cv,!is.na(mod1.cv$Pop_dens_200m.s))

tryx <-svm(res.m1~ Elev_200m.s + Pop_dens_200m.s+Road_dis_200m.s +sp_ndvi_L_200m.s+
             Dis_Railways_200m.s +P_In_200m.s+P_Ur_200m.s ,
           type="nu-regression",cross=10,data=mod1.cv)

mod1.cv$predsvmloc <- predict(object=tryx)

## reg
mod1.cv$pred.m1.bothsvm <- mod1.cv$pred.m1 + mod1.cv$predsvmloc
res[res$type=="PM25", 'm1cvloc.R2'] <- print(summary(lm(PM25~pred.m1.bothsvm,data=mod1.cv))$r.squared)
res[res$type=="PM25", 'm1cvloc.Ise'] <-print(summary(lm(PM25~pred.m1.bothsvm,data=mod1.cv))$coef[1,2])
res[res$type=="PM25", 'm1cvloc.slope'] <-print(summary(lm(PM25~pred.m1.bothsvm,data=mod1.cv))$coef[2,1])
res[res$type=="PM25", 'm1cvloc.slopese'] <-print(summary(lm(PM25~pred.m1.bothsvm,data=mod1.cv))$coef[2,2])

#RMSPE
m1.fit.all.loc<-lm(PM25~pred.m1.bothsvm,data=mod1.cv)
res[res$type=="PM25", 'm1cvloc.rmspe'] <- print(rmse(residuals(m1.fit.all.loc)))

# spatial cv.R2 & RMSPE after local stage
spatialall.cv.loc<-mod1.cv %>%
  dplyr::group_by(stn) %>%
  dplyr::summarise(barpm = mean(PM25, na.rm=TRUE), barpred = mean(pred.m1.bothsvm, na.rm=TRUE)) 
m1.fit.all.cv.loc.s <- lm(barpm ~ barpred, data=spatialall.cv.loc)
res[res$type=="PM25", 'm1cvloc.R2.space'] <-  print(summary(lm(barpm ~ barpred, data=spatialall.cv.loc))$r.squared)
res[res$type=="PM25", 'm1cvloc.R2.time'] <- print(rmse(residuals(m1.fit.all.cv.loc.s)))

# temporal cv.R2 & RMSPE after local stage
tempoall.loc.cv<-left_join(mod1.cv,spatialall.cv.loc)
tempoall.loc.cv$delpm <-tempoall.loc.cv$PM25-tempoall.loc.cv$barpm
tempoall.loc.cv$delpred <-tempoall.loc.cv$pred.m1.both-tempoall.loc.cv$barpred
mod_temporal.loc.cv <- lm(delpm ~ delpred, data=tempoall.loc.cv)
res[res$type=="PM25", 'm1cv.R2.time'] <-  print(summary(lm(delpm ~ delpred, data=tempoall.loc.cv))$r.squared)

# Save AQUA mod1.cv
saveRDS(mod1.cv,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_CV/mod1.AQ.2003_2015.PM25_Daily_clean_CV.rds")
# Save AQUA result
saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2012/results.mod1.AQ.2003_2015.PM25_Daily_clean.rds")

### Run mod 2 

mod2 <- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod2/AQUA/mod2.AQ.2003_2015.PM25_Daily_Re.rds")
mod2=as.data.table(mod2)

# Data filtering
# mod2<-filter(mod2,RelAZ < 90)
mod2<-filter(mod2,UN < 0.04 & UN > 0)
mod2<-filter(mod2,aod_047 < 3.5)
mod2$c=substr(mod2$day,1,4)

# Predict PM for mod2 from formula craeted in the calibration stage
mod2$pred.m2 = predict(object=mod1fit,newdata=mod2,allow.new.levels=TRUE,re.form=NULL)
head(mod2)
summary(mod2$pred.m2)
mod2 <- mod2[,c("day","aodid","PM25_D_mean","lon_aod","lat_aod","pred.m2","aod_047","PM25_D_closest","PM25_IDW","c")]

gc()

# Delete implossible values
mod2<-filter(mod2,!is.na(pred.m2))
mod2<-filter(mod2,pred.m2 > 0)
summary(mod2$pred.m2)

# Add bimons
mod2$m = as.numeric(format(mod2$day, "%m")) 
mod2$bimon = (mod2$m + 1) %/% 2
mod2$bimon=paste(mod2$bimon,mod2$c,sep=".")
gc()
names(mod2)

saveRDS(mod2,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod2_pred/mod2.AQ_Daily_clean.2003_2015.PM.predm2.rds")
#keep(mod2,res,rmse,splitdf, sure=TRUE) 

#check spatial patterns by plotting a map in mod2
out <-mod2 %>%
  dplyr::group_by(aodid) %>%
  dplyr::summarise(x=mean(lon_aod, na.rm=TRUE), y =mean(lat_aod, na.rm=TRUE), predm2=mean(pred.m2, na.rm=TRUE), aodm=mean(aod_047))
out<-na.omit(out)
out=as.data.table(out)

mod2=as.data.table(mod2)

# Plot the PM2.5 Prediction map
out=out[out$predm2<40,]
ggplot() +  geom_point(data=out, aes(x=x, y=y,color=predm2),alpha=1, size=1)+
  scale_color_gradient("PM25 prediction",low="lightblue2",high="red")+ggtitle("PM25 prediction for 2003-2015")
write.csv(out,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/Temp/tst_2003_2015_with_no2_predictor.csv")

mod2.fit.s2 = lme(pred.m2 ~ PM25_IDW,random = list(aodid= ~1 + PM25_IDW),control=lmeControl(opt = "optim"), data= mod2 )
#correlate to see everything from mod2 and the mpm works
mod2$pred.m2.s2 = predict(mod2.fit.s2)
mod2$resid = residuals(mod2.fit.s2)
#check R2 
print(summary(lm(pred.m2~pred.m2.s2,data=mod2))$r.squared)

#run lme regression, this *should* include the thin plate spline yet will not run (computational limitations) thus we break it down into 2 components  

mod2= arrange(mod2,day, aodid)
mod2=as.data.table(mod2)

library(foreach)
library(doMC)
registerDoMC(cores=30)

mod2_split = split(mod2, mod2$bimon) #split the data to bi-monthly periods, in each year (overall 78 periods)

# for each period produce a fit of a gam model (interpolation of residualds) 

fits = foreach(i = 1:length(mod2_split)) %dopar% {
  
  gam(resid ~ s(lon_aod,lat_aod),  data= mod2_split[[i]] )
  
}

# use the fits (from the previous stage) for removing the residuals from mod2 prediction

final = foreach(i = 1:length(mod2_split)) %dopar% {
  
  cbind(mod2_split[[i]], pred.m2.s3 = (mod2_split[[i]]$pred.m2.s2 - fits[[i]]$fitted))
  
}

final = do.call(rbind, final)
## sort 
final= arrange(final,day, aodid)

## rerun the lme on the predictions including the spatial spline (smooth)
mod2.fit.s4 <- lme(pred.m2.s3 ~PM25_IDW ,random = list(aodid= ~1 + PM25_IDW),control=lmeControl(opt = "optim"),data= final)
# mod2.fit.s4 <- lme(pred.m2.s3 ~PM25_IDW ,random = list(aodid= ~1 + PM25_IDW+m ),control=lmeControl(opt = "optim"),data= final)
mod2$pred.m2.s4 = predict(mod2.fit.s4)
## check correlations
print(summary(lm(pred.m2 ~ pred.m2.s4,data=mod2))$r.squared) 

## real mod3 starts here :)
#mod 3 (5-8 h)
mod3 <- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3/AQUA/mod3.AQ.2003_2015.PM25_Daily.rds")
mod3$c=substr(mod3$day,1,4)
mod3 <- dplyr::select(mod3,day,aodid,PM25_D_mean,PM25_IDW,lon_aod,lat_aod,RelAZ,c)
mod3[, m := as.numeric(format(day, "%m")) ]
mod3$bimon = (mod3$m + 1) %/% 2
mod3$bimon = paste(mod3$bimon,mod3$c,sep=".")
summary(mod3)

#generate m.3 mix model predictions 

mod3= arrange(mod3,day, aodid)
mod3$pred.m3.mix <-  predict(mod2.fit.s4,mod3)
summary(mod3$pred.m3.mix)

#create unique grid
ugrid <-mod3 %>%
  dplyr::group_by(aodid) %>%
  dplyr::summarise(lon_aod = mean(lon_aod, na.rm=TRUE),  lat_aod = mean(lat_aod, na.rm=TRUE)) 

ugrid<-as.data.table(ugrid)

# Parallel process 

library(foreach)
library(doMC)
registerDoMC(cores=30)

mod3_split = split(mod3, mod3$bimon) #split back into bimons to include the gam prediction in final prediction        
uniq_gid_bimon <- ugrid

final_mod3 = foreach(i = 1:length(mod3_split)) %dopar% {
  uniq_gid_bimon$gpred <- predict.gam(fits[[i]],uniq_gid_bimon)
  setkey(uniq_gid_bimon,aodid)
  mod3_split[[i]]=as.data.table(mod3_split[[i]])
  setkey(mod3_split[[i]],aodid)
  merge(mod3_split[[i]], uniq_gid_bimon[,list(aodid,gpred)], all.x = T)
}

final_mod3 = do.call(rbind, final_mod3)
## sort 
final_mod3= arrange(final_mod3,day, aodid)

mod3=final_mod3
rm(final_mod3)
gc()

mod3=as.data.table(mod3)
# create pred.m3
mod3$pred.m3 <-mod3$pred.m3.mix+mod3$gpred
hist(mod3$pred.m3)
summary(mod3$pred.m3)

# filter unreasonable values
mod3=filter(mod3, pred.m3 >0)
mod3=as.data.table(mod3)

saveRDS(mod3,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod3_pred/mod3.pred.AQ_Daily_clean.2003_2015.rds")

keep(data.m3,mod3,res,rmse, sure=TRUE) 
gc()

#calculate stage 3 R2- CV ten folds approach will take 6 weeks...we don't currently do CV for stage 3.
mod1 <-readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod1_pred/mod1_Daily.AQ.2003_2015.PM25_clean.predm1.rds")
mod1$aodid=paste(formatC(round(mod1$lon_aod,3),format='f',3),formatC(round(mod1$lat_aod,3),format='f',3),sep="-")

mod1<-mod1[,c("aodid","day","PM25","pred.m1","stn"),with=FALSE]
#R2.m3
setkey(mod3,day,aodid)
setkey(mod1,day,aodid)
mod1 <- merge(mod1,mod3[, list(day,aodid,pred.m3)], all.x = T)
m3.fit.all<- summary(lm(PM25~pred.m3,data=mod1))
res[res$type=="PM25", 'm3.R2'] <- print(summary(lm(PM25~pred.m3,data=mod1))$r.squared)    
res[res$type=="PM25", 'm3.I'] <-print(summary(lm(PM25~pred.m3,data=mod1))$coef[1,1])
res[res$type=="PM25", 'm3.Ise'] <-print(summary(lm(PM25~pred.m3,data=mod1))$coef[1,2])
res[res$type=="PM25", 'm3.slope'] <-print(summary(lm(PM25~pred.m3,data=mod1))$coef[2,1])
res[res$type=="PM25", 'm3.slopese'] <-print(summary(lm(PM25~pred.m3,data=mod1))$coef[2,2])
#RMSPE
res[res$type=="PM25", 'm3.rmspe'] <- print(rmse(residuals(m3.fit.all)))

#spatial
###to check
spatialall<-mod1 %>%
  dplyr::group_by(stn) %>%
  dplyr::summarise(barpm = mean(PM25, na.rm=TRUE), barpred = mean(pred.m3, na.rm=TRUE)) 
m1.fit.all.spat<- lm(barpm ~ barpred, data=spatialall)
res[res$type=="PM25", 'm3.R2.space'] <-  print(summary(lm(barpm ~ barpred, data=spatialall))$r.squared)
res[res$type=="PM25", 'm3.rmspe.space'] <- print(rmse(residuals(m1.fit.all.spat)))

#temporal
tempoall<-left_join(mod1,spatialall)
tempoall$delpm <-tempoall$PM25-tempoall$barpm
tempoall$delpred <-tempoall$pred.m3-tempoall$barpred
mod_temporal <- lm(delpm ~ delpred, data=tempoall)
res[res$type=="PM25", 'm3.R2.time'] <-  print(summary(lm(delpm ~ delpred, data=tempoall))$r.squared)

saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2003_2015/results.mod1.AQ.2003_2015.PM25_Daily_clean.rds")

#create final prediction data set for use in health outcome studies

#import mod2
mod2<- readRDS( "/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/mod2_pred/mod2.AQ_Daily_clean.2003_2015.PM.predm2.rds")
mod2=as.data.table(mod2)
mod2<-mod2[,c("aodid","day","pred.m2"),with=FALSE]

#----------------> store the best available
mod3best <- mod3[, list(aodid, lon_aod, lat_aod, day, pred.m3)]
setkey(mod3best, day, aodid)
setkey(mod2, day, aodid)
mod3best <- merge(mod3best, mod2[,list(aodid, day, pred.m2)], all.x = T)
setkey(mod1,day,aodid)
mod3best <- merge(mod3best, mod1[,list(aodid,day,pred.m1,PM25)], all.x = T,allow.cartesian = T)
mod3best[,bestpred := pred.m3]
mod3best[!is.na(pred.m2),bestpred := pred.m2]
mod3best[!is.na(pred.m1),bestpred := pred.m1]
summary(mod3best$bestpred)
mod3best[bestpred < 0 , bestpred  := 0.5]
mod3best<-dplyr::select(mod3best,day,aodid,lon_aod,lat_aod,bestpred)
#save
saveRDS(mod3best,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/Best_pred/bestpred.AQ_Daily_clean.2003_2015.PM25.rds")

#save res
saveRDS(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2003_2015/results.mod1.AQ.2003_2015.PM25_Daily_clean.rds")
write.csv(res,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/summary_table/2003_2015/results.mod1.AQ.2003_2015.PM25_Daily_clean.csv")

mod3best<-filter(mod3best,!is.na(bestpred))

#save for plotting in QGIS
out <- mod3best %>% group_by(aodid) %>%
  summarise(lon=mean(lon_aod, na.rm=TRUE), lat =mean(lat_aod, na.rm=TRUE), bestpred=mean(bestpred, na.rm=TRUE))
out<-na.omit(out)
write.csv(out,"/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/Best_pred/map.bestpred.AQ_Daily_clean.2003_2015.PM25.csv")

# Plot the PM2.5 Prediction map
tiff(width = 800, height = 960,res=150,filename="/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/Maps/mod3bet_PM25_pred_2003_2015.tiff")
ggplot() +  geom_point(data=out, aes(x=lon, y=lat,color=bestpred),alpha=1, size=1)+
  scale_color_gradient("PM2.5 estimation",low="lightblue",high="red")+ ggtitle("Average PM2.5 estimation for 2003-2015")
dev.off()
############################## Explore spatial patterns per year

mod3best$year=year(mod3best$day)

years=unique(mod3best$year)

for (i in 1:length(years))
{
  mod3best_s=filter(mod3best, year==years[i]) # filter the data per spesific year
  # average the data annualy
  out <- mod3best_s %>% group_by(aodid) %>%
    summarise(lon=mean(lon_aod, na.rm=TRUE), lat =mean(lat_aod, na.rm=TRUE), bestpred=mean(bestpred, na.rm=TRUE))
  out<-na.omit(out)
  
  # plot the data
  mypath <- file.path("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/pm25_prediction_graphs",paste("pred_",years[i] , ".jpeg", sep = ""))
  g = ggplot() +  geom_point(data=out, aes(x=lon, y=lat,color=bestpred),alpha=1, size=1)+
    scale_color_gradient("PM10 prediction",low="lightgreen",high="red")+ ggtitle(paste("PM2.5 prediction for", years[i], sep = " "))
  ggsave(filename = mypath, plot = g, width = 5, height = 6)
}


# Plotting time series of certain stations

mod3best=readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/Best_pred/bestpred.AQ_Daily_clean.2003_2015.PM25.rds")
setkey(mod1,day,aodid)
setkey(mod3,day,aodid)
mod1 <- merge(mod1,mod3best[, list(day,aodid,bestpred)], all.x = T)
stn=unique(mod1$stn)
print(summary(lm(PM25~bestpred,data=mod1))$r.squared)

for (i in 1:length(stn))
{
  mod1_f=filter(mod1, stn== stn[i])
  mypath <- file.path("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/results/","Graphs_by_station",paste("PM2.5_", stn[i], ".png", sep = ""))
  png(file=mypath)
  plot(mod1_f$day,mod1_f$PM25,type= "p",ylim=c(0,100),main=paste(stn[i]),cex=.5, ylab="PM2.5 concentration",xlab="date")
  points(mod1_f$day,mod1_f$bestpred,type= "p",ylim=c(0,100),main=paste(stn[i]),cex=.5,col="red")
  legend((max(mod1_f$day)-2*365),90,c("Obs","Pred"),pch=c(1,1),col=c("black","red"))
  dev.off()
}

library(ggplot2)

ggplot(mod1, aes(x = bestpred, y = PM25)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

fit <- lm(PM25 ~ bestpred, data = mod1)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       " P = ",signif(summary(fit)$coef[2,4], 4)))
}

ggplotRegression(fit1)
