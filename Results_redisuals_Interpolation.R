## This code plots the results of Residuals intepolation improvment

# Load needed code
source("/media/qnap/Data/code/R_functions/rmspe.r")

# Load resutls of residuals Inrepolation
res <- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/resid_check/PM25_cv.rds")

# Station 1
plot(res[[1]]$PM25~res[[1]]$pred.m3.mix,ylab="Observed PM2.5",xlab="Predicted PM2.5", main=res[[1]]$stn[1])
points(res[[1]]$PM25~res[[1]]$pred.m3.mix.I, col="red")

# Compute residuals
res[[1]]$resid.mix <- res[[1]]$PM25-res[[1]]$pred.m3.mix
res[[1]]$resid.mix.I <- res[[1]]$PM25-res[[1]]$pred.m3.mix.I

# Plot station 1 results
plot(res[[1]]$PM25~res[[1]]$pred.m3.mix,ylab="Observed PM2.5",xlab="Predicted PM2.5", main=res[[1]]$stn[1],col.lab=rgb(0,0.5,0))
points(res[[1]]$PM25~res[[1]]$pred.m3.mix.I, col="red")
graphics::text(x=55,y=350,labels=paste("RMSE Mixed Model=", round(rmse(res[[1]]$resid.mix),1)))
graphics::text(x=55,y=300,labels=paste("RMSE Improved =", round(rmse(res[[1]]$resid.mix.I),1)))

# Station 2
setnames(res[[2]],"pred.m3.mix.y","pred.m3.mix")

plot(res[[2]]$PM25~res[[2]]$pred.m3.mix,ylab="Observed PM2.5",xlab="Predicted PM2.5", main=res[[2]]$stn[2],col.lab=rgb(0,0.5,0))
points(res[[2]]$PM25~res[[2]]$pred.m3.mix.I, col="red")
graphics::text(x=55,y=350,labels=paste("RMSE Mixed Model=", round(rmse(res[[2]]$resid.mix),2)))
graphics::text(x=55,y=300,labels=paste("RMSE Improved =", round(rmse(res[[2]]$resid.mix.I),2)))

