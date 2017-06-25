## This code plots the results of Residuals intepolation improvment

# Load needed code
source("/media/qnap/Data/code/R_functions/rmspe.r")

# Load resutls of residuals Inrepolation
res <- readRDS("/media/qnap/Projects/P028.IL.Israel.MAIAC.PM.V2/work/RDS_files/resid_check/PM25_cv.rds")

# correct the names 
for (i in 1:length(res))
{
  if (any(colnames(res[[i]])=="pred.m3.mix.x"))
  {
    setnames(res[[i]],"pred.m3.mix.x","pred.m3.mix")
  }  
}


# An empty data frame for saving the performance results
df <- data.frame(stn=character(),RMSE.MIX=numeric(0), RMSE.MIX.I=numeric(0), R2.MIX=numeric(0), R2.MIX.I=numeric(0),stringsAsFactors=FALSE)

for (i in 1:length(res))
{
  # Compute residuals
  res[[i]]$resid.mix <- res[[i]]$PM25-res[[i]]$pred.m3.mix
  res[[i]]$resid.mix.I <- res[[i]]$PM25-res[[i]]$pred.m3.mix.I
  
  # Leave only days where residuals correction was applied
  res[[i]]=res[[i]][res[[i]]$resid_IDW !=0.0000000]
  # Calculate RMSE and R2
  df[i,"RMSE.MIX"] <- round(rmse(abs(res[[i]]$resid.mix)),2)
  df[i,"RMSE.MIX.I"] <- round(rmse(abs(res[[i]]$resid.mix.I)),2)
  df[i,"R2.MIX"] <- summary(lm(PM25~pred.m3.mix,data=res[[i]]))$r.squared
  df[i,"R2.MIX.I"] <- summary(lm(PM25~pred.m3.mix.I,data=res[[i]]))$r.squared
  df[i,"stn"] <- as.character(res[[i]]$stn[1])
  
  # Plot the results

  # plot(res[[2]]$PM25~res[[2]]$pred.m3.mix,ylab="Observed PM2.5",xlab="Predicted PM2.5", main=res[[2]]$stn[2],col.lab=rgb(0,0.5,0))
  # points(res[[2]]$PM25~res[[2]]$pred.m3.mix.I, col="red")
  # graphics::text(x=55,y=350,labels=paste("RMSE Mixed Model=", round(rmse(res[[2]]$resid.mix),2)))
  # graphics::text(x=55,y=300,labels=paste("RMSE Improved =", round(rmse(res[[2]]$resid.mix.I),2)))
  
}
