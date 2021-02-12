library("signal")
library(Metrics)
library(data.table)
library(ggplot2)
###
library(REddyProc) 
library(dplyr)
library(bigleaf) #umolCO2.to.gc

#before gapfilling, reogranize the original dataset, 
#like the variables saved in usmms.amerif.reorganized.csv
# Half-Hourly data set
#df_usmms_all <- read.csv("C:/zIUwork/research1/r1datasets/ameriflux/amerifluxHH/AMF_US-MMS_BASE_HR_17-5.csv",header = TRUE)
#df_usmms_9919 <- df_usmms_all[1:184080,] #1999-2019
#df_usmms_9919[df_usmms_9919 == -9999] <- NA
#df_usmms_9919$Rg_cal_swlw <- (df_usmms_9919$SW_IN_1_1_1 - df_usmms_9919$SW_OUT_1_1_1) - (df_usmms_9919$LW_OUT_1_1_1 - df_usmms_9919$LW_IN_1_1_1)
#if VPD is not provided
#df_usmms_9919$VPD <- fCalcVPDfromRHandTair(df_usmms_9919$RH_1_1_1, df_usmms_9919$TA_1_1_1)
# NEE umolCO2 m-2 s-1, Rg: (W m-2): Net radiation, Tair (deg C): Air temperature, VPD hPa, Ustar (m s-1): Friction velocity
# add year
#df_usmms_9919$Year <- as.numeric(df_usmms_9919$TIMESTAMP_START) %/% 100000000
# select the columns will be used in the study, VPD and Rg provided is before 2014,thus,we need calculate
#df_select1 <- data.frame(Year = df_usmms_9919$Year, 
#                         NEE = df_usmms_9919$FC_1_1_1, LE= df_usmms_9919$LE_1_1_1, H = df_usmms_9919$H_1_1_1, #H sensible heat 
#                         Rg = df_usmms_9919$Rg_cal_swlw, Tair = df_usmms_9919$TA_1_1_1, Tsoil = df_usmms_9919$TS_2_1_1, 
#                         rH = df_usmms_9919$RH_1_1_1, VPD = df_usmms_9919$VPD, Ustar = df_usmms_9919$USTAR_1_1_1 ,
#                         SWC = df_usmms_9919$SWC_PI_1, WS = df_usmms_9919$WS_1_1_1) #soil water content, water table depth, wind speed(m s-1) 

##stpe1: Gapfilling
df_new <- read.csv(file = "C:\\ztest\\EFI-Qing\\usmms.amerif.reorganized.csv", header = TRUE)
#make new dataframe
df_usmms <- data.frame(Year = df_new$Year, DoY = df_new$DoY, Hour = df_new$Hour,
                       NEE = df_new$NEE, LE= df_new$LE, H = df_new$H,SWC = df_new$SWC,
                       Rg = df_new$Rg, Tair = df_new$Tair, Tsoil = df_new$Tsoil, 
                       rH = df_new$rH, VPD = df_new$VPD, Ustar = df_new$Ustar,WS = df_new$WS)#(m s-1): Friction velocity
#length of na for each variable
maxbad <- nrow(df_usmms) / 2
n_Rg <- length(which(is.na(df_usmms$Rg)))
n_Tair <- length(which(is.na(df_usmms$Tair)))
n_VPD <- length(which(is.na(df_usmms$VPD)))
n_rH <- length(which(is.na(df_usmms$rH))) #relative humidity
n_LE <- length(which(is.na(df_usmms$LE))) 
n_H <- length(which(is.na(df_usmms$H)))  #sensible heat
n_Tsoil <- length(which(is.na(df_usmms$Tsoil))) 
n_SWC <- length(which(is.na(df_usmms$SWC))) 
n_WS <- length(which(is.na(df_usmms$WS)))
# Create EddyProc object
EddyDataWithPosix.F <- fConvertTimeToPosix(df_usmms, 'YDH', Year.s='Year', Day.s='DoY', Hour.s='Hour')
EddySetups.C <- sEddyProc$new('US-MMS',EddyDataWithPosix.F, c('NEE','LE','H','SWC','Rg','Tair','Tsoil','rH','VPD','Ustar','WS'))
EddySetups.C$sSetLocationInfo(Lat_deg.n=39.3232, Long_deg.n=-86.4131, TimeZone_h.n=-6)  #Location of DE-Tharandt
(uStarTh <- EddySetups.C$sEstUstarThreshold()$uStarTh)
# Gapfilling sing MDS method
if (n_Rg > 0 && n_Rg < maxbad) {
  EddySetups.C$sMDSGapFill("Rg", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_Tair > 0 && n_Tair < maxbad) {
  EddySetups.C$sMDSGapFill("Tair", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_VPD > 0 && n_VPD < maxbad) {
  EddySetups.C$sMDSGapFill("VPD", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_LE > 0 && n_LE < maxbad) {
  EddySetups.C$sMDSGapFill("LE", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_H > 0 && n_H < maxbad) {
  EddySetups.C$sMDSGapFill("H", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_Tsoil > 0 && n_Tsoil < maxbad) {
  EddySetups.C$sMDSGapFill("Tsoil", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_SWC > 0 && n_SWC < maxbad) {
  EddySetups.C$sMDSGapFill("SWC", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
if (n_WS > 0 && n_WS < maxbad) {
  EddySetups.C$sMDSGapFill("WS", FillAll.b = FALSE, V1.s = "Rg", V2.s = "VPD", V3.s = "Tair",
                           Verbose.b = FALSE)
}
colnames(EddySetups.C$sExportResults()) # Note the suffix in output columns
FilledEddyData.F <- EddySetups.C$sExportResults()
CombinedData.F <- cbind(df_usmms, FilledEddyData.F)

##Get hourly data from half hourly results
evenrow <- seq(1,nrow(CombinedData.F)-1,by = 2)#HR hour
df_HR <- CombinedData.F[evenrow,]
## Calculate daily NEE and GPP from hourly result
df_DD <- as.data.table(df_HR)[, lapply(.SD, function(x) mean(x, na.rm=TRUE)), by = c("Year","DoY")]
setcolorder(df_DD, colnames(df_HR))
df_DD$Hour <- 1
write.csv(df_DD, file = "C:\\ztest\\EFI-Qing\\gapfilled.daily.usmms.csv", row.names = FALSE)#HH half hour
