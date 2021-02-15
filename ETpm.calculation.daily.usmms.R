library("signal")
library(Metrics)
library(data.table)
library(ggplot2)
###
library(REddyProc) 
library(dplyr)
library(bigleaf) #umolCO2.to.gc

##step2: Calculate daily ET
#find latitude and elevation
flux_lat <- 39.3232
Elev <- 275
##1.0 add canopy height and measured height for wind speed and humidity
zm <- 46#(measurement height for wind speed and humidity)
h <- 27#canopy height

##1.2 make the used data frame
df_flux <- read.csv("C:\\ztest\\EFI-Qing\\gapfilled.daily.usmms.csv",header = TRUE)
df_flux$Lat <- rep(flux_lat,nrow(df_flux))
df_PM <- data.frame(Lat = df_flux$Lat, Year = df_flux$Year, DOY = df_flux$DoY) 
#df_PM$TIMESTAMP_START <- df_flux$TIMESTAMP_START
df_PM$Tem <- df_flux$Tair_f #temperature in degree C
#df_PM$Pre <- df_flux$Pre_HR * 24 * 0.0394 #Precipitation inch/day
df_PM$LE <- df_flux$LE_f # w/m2 Latent heat flux, gapfilled using MDS method
df_PM$LE[df_PM$LE <0] <- 0 #PET is zero
df_PM$ustar <- df_flux$Ustar #Friction velocity
df_PM$Hs <- df_flux$H_f #Sensible heat flux, gapfilled using MDS method
df_PM$U <- df_flux$WS_f #Wind speed, consolidated from WS and WS_ERA

#Rn is the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl)
df_PM$Rn <- df_flux$Rg_f  #net radiation
df_PM$VPD <- df_flux$VPD_f * 0.1 #from hPa to Kpa #Vapor Pressure Deficit consolidated from VPD_F_MDS and VPD_ERA
#add lai needed 
df_PM$SWC <- df_flux$SWC_f #soil water content
df_PM$lai <- df_flux$lai
df_PM$lai_bise <- df_flux$lai_bise
##1.3 parameters
zd <- 2/3 * h #zd is zero plane displacement height 
zo <- 0.1 * h #roughness length governing momentum transfer and governing transfer of heat and vapour 

##1.4 Calculate ET
#define a variable, 'avgvar', that is 1 if tower data are half-hourly, 2 if data are hourly. 
#avgvar <- 48
#Convert LE to ET
Lv <- 2.500 * 10^6 - 2.386 *10^3*df_PM$Tem  #J/kg   #this is the latent heat of vaporization, Allen 1998 2.45 MJ/m^3
#ET <- df_PM$LE / Lv* 1.8 * avgvar  #the 1.8 is the number of seconds in a half-hour/1000
ETlv_mmday <- 86400 * df_PM$LE / Lv #kg/m2/s = 86400mm/day ET calculated directly from LE
#ETlv_inweek <- ETlv_mmday * 7 * 0.03937 #inch/week, 7 days a week
df_PM$ETlv_mmday <- ETlv_mmday
##1.5 define some constants necessary for the Penman-Monteith Equation, a few of which depend on temperature
cp <- 1006 #specific heat capacity of dry air, J/kg/°C.  Cp = 1.013 10-3 [MJ kg-1 °C-1]
rho_w <- 1000 #density of water, kg/m^3
k <- 0.4 #von Karman constant, unitless
rho_a <- 101.3*10^3/(287.058) /(df_PM$Tem+273.15) #density of dry air kg/m^3  #287.058 ideal gas law: dry air density = P/ (specific gas constant * Tk) 
#S <- 2508/(df_PM$Tem+237.3)^2 * exp(17.3*df_PM$Tem/(df_PM$Tem+237)) #slope of saturation water vapor function, #kPa/KSTst
delta <-  4098*0.6108*exp(17.27*df_PM$Tem/(df_PM$Tem+237.3))/(df_PM$Tem+237.3)^2 #slope of saturation water vapor function, #kPa/K 
pres <- 101.3 * (((293-0.0065*Elev)/293)^5.26)#air pressure, kPa
gamma <- cp * pres/(0.622*Lv)  #psychrometric constant, kPa/K,cp: 1.013 10-3 [MJ kg-1 °C-1],
Ta_K <- df_PM$Tem + 273.15 #air temp in kelvin
Rho <- (1.3079 - 0.0045*df_PM$Tem)  #density of dry air kg m-3
Cp <- 1005 #specific heat capacity of air, J/kg/°C
##1.6 correction parameter for ga
#find the aerodynamic conductance 
#OL is Monin-Obhukov lengh #stab is the atmospheric stability facture
#consult SI of Novick et al. 2016 (Nature Climate Change) for details
eps <- 0.1^8
OL <- -Rho * (df_PM$ustar + eps)^3 /(0.4*9.81*((df_PM$Hs + eps)/(Cp*Ta_K))) #the ep is a very small number, #which prevents ustart and H from being exactly zero, which mucks up the calculations
stab <- (zm-2/3*zd) / OL
psiH <- 6*log(1+stab) #atmospheric stability correction
if (length(psiH) == length(stab)) {
  psiH[which(stab <0)] <- -2*log((1+(1-16*stab[which(stab<0)])^0.5)/2)
} else {
  print ("check the dataset")
}
##1.7 ga (m/s), Gs (m/s) calculation
#atmospheric stability correction
ga <- df_PM$U * k^2 /((log((zm-zd)/zo) + psiH)^2) #m/s
df_PM$ga <- ga
#df_PM$ga[df_PM$ga > 1] <- NA
#df_PM$ga[df_PM$ga > 1] <- NA
#Finally, find the reference surface conductance by inverting the penman monteith equation
Gs <- gamma*df_PM$LE*ga / ( delta*df_PM$Rn + rho_a*cp*ga*df_PM$VPD - df_PM$LE*(delta+gamma) )                     
df_PM$Gs <- Gs
#df_PM$Gs[df_PM$Gs > 0.1] <- NA
#df_PM$Gs[df_PM$Gs < 0] <- NA
#Gs m/s
#gamma kPa/K
#LE W/m2 
#ga m/s
#delta kPa/K
#Rn W/m2
#rho_a kg/m^3
#cp j/kg/K  1j/kg/K = 1j/kg/C
#VPD kpa
#Lv J/kg
#w/m^2 = 1 J/m^2/s= 0.0864*10^6/(24*60*60) J m-2 s-1
##1.8 Peman-Monteith PET calculation
dfref <- df_PM[abs(df_PM$VPD - 1) < 0.05,]
if (nrow(dfref) == 0) {
  print("check condition setting,no data at VPD~1kpa")
} else {
  Gsref <- mean(dfref$Gs,na.rm = TRUE)
}
df_PM$Gsref <- Gsref
##1.10 calculate PM-ET
ETpm_mmday <- 86400 * (delta*df_PM$Rn + rho_a*cp*ga*df_PM$VPD) / (( delta + gamma*(1+(ga/Gs))) * Lv)
df_PM$ETpm_mmday <- ETpm_mmday
write.csv(df_PM, file = "C:\\ztest\\EFI-Qing\\ouputs.daily.usmms.csv", row.names = FALSE)#HH half hour


#plot(df_PM$ETlv_mmday,df_PM$ETpm_mmday)
