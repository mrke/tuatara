
# R Implementation of an integration of the microclimate model of Warren Porter's Niche Mapper system 
# Michael Kearney November 2013

# This version uses the National Institute of Water and Atmospheric Research (NIWA) daily 5km climate
# layers for New Zealand for air temperature, relative humidity, rainfall, wind speed (1997 onwards), cloud cover
# and soil moisture  Cloud cover is  based on daily solar layers relative to clear sky calculations from NicheMapR).
# It also uses a global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html
# Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
# Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
# by choosing the option 'rungads<-1' 

# required R packages
# raster
# sp
# ncdf
# XML
# dismo
# chron
# rgdal
# zoo
# RODBC

spatial<-"/Spatial_Data/Climate/New Zealand/weather/" # place where climate input files are kept
mac<-0

############## location and climatic data  ###################################
sitemethod <- 0 # 0=specified single site long/lat, 1=place name search using geodis (needs internet)
longlat<-c(173.82,-40.823) # Stephens Island: c(173.82,-40.823)
loc <- "Arthurs Pass, New Zealand" # type in a location here, used if option 1 is chosen above
terrain<-0 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-0 # include soil data for New Zealand (1) or not (0)?
snowmodel<-0 # run snow version? (slower!)
ystart <- 1997# start year for weather generator calibration dataset or AWAP database
yfinish <- 2007# end year for weather generator calibration dataset
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)

############# microclimate model parameters ################################
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
# Next for parameters are segmented velocity profiles due to bushes, rocks etc. on the surface, IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO!
Z01 <- 0. # Top (1st) segment roughness height(m)
Z02 <- 0. # 2nd segment roughness height(m)
ZH1 <- 0. # Top of (1st) segment, height above surface(m)
ZH2 <- 0. # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 1.25 # Integrator error for soil temperature calculations
DEP <- c(0.,1.5,  3.5, 5.,  10,  15,  30.,  60.,  100.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2560. # soil minerals density (kg/m3)
SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
BulkDensity <- 2560. # soil bulk density (kg/m3)
cap<-0 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- 22 # clay content for matric potential calculations (%)
SoilMoist <- 0 # fractional soil moisture (decimal %)
rainmult<-1 # rain multiplier for surface soil moisture (use to induce runoff), proportion
runmoist<-0 # run soil moisture model (0=no, 1=yes)?
SoilMoist_Init<-rep(0.0,10) # initial soil water content, m3/m3
evenrain<-1 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (2)
maxpool<-10 # max depth for water pooling on the surface, mm (to account for runoff)
soiltype<-8
CampNormTbl9_1<-read.csv('microclimates/CampNormTbl9_1.csv')
fieldcap<-CampNormTbl9_1[soiltype,7] # field capacity, mm
wilting<-CampNormTbl9_1[soiltype,8]  # use value from digital atlas of Australian soils # wilting point, mm
PE<-rep(CampNormTbl9_1[soiltype,4],19)
KS<-rep(CampNormTbl9_1[soiltype,6],19)
BB<-rep(CampNormTbl9_1[soiltype,5],19) 
BD<-rep(2.640,19) # Mg/m3, soil bulk density for soil moisture calcs
L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000
LAI<-0.1 # leaf area index, used to partition traspiration/evaporation from PET
REFL<-0.2 # soil reflectance (decimal %)
slope<-0. # slope (degrees, range 0-90)
aspect<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
PCTWET<-0. # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise          												
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-80. # maximum available shade (%)
runshade<-1. # run the model twice, once for each shade level (1) or just for the first shade level (0)?
manualshade<-1 # if using soildata, which includes shade, this will override the data from the database and force max shade to be the number specified above
Usrhyt <- 5# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
rainwet<-1.5 # mm rain that causes soil to become 90% wet
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.325 # snow density (mg/m3)
snowmelt<-1 # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1.2 # undercatch multipier for converting rainfall to snow
rainmelt<-0.013 # paramter in equation that melts snow with rainfall as a function of air temp.
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value
write_input<-0 # write csv files of final input to working directory? 1=yes, 0=no.

# run the model
niche<-list(mac=mac,L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
source('microclimates/NicheMapR_setup_micro_NZdaily.R')
nicheout<-NicheMapR(niche)

tzone<-paste("Etc/GMT-12",sep="") # doing it this way ignores daylight savings!
dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
ndays<-length(dates2)
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours") 
leap1<-(format(dates2, "%m/%d")!= "02/29") # used for removing leap years
leap2<-(format(dates, "%m/%d")!= "02/29") # used for removing leap years

# get output
dim<-nicheout$dim
metout<-as.data.frame(nicheout$metout[1:(dim*24),]) # above ground microclimatic conditions, min shade
shadmet<-as.data.frame(nicheout$shadmet[1:(dim*24),]) # above ground microclimatic conditions, max shade
soil<-as.data.frame(nicheout$soil[1:(dim*24),]) # soil temperatures, minimum shade
shadsoil<-as.data.frame(nicheout$shadsoil[1:(dim*24),]) # soil temperatures, maximum shade
soilmoist<-as.data.frame(nicheout$soilmoist[1:(dim*24),]) # soil water content, minimum shade
shadmoist<-as.data.frame(nicheout$shadmoist[1:(dim*24),]) # soil water content, maximum shade
humid<-as.data.frame(nicheout$humid[1:(dim*24),]) # soil humidity, minimum shade
shadhumid<-as.data.frame(nicheout$shadhumid[1:(dim*24),]) # soil humidity, maximum shade
soilpot<-as.data.frame(nicheout$soilpot[1:(dim*24),]) # soil water potential, minimum shade
shadpot<-as.data.frame(nicheout$shadpot[1:(dim*24),]) # soil water potential, maximum shade
rainfall<-as.data.frame(nicheout$RAINFALL)
MAXSHADES<-as.data.frame(nicheout$MAXSHADES)
elev<-as.numeric(nicheout$ALTT)
REFL<-as.numeric(nicheout$REFL)
longlat<-as.matrix(nicheout$longlat)
ectoin<-rbind(elev,REFL,longlat,0,0,1990,1990+nyears-1)

# removing leap years because at the moment the ectotherm model doesn't handle them
metout<-metout[leap2,]
shadmet<-shadmet[leap2,]
soil<-soil[leap2,]
shadsoil<-shadsoil[leap2,]
soilmoist<-soilmoist[leap2,]
shadmoist<-shadmoist[leap2,]
humid<-humid[leap2,]
shadhumid<-shadhumid[leap2,]
soilpot<-soilpot[leap2,]
shadpot<-shadpot[leap2,]
rainfall<-rainfall[leap1,]

# write ouput for ectotherm model
write.csv(metout,'microclimates/micro output/metout.csv')
write.csv(shadmet,'microclimates/micro output/shadmet.csv')
write.csv(soil,'microclimates/micro output/soil.csv')
write.csv(shadsoil,'microclimates/micro output/shadsoil.csv')
write.csv(soilmoist,'microclimates/micro output/soilmoist.csv')
write.csv(shadmoist,'microclimates/micro output/shadmoist.csv')
write.csv(soilpot,'microclimates/micro output/soilpot.csv')
write.csv(shadpot,'microclimates/micro output/shadpot.csv')
write.csv(humid,'microclimates/micro output/humid.csv')
write.csv(shadhumid,'microclimates/micro output/shadhumid.csv')
write.csv(rainfall,'microclimates/micro output/rainfall.csv')
write.csv(ectoin,'microclimates/micro output/ectoin.csv')
write.csv(DEP,'microclimates/micro output/DEP.csv')
write.csv(MAXSHADES,'microclimates/micro output/MAXSHADES.csv')

metout<-cbind(dates[leap2],metout)
shadmet<-cbind(dates[leap2],shadmet)
soil<-cbind(dates[leap2],soil)
shadsoil<-cbind(dates[leap2],shadsoil)
soilmoist<-cbind(dates[leap2],soilmoist)
shadmoist<-cbind(dates[leap2],shadmoist)
humid<-cbind(dates[leap2],humid)
shadhumid<-cbind(dates[leap2],shadhumid)
soilpot<-cbind(dates[leap2],soilpot)
shadpot<-cbind(dates[leap2],shadpot)

rainfall<-as.data.frame(cbind(dates2[leap1],rainfall))
colnames(rainfall)<-c('dates','rainfall')
colnames(metout)[1]<-"dates"
colnames(shadmet)[1]<-"dates"
colnames(soil)[1]<-"dates"
colnames(shadsoil)[1]<-"dates"
colnames(soilmoist)[1]<-"dates"
colnames(shadmoist)[1]<-"dates"
colnames(humid)[1]<-"dates"
colnames(shadhumid)[1]<-"dates"
colnames(soilpot)[1]<-"dates"
colnames(shadpot)[1]<-"dates"

dstart<-as.POSIXct(as.Date(paste('01/01/',ystart,sep=""), "%d/%m/%Y"))-3600*11
dfinish<-as.POSIXct(as.Date(paste('31/12/',yfinish,sep=""), "%d/%m/%Y"))-3600*10
plotsoil<-subset(soil,  soil$dates > dstart & soil$dates < dfinish )
plotmetout<-subset(metout,  metout$dates > dstart & metout$dates < dfinish )

plot(plotsoil$dates, plotsoil[,4],type='l',col = "red",lty=1,ylim = c(-10,80),ylab='temperature (C)',xlab='date')
points(plotsoil$dates, plotsoil[,5],type='l',col = 3,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,6],type='l',col = 4,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,7],type='l',col = 5,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,8],type='l',col = 6,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,9],type='l',col = 7,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,10],type='l',col = 8,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,11],type='l',col = 9,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,12],type='l',col = 10,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,13],type='l',col = 11,lty=1,ylim = c(0,50),ylab='relative humdity (%)',xlab='date')

