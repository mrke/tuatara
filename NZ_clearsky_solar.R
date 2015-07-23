
getsolar<-function(longlat){
  ######################### times and location info #######################################################
  mac<-0 # choose mac (1) or pc (0)
  julnum<-365 # number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
  julday<-seq(1,365) # middle day of each month
  idayst <- 1 # start month
  ida<-julnum # end month
  if(julnum<365){
    microdaily<-0 # run microclimate model as normal, where each day is iterated 3 times starting with the initial condition of uniform soil temp at mean monthly temperature
  }else{
    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day
  }
  HEMIS <- ifelse(longlat[2]<0,2.,1.) # chose hemisphere based on latitude
  ALAT <- abs(trunc(longlat[2])) # degrees latitude
  AMINUT <- (abs(longlat[2])-ALAT)*60 # minutes latitude
  ALONG <- abs(trunc(longlat[1])) # degrees longitude
  ALMINT <- (abs(longlat[1])-ALONG)*60 # minutes latitude
  ALREF <- ALONG # reference longitude for time zone
  #########################################################################################################
  
  ############################### microclimate model parameters ###########################################
  EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
  RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
  # Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
  #IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter RUF)
  Z01 <- 0. # Top (1st) segment roughness height(m)
  Z02 <- 0. # 2nd segment roughness height(m)
  ZH1 <- 0. # Top of (1st) segment, height above surface(m)
  ZH2 <- 0. # 2nd segment, height above surface(m)  
  SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
  ERR <- 2.0 # Integrator error for soil temperature calculations
  DEP <- c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
  Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
  Density <- 2560. # soil minerals density (kg/m3)
  SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
  BulkDensity <- 1300 # soil bulk density (kg/m3)
  SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
  Clay <- 20 # clay content for matric potential calculations (%)
  SoilMoist <- 0 # fractional soil moisture (decimal %)
  REFL<-0.10 # soil reflectance (decimal %)
  ALTT<-226 # altitude (m)
  slope<-0. # slope (degrees, range 0-90)
  azmuth<-180. # aspect (degrees, 0 = North, range 0-360)
  hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
  VIEWF <- 1-sum(sin(hori*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
  PCTWET<-0 # percentage of surface area acting as a free water surface (%)
  SNOW <- rep(0,julnum) # indicates if snow is on the surface (1 is yes, 0 is no), will remove this ultimately
  CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
  TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise        													
  TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
  minshade<-0. # minimum available shade (%)
  maxshade<-90. # maximum available shade (%)
  runshade<-0 # run the model twice, once for each shade level (1) or just for the first shade level (0)?
  runmoist<-0 # run soil moisture model (0=no, 1=yes)?
  Usrhyt <- 1# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
  # Aerosol profile
  # the original profile from Elterman, L. 1970. Vertical-attenuation model with eight surface meteorological ranges 2 to 13 kilometers. U. S. Airforce Cambridge Research Laboratory, Bedford, Mass.
  #TAI<-c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
  # the values extracted from GADS for Madison
  x <- rbind(longlat) # get long/lat in a form usable by the geocode and extract function
  
  maindir<-getwd()
  setwd('../micro_global/')
  ####### get solar attenuation due to aerosols with program GADS #####################
  lat5s<-seq(-90,90,5) #lat range for GADS
  lon5s<-seq(-180,175,5) #long range for GADS
  lat5<-(45)# lat5s[which.min(abs(lat5s-x[2]))]
  lon5<-(-90)#lon5s[which.min(abs(lon5s-x[1]))]
  relhum<-1.
  season<-0.
  gadin<-list(lat5=lat5,lon5=lon5,relhum=relhum,season=season)
  source('gads/gads.R')
  gadout<-gads(gadin)
  optdep.summer<-as.data.frame(gadout$optdep)
  season<-1.
  gadin<-list(lat5=lat5,lon5=lon5,relhum=relhum,season=season)
  gadout<-gads(gadin)
  optdep.winter<-as.data.frame(gadout$optdep)
  optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
  optdep<-as.data.frame(optdep)
  colnames(optdep)<-c("LAMBDA","OPTDEPTH")
  a<-lm(OPTDEPTH~poly(LAMBDA, 6, raw=TRUE),data=optdep)
  LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
  TAI<-predict(a,data.frame(LAMBDA))
  setwd(maindir) #getting back to working directory
  
  
  ###################################################################################################################
  
  ######################### Time varying environmental data ##########################
  TMAXX<-rep(12,365) # maximum air temperatures (deg C)
  TMINN<-rep(12,365) # minimum air temperatures (deg C)
  RAINFALL<-rep(0,365) # monthly mean rainfall (mm)
  CCMAXX<-rep(0,365) # max cloud cover (%)
  CCMINN<-rep(0,365) # min cloud cover (%)
  WNMAXX<-rep(1,365) # max wind speed (m/s)
  WNMINN<-rep(1,365) # min wind speed (m/s)
  RHMAXX<-rep(50,365) # max relative humidity (%)
  RHMINN<-rep(50,365) # min relative humidity (%)
  tannul<-mean(c(TMAXX,TMINN)) # annual mean temperature for getting monthly deep soil temperature (deg C)
  tannulrun<-rep(tannul,julnum) # monthly deep soil temperature (2m) (deg C)
  SoilMoist<-rep(0,365) # soil moisture (decimal %, 1 means saturated)
  SoilMoist_Init<-rep(0.2,10) # initial soil water content, m3/m3
  # creating the arrays of environmental variables that are assumed not to change with month for this simulation 
  MAXSHADES <- rep(maxshade,julnum) # daily max shade (%)
  MINSHADES <- rep(minshade,julnum) # daily min shade (%)
  SLES <- rep(SLE,julnum) # ground emissivities
  REFLS<-rep(REFL,julnum) # soil reflectances
  PCTWET<-rep(PCTWET,julnum) # soil wetness
  ####################################################################################
  
  ################ soil properties  ################################################## 
  # set up a profile of soil properites with depth for each day to be run
  Intrvls <-(1:julnum) # user-supplied last Julian day in each time interval sequence
  Numint <- julnum  # number of time intervals
  Numtyps <- 2 # number of soil types
  Nodes <- matrix(data = 0, nrow = 10, ncol = 7300) # array of all possible soil nodes for max time span of 20 years
  Nodes[1,1:julnum]<-3 # deepest node for first substrate type
  Nodes[2,1:julnum]<-9 # deepest node for second substrate type
  #SoilMoist<-rep(SoilMoist,timeinterval) # soil moisture
  Density<-Density/1000 # density of minerals - convert to Mg/m3
  BulkDensity<-BulkDensity/1000 # density of minerals - convert to Mg/m3
  moists2<-matrix(nrow=10, ncol = julnum, data=0) # set up an empty vector for soil moisture values through time
  
  moists2[1:10,]<-SoilMoist_Init
  moists<-moists2 # final soil moisture vector
  
  # now make the soil properties matrix
  # columns are: 
  #1) bulk density (Mg/m3)
  #2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
  #3) clay content (%)
  #4) thermal conductivity (W/mK)
  #5) specific heat capacity (J/kg-K)
  #6) mineral density (Mg/m3)
  soilprops<-matrix(data = 0, nrow = 10, ncol = 6) # create an empty soil properties matrix
  soilprops[1,1]<-BulkDensity # insert soil bulk density to profile 1
  soilprops[2,1]<-BulkDensity # insert soil bulk density to profile 2
  soilprops[1,2]<-SatWater # insert saturated water content to profile 1
  soilprops[2,2]<-SatWater # insert saturated water content to profile 2
  soilprops[1,3]<-Clay     # insert percent clay to profile 1
  soilprops[2,3]<-Clay     # insertpercent clay to profile 2
  soilprops[1,4]<-Thcond # insert thermal conductivity to profile 1
  soilprops[2,4]<-Thcond # insert thermal conductivity to profile 2
  soilprops[1,5]<-SpecHeat # insert specific heat to profile 1
  soilprops[2,5]<-SpecHeat # insert specific heat to profile 2
  soilprops[1,6]<-Density # insert mineral density to profile 1
  soilprops[2,6]<-Density # insert mineral density to profile 2
  soilinit<-rep(tannul,length(DEP)) # make iniital soil temps equal to mean annual
  #########################################################################################  
  
  #  soil moisture parameters for sand (Table 9.1 in Campbell and Norman, 1995)
  PE<-rep(0.7,19) #air entry potential J/kg 
  KS<-rep(0.0058,19) #saturated conductivity, kg s/m3
  BB<-rep(1.7,19) #soil 'b' parameter
  BD<-rep(1.3,19) # soil bulk density, Mg/m3
  L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000
  LAI<-0.1 # leaf area index, used to partition traspiration/evaporation from PET
  rainmult<-1 # rainfall multiplier to impose catchment
  maxpool<-10 # max depth for water pooling on the surface, mm (to account for runoff)
  evenrain<-0 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (2)
  
  ####ignore these for now, they are currently needed as input but are only for the snow version ##########
  snowtemp<--100.5 # temperature at which precipitation falls as snow (used for snow model)
  snowdens<-0.4 # snow density (mg/m3)
  snowmelt<-1. # proportion of calculated snowmelt that doesn't refreeze
  undercatch<-1. # undercatch multipier for converting rainfall to snow
  rainmelt<-0.016 # paramter in equation that melts snow with rainfall as a function of air temp
  #########################################################################################################  
  
  # intertidal simulation input vector
  tides<-matrix(data = 0., nrow = 24*7300, ncol = 3) # make an empty matrix
  
  # microclimate input parameters list
  microinput<-c(julnum,RUF,ERR,Usrhyt,Numtyps,Numint,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain)
  
  # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
  micro<-list(microinput=microinput,tides=tides,julday=julday,SLES=SLES,DEP=DEP,Intrvls=Intrvls,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX,TMINN=TMINN,RHMAXX=RHMAXX,RHMINN=RHMINN,CCMAXX=CCMAXX,CCMINN=CCMINN,WNMAXX=WNMAXX,WNMINN=WNMINN,SNOW=SNOW,REFLS=REFLS,PCTWET=PCTWET,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists,RAINFALL=RAINFALL,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,L=L,LAI=LAI)
  setwd('../microclimate/')
  if(mac==1){
    source('microrun_mac.R') # Fortran wrapper for the microclimate model  
  }else{
    source('microrun.R') # Fortran wrapper for the microclimate model  
  }
  microut<-microclimate(micro) # run the model in Fortran
  
  solar<-as.data.frame(microut$metout[1:(julnum*24),13]) # retrieve above ground microclimatic conditions, min shade
  julday<-as.data.frame(microut$metout[1:(julnum*24),1]) # retrieve above ground microclimatic conditions, min shade
  result<-cbind(julday,solar)
  colnames(result)<-c('julday','solr')
  
  
  
  int_solar<-function(sol){
    hours<-seq(0,86400-3600,3600)  
    AUC = trapz(hours,sol)  
    return(AUC)
  }
  
  integ<-aggregate(result$solr,by=list(result$julday),FUN=int_solar)/1000000
  
  clear<-integ$x
  setwd(maindir) #getting back to working directory
  
  return(clear)
}

library(raster)
library(ncdf)
library(pracma)
DEP <- c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature

res<-0.05 # raster resolution
directories<-list.dirs('c:/Spatial_Data/Climate/New Zealand/csv_files/') # get directories (years) to loop through
directories<-directories[-1] # remove first directory
files<-list.files(directories[1]) # get all file names in current directory (year)
alldata<-read.csv(paste(directories[1],'/',files[1],sep=''),header=TRUE)

lat1<-min(alldata$Lat.)-res/2 # min latitude
lat2<-max(alldata$Lat.)+res/2 # max latitude
lon1<-min(alldata$Longt.)-res/2 # min longitude
lon2<-max(alldata$Longt.)+res/2 # max longitude
quadwid<-(lon2-lon1)/res
quadlen<-(lat2-lat1)/res

xx<-cbind(alldata[,3],alldata[,2]) # list of co-ordinates

for(i in 171:nrow(xx)){
  longlat<-xx[i,]
  clear_sky<-cbind(longlat[1],longlat[2],seq(1,365),getsolar(longlat))
  write.table(clear_sky,file = 'weather/clearsky.csv',row.names = F,col.names = F,append=TRUE,sep=',')
  cat(i,'\n')
}

gridout <- raster(ncol=quadwid, nrow=quadlen, xmn=lon1, xmx=lon2, ymn=lat1, ymx=lat2)

alldata<-read.csv('weather/clearsky.csv',head=FALSE)
colnames(alldata)<-c('long','lat','day','solarMJ')
for(i in 1:365){
  data<-subset(alldata, day==i)
  if(i==1){
    x<-cbind(data[,1],data[,2])
  }
  vals <- cbind(data[,1],data[,2],data[,4]) # list of coordinates and variables
  grid <- stack(rasterize(x, gridout, vals[,3])) # grid of all variables
  if(i==1){
    s<-grid
  }else{
    s<-stack(s,grid)
  }
  cat(i,'\n')
}
# write the results to a file
writeRaster(s, filename='clearsky.nc', overwrite=TRUE)
data<-open.ncdf( 'clearsky.nc', write=TRUE, readunlim=TRUE)
put.var.ncdf( data, varid='z', vals=seq(1,365))
close.ncdf(data) 



################# code to get cloud cover from clear sky data and observed data ###########
longlat<-c(174.3,-36)
year<-2013


solar_grid<-stack(paste('weather/',year,'_Rad.nc',sep=""))
clear_grid<-stack('weather/clearsky.nc')
x<-rbind(longlat)
solar_seq<-as.numeric(extract(solar_grid,x))
clear<-as.numeric(extract(clear_grid,x))
if(length(solar_seq)==366){ # add day for leap year if needed
 clear<-c(clear[1:59],clear[59],clear[60:365])
}

plot(clear,ylim=c(0,35))
points(solar_seq,col='red')

cloud<-(1-solar_seq/clear)*100  

cloud[cloud<0]<-0
cloud[cloud>100]<-100
plot(cloud,ylim=c(0,100))

doy<-1
for(doy in 1:10){
plot(((1-solar_grid[[doy]]/clear_grid[[doy]])*100),zlim=c(0,100),main=paste('day ',doy,' of ',year,sep=""))
}