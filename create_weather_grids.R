library(raster)
library(ncdf)

variables<-c('MSLP','Pet','Rain','RH','SoilM','ETmp','Rad','Tmax','Tmin','VP','Wind')
res<-0.05 # raster resolution
tzone<-paste("Etc/GMT+",12,sep="") # doing it this way ignores daylight savings!
directories<-list.dirs('c:/Spatial_Data/Climate/New Zealand/csv_files/') # get directories (years) to loop through
directories<-directories[-1] # remove first directory

for(i in 1:length(directories)){
  files<-list.files(directories[i]) # get all file names in current directory (year)
  
  
  for(j in 1:length(files)){
    
    alldata<-read.csv(paste(directories[i],'/',files[j],sep=''),header=TRUE)
    
    # get dimensions to create an empty grid - only need to do this once
    if(i==1 & j==1){
      lat1<-min(alldata$Lat.)-res/2 # min latitude
      lat2<-max(alldata$Lat.)+res/2 # max latitude
      lon1<-min(alldata$Longt.)-res/2 # min longitude
      lon2<-max(alldata$Longt.)+res/2 # max longitude
      quadwid<-(lon2-lon1)/res
      quadlen<-(lat2-lat1)/res
    }
    
    gridout <- raster(ncol=quadwid, nrow=quadlen, xmn=lon1, xmx=lon2, ymn=lat1, ymx=lat2)
    
    x<-cbind(alldata[,3],alldata[,2]) # list of co-ordinates
    vals <- cbind(alldata[,3],alldata[,2],alldata[,(5:15)]) # list of coordinates and variables
    grid <- stack(rasterize(x, gridout, vals[,3:13])) # grid of all variables
    
    # stacking all days together
    if(j==1){
      MSLP<-grid[[1]]
      Pet<-grid[[2]]
      Rain<-grid[[3]]
      RH<-grid[[4]]
      SoilM<-grid[[5]]
      ETmp<-grid[[6]]
      Rad<-grid[[7]]
      Tmax<-grid[[8]]
      Tmin<-grid[[9]]
      VP<-grid[[10]]
      if(as.numeric(substr(files[j],1,4))>1996){ # only have wind data from 1997
        Wind<-grid[[11]]
      }
    }else{
      MSLP<-stack(MSLP,grid[[1]])
      Pet<-stack(Pet,grid[[2]])
      Rain<-stack(Rain,grid[[3]])
      RH<-stack(RH,grid[[4]])
      SoilM<-stack(SoilM,grid[[5]])
      ETmp<-stack(ETmp,grid[[6]])
      Rad<-stack(Rad,grid[[7]])
      Tmax<-stack(TMax,grid[[8]])
      Tmin<-stack(Tmin,grid[[9]])
      VP<-stack(VP,grid[[10]])
      if(as.numeric(substr(files[j],1,4))>1996){
        Wind<-stack(Wind,grid[[11]])
      }
    }
    cat(j,'\n')
  } # end loop through days
  
  dates<-seq(ISOdate(as.numeric(substr(files[j],1,4)),1,1,tz=tzone)-3600*12, ISOdate(as.numeric(substr(files[j],1,4))+1,1,1,tz=tzone)-3600*13, by="days")   
  
  # create a stack for writing output
  if(as.numeric(substr(files[j],1,4))>1996){
    rasters<-stack(MSLP,Pet,Rain,RH,SoilM,ETmp,Rad,TMax,Tmin,VP,Wind)
    a<-14
  }else{
    rasters<-stack(MSLP,Pet,Rain,RH,SoilM,ETmp,Rad,TMax,Tmin,VP)
    a<-13
  }
  filenames<-paste('weather/',substr(files[j],1,4),"_",variables,".nc",sep="")
  
  # write rasters, re-open and add in dates, then close again
  for(m in 1:a){ 
    writeRaster(rasters[[m]], filename=filenames[m], overwrite=TRUE)
    data<-open.ncdf( filenames[m], write=TRUE, readunlim=TRUE)
    put.var.ncdf( data, varid='z', vals=dates)
    close.ncdf(data)  
  }
  
  
} # end loop through years