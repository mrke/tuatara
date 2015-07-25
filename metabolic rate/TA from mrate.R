# read mrate (ml O2/g) vs. temperature (deg C) data from Yuni and compute Arrhenius temperature
data<-read.csv('metabolic rate/mrate_Tb.csv')
with(data,plot(mrate~Tb))

data$invT<-1/(data$Tb+273) # 1/K
data$lnmrate<-log(data$mrate)
with(data,plot(lnmrate~invT))

TA.data<-data[2:7,] # remove lowest temperature
with(TA.data,plot(lnmrate~invT))
arrhen<-lm(TA.data$lnmrate~TA.data$invT) # get Arrhenius temperature
TA<-arrhen$coefficients[2]
list(arrhen)


