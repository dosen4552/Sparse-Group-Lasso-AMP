### This code plots alpha range for different delta
source("SGL_StateEvo_Calibration.R")

delta=0.2;g_range=seq(0,1,0.01);upper=lower=0*g_range
for(i in 1:length(g_range)){
  result=Arange(g_range[i],delta,length_precision = 10000)
  upper[i]=result$amax
  lower[i]=result$amin
}
plot(g_range,upper,type='l')
lines(g_range,lower)

delta=0.5;g_range=seq(0,1,0.01);upper=lower=0*g_range
for(i in 1:length(g_range)){
  result=Arange(g_range[i],delta,length_precision = 10000)
  upper[i]=result$amax
  lower[i]=result$amin
}
plot(g_range,upper,type='l')
lines(g_range,lower)

delta=0.8;g_range=seq(0,1,0.01);upper=lower=0*g_range
for(i in 1:length(g_range)){
  result=Arange(g_range[i],delta,length_precision = 10000)
  upper[i]=result$amax
  lower[i]=result$amin
}
par(mar=c(4,5,1,1))
plot(g_range,upper,type='l',lwd=2,xlim=c(0,1),ylim=c(0,20),xaxs='i',
     xlab=expression(gamma),ylab=expression(alpha),cex.lab=2,cex.axis=1.5)
lines(g_range,lower,lwd=2)
polygon(c(g_range,rev(g_range)[2:length(g_range)]),c(lower,rev(upper)[2:length(g_range)]), 
        border=NA,col=2,density=10)
