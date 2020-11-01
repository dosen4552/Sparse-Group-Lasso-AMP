### This code plots various (group_info,gamma) combination to study their effect on performance

source('SGL_StateEvo_Calibration.R')
library(glmnet)
library(SGL)

############ MSE and SE against dimension for fix gamma and fix perfect group info
set.seed(1)
d=0.25;e=0.5
dim_seq=c(100,200,300,400,600,800,1000)
record1=record2=matrix(0,100,length(dim_seq));sigma=1
for(i in 1:nrow(record1)){
  start=Sys.time() # 1min per iteration
  for(j in 1:length(dim_seq)){
    pp=dim_seq[j]
    X=matrix(rnorm(pp*pp*d,sd = 1/sqrt(pp*d)),pp*d,pp)
    truth=c(rep(1,pp*e),rep(0,pp*(1-e)))
    y=X%*%truth+rnorm(pp*d)*sigma
    alpha=1;g_y=c(rep(1,pp*e),rep(2,pp*(1-e)))
    lambda=alpha_to_lambda(alpha,0.5,g_y,truth,d,sigma=sigma)
    sgl=SGL(data=list(x=X, y=y),index=g_y,lambdas = lambda/pp/d,standardize = FALSE,alpha=0.5)
    record1[i,j]=mean((sgl$beta-truth)^2)
    tau=alpha_to_tau(alpha,0.5,g_y,truth,d,sigma=sigma)
    record2[i,j]=d*(tau^2-sigma^2)
  }
  print(Sys.time()-start)
}

par(mar=c(5,5,1,1))
plot(dim_seq,colMeans(record1),type='o',lwd=2,ylim=c(0.35,0.55),xlab='Dimension p',ylab='MSE',
     cex.lab=1.5,cex.axis=1.5,pch=20)
lines(dim_seq,colMeans(record2),col=2,lwd=2,type='o',pch=20)
legend('topright',c('Empirical MSE','Estimated MSE'),lty=1,col=c(1,2),lwd=2,cex=1.5)
std=apply(record1,2,sd)
arrows(dim_seq, colMeans(record1)-std, dim_seq, colMeans(record1)+std, length=0.1, angle=90,code=3, lwd = 1.5)

### add legend for previous plot
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("bottom", c(expression(gamma ==0.1~' '),expression(gamma ==0.3 ),
                   expression(gamma ==0.5),expression(gamma ==0.7),
                   expression(gamma ==0.9)),lty=1,col=1:6,lwd=3,cex=1,ncol=5)


############ SE/Power/FDR for different gamma against lambda, fix perfect group info
set.seed(1)
p=400;d=0.25;e=0.5
X=matrix(rnorm(p*p*d),p*d,p)
truth=c(rep(1,p*e),rep(0,p*(1-e)))
y=X%*%truth # plus possible noise
g_y=c(rep(1,p*e),rep(2,p*(1-e)))

lambda_seq=c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,1.3,1.6,2,3,4)
recordMSE=recordPower=recordFDR=matrix(0,5,length(lambda_seq))
for(i in 1:nrow(recordMSE)){
  start=Sys.time() # 9min per gamma
  for(j in 1:length(lambda_seq)){
    lambda=lambda_seq[j]
    alpha=lambda_to_alpha(lambda,0.2*i-0.1,g_y,truth,d)
    tau=alpha_to_tau(alpha,0.2*i-0.1,g_y,truth,d)
    recordMSE[i,j]=d*tau^2
    
    power=fdr=0
    iter=100
    for(k in 1:iter){
      sminimizer=eta(truth/tau+rnorm(p),0.2*i-0.1,alpha,g_y)
      power=power+sum((sminimizer!=0)*(truth!=0))/max(sum(truth!=0),1)/iter
      fdr=fdr+sum((sminimizer!=0)*(truth==0))/max(sum(sminimizer!=0),1)/iter
    }
    recordPower[i,j]=power
    recordFDR[i,j]=fdr
    print(d*tau^2)
  }
  print(Sys.time()-start)
}

par(mar=c(4,5,1,1))
plot(lambda_seq,recordMSE[1,],type='o',lwd=2,ylim=c(0.3,0.51),xlab=expression(lambda),ylab='MSE',
     cex.lab=2.5,cex.axis=2,xaxs='i',pch=20)
for (col in 2:nrow(recordMSE)){
  lines(lambda_seq,recordMSE[col,],col=col,lwd=2,type='o',pch=20)
}



plot(lambda_seq,recordPower[1,],type='o',lwd=2,ylim=c(0,1),xlab=expression(lambda),ylab='Power',
     cex.lab=2.5,cex.axis=2,xaxs='i',yaxs='i',pch=20)
for (col in 2:nrow(recordPower)){
  lines(lambda_seq,recordPower[col,],col=col,lwd=2,type='o',pch=20)
}


plot(lambda_seq,recordFDR[1,],type='o',lwd=2,ylim=c(0,0.51),xlab=expression(lambda),ylab='FDR',
     cex.lab=2.5,cex.axis=2,xaxs='i',yaxs='i',pch=20)
for (col in 2:nrow(recordFDR)){
  lines(lambda_seq,recordFDR[col,],col=col,lwd=2,type='o',pch=20)
}


############# SE for different gamma against lambda, fix bad group info
g_y=rep(1,p)
record=matrix(0,5,length(lambda_seq))
for(i in 1:nrow(record)){
  start=Sys.time() # 5min per gamma
  for(j in 1:length(lambda_seq)){
    lambda=lambda_seq[j]
    alpha=lambda_to_alpha(lambda,0.2*i-0.1,g_y,truth,d)
    tau=alpha_to_tau(alpha,0.2*i-0.1,g_y,truth,d)
    record[i,j]=d*tau^2
    print(d*tau^2)
  }
  print(Sys.time()-start)
}
plot(lambda_seq,record[1,],type='o',lwd=2,ylim=c(0.37,0.51),xlab=expression(lambda),ylab='MSE',
     cex.lab=2.5,cex.axis=2,xaxs='i',yaxs='i',pch=20)
for (col in 2:nrow(record)){
  lines(lambda_seq,record[col,],col=col,lwd=2,type='o',pch=20)
}


################# fix gamma=0.5, for each correctness, plot MSE against lambda
set.seed(1)
p=400
X=matrix(rnorm(p*p*d),p*d,p)
truth=c(rep(1,p*e),rep(0,p*(1-e)))
y=X%*%truth # plus possible noise

lambda_seq=c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,1.3,1.6,2,3,4,5)
proportion=c(0,0.25,0.5,0.75,1)
recordCorrect=matrix(0,length(proportion),length(lambda_seq))
for(i in 1:nrow(recordCorrect)){
  g_y=c(rep(1,p*e*proportion[i]),rep(2,p-p*e*proportion[i]))
  start=Sys.time() # 4min per i
  for(j in 1:length(lambda_seq)){
    lambda=lambda_seq[j]
    alpha=lambda_to_alpha(lambda,0.5,g_y,truth,d)
    tau=alpha_to_tau(alpha,0.5,g_y,truth,d)
    recordCorrect[i,j]=d*tau^2
    print(d*tau^2)
  }
  print(Sys.time()-start)
}
plot(lambda_seq,recordCorrect[1,],type='o',lwd=2,ylim=c(0.34,0.51),xlab=expression(lambda),
     ylab='MSE',cex.lab=1.5,cex.axis=1.5,xaxs='i',yaxs='i',pch=20)
for (col in 2:nrow(recordCorrect)){
  lines(lambda_seq,recordCorrect[col,],col=col,lwd=2,type='o',pch=20)
}
legend('bottomright',c(as.expression(bquote('Correctness' ==0)),
                       as.expression(bquote('Correctness' ==0.25)),
                       as.expression(bquote('Correctness' ==0.5)),
                       as.expression(bquote('Correctness' ==0.75)),
                       as.expression(bquote('Correctness' ==1))),
       lty=1,col=1:nrow(recordCorrect),lwd=2,cex=1.5)









