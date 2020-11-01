#### This code runs multiple simulations of ISTA/FISTA/AMP/blockwise subgradient
#### comparing these optimizers in time, loss, MSE, number of iterations

library(SGL)
library(latex2exp)
source('SGL_StateEvo_Calibration.R')
source('SGL_optimizers.R')
source('SGL_proximal.R')

##############################Initialization (always run)######################################
iter = 500
max_iter = 1000 # to get true minimizer
gamma = 0.5
n = 2000 # number of observations
N = 4000  # number of features
delta = n / N  # ratio
g_y=rep(1,N) # group information

##################Gaussian (and others) iid design ################################
set.seed(1)
X0 <- c(rep(5,N*0.2), rep(0, N * 0.8))
A <- matrix( rnorm(n*N,mean=0,sd=sqrt(1/n)), n, N)   # Gaussain Matrix with mean 0 variance 1/n
#A=matrix(rexp(n*N,rate = sqrt(n))-1/sqrt(n),n,N)  #  iid Exponential design matrix 
#A=matrix((rbinom(n*N,1,0.5)*2-1)/sqrt(N*delta),n,N)  # iid Binomial design matrix 
sigma <- 0
y = A %*% X0+sigma*rnorm(n)

### Given alpha
alpha = 1.24
lambda=alpha_to_lambda(alpha,gamma,g_y,X0,delta,sigma = sigma) # lambda=2
print(paste0('Under current prior, for alpha= ',alpha, '; lambda= ',lambda))

sss=Sys.time() # 3min
# run a full result as if true minimizer by FISTA (ISTA for long time also works)
ans = FISTA_SGL(A,y,g_y,gamma,lambda,max_iter = max_iter)
# get per iteration result until max_iter for ISTA/FISTA/AMP
FISTA_list = list(grad=ans$grad[1:iter],time=ans$time[1:iter])
ISTA_list = ISTA_SGL(A,y,g_y,gamma,lambda,max_iter = iter)
AMP_list = AMP_SGL_fromAlpha_list(A,y,g_y,gamma,alpha,X0,max_iter = iter)
Sys.time()-sss

### get blockwise subgradient descent result; this is not per iteration so select some iterations
iters_try=c(1:5,seq(6,10,2),15,20,30,40,50,100)
# initialize
cost_Subgradient_list=time_Subgradient=mse_Subgradient_list=iters_try*0
# total 15min, directly save MSE and cost
for (i in 1:length(iters_try)){
  time_start=Sys.time()
  subgrad=SGL(data=list(x=A, y=y),index=g_y,lambdas = lambda/n,standardize = FALSE,alpha=gamma,maxit = iters_try[i])$beta
  mse_Subgradient_list[i]=mean((subgrad-ans$grad[[max_iter]])^2)
  cost_Subgradient_list[i]=full_cost(y,A,subgrad,gamma,lambda,g_y)
  print(Sys.time()-time_start)
  time_Subgradient[i]=difftime(Sys.time(),time_start,units='secs')
  print(paste0('For ',iters_try[i],' iterations, MSE is ',mse_Subgradient_list[i],
               ' cost is ',cost_Subgradient_list[i]))
}
# save result
mse_Subgradient_list=c(mean(ans$grad[[max_iter]]^2),mse_Subgradient_list)

######################################################
### calculate per iteration MSE/cost for ISTA, FISTA, AMP

# initialize
mse_ISTA_list =mse_FISTA_list =mse_AMP_list = rep(0,iter)
cost_ISTA_list =cost_FISTA_list =cost_AMP_list = rep(0,iter)
# compute
for(i in 1:iter){
  mse_ISTA_list[i] =  mean((ISTA_list$grad[[i]] - ans$grad[[max_iter]])^2 )
  mse_FISTA_list[i] = mean((FISTA_list$grad[[i]] - ans$grad[[max_iter]])^2 )
  mse_AMP_list[i] = mean((AMP_list$grad[[i]] - ans$grad[[max_iter]])^2 )
  cost_ISTA_list[i] =  full_cost(y,A,ISTA_list$grad[[i]],gamma,lambda,g_y)
  cost_FISTA_list[i] = full_cost(y,A,FISTA_list$grad[[i]],gamma,lambda,g_y)
  cost_AMP_list[i] = full_cost(y,A,AMP_list$grad[[i]],gamma,lambda,g_y)
}

################# MSE plot for ISTA/FISTA/AMP/blockwise
par(mar=c(4,5,1,1))
plot((1:iter), mse_ISTA_list[1:iter],lwd = 3, type = 'l', col = 'blue',xlim = c(0,300),lty = 3,xlab = 'Iteration',ylab = 'Optimization error',cex.lab = 2, cex.axis = 1.5)
lines((1:iter), mse_FISTA_list[1:iter], lwd = 3, type = 'l', col = 'green',lty = 2)
lines((1:iter), mse_AMP_list[1:iter],lwd = 2 ,type = 'l', col = 'red',lty = 1)
lines(c(1,1+iters_try,300), c(mse_Subgradient_list,0),lwd = 2 ,type = 'l', col = 'black',lty = 4)
legend('topright', legend=c("AMP", "FISTA","ISTA",'Blockwise'),
       col=c("red", "green","blue",'black'), lty=1:4,lwd=2,cex=1.5)


################# cost(loss) plot for ISTA/FISTA/AMP/blockwise
par(mar=c(5,5,2,1))
plot(ISTA_list$time,cost_ISTA_list,lwd = 3, type = 'l', col = 'blue',xlim = c(0,38),lty = 3,
     xlab = 'Time in sec',ylab = 'Loss',xaxs='i',cex.lab = 2, cex.axis = 1.5)
lines(FISTA_list$time,cost_FISTA_list, lwd = 3, type = 'l', col = 'green',lty = 2)
lines(AMP_list$time,cost_AMP_list,lwd = 2 ,type = 'l', col = 'red',lty = 1)
lines(time_Subgradient,cost_Subgradient_list,lwd = 3 ,type = 'l', col = 'black',lty = 4)
legend('topright', legend=c("AMP", "FISTA","ISTA",'Blockwise'),
       col=c("red", "green","blue",'black'), lty=1:4,lwd=3,cex=1.5)


###############################
########## table plot comparing ISTA/FISTA/AMP

# setting
set.seed(1)
max_iter = 1000

n=2000
N=4000
X0 <- c(rep(5,N*0.1), rep(0, N * 0.9))
A <- matrix( rnorm(n*N,mean=0,sd=sqrt(1/n)), n, N)   # Gaussain Matrix with mean 0 variance 1/n
sigma <- 0
y = A %*% X0+sigma*rnorm(n)

# Given alpha
alpha = 1.2
lambda=alpha_to_lambda(alpha,gamma,g_y,X0,delta,sigma = sigma) # lambda=1
print(paste0('Under current prior, for alpha= ',alpha, '; lambda= ',lambda))

# run a full result as if true minimizer by FISTA (ISTA for long time also works)
ans = FISTA_SGL(A,y,g_y,gamma,lambda,max_iter = max_iter)
# get per iteration result until max_iter for ISTA/FISTA/AMP
FISTA_list = list(grad=ans$grad[1:500],time=ans$time[1:500])
ISTA_list = ISTA_SGL(A,y,g_y,gamma,lambda,max_iter = 3000)
AMP_list = AMP_SGL_fromAlpha_list(A,y,g_y,gamma,alpha,X0,max_iter = 200)

# calculate per iteration MSE
mse_ISTA_list =mse_FISTA_list =mse_AMP_list = vector()

for(i in 1:200){ mse_AMP_list[i] =  mean((AMP_list$grad[[i]] - ans$grad[[max_iter]])^2 )}
for(i in 1:500){ mse_FISTA_list[i] =  mean((FISTA_list$grad[[i]] - ans$grad[[max_iter]])^2 )}
for(i in 1:3000){ mse_ISTA_list[i] =  mean((ISTA_list$grad[[i]] - ans$grad[[max_iter]])^2 )}

for(i in 2:5){
  print(paste0('ISTA error 1e-',i,' is ',which(mse_ISTA_list<10^(-i))[1]))
  print(paste0('FISTA error 1e-',i,' is ',which(mse_FISTA_list<10^(-i))[1]))
  print(paste0('AMP error 1e-',i,' is ',which(mse_AMP_list<10^(-i))[1]))
}


#####################
######## log loss plot comparing ISTA/FISTA/blockwise
iter = 2000
max_iter = 2000 # to get true minimizer

sss=Sys.time() # 6min
# run a full result as if true minimizer by FISTA (ISTA for long time also works)
ans = FISTA_SGL(A,y,g_y,gamma,lambda,max_iter = max_iter)
# get per iteration result until max_iter for ISTA/FISTA only
FISTA_list = list(grad=ans$grad[1:iter],time=ans$time[1:iter])
ISTA_list = ISTA_SGL(A,y,g_y,gamma,lambda,max_iter = iter)
Sys.time()-sss

# calculate per iteration MSE/cost for ISTA/FISTA only
mse_ISTA_list =mse_FISTA_list = rep(0,iter)
cost_ISTA_list =cost_FISTA_list = rep(0,iter)
for(i in 1:iter){
  mse_ISTA_list[i] =  mean((ISTA_list$grad[[i]] - ans$grad[[max_iter]])^2 )
  mse_FISTA_list[i] = mean((FISTA_list$grad[[i]] - ans$grad[[max_iter]])^2 )
  cost_ISTA_list[i] =  full_cost(y,A,ISTA_list$grad[[i]],gamma,lambda,g_y)
  cost_FISTA_list[i] = full_cost(y,A,FISTA_list$grad[[i]],gamma,lambda,g_y)
}


# get per iteration result until max_iter for subgradient
iters_try=c(1:5,seq(6,10,2),15,20,30,40,50,75,100,125,150)
cost_Subgradient_list=time_Subgradient=mse_Subgradient_list=iters_try*0
# total 1 hour
for (i in 1:length(iters_try)){
  time_start=Sys.time()
  subgrad=SGL(data=list(x=A, y=y),index=g_y,lambdas = lambda/n,standardize = FALSE,alpha=gamma,maxit = iters_try[i],thresh=1e-5)$beta
  mse_Subgradient_list[i]=mean((subgrad-ans$grad[[max_iter]])^2)
  cost_Subgradient_list[i]=full_cost(y,A,subgrad,gamma,lambda,g_y)
  print(Sys.time()-time_start)
  time_Subgradient[i]=difftime(Sys.time(),time_start,units='secs')
  print(c(iters_try[i],mse_Subgradient_list[i],cost_Subgradient_list[i]))
}

mse_Subgradient_list=c(mean(ans$grad[[max_iter]]^2),mse_Subgradient_list)

### Three wall-clock-time-related plots
# log error vs log time
par(mar=c(5,5,1,1))
final_cost=full_cost(y,A,ans$grad[[max_iter]],gamma,lambda,g_y)
plot(log(ISTA_list$time),log(cost_ISTA_list-final_cost),lwd = 3, type = 'l', col = 'blue',
     xlim = c(-2,7),ylim=c(-7,8),lty = 3,xlab = 'Log time in sec',ylab = 'Log error',xaxs='i',
     cex.lab = 2, cex.axis = 1.5)
lines(log(FISTA_list$time),log(cost_FISTA_list-final_cost), lwd = 3, type = 'l', col = 'green',lty = 1)
lines(log(time_Subgradient),log(cost_Subgradient_list-final_cost),lwd = 3 ,type = 'l', col = 'black',lty = 4)
legend('bottomleft', legend=c("FISTA","ISTA",'Blockwise'),
       col=c("green","blue",'black'), lty=c(1,3,4),lwd=3,cex=1.5)

# time vs iter
plot(1:iter,FISTA_list$time, lwd = 3, type = 'l', col = 'green',lty = 1,
     xlim=c(1,100),ylim=c(0,500),xlab = 'Iteration',ylab = 'Time in sec',
     xaxs='i',cex.lab = 2, cex.axis = 1.5)
lines(1:iter,ISTA_list$time,lwd = 2, type = 'l', col = 'blue',lty = 3)
lines(iters_try,time_Subgradient,lwd = 3 ,type = 'l', col = 'black',lty = 4)
legend('topleft', legend=c("FISTA","ISTA",'Blockwise'),
       col=c("green","blue",'black'), lty=c(1,3,4),lwd=3,cex=1.5)

# log time vs log iter
plot(log(1:iter),log(FISTA_list$time), lwd = 3, type = 'l', col = 'green',lty = 1,
     xlim=c(0,7),ylim=c(-2,7),xlab = 'Log iteration',ylab = 'Log time in sec',
     xaxs='i',cex.lab = 2, cex.axis = 1.5)
lines(log(1:iter),log(ISTA_list$time),lwd = 3, type = 'l', col = 'blue',lty = 3)
lines(log(iters_try),log(time_Subgradient),lwd = 3 ,type = 'l', col = 'black',lty = 4)
legend('topleft', legend=c("FISTA","ISTA",'Blockwise'),
       col=c("green","blue",'black'), lty=c(1,3,4),lwd=3,cex=1.5)


####################### reset settings, showing state evolution characterizes solution
start=Sys.time() # 15min
set.seed(1)
iter = 100
gamma = 0.5
n = 1000 # number of observations
N = 4000  # number of features
delta = n / N  # ratio
g_y=c(rep(1,N*0.5),rep(2,N*0.5)) # perfect group information

X0 <- c(rep(1,N*0.5), rep(0, N * 0.5))
A <- matrix( rnorm(n*N,mean=0,sd=sqrt(1/n)), n, N)   # Gaussain Matrix with mean 0 variance 1/n
sigma <- 1
y = A %*% X0+sigma*rnorm(n)

## modify alpha_to_lambda for speed
# simultaneously output tau and lambda instead use separate functions
alpha_to_lambda2=function(alpha,gamma,group,prior,delta,max_iter=100,second_iter=100,sigma=0){
  tau=alpha_to_tau(alpha,gamma,group,prior,delta,max_iter=max_iter,sigma=sigma)
  
  p=length(prior)
  E=0 # expectation term in calibration
  for (t in 1:second_iter){
    prox_prime=eta_prime(prior+tau*rnorm(p),gamma,alpha*tau,group)
    E=E+mean(prox_prime)/delta/second_iter
  }
  return(c((1-E)*alpha*tau,tau))
}

# try different alpha (equivalently lambda) and check whether empirical and estimated MSE similar
MSE_record=c()
seMSE_record=c()
alpha_seq=c(seq(0.96,1.44,0.02))
for(alpha in alpha_seq){
  result=alpha_to_lambda2(alpha,gamma,g_y,X0,delta,sigma=sigma,max_iter = 200,second_iter = 1000)
  lambda=result[1]
  tau=result[2]
  print(lambda)
  grad=FISTA_SGL(A,y,g_y,gamma,lambda,max_iter=iter)$grad[[iter]]
  se=eta(X0+tau*rnorm(N),gamma,alpha*tau,g_y)
  MSE_record=c(MSE_record,mean((grad-X0)^2))
  seMSE_record=c(seMSE_record,mean((se-X0)^2))
  print(Sys.time()-start)
}
# MSE plot: empirical and estimated agree
par(mar=c(4,0,1,1),pty="s")
plot(MSE_record,seMSE_record,xlim=c(0.38,0.52),ylim=c(0.38,0.52),xlab='Empirical MSE',
     ylab='Estimated MSE',cex.lab=2.5,cex.axis=1.5,xaxs='i',yaxs='i',cex=2,pch=20)
abline(0,1,lwd=3,col=8)

# plots of quantiles between solution and characterization for fixed alpha (equivalently lambda)
par(mar=c(4,0,1,1),pty="s")
set.seed(2)
alpha=1.13
result=alpha_to_lambda2(alpha,gamma,g_y,X0,delta,sigma=sigma,max_iter = 500,second_iter = 1000)
lambda=result[1] #lam=1
tau=result[2]
print(paste0('Under current prior, for alpha= ',alpha,'; tau= ',tau,', lambda= ',lambda))

# get true minimizer
grad=FISTA_SGL(A,y,g_y,gamma,lambda,max_iter=100)$grad[[100]]
# get charaterization
se=eta(X0+tau*rnorm(N),gamma,alpha*tau,g_y)

plot(sort(grad),sort(se),xlim=c(-1.1,1.5),ylim=c(-1.1,1.5),xlab=TeX('$\\;\\hat{\\beta}$'),
     ylab=TeX('$\\eta_{\\gamma}(\\Pi{+}\\tau Z;\\alpha\\tau)$'),cex.lab=2,xaxs='i',yaxs='i',
     cex=2,pch=20,xaxt='n',cex.axis=1.5)
axis(1, at=seq(-1,1.5,0.5), labels= seq(-1,1.5,0.5),cex.axis=1.5)  
abline(0,1,lwd=3,col=8)
