### This code defines state evolution and calibration for AMP

source("SGL_proximal.R")

# compute E[eta_{soft}^2(Z;t)]
ETA2 <- function(t){
  f <- (1 + t^2)* pnorm(-t) - t * dnorm(t)
  return(2*f+1e-10)
}

################################## alpha range
judge=function(t,gamma,delta){return(ETA2(gamma*t)-2*(1-gamma)*t*sqrt(ETA2(gamma*t))+((1-gamma)*t)^2-delta)}

Arange= function(gamma,delta,tol=1e-8,length_precision=100000){
  amin=0
  
  if(gamma==1){
    # find a point which is in A
    bisecmax=1
    while(judge(bisecmax,gamma,delta)>0){bisecmax=2*bisecmax}
    # bisection
    while ((bisecmax-amin)>tol){
      bisecmid=(bisecmax+amin)/2
      if (judge(bisecmid,gamma,delta)>0){amin=bisecmid}
      else{bisecmax=bisecmid}
    }
    amax=NA
  }else{
    # find a point which is not in A
    amax=1
    while(judge(amax,gamma,delta)<0){amax=2*amax}
    range=seq(amin,amax,length.out = length_precision)
    negative=which(judge(range,gamma,delta)<0)
    amin=range[min(negative)]
    amax=range[max(negative)]
  }
  return(list(amin=amin,amax=amax))
}


################################## state evolution equation
F= function(tau,alpha,gamma,group,prior,delta,iteration=100,sigma=0){
  p=length(prior)
  result=0
  for (i in 1:iteration){
    result=result+mean((eta(prior+tau*rnorm(p),gamma,lambda=alpha*tau,group)-prior)^2)/iteration
  }
  return(sigma^2+result/delta)
}


################################## alpha to tau calibration
alpha_to_tau=function(alpha,gamma,group,prior,delta,max_iter=100,sigma=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha,gamma,group,prior,delta,sigma=sigma))
    record_tau[t]=tau #record each tau
  }
  return(mean(record_tau))
}

################################## alpha to lambda calibration
alpha_to_lambda=function(alpha,gamma,group,prior,delta,max_iter=100,second_iter=100,sigma=0){
  tau=alpha_to_tau(alpha,gamma,group,prior,delta,max_iter=max_iter,sigma=sigma)
  
  p=length(prior)
  E=0 # expectation term in calibration
  for (t in 1:second_iter){
    prox_prime=eta_prime(prior+tau*rnorm(p),gamma,alpha*tau,group)
    E=E+mean(prox_prime)/delta/second_iter
  }
  return((1-E)*alpha*tau)
}
################################## lambda to alpha calibration (implicitly)
lambda_to_alpha=function(lambda,gamma,group,prior,delta,tol=1e-2,sigma=0){
  
  temp=Arange(gamma,delta)
  alpha1=temp$amin
  alpha2=temp$amax
  
  ### bisection to find the alpha_seq which is parallel to lambda_seq
  while ((alpha2-alpha1)>tol){
    middle_alpha=(alpha1+alpha2)/2
    middle_lambda=alpha_to_lambda(middle_alpha,gamma=gamma,group=group,prior=prior,delta=delta,sigma=sigma)
    if (middle_lambda>lambda){
      alpha2=middle_alpha
    }else if (middle_lambda<lambda){
      alpha1=middle_alpha
    }else{
      break
    }
  }
  return(middle_alpha)
}
