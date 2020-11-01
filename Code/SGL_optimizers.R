### This code defines ISTA/FISTA/AMP for SGL, especially ISTA/FISTA take in lambda as penalty,
### where AMP takes in alpha, which corresponds to lambda via calibration.
### The output is per iteration iterate and time consumed.

# ISTA for SGL, i.e. proximal gradient descent
ISTA_SGL = function(A ,y, g_y , gamma , lambda,initial=rep(0,ncol(A)), max_iter = 1000){
  x = initial
  s = 1/(2*norm(t(A) %*% A,type = "F"))
  time_list= rep(0,max_iter)
  x_list =list()
  time_start=Sys.time()
  for(i in 1:max_iter){
    x_list[[i]] = x
    x = eta(x - s * t(A) %*% (A %*% x - y) ,gamma,lambda * s ,g_y)
    time_list[i]=difftime(Sys.time(),time_start,units='secs')
  }
  return(list(grad=x_list,time=time_list))
}

# FISTA for SGL, i.e. Nesterov-accelarated proximal gradient descent
FISTA_SGL = function(A,y, g_y ,  gamma , lambda, initial=rep(0,ncol(A)),max_iter = 1000){
  x = initial
  k = initial
  a = 1
  s = 1/(2*norm(t(A) %*% A,type = "F"))
  time_list=rep(0,max_iter)
  x_list = list()
  time_start=Sys.time()
  for(i in 1:max_iter){
    x_list[[i]] = x
    a_t = (1 + sqrt(1 + 4 * a^2))/2
    b = (1 - a)/a_t
    k_t = eta(x - s * t(A) %*% (A %*% x - y) ,gamma,lambda * s ,g_y)
    x = (1 - b) * k_t + b * k
    k = k_t
    a = a_t
    time_list[i]=difftime(Sys.time(),time_start,units='secs')
  }
  return(list(grad=x_list,time=time_list))
}

# AMP for SGL
AMP_SGL_fromAlpha_list = function(A, y,g_y, gamma,alpha, prior, sigma = 0,initial=rep(0,ncol(A)), max_iter = 1000){
  n = nrow(A)
  N = ncol(A)
  delta = n / N
  X0 = prior
  x <- initial # zero array
  z <- rep(0,n) # zero array
  tau<- sqrt(sigma^2 + mean(X0^2)/delta)
  time_list=rep(0,max_iter)
  x_list = list()
  time_start=Sys.time()
  for(i in 1:max_iter){
    x_list[[i]] = x
    obs=t(A) %*% z + x
    x = eta(obs, gamma ,alpha * tau,g_y)
    aver = eta_prime(obs, gamma ,alpha * tau, g_y )
    z <- y - A %*% x + (1 /delta) *z*mean(aver)
    second_moment <- (eta(X0 + tau * rnorm(N),gamma,alpha * tau, g_y ) - X0)^2
    tau <- sqrt(sigma^2 + (1/delta)*mean(second_moment))
    time_list[i]=difftime(Sys.time(),time_start,units='secs')
  }
  return(list(grad=x_list,time=time_list))
}