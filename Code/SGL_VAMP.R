### This code develops SGL VAMP and compare with FISTA/ISTA

library(MASS)
source('SGL_proximal.R')
source('SGL_optimizers.R')

# full cost in single group case
cost_compare = function(A,x,y,gamma,lambda){
  G = (1/2) * sum((A %*% x - y)^2)
  ans =  (1 - gamma) * lambda * sqrt(length(x))*sqrt(sum(x^2))
  return(G + ans + gamma * lambda * sum(abs(x)))
}

#### SGL VAMP on Boston dataset #################
set.seed(0)
Boston <- na.omit(Boston) 
A=as.matrix(Boston[1:506,1:13])
y=as.vector(Boston[1:506,14])
p=ncol(A);n=nrow(A);delta=n/p
lambda = 5
g_y = rep(1,p)
gamma = 0.5




f_prox=function(A, B,XtX,Xty,iter=5){
  ### A is rho in paper, B is u^t, X is A
  R=chol(XtX + A * diag(nrow(XtX)))
  mean = backsolve(R,forwardsolve(t(R),Xty + B))
  var = mean(diag(chol2inv(R)))
  return(list("mean" = mean, "var" = var))
}

######### z^t and sigma_z^t
g_prox=function(A, B,alpha,gamma,g_y){
  # A is 1/[sigma_x/(1-sigma_x*rho)] in z^t, alpha is lambda in z^t, B is x^t/sigma_x-u^t
  prox_solution=eta(B,gamma,alpha,g_y)
  mean = prox_solution / A
  var = eta_prime(B/A, gamma, alpha/A,g_y) / A
  return(list("mean" = mean, "var" = var))
}




# Initialization, note VAMP converges in one iteration
VAMP_iter = 1
x=rep(0,p);z=y;
rho=1;u=rep(0,p);damp=0.1
XtX=t(A)%*%A
Xty=t(A)%*%y
time_list_VAMP = rep(0,3)
loss_VAMP = rep(0,3)

### VAMP recursion
s=Sys.time()
for (i in 1:VAMP_iter){
  x_relavant=f_prox(rho,u,XtX,Xty)
  x=x_relavant$mean
  sigmax=x_relavant$var
  z_relavant=g_prox(1 / sigmax - rho, x/sigmax - u,1.5,0.5,g_y)
  z=z_relavant$mean
  sigmaz=z_relavant$var
  rho = rho + (1 - damp) * (1 /sigmaz - 1 /sigmax)
  u = u + (1 - damp) * (z /sigmaz - x/sigmax)
}
# save time and loss
time_list_VAMP[2] = Sys.time()-s
time_list_VAMP[3] = 100
loss_VAMP[1] = cost_compare(A,rep(0,p),y,gamma,lambda )
loss_VAMP[2] = cost_compare(A,x,y,gamma,lambda )
loss_VAMP[3] = cost_compare(A,x,y,gamma,lambda )
print(Sys.time()-s)

### FISTA 
FISTA_iter = 10000
list_FISTA = FISTA_SGL(A,y,g_y,gamma,lambda,initial=rep(0,ncol(A)), max_iter = FISTA_iter)
minimizer_FISTA = list_FISTA$grad[[FISTA_iter]]
time_list_FISTA = list_FISTA$time
for(i in 2:FISTA_iter){
  time_list_FISTA[[i]] = time_list_FISTA[[i-1]] + time_list_FISTA[[i]]
}
loss_FISTA = rep(0,FISTA_iter)
for(i in 1:FISTA_iter){
  loss_FISTA[i] = cost_compare(A,list_FISTA$grad[[i]],y,gamma,lambda)
}

### ISTA
ISTA_iter = 10000
list_ISTA = ISTA_SGL(A,y,g_y,gamma,lambda,initial=rep(0,ncol(A)), max_iter = ISTA_iter)
minimizer_ISTA = list_ISTA$grad[[ISTA_iter]]
time_list_ISTA = list_ISTA$time
for(i in 2:FISTA_iter){
  time_list_ISTA[[i]] = time_list_ISTA[[i-1]] + time_list_ISTA[[i]]
}
loss_ISTA = rep(0,ISTA_iter)
for(i in 1:ISTA_iter){
  loss_ISTA[i] = cost_compare(A,list_ISTA$grad[[i]],y,gamma,lambda)
}


### Plot loss against time
par(mar=c(4,5,1,1))
plot(time_list_VAMP,loss_VAMP,lwd = 3, type = 'l', col = 'red',xlim = c(0,2),lty = 1,
     xlab = 'Time in sec',ylab = 'Loss',xaxs='i',cex.lab = 1.5, cex.axis = 1.5)
lines(time_list_FISTA,loss_FISTA, lwd = 3, type = 'l', col = 'green',lty = 2)
lines(time_list_ISTA,loss_ISTA,lwd = 3 ,type = 'l', col = 'blue',lty = 3)
legend('topright', legend=c("VAMP", "FISTA","ISTA"),
       col=c("red", "green","blue"), lty=1:3,lwd=2,cex=1.5)



### Loss for FISTA/ISTA after 10000 iterations
print(cost_compare(A,x,y,gamma,lambda))
print(cost_compare(A,minimizer_FISTA,y,gamma,lambda ))
print(cost_compare(A,minimizer_ISTA,y,gamma,lambda ))









