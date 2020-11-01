### This code contains ingredients of SGL proximal operation

# objective function of SGL proximal operator for single group
cost <- function(s,beta, gamma, lambda){
  G = (1/2) * sum((s - beta)^2)
  ans =  (1 - gamma) * lambda * sqrt(length(beta))*sqrt(sum(beta^2))
  return(G + ans + gamma * lambda * sum(abs(beta)))
}

# full objective function of SGL proximal operator
full_cost <- function(y,X,beta, gamma, lambda,g_y){
  G = (1/2) * sum((y - X%*%beta)^2)
  for(i in unique(g_y)){G = G+(1 - gamma) * lambda * sqrt(sum(g_y==i))*sqrt(sum(beta[g_y==i]^2))}
  return(G + gamma * lambda * sum(abs(beta)))
}

# soft-thresholding
eta_s <- function(x, theta){
  y <- rep(0, length(x))
  for(i in 1:length(x)){if(abs(x[i])> theta){y[i] = x[i] - sign(x[i])*theta}}
  return(y)
}

# SGL proximal operator
eta <- function(y, gamma, lambda, g_y){
  ans <- rep(0, length(y))
  for(i in unique(g_y)){
    X=eta_s(y[g_y == i],gamma*lambda)
    if(mean(X^2) != 0 ){
      ans[g_y == i] <- X* (1 - (1 - gamma)*lambda/sqrt(mean(X^2)))
    }
    if(cost(y[g_y == i],  ans[g_y == i] ,gamma, lambda) >= 1/2 * sum(y[g_y == i]^2)  ){
      ans[g_y == i] = 0
    }
  }
  return(ans)
}

# SGL proximal operator derivative
eta_prime = function(y, gamma, lambda, g_y ){
  bns <- rep(0, length(y))
  ans <- rep(0, length(y))
  for(i in unique(g_y)){
    X=eta_s(y[g_y == i],gamma*lambda)
    if(mean(X^2) != 0){
      bns[g_y == i] <- X* (1 - (1 - gamma)*lambda/sqrt(mean(X^2)))
      ans[g_y == i] = (X!=0) * (1 - (1 - gamma) * lambda/ sqrt(mean(X^2))* (1 - X^2/ sum(X^2 )))
    }
    if(cost(y[g_y == i],  bns[g_y == i] ,gamma, lambda) >= 1/2 * sum(y[g_y == i]^2)  ){
      ans[g_y == i] = 0
    }
  }
  return(ans)
}