conjVAR <- function(Y,X,var.ar1,p=12,ihor=36){
  # some dimensions
  cons <- TRUE
  k <- ncol(X)
  n <- ncol(Y)
  T <- nrow(Y)
  
  # OLS estimate
  A_OLS <- solve(crossprod(X))%*%crossprod(X,Y)
  Yfit <- X %*% A_OLS # fitted values
  epsilon <- Y - X %*% A_OLS
  
  Sigma_OLS <- crossprod(epsilon)/(T-k)
  
  delta <- var.ar1
  theta <- 0.1
  theta_alpha <- 0.1
  omega <- 1
  
  # univariate models
  s_arp <- rep(NA,n)
  for(nn in 1:n){
    s_arp[nn] <- sqrt(arima(Y[,nn],order=c(p,0,0))$sigma2)
  }
  
  # Minnesota dummies
  X_ <- cbind(rbind(kronecker(diag(1:p),diag(s_arp))/theta,matrix(0,n+1,n*p)),
              c(rep(0,n*(p+1)),theta_alpha))
  Y_ <- rbind(diag(delta*s_arp)/theta,
              matrix(0,n*(p-1),n),
              diag(s_arp),
              rep(0,n))
  
  V_ <- solve(crossprod(X_))
  A_ <- V_ %*% crossprod(X_,Y_)
  
  # augmented data
  YY <- rbind(Y,Y_)
  XX <- rbind(X,X_)
  
  S0 <- crossprod(Y_ - X_%*%A_) # prior scaling x prior DoF
  s0 <- nrow(Y_)*omega + n + 2  # prior DoF (n + 2 ensure that prior is proper)
  
  # Compute posterior moments
  V1 <- solve(crossprod(XX))
  V1.chol <- t(chol(V1))  # cholesky of V1
  
  A1 <- V1 %*% crossprod(XX,YY)
  a1 <- c(A1) # vectorize coefficients
  
  S1 <- crossprod(YY-XX%*%A1)
  s1 <- T + s0
  
  # initialize companion matrix
  M <- matrix(0,k-cons,k-cons)
  M[(n+1):nrow(M),1:(n*(p-1))] <- diag(n*(p-1))
  
  nsave <- 1000
  nburn <- 0 # no burn-in needed (conjugate)
  ntot <- nsave + nburn
  
  irf_store <- array(NA,dim=c(nsave,ihor,n))
  pb <- txtProgressBar(min = 1, max = ntot, style = 3)
  for(irep in 1:ntot){
    # sample from inverse wishart
    Sigmai <- matrix(rWishart(1,s1,solve(S1)),n,n)
    Sigma_draw <- solve(Sigmai)
    
    a_shks <- matrix(rnorm(k*n),k,n) #define a K x M matrix of standard normal errors
    a_shks <- as.vector(V1.chol %*% a_shks %*% chol(Sigma_draw))
    a_draw <- a1 + a_shks
    A_draw <- matrix(a_draw,k,n)

    if(irep > nburn){
      irf <- matrix(NA,ihor,ncol(X)-cons)
      M[1:n,] <- t(A_draw[-k,])
      delta0 <- c(t(chol(Sigma_draw))[,1],rep(0,n*(p-1)))
      for(hh in 1:ihor){
        if(hh == 1){
          irf[hh,] <- delta0
        }else{
          irf[hh,] <- M %*% irf[hh-1,]
        }
      }
      irf_store[irep-nburn,,] <- irf[,1:n]
    }
    setTxtProgressBar(pb, irep)
  }
  return(irf_store)
}