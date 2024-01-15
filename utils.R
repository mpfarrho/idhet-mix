# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 2 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

get_jointlik_lev <- function(ht_lags,ytt,svpara,gamma,beta,Omega,it,mean.exp=TRUE){
  require(mvnfast)
  
  ht_y <- ht_lags[,1]
  ht_x <- ht_lags[,2]
  
  yy <- c(ytt[it,],ht_y[it])
  if(mean.exp){
    mm <- c(gamma * exp(ht_y[it]),svpara[1] * ht_x[it])
  }else{
    mm <- c(gamma * ht_y[it],svpara[1] * ht_x[it])
  }
  vv_11 <- (exp(ht_y[it]) * (beta %*% t(beta))) + Omega
  vv_12 <- sqrt(svpara[2]) * svpara[3] * exp(ht_y[it]/2) * beta
  vv <- rbind(cbind(vv_11,vv_12),
              cbind(t(vv_12),svpara[2]))
  return(mvnfast::dmvn(X = yy, mu = mm, sigma = vv, log=TRUE))
}

# -----------------------------------------------------------------------------------------------
# estimate sparse finite mixture factor model
sfm.fm <- function(y,d,x=NA,
                   nburn=1000,nsave=1000,thinfac=1,
                   mix=TRUE,J=20,sample.e0=TRUE,
                   Q=2,sv=TRUE,svm=FALSE,svf=FALSE,lev=FALSE,
                   quiet=FALSE,checkmode=FALSE){
  require(Matrix)
  require(stochvol)
  require(TruncatedNormal)
  require(invgamma)
  require(mvnfast)
  if(svm) stop("Not available.")
  
  nthin <- round(nsave/thinfac)
  ntot <- nsave+nburn
  thin.set <- floor(seq(nburn+1,ntot,length.out=nthin))
  in.thin <- 0
  
  # split into non-ruling and ruling periods
  Y_NR <- y[d==0,]
  Y_R <- y[d==1,]
  X_NR <- x[d==0,,drop=FALSE]
  X_R <- x[d==1,,drop=FALSE]
  
  Y_RNR <- rbind(Y_R,Y_NR)
  X_RNR <- rbind(X_R,X_NR)

  # dimensions
  T <- nrow(y)
  T_R <- sum(d)
  T_NR <- T-T_R
  n <- ncol(y)
  
  # -----------------------------------------------------------------------------------------
  # priors
  q <- 1   # number of ruling factors
  Qf <- q + Q # total number of factors
  tau <- 1 # lambda_j ~ N(0,tau)
  
  norm.fac <- TRUE # normalize factor in each sweep
  id.signs <- TRUE # identify sign of always active factors
  mean.exp <- FALSE # volatility or log-volatility in the mean

  # theta_j ~ G^-1(c0,d0) 
  c0 <- 0.1
  d0 <- 0.1 
  
  # prior for the dirichlet weights, e0 ~ G(a_gam,J*a_gam)
  a_gam <- 10
  b_gam <- J*a_gam
  c_prop <- 1
  accept <- 0
  e0 <- 1/J # initial value
  
  # auxiliary objects for mixture
  sig2_draw <- rep(1,J)
  si_draw <- matrix(0,T,J)
  eta_draw <- rep(1,J)/J
  omega <- matrix(1/J,T,J)
  
  i_T <- matrix(1,T,1)
  i_n <- matrix(1,n,1)
  I_Tn <- matrix(1,T,n)
  I_J <- diag(J)
  I_T <- Matrix::Matrix(0,T,T,sparse=TRUE)
  I_TNR <- Matrix::Matrix(0,T_NR,T_NR,sparse=TRUE)
  diag(I_T) <- diag(I_TNR) <- 1
  
  for(tt in 1:T){
    si_draw[tt,] <- I_J[sample(1:J,1,prob=omega[tt,],replace=T),]
  }
  lik <- matrix(NA,T,J) # to hold the likelihoods over time
  
  # -----------------------------------------------------------------------------------------
  # initialize factors
  f_draw <- matrix(rnorm(T_R*q),T_R,q) # factor that is active on ruling dates
  F_tmp <- matrix(rnorm(T*Q),T,Q)     # factor that is always active
  F_draw <- cbind(matrix(c(as.numeric(f_draw),rep(0,T_NR)),T,q),F_tmp) # combine ruling factor and always active factor
  F_draw_t <- F_draw*0
  
  # initialize loadings
  lambda <- solve(crossprod(F_draw))%*%crossprod(F_draw,rbind(Y_R,Y_NR))
  ft_lambda <- F_draw %*% lambda
  
  if(sv) fHt <- matrix(rnorm(T_R,0,0.1),T_R,1) else fHt <- matrix(0,T_R,1)
  if(svf){
    svf_draw <- svf_latent <- svf_priors <- list()
    svf_para <- matrix(NA,Q,4)
    for(qq in 1:Q){
      svf_draw[[qq]] <- list(mu = 0, phi = 0.99, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      svf_latent[[qq]] <- rep(0,T)
      
      svf_priors[[qq]] <- specify_priors(
        mu = sv_normal(mean=0, sd=1), # identify scale of factors by normalizing unconditional variance to one
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.1)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
  }
  FHt <- matrix(0,T,Q) # sv paths of always on factor
  
  # priors for the log vola process
  a_phi <- 3
  b_phi <- 3
  a_sig <- 3
  b_sig <- 0.3
  a_rho <- 5
  b_rho <- 5
  
  sv_para <- sv_para_prop <- c(0.95,0.1,0)
  scale_sv_para <- c(0.01,0.1,0.01)
  acc_sv_para <- 0
  
  diagFhti <- Matrix(0,T_NR*Q,T_NR*Q,sparse=TRUE)
  diag(diagFhti) <- 1
  
  # check whether exogenous controls are present
  n_x <- ncol(X_RNR)
  x_delta <- 0
  
  # initialize exogenous parameters
  delta <- matrix(0,n,n_x) # coefficient on controls
  delta_V0 <- rep(1,n_x) # prior variance
  delta_V0i <- 1/delta_V0
  
  # -----------------------------------------------------------------------------------------
  # storage
  F_store <- array(NA,c(nthin,T,Qf))
  FHt_store <- array(NA,c(nthin,T,Q))
  svf_store <- array(NA,c(nthin,Q,3))
  
  f_store <- array(NA,c(nthin,T_R))
  fht_store <- array(NA,c(nthin,T_R))
  svp_store <- array(NA,c(nthin,3))
  
  lambda_store <- array(NA,c(nthin,n,Qf))
  irf_store <- array(NA,c(nthin,n,2))
  delta_store <- array(NA,c(nthin,n,n_x))
  
  si_store <- array(NA,c(nthin,T,J))
  si_NE_store <- array(NA,c(nthin,T))
  sig2_store <- array(NA,c(nthin,J))
  sig2t_store <- array(NA,c(nthin,T))
  nJ_store <- array(NA,c(nthin,1))
  
  # -----------------------------------------------------------------------------------------
  # start sampling loop
  if(!quiet) pb <- txtProgressBar(min = 1, max = ntot, style = 3)
  
  for(irep in 1:ntot){
    eps <- Y_RNR - x_delta - ft_lambda # residuals after controlling for everything else
    if(mix){
      # Step 1a: Sample clustered variances
      for(j in 1:J){
        id_j <- si_draw[,j]==1
        T_j <- sum(si_draw[,j])
        if(T_j>0){
          eps_j <- as.numeric(eps[id_j,])
          c0p <- c0 + NROW(eps_j)/2
          d0p <- as.numeric(d0 + crossprod(eps_j)/2)
          sig2_draw[j] <- 1/rgamma(1,c0p,d0p)
        }else{
          sig2_prior <- 1/rgamma(1,c0,d0)
          sig2_prior[sig2_prior>100] <- 100
          sig2_draw[j] <- sig2_prior
        }
      }
      
      # Step 1b: Sample state allocation
      for(j in 1:J) lik[,j] <- dnorm(eps,mean=0,sd=sqrt(sig2_draw[j]),log=TRUE) %*% i_n
      omega <- t(matrix(eta_draw,J,T))*exp(lik)
      omega[omega<1e-80] <- 1e-80
      omega <- omega/apply(omega,1,sum) # normalize for probability
      
      si_draw[1:T_R,] <- 0
      for(tt in (T_R+1):T){
        si_draw[tt,] <- I_J[sample(1:J,1,prob=omega[tt,],replace=TRUE),]
      }
      nk <- apply(si_draw[(T_R+1):T,],2,sum)
      eta_draw <- as.numeric(bayesm::rdirichlet(e0+nk))
      
      sl.NE <- which(sig2_draw==min(sig2_draw[nk>0]))
      si_draw[1:T_R,sl.NE] <- 1 # assign rulings to minimum-variance regime
      sig2_t <- apply(t(sig2_draw*t(si_draw)),1,sum)
      
      if(sample.e0){
        e0_j <- e0 # current value
        le0_p <- log(e0_j) + rnorm(1,0,c_prop) # log-proposal
        e0_p <- exp(le0_p) # proposal
        
        # likelihood ratio MH
        lp1 <- (e0_p - e0_j) * sum(log(omega[-c(1:T_R),])) + lgamma(J * e0_p) - lgamma(J * e0_j) - J *
          (lgamma(e0_p) - lgamma(e0_j)) + (a_gam - 1) * (log(e0_p) - log(e0_j)) - (e0_p - e0_j) *
          b_gam + log(e0_p) - log(e0_j)
        p1 <- min(exp(lp1), 1)
        p2 <- runif(1)
        
        if (p2 <= p1) {
          e0_j <- e0_p
          e0 <- e0_j
          accept <- accept+1
        }
        if (irep<(0.5*nburn)){
          if (accept/irep < 0.2) c_prop <- 0.99*c_prop
          if (accept/irep > 0.4) c_prop <- 1.01*c_prop
        }
      }
    }else{
      eps_j <- as.numeric(eps)
      c0p <- c0 + NROW(eps_j)/2
      d0p <- as.numeric(d0 + crossprod(eps_j)/2)
      sig2_draw[] <- 1/rgamma(1,c0p,d0p)
      sig2_t <- rep(sig2_draw[1],T)
      
      si_draw[] <- 0
      si_draw[,1] <- 1
    }
    
    # Step 2: In case there are exogenous controls, sample the parameters
    normalizer <- 1/sqrt(sig2_t)
    YY <- (Y_RNR - ft_lambda)*normalizer
    
    # loop through dependent variables
    for(nn in 1:n){
      yy_nn <- YY[,nn]
      xx_nn <- X_RNR*normalizer
      
      Vp <- solve(crossprod(xx_nn) + delta_V0i * diag(n_x))
      bp <- Vp %*% crossprod(xx_nn,yy_nn)
      
      delta[nn,] <- bp + t(chol(Vp)) %*% rnorm(n_x,0,1)
    }
    x_delta <- X_RNR %*% t(delta) # fitted values
    
    # Step 3: Sample the factor loadings
    YY <- (Y_RNR - x_delta)/sqrt(sig2_t)
    XX <- F_draw/sqrt(sig2_t)
    
    Vp <- solve(crossprod(XX) + diag(Qf)/tau)
    Vp.chol <- t(chol(Vp))
    for(nn in 1:n){
      bp <- Vp %*% crossprod(XX,YY[,nn])
      lambda[,nn] <- bp + Vp.chol%*%rnorm(Qf,0,1)
    }
    
    # Step 4: Sample the factors
    YY <- Y_RNR - x_delta
    normalizer <- 1/sqrt(sig2_t)
    
    # Step 4a: Factor matrix during ruling events (all factors active)
    for(tt in 1:T_R){
      yt <- YY[tt,]*normalizer[tt]
      Lt <- t(lambda)*normalizer[tt]
      
      ft_sig2 <- solve(crossprod(Lt) + diag(Qf)/c(exp(fHt[tt,]),exp(FHt[tt,])))
      ft_mu <- ft_sig2%*%crossprod(Lt,yt)
      ft_sig <- t(chol(ft_sig2))
      
      f_tt <- F_draw[tt,] <- ft_mu + ft_sig %*% rnorm(Qf)
      f_draw[tt,] <- f_tt[1:q]
    }
    if(norm.fac) f_draw <- (f_draw - mean(f_draw))/sd(f_draw)
    
    # Step 4b: Sample always active factor (only always active factors appear)
    YYNR <- c(t(YY[-c(1:T_R),]*normalizer[-c(1:T_R)]))
    XX <- Matrix::kronecker(I_TNR*normalizer[-c(1:T_R)],t(lambda[-c(1:q),]))
    
    VF_post <- Matrix::solve(diagFhti + Matrix::crossprod(XX))
    F_post <- VF_post %*% (Matrix::crossprod(XX,YYNR))
    f_draw_TNR <- F_post + t(Matrix::chol(VF_post))%*%rnorm(T_NR*Q)
    F_draw[(T_R+1):T,(q+1):Qf] <- matrix(f_draw_TNR,ncol=Q,byrow=TRUE)
    if(svf & norm.fac) F_draw[,(q+1):Qf] <- apply(F_draw[,(q+1):Qf],2,function(x) (x-mean(x)) / sd(x))
    if(norm.fac) F_draw[1:T_R,1] <- (F_draw[1:T_R,1] - mean(F_draw[1:T_R,1]))/sd(F_draw[1:T_R,1])
    ft_lambda <- F_draw %*% lambda # compute new fitted values for the factors
    
    # Step 5: Sample the SV process
    if(sv){
      eps <- Y_R - x_delta[1:T_R,]
      gamma <- rep(0,n)
      sv_draw <- SV_JPR(mu_h0=0,sig_h0=0.2,c_h=0,sv_para=sv_para,
                        Em=eps,hv=fHt,Hv=FHt[1:T_R,],L=t(lambda),gamma=gamma,Sigma=min(sig2_draw),
                        Qf=Qf,Q=Q,TT=T_R,n=n,scale_h=0.1,mean.exp=mean.exp)
      fHt[] <- sv_draw$hv
      hraw <- c(sv_draw$h0,sv_draw$hv)
      
      # sampling the SV parameters (potentially with leverage)
      ht_lags <- embed(hraw,2)
      lik_t <- matrix(NA,T_R,2)

      sv_para_prop[1] <- TruncatedNormal::rtmvnorm(1,sv_para[1],sigma=scale_sv_para[1],lb=-1,ub=1)
      sv_para_prop[2] <- exp(rnorm(1,0,sqrt(scale_sv_para[2])))*sv_para[2]
      if(lev){
        sv_para_prop[3] <- TruncatedNormal::rtmvnorm(1,sv_para[3],sigma=scale_sv_para[3],lb=-1,ub=1)
      }else{
        sv_para_prop[3] <- 0
      }

      # evaluate likelihood
      Sigma <- (diag(n) * min(sig2_draw))
      eps <- Y_R - x_delta[1:T_R,]
      gamma <- rep(0,n)
      
      for(tt in 1:T_R){
        Omega <- Sigma + (t(lambda[-1,,drop=FALSE]) %*% (diag(Q)*exp(FHt[tt,])) %*% lambda[-1,,drop=FALSE])
        lik_t[tt,1] <- get_jointlik_lev(ht_lags=ht_lags, ytt = eps, svpara = sv_para_prop, gamma=gamma, beta=t(lambda[1,,drop=FALSE]), Omega=Omega, it=tt, mean.exp=mean.exp)
        lik_t[tt,2] <- get_jointlik_lev(ht_lags=ht_lags, ytt = eps, svpara = sv_para, gamma=gamma, beta=t(lambda[1,,drop=FALSE]), Omega=Omega, it=tt, mean.exp=mean.exp)
      }
      lik_sv <- apply(lik_t,2,sum)

      if(lev){
        pri <- c(dbeta((sv_para_prop[1]+1)/2,a_phi,b_phi,log=TRUE) + dinvgamma(sv_para_prop[2],a_sig,b_sig,log=TRUE) + dbeta((sv_para_prop[3]+1)/2,a_rho,b_rho,log=TRUE),
                 dbeta((sv_para[1]+1)/2,a_phi,b_phi,log=TRUE) +      dinvgamma(sv_para[2],a_sig,b_sig,log=TRUE) +      dbeta((sv_para[3]+1)/2,a_rho,b_rho,log=TRUE))
        cor <- c(  TruncatedNormal::dtmvnorm(sv_para[1],mu = sv_para_prop[1],sigma = scale_sv_para[1],lb=-1,ub=1,B=100,log=TRUE) +
                   dlnorm(sv_para[2],log(sv_para_prop[2]),sqrt(scale_sv_para[2]),log=TRUE) +
                   TruncatedNormal::dtmvnorm(sv_para[3],mu = sv_para_prop[3],sigma = scale_sv_para[3],lb=-1,ub=1,B=100,log=TRUE),
                   TruncatedNormal::dtmvnorm(sv_para_prop[1],mu = sv_para[1],sigma = scale_sv_para[1],lb=-1,ub=1,B=100,log=TRUE) +
                   dlnorm(sv_para_prop[2],log(sv_para[2]),sqrt(scale_sv_para[2]),log=TRUE) +
                   TruncatedNormal::dtmvnorm(sv_para_prop[3],mu = sv_para[3],sigma = scale_sv_para[3],lb=-1,ub=1,B=100,log=TRUE))
      }else{
        pri <- c(dbeta((sv_para_prop[1]+1)/2,a_phi,b_phi,log=TRUE) + dinvgamma(sv_para_prop[2],a_sig,b_sig,log=TRUE),
                 dbeta((sv_para[1]+1)/2,a_phi,b_phi,log=TRUE) +      dinvgamma(sv_para[2],a_sig,b_sig,log=TRUE))
        cor <- c(TruncatedNormal::dtmvnorm(sv_para[1],mu = sv_para_prop[1],sigma = scale_sv_para[1],lb=-1,ub=1,B=100,log=TRUE) +
                 dlnorm(sv_para[2],log(sv_para_prop[2]),sqrt(scale_sv_para[2]),log=TRUE),
                 TruncatedNormal::dtmvnorm(sv_para_prop[1],mu = sv_para[1],sigma = scale_sv_para[1],lb=-1,ub=1,B=100,log=TRUE) +
                 dlnorm(sv_para_prop[2],log(sv_para[2]),sqrt(scale_sv_para[2]),log=TRUE))
      }

      # accept/reject and tuning of proposals
      post <- lik_sv + pri + cor
      u <- log(runif(1,0,1))<(post[1]-post[2])
      sv_para <- u*sv_para_prop + (1-u)*sv_para
      acc_sv_para <- acc_sv_para + u

      if(irep < (nburn/2)){
        if((acc_sv_para/irep) < 0.2){
          scale_sv_para <- scale_sv_para*0.99
        }else if((acc_sv_para/irep) > 0.4){
          scale_sv_para <- scale_sv_para*1.01
        }
      }
    }
    
    # sample sv of always-active factors
    if(svf){
      F_draw_t[d==1,] <- F_draw[1:T_R,]
      F_draw_t[d==0,] <- F_draw[(T_R+1):T,]
      
      for(qq in (q+1):Qf){
        svfdraw_qq <- svsample_fast_cpp(F_draw_t[,qq],startpara=svf_draw[[qq-q]],startlatent=svf_latent[[qq-q]],priorspec=svf_priors[[qq-q]])
        svf_draw[[qq-q]][c("mu", "phi", "sigma")] <- as.list(svfdraw_qq$para[, c("mu", "phi", "sigma")])
        svf_para[qq-q,] <- svfdraw_qq$para[,c("mu","phi","sigma","nu")]
        FHt_tmp <- svf_latent[[qq-q]] <- svfdraw_qq$latent
        
        FHt[1:T_R,qq-q] <- FHt_tmp[d==1]
        FHt[(T_R+1):T,qq-q] <- FHt_tmp[d==0]
      }
      FHt[FHt>3] <- 3
      FHt[FHt<(-3)] <- (-3)
      diag(diagFhti) <- c(t(exp(-FHt[(T_R+1):T,])))
    }else{
      FHt[] <- 0
      diag(diagFhti) <- 1
    }
    
    # identification
    if(id.signs & (irep < (nburn/2))){
      # first moment
      Flvl_sign <- sign(F_draw[1,1]) # first ruling normalized to be positive

      F_draw[,1] <- F_draw[,1]*Flvl_sign
      f_draw <- f_draw*Flvl_sign
      lambda[1,] <- lambda[1,]*Flvl_sign
    }
    
    # storage etc.
    if(irep %in% thin.set){
      in.thin <- in.thin+1
      F_store[in.thin,,] <- F_draw
      f_store[in.thin,] <- f_draw
      
      fht_store[in.thin,] <- exp(fHt)
      FHt_store[in.thin,,] <- exp(FHt)
      
      lambda_store[in.thin,,] <- t(lambda)
      delta_store[in.thin,,] <- delta
      
      # storage of irfs
      irf_store[in.thin,,1] <- as.numeric(lambda[1,])
      if(sv) svp_store[in.thin,] <- sv_para
      if(svf) svf_store[in.thin,,] <- svf_para[,1:3]
      
      if(mix) si_NE_store[in.thin,] <- si_draw[,sl.NE] else si_NE_store[in.thin,] <- 1
      if(mix) nJ_store[in.thin,1] <- sum(nk>0)
      
      si_store[in.thin,,] <- si_draw
      sig2_store[in.thin,] <- sig2_draw
      sig2t_store[in.thin,] <- sig2_t
    }
    if(!quiet) setTxtProgressBar(pb, irep)
  }
  
  F_store_T <- array(0,c(nthin,T))
  si_NE_store_T <- array(0,c(nthin,T))
  sig2t_store_T <- array(NA,c(nthin,T))
  
  cov_store <- array(NA,c(nthin,T_R,n,n,2))
  covAA_store <- covAAnoME_store <- array(NA,c(nthin,n,n))
  var_store <- array(NA,c(nthin,T_R,n))
  var2_store <- array(NA,c(nthin,T_R,n,3))
  
  for(irep in 1:nthin){
    f_draw <- F_store[irep,,1]
    si_NE <- si_NE_store[irep,]
    sig2t <- sig2t_store[irep,]
    
    F_store_T[irep,d==1] <- f_draw[1:T_R]
    si_NE_store_T[irep,d==0] <- si_NE[(T_R+1):length(si_NE)]
    sig2t_store_T[irep,d==0] <- sig2t[(T_R+1):length(sig2t)]
    sig2t_store_T[irep,d==1] <- sig2t[1:T_R]
    
    lambda_aa <- lambda <- lambda_store[irep,,]
    lambda_aa[,1] <- 0
    for(tt in 1:T_R){
      if(sv){
        diagHt <- diag(Qf) * c(fht_store[irep,tt],FHt_store[irep,tt,])
        eHt <- fht_store[irep,tt]
      }else{
        diagHt <- diag(Qf)
        eHt <- 1
      }
      ehLL <- eHt * (lambda[,1] %*% t(lambda[,1]))
      ehLLA <- lambda[,-1] %*% diagHt[-1,-1] %*% t(lambda[,-1])
      ehLLME <- lambda %*% diagHt %*% t(lambda) + (diag(n) * sig2t[tt])
      
      var_store[irep,tt,] <- diag(ehLL)/diag(ehLLME)
      cov_store[irep,tt,,,1] <- ehLL
      cov_store[irep,tt,,,2] <- ehLLME
      
      # variance shares by ruling
      eh_normalizer <- diag(ehLL) + diag(ehLLA) + rep(sig2t[tt],n)
      var2_store[irep,tt,,1] <- diag(ehLL)/eh_normalizer
      var2_store[irep,tt,,2] <- diag(ehLLA)/eh_normalizer
      var2_store[irep,tt,,3] <- rep(sig2t[tt],n)/eh_normalizer
    }
    
    covaa_run <- matrix(0,n,n)
    tt_count <- 0
    for(tt in (T_R+1):length(sig2t)){
      tt_count <- tt_count + 1
      LLaa <- lambda_aa %*% (diag(Qf)*c(rep(1,q),FHt_store[irep,tt,])) %*% t(lambda_aa)
      covaa_run <- covaa_run + (LLaa + (diag(n) * sig2t[tt]))
    }
    covAA_store[irep,,] <- covaa_run/tt_count
    covAAnoME_store[irep,,] <- LLaa
  }
  
  cov_tavg_store <- apply(cov_store,c(1,3,4,5),mean)
  var2_agg <- apply(var2_store,c(1,2,4),mean)
  var2_ciss <- apply(var2_store[,,grepl("CISS_",colnames(y)),],c(1,2,4),mean)
  
  var2_ag_post <- apply(var2_agg,c(2,3),median) / apply(apply(var2_agg,c(2,3),median),1,sum)
  var2_ciss_post <- apply(var2_ciss,c(2,3),median) / apply(apply(var2_ciss,c(2,3),median),1,sum)
  
  out <- list("F"=F_store,"Ft"=F_store_T,"ht"=fht_store,
              "lambda"=lambda_store,"delta"=delta_store,
              "irf"=irf_store,
              "si"=si_NE_store_T,"nJ"=nJ_store,
              "sig2"=sig2t_store,"sig2t"=sig2t_store_T,
              "varexpl"=var_store,"cov"=cov_tavg_store,"cov_AA"=covAA_store,"cov_AAj"=covAAnoME_store,
              "varexpl_ag"=var2_ag_post,"varexpl_ciss"=var2_ciss_post,
              "svpara"=svp_store)
}

# -----------------------------------------------------------------------------------------------
# JPR-algorithm
SV_JPR <- function(mu_h0,sig_h0,c_h=0,sv_para,
                   Em,hv,Hv,L,gamma,Sigma,
                   Qf,Q,TT,n,scale_h=0.1,mean.exp=TRUE){
  require(mvtnorm)
  phi_h <- sv_para[1]
  sig_h <- sv_para[2]
  rho <- sv_para[3]
  
  for (it in 0:TT){
    if (it==0){
      ht_sig <- 1/(1/sig_h0 + phi_h^2/sig_h)
      ht_mu <- ht_sig*(mu_h0/sig_h0 + phi_h*(hv[1]-c_h)/sig_h)
      h0d <- ht_mu + sqrt(ht_sig) * rnorm(1)
    }else{
      # define prior mean and variance t-by-t
      if (it==1){
        h_mut <- ((1-phi_h)*c_h + phi_h*(h0d + hv[2]))/(1+phi_h^2)
        h_sig <- sig_h/(1+phi_h^2)
        ht1 <- h0d
      }else if (it==TT){
        ht1 <- hv[it-1]
        h_mut <- c_h + phi_h*ht1
        h_sig <- sig_h
      }else{
        ht1 <- hv[it-1]
        h_mut <- ((1-phi_h)*c_h + phi_h*(ht1 + hv[it+1]))/(1+phi_h^2)
        h_sig <- sig_h/(1+phi_h^2)
      }
      
      h_old <- hv[it]                            # old draw
      h_prop <- h_old + rnorm(1,0,sqrt(scale_h)) # proposal
      
      # adjust mean
      if(mean.exp){
        e.0 <- Em[it,]-gamma*exp(h_prop)
        e.1 <- Em[it,]-gamma*exp(h_old)
      }else{
        e.0 <- Em[it,]-gamma*h_prop
        e.1 <- Em[it,]-gamma*h_old
      }
      
      # compute covariance matrix and means
      bbt <- L[,1,drop=FALSE] %*% t(L[,1,drop=FALSE])
      Omega <- (diag(n) * Sigma) + (L[,-1,drop=FALSE] %*% (diag(Q)*exp(Hv[it,])) %*% t(L[,-1,drop=FALSE]))
      
      m.0 <- rho/sqrt(sig_h) * (h_prop - phi_h * ht1) * exp(h_prop/2) * L[,1,drop=FALSE]
      m.1 <- rho/sqrt(sig_h) * (h_old - phi_h * ht1) * exp(h_old/2) * L[,1,drop=FALSE]
      v.0 <- (1-rho^2) * (exp(h_prop) * bbt) + Omega
      v.1 <- (1-rho^2) * (exp(h_old)  * bbt) + Omega
      
      alphanum <- mvnfast::dmvn(e.0,mu=m.0,sigma=v.0,log=TRUE) + dnorm(h_prop,h_mut,sqrt(h_sig),log=TRUE)
      alphaden <- mvnfast::dmvn(e.1,mu=m.1,sigma=v.1,log=TRUE) + dnorm(h_old,h_mut,sqrt(h_sig),log=TRUE)
      
      u <- log(runif(1,0,1))<(alphanum-alphaden)
      hv[it] <- u*h_prop+(1-u)*h_old
    }
  }
  return(list(h0=h0d,hv=hv))
}

# draw conjugate regression coefs
get.svpara <- function(y,X,b0,V0,a0,a1){
  n <- length(y)
  k <- NCOL(X)
  
  par1 <- (a0+n)/2
  var <- solve(crossprod(X)+solve(V0))
  mean <- matrix(var%*%(crossprod(X,y)+crossprod(solve(V0),b0)),k,1)
  par2 <- a0*a1 + sum((y-crossprod(t(X),mean))^2)
  par2 <- (par2 + crossprod(t(crossprod(mean-b0,V0)),mean-b0))/2
  
  sig2 <- 1/rgamma(1,par1,par2)
  var <- var*sig2
  mean <- mean + crossprod(t(chol(var)),rnorm(k))
  
  return(c(mean,sig2))
}
