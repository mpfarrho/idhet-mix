run <- commandArgs(trailingOnly = TRUE)

require(dplyr)
require(tidyr)

source("utils.R")
load("data/dailydata.rda")
load("data/mpu_indexes.rda")
df <- left_join(df,mpu_dd)

source("spec_grid.R") # specify IRF grid

mix <- grid.full[run,"mix"]
pub <- grid.full[run,"pub"]
exo <- grid.full[run,"exo"]
J <- grid.full[run,"J"]
Q <- grid.full[run,"Q"]
sv <- grid.full[run,"sv"]
svm <- grid.full[run,"svm"]
svf <- grid.full[run,"svf"]
lev <- grid.full[run,"lev"]
h <- grid.full[run,"h"]

dirn <- "irf_mcmc"
dir.create(dirn,showWarnings = FALSE)
dir.create("logs",showWarnings = FALSE)

spec <- paste0(ifelse(mix,"mix","nomix"),"_",
               ifelse(pub,"pub","decision"),"_",
               ifelse(exo,"exo","noexo"),"_",
               paste0("J",J),
               paste0("Q",Q),
               ifelse(sv,"sv",""),
               ifelse(svm,"-m",""),"_",
               ifelse(lev,"lev","iid"),"_",
               ifelse(svf,"facsv","fachom"),"_",
               paste0("h",formatC(h,width=2,flag="0")))
fnam <- paste0(dirn,"/",spec,".rda")
sl.excl <- c("gold_usd_d","oil_brent_d",
             "bofaml_us_hyld_oas_d","bofaml_us_aa_oas_d","bofaml_us_bbb_oas_d",
             "gs1_d","gs10_d","usbkeven5y_d","vixcls_d","nasdaqcom_d")

# ------------------------------------------------------------------------
# construct lag structure
if(exo){
  load("data/exoshocks.rda")
  df_exo <- left_join(select(df,"date"),exoshocks) %>% as.data.frame()
  df_exo[is.na(df_exo)] <- 0
  Xexo <- matrix(as.matrix(df_exo[,-1]),ncol=ncol(df_exo)-1)
}

Xraw <- matrix(as.matrix(df[,-1]),ncol=ncol(df)-1)
colnames(Xraw) <- colnames(df)[-1]

shock <- Xraw[,1:4]
endog <- Xraw[,-c(1:4)]
endog <- endog[,!(colnames(endog) %in% sl.excl)]

# lag endogenous variables forward
Xtmp <- cbind(shock,mlag(endog,lag=p))
yh <- endog[(p+h+1):nrow(endog),]
Xh <- Xtmp[(p+1):(nrow(Xtmp)-h),]
Dyth <- yh-Xh[,-c(1:4)] # deterministic RW behavior for financial series

# ruling dates wrt. time t
if(pub) dth <- Xh[,"ruling_pub"] else dth <- Xh[,"ruling"]

# normalize dependent variable
Dyth_mu <- apply(Dyth,2,mean)
Dyth_sd <- apply(Dyth,2,sd)
Dyth <- apply(Dyth,2,function(x) (x-mean(x))/sd(x))

# other controls or exogenous shocks etc.
if(exo){
  xth <- cbind(Xexo[(p+1):(nrow(Xexo)-h),],1)
}else{
  xth <- matrix(1,nrow(Dyth),1)
}

# estimate model
message(fnam)
if(!file.exists(fnam)){
  checkmode <- TRUE
  est <- try(sfm.fm(y=Dyth,d=dth,x=xth,
                nburn=nburn,nsave=nsave,thinfac=thinfac,
                mix=mix,J=J,sample.e0=TRUE,
                Q=Q,sv=sv,svm=svm,lev=lev,svf=svf,
                quiet=FALSE,checkmode=checkmode),silent=TRUE)
  
  if(is(est,"try-error")){
    cat(est[1],file=paste0("logs/",spec,".txt"))
  }
  
  est$stmu <- Dyth_mu
  est$stsd <- Dyth_sd
  
  # summarize some of the posteriors
  est$F <- est$F[,1:sum(dth),1]
  est$Ft <- NULL
  est$sig2 <- NULL
  
  est$varexpl <- apply(est$varexpl,c(1,3),mean,na.rm=TRUE)
  est$cov <- apply(est$cov,c(2,3,4),median,na.rm=TRUE)
  est$cov_AA <- apply(est$cov_AA,c(2,3),median,na.rm=TRUE)
  est$cov_AAj <- apply(est$cov_AAj,c(2,3),median,na.rm=TRUE)
  
  save(est,file=fnam)
}


