library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)

source("utils.R")
source("spec_grid.R")
load("data/dailydata.rda")

sl.excl <- sl.excl[-c(2)]

# colors for IRFs
col_dark1 <- "#7E9ED3"
col_lght1 <- "#D8E2F0"
col_dark2 <- "#C48080"
col_lght2 <- "#F0D8D8"

col_grey1 <- "grey60"
red <- "#D63935"

middleblue2 <- "#75A1B7"
font.type <- "serif"

save.loc <- "plots"
dir.create(save.loc,showWarnings = FALSE)

# ------------------------------------------------------------------------
# specify IRF grid
grid.hor <- 0
grid.mix <- c(TRUE) #,FALSE
grid.pub <- c(TRUE)
grid.J <- c(30)
grid.sv <- c(TRUE)
grid.svm <- c(FALSE)
grid.svf <- c(TRUE)
grid.lev <- c(TRUE)
grid.exo <- c(TRUE)

grid.full <- expand.grid("mix"=grid.mix,
                         "pub"=grid.pub,
                         "J"=grid.J,
                         "Q"=grid.Q,
                         "h"=0:grid.hor,
                         "sv"=grid.sv,
                         "svm"=grid.svm,
                         "svf"=grid.svf,
                         "lev"=grid.lev,
                         "exo"=grid.exo)
grid.full <- grid.full[!(grid.full$sv==FALSE & grid.full$svm!=FALSE),]
rownames(grid.full) <- 1:nrow(grid.full)

# setup
p <- 1 # lags
nsave <- 9000
nburn <- 3000
thinfac <- 3

run <- 1
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

# ----------------------------------------------------------------------------------------
# VAR
source("utils_conjVAR.R")
load("data/exoshocks.rda")
load("data/FACTOR_mix_pub_exo_J30Q2sv_lev_facsv.rda")
shock_df <- F_ls[[1]][,c("date","h0")] %>%
  separate(col=date,into=c("Y","M","D"),sep="-") %>%
  mutate(D="01") %>%
  mutate("DATE"=paste0(Y,"-",formatC(M,width=2,flag="0"),"-",D)) %>%
  select(!c(Y,M,D))
colnames(shock_df)[1] <- "shock"

exoshock_df <- exoshocks %>%
  separate(col=date,into=c("Y","M","D"),sep="-") %>%
  mutate(D="01") %>%
  group_by(Y,M) %>%
  summarize(mptg = sum(mptg), mpfg = sum(mpfg), mpqe = sum(mpqe), oil = sum(oil)) %>%
  mutate("DATE"=paste0(Y,"-",formatC(M,width=2,flag="0"),"-01")) %>% ungroup() %>%
  select(DATE,mptg,mpfg,mpqe,oil)

load("data/VARdata.rda")

var.tfm <- c("shock"=0,"euribor_6m"=0,"eagovb_2y"=0,"estoxx50"=4,"ip"=3,"hicp"=3,"ciss"=0,"bbspread"=0,"mpu"=1)
var.ar1 <- c(0,1,1,0,1,1,1,1,0)
var.scl <- c(1,100,100,1,1,1,100,100,1)
var.lab <- c("raw","BPs","BPs","%","%","%","100 x raw","BPs","raw")
var.cln <- c("Ruling Shock","OIS 6m","Bond 2y","EuroStoxx 50","IP","HICP","CISS","AA OAS","MPU")

var.df <- data.frame("code"=names(var.tfm),
                     "scale"=var.lab,
                     "label"=var.cln)

Yraw <- ea_ts[,names(var.tfm)]
for(i in 1:ncol(Yraw)){
  if(var.tfm[i]==0){
    next
  }else if(var.tfm[i]==1){
    Yraw[,i] <- log(Yraw[,i])
  }else if(var.tfm[i]==2){
    Yraw[,i] <- c(NA,diff(1200*log(Yraw[,i])))
  }else if(var.tfm[i]==3){
    Yraw[,i] <- c(rep(NA,12),diff(100*log(Yraw[,i]),12))
  }else if(var.tfm[i]==4){
    Yraw[,i] <- c(NA,diff(100*log(Yraw[,i])))
  }
}

Yraw <- window(Yraw,start=c(2000,1))
Ysd <- c(1,apply(Yraw[,-1],2,sd))
Yraw[,-1] <- apply(Yraw[,-1],2,function(x) (x-mean(x))/sd(x))

p <- 12
ihor <- 25

Ylags <- embed(Yraw,p+1)
Y <- Ylags[,1:ncol(Yraw)]
X <- cbind(Ylags[,-c(1:ncol(Yraw))],1)

est_VAR <- conjVAR(Y=Y,X=X,var.ar1=var.ar1,p=p,ihor=ihor)
irf_resc <- est_VAR*NA
for(irep in 1:dim(est_VAR)[[1]]){
  irf_resc[irep,,] <- t(var.scl * Ysd * t(est_VAR[irep,,]))
}

# irf_post <- apply(est_VAR,c(2,3),quantile,probs=c(0.05,0.16,0.5,0.84,0.95))
irf_post <- apply(irf_resc,c(2,3),quantile,probs=c(0.05,0.16,0.5,0.84,0.95))
dimnames(irf_post) <- list(paste0("p",c(0.05,0.16,0.5,0.84,0.95)*100),
                           paste0("h",0:(ihor-1)),
                           colnames(Yraw))

irf_melt <- melt(irf_post) %>%
  pivot_wider(names_from = Var1, values_from = value) %>%
  left_join(var.df,by=c("Var3"="code")) %>%
  mutate(id = paste0(label," [",scale,"]")) %>%
  rename("horizon"="Var2","variable"="Var3") %>%
  mutate(horizon = as.character(horizon))

toplot <- c("ciss","euribor_6m","eagovb_2y","bbspread","estoxx50","hicp","ip","mpu")
pp_ls <- list()
for(i in toplot){
  tmp <- subset(irf_melt,variable %in% i)
  ycode <- unique(tmp$variable)
  ylab <- unique(tmp$id)
  pp_ls[[i]] <- tmp %>% ggplot() +
    geom_ribbon(aes(ymin=p5,
                    ymax=p95,
                    x=as.numeric(gsub("h","",horizon))),linewidth=0,
                fill=col_dark1,alpha=0.4) +
    geom_ribbon(aes(ymin=p16,
                    ymax=p84,
                    x=as.numeric(gsub("h","",horizon))),linewidth=0,
                fill=col_dark1,alpha=0.6) +
    geom_line(aes(y=p50,
                  x=as.numeric(gsub("h","",horizon))),color="black") +
    geom_vline(xintercept = 2, linetype = "dashed") +
    scale_x_continuous(breaks = c(0,6,12,18,24)) +
    
    geom_hline(yintercept=0,linewidth=1) +
    xlab("Horizon (months)") + ylab(ylab) +
    scale_size_manual(values=c(0.3,0.5,0.7,0.5,0.3)) +
    coord_cartesian(expand=FALSE,clip="off") +
    theme_cowplot() +
    theme(text = element_text(family = font.type),
          legend.position="none",axis.line.x=element_blank(),
          panel.grid.major = element_line(linewidth=0.3,color="grey90"))
  pdf(file=paste0(save.loc,"/",ycode,".pdf"),width=3,height=3)
  print(pp_ls[[i]])
  dev.off()
}

pdf(file=paste0(save.loc,"/irfm_all.pdf"),width=12,height=4.5)
print(plot_grid(plotlist = pp_ls,align="hv",ncol=4))
dev.off()

# ---------------------------------------------------------------------------------------------
# robustness wrt rulings
irf_ls <- list()
for(case in 0:3){
  Yraw_tmp <- Yraw
  rulings <- which(Yraw_tmp[,1]!=0)
  if(case == 0){
    Yraw_tmp <- Yraw
    rob <- "default"
  }else if(case == 1){
    Yraw_tmp[min(rulings),1] <- 0
    rob <- "nofirst"
  }else if(case == 2){
    Yraw_tmp[max(rulings),1] <- 0
    rob <- "nolast"
  }else if(case == 3){
    Yraw_tmp[c(min(rulings),max(rulings)),1] <- 0
    rob <- "nolarge"
  }
  
  Ylags <- embed(Yraw_tmp,p+1)
  Y <- Ylags[,1:ncol(Yraw)]
  X <- cbind(Ylags[,-c(1:ncol(Yraw))],1)
  
  est_VAR <- conjVAR(Y=Y,X=X,var.ar1=var.ar1,p=p,ihor=ihor)
  irf_resc <- est_VAR*NA
  for(irep in 1:dim(est_VAR)[[1]]){
    irf_resc[irep,,] <- t(var.scl * Ysd * t(est_VAR[irep,,]))
  }
  
  irf_post <- apply(irf_resc,c(2,3),quantile,probs=c(0.16,0.84))
  dimnames(irf_post) <- list(paste0("p",c(0.16,0.84)*100),
                             paste0("h",0:(ihor-1)),
                             colnames(Yraw))
  
  irf_ls[[rob]] <- melt(irf_post) %>%
    pivot_wider(names_from = Var1, values_from = value) %>%
    left_join(var.df,by=c("Var3"="code")) %>%
    mutate(id = paste0(label," [",scale,"]")) %>%
    rename("horizon"="Var2","var"="Var3") %>%
    mutate(horizon = as.character(horizon)) %>%
    mutate(robustness = rob)
}

irf_rob <- do.call("bind_rows",irf_ls)
irf_wide <- melt(irf_rob,id=c("horizon","var","scale","label","id","robustness"))

toplot <- c("ciss","euribor_6m","eagovb_2y","bbspread","estoxx50","hicp","ip","mpu")
pp_ls <- list()
for(i in toplot){
  tmp <- subset(irf_rob,var %in% i)
  tmp_wide <- subset(irf_wide, var %in% i)
  
  ycode <- unique(tmp$var)
  ylab <- unique(tmp$id)
  
  pp_ls[[i]] <- tmp %>% subset(robustness == "default") %>% ggplot() +
    geom_ribbon(aes(ymin=p16,
                    ymax=p84,
                    x=as.numeric(gsub("h","",horizon))),linewidth=0,
                fill="grey70",alpha=0.6) +
    geom_line(aes(y=value,
                  x=as.numeric(gsub("h","",horizon)),linewidth = variable),
              data=subset(tmp_wide,robustness == "nofirst"),color="#D63935") +
    geom_line(aes(y=value,
                  x=as.numeric(gsub("h","",horizon)),linewidth = variable),
              data=subset(tmp_wide,robustness == "nolast"),color="#1f78b4") +
    geom_line(aes(y=value,
                  x=as.numeric(gsub("h","",horizon)),linewidth = variable),
              data=subset(tmp_wide,robustness == "nolarge"),color="#33a02c") +
    
    geom_vline(xintercept = 2, linetype = "dashed") +
    scale_x_continuous(breaks = c(0,6,12,18,24)) +
    scale_linewidth_manual(values = c(1,1)) +
    
    geom_hline(yintercept=0,linewidth=1) +
    xlab("Horizon (months)") + ylab(ylab) +
    coord_cartesian(expand=FALSE,clip="off") +
    theme_cowplot() +
    theme(text = element_text(family = font.type),
          legend.position="none",axis.line.x=element_blank(),
          panel.grid.major = element_line(linewidth=0.3,color="grey90"))
  
  pdf(file=paste0(save.loc,"/rob_",ycode,".pdf"),width=3,height=3)
  print(pp_ls[[i]])
  dev.off()
}

pdf(file=paste0(save.loc,"/irfm_rob_all.pdf"),width=12,height=4.5)
print(plot_grid(plotlist = pp_ls,align="hv",ncol=4))
dev.off()


