rm(list=ls())
gc()

library(dplyr); library(tidyr); library(reshape2) # shaping data
library(ggplot2); library(cowplot) # plotting
library(invgamma) # other

# ------------------------------------------------------------------------
# select to plot
source("spec_grid.R")
scale.twice <- FALSE
font.type <- "serif"

mix <- TRUE
pub <- TRUE
J <- 30
Q <- 2
sv <- TRUE
svm <- FALSE
svf <- TRUE
lev <- TRUE
exo <- TRUE

# colors for IRFs
col_dark1 <- "#7E9ED3"
col_lght1 <- "#D8E2F0"
col_dark2 <- "#C48080"
col_lght2 <- "#F0D8D8"

col_grey1 <- "grey60"
red <- "#D63935"
middleblue2 <- "#75A1B7"

# some other stuff
load("alldata.rda")
load("mpu_indexes.rda")
df <- left_join(df,mpu_dd)

ruling.dates <- ruling.dates.all[-((nrow(ruling.dates.all)-2):nrow(ruling.dates.all)),]
ruling.dates$case <- factor(ruling.dates$case,levels=c("Pringle","Gauweiler","Weiss"))
ruling.dates <- ruling.dates[ruling.dates$ruling_pub!=0,] # check here to select pub/ruling event!
add.dates <- data.frame("date"=c("2012-07-26","2015-08-25","2016-06-24","2016-11-09","2020-02-28"),
                        "event"=c("'Whatever it takes'","Immigration crisis","Brexit referendum","Trump election","Covid-19"))
var.labs <- var.labs[!var.labs$code %in% sl.excl,]
var.labs$clean <- gsub(" \\(EA\\)","",var.labs$clean)

var.order <- var.labs$clean
var.order <- var.order[c(11:20,32:41,10,21:31,4:8,1:3,9,42)]

var.order <- c(var.order,"MPU")
var.labs <- rbind(var.labs,c("mpu","MPU"))

# ------------------------------------------------------------------------
# load model grid
save.dir <- "collect_mcmc"
spec0 <- paste0(ifelse(mix,"mix","nomix"),"_",
                ifelse(pub,"pub","decision"),"_",
                ifelse(exo,"exo","noexo"),"_",
                paste0("J",J),
                paste0("Q",Q),
                ifelse(sv,"sv",""),
                ifelse(svm,"-m",""),"_",
                ifelse(lev,"lev","iid"),"_",
                ifelse(svf,"facsv","fachom"))
snam <- paste0(save.dir,"/",spec0,".rda")
load(snam)

irf_store <- save_ls$irf_store
irf_rs_store <- save_ls$irf_rs_store
F_store <- save_ls$F_store
var_store <- save_ls$var_store
svpara_store <- save_ls$svpara_store
cov_store <- save_ls$cov_store
nJ_store <- save_ls$nJ_store
si_store <- save_ls$si_store
sig2_store <- save_ls$sig2_store

F_vola <- apply(log(F_store[,,,"vola"]),c(2,3),median)

# ------------------------------------------------------------------------
# get variable names and decide on scaling
varnames <- dimnames(irf_rs_store)[[2]]
ciss_sd <- sd(apply(df[,grepl("CISS",colnames(df))],1,mean))
scaleunit <- c(rep(100,14),rep(1,23),rep(100,10),rep(1,1),1)
scaleinfo <- data.frame("variable"=varnames,"scale"=scaleunit)
scaleinfo <- scaleinfo[!scaleinfo$variable %in% sl.excl,]

scaleinfo$label <- "BPs"
scaleinfo$label[scaleinfo$scale==1] <- "%"
scaleinfo$label[grepl("CISS",scaleinfo$variable)] <- "raw"
scaleinfo$label[scaleinfo$variable=="mpu"] <- "raw"

# table for data
tab.data <- left_join(var.labs,scaleinfo,by=c("code"="variable")) %>% select(clean,code,label)
colnames(tab.data) <- c("Label","Code","Scale")
tab.data$Label <- gsub(" \\(EA\\)","",tab.data$Label)

# tab.data <- data.frame(tab.data,
#                        "Source"=c(rep("FRED",3),
#                                   rep("ECB SDW",11),
#                                   rep("FRED",1),
#                                   rep("ECB SDW",11),
#                                   rep("Google Finance",11),
#                                   rep("Macrobond",10),
#                                   rep("FRED",1))) %>% arrange(Source)
# tab.data$Code <- gsub("_d","",tab.data$Code)
# print(xtable::xtable(tab.data),include.rownames = FALSE)

# groups for scales on y-axis
scaleinfo$group <- ""
scaleinfo$group[grepl("CISS",scaleinfo$variable)] <- "CISS"
scaleinfo$group[grepl("spread_de",scaleinfo$variable)] <- "SpreadsCore"
scaleinfo$group[grepl(paste0(c("es10","it10","pt10","ie10"),collapse="|"),scaleinfo$variable)] <- "SpreadPeri"
scaleinfo$group[grepl(paste0(c("gr10"),collapse="|"),scaleinfo$variable)] <- "SpreadGR"
scaleinfo$group[grepl("eureon",scaleinfo$variable)] <- "Rates"
scaleinfo$group[grepl("bofaml",scaleinfo$variable)] <- "Spreadrates"
scaleinfo$group[grepl("stockmarket",scaleinfo$variable) | grepl("ESTOXX50",scaleinfo$variable)] <- "Stocks"
scaleinfo$group[grepl(paste0(c("DEXUSEU","gold_eur_d"),collapse="|"),scaleinfo$variable)] <- "Other"
scaleinfo$group[grepl(paste0(c("mpu"),collapse="|"),scaleinfo$variable)] <- "MPU"

# ------------------------------------------------------------------------
# posteriors
F_post <- apply(F_store,2:4,quantile,probs=conf,na.rm=TRUE)
irf_post <- apply(irf_rs_store,2:4,quantile,probs=conf,na.rm=TRUE); rm(irf_rs_store)
irf_og_post <- apply(irf_store,2:4,quantile,probs=conf,na.rm=TRUE); rm(irf_store)

sig2_post <- apply(sig2_store,2:3,quantile,probs=conf,na.rm=TRUE)
si_post <- apply(si_store,2:3,median,na.rm=TRUE)
var_post <- apply(var_store,2:3,quantile,probs=conf,na.rm=TRUE)
nJ_post <- apply(nJ_store,2,mean,na.rm=TRUE)

dimnames(F_post)[[1]] <- dimnames(irf_post)[[1]] <- dimnames(irf_og_post)[[1]] <- dimnames(sig2_post)[[1]] <- dimnames(var_post)[[1]] <- paste0("p",100*conf)

# reshape everything for easier handling
irf_melt <- melt(irf_post) %>% 
  rename("conf"="Var1",
         "variable"="Var2",
         "horizon"="Var3",
         "moment"="Var4")

irf_og_melt <- melt(irf_og_post) %>% 
  rename("conf"="Var1",
         "variable"="Var2",
         "horizon"="Var3",
         "moment"="Var4")
rm(irf_post,irf_og_post)

cov_melt <- melt(cov_store) %>%
# cov_melt <- melt(cor_store) %>%
  rename("horizon"="Var1",
         "variable_x"="Var2",
         "variable_y"="Var3",
         "type"="Var4")

F_melt <- melt(F_post) %>% rename("conf"="Var1",
                                  "date"="Var2",
                                  "horizon"="Var3",
                                  "moment"="Var4") %>%
  pivot_wider(names_from = conf, values_from = value)

F_full <- melt(F_store) %>% rename("mcmc"="Var1",
                                   "date"="Var2",
                                   "horizon"="Var3",
                                   "moment"="Var4")
rm(F_store)

sv_melt <- melt(svpara_store) %>% rename("mcmc"="Var1",
                                         "parameter"="Var2",
                                         "horizon"="Var3",)

sig2_melt <- melt(sig2_post) %>% 
  rename("conf"="Var1",
         "date"="Var2",
         "horizon"="Var3")
si_melt <- melt(si_post) %>% 
  rename("date"="Var1",
         "horizon"="Var2")
var_melt <- melt(var_post) %>% 
  rename("conf"="Var1",
         "variable"="Var2",
         "horizon"="Var3")

spec <- paste0(ifelse(mix,"mix","nomix"),"_",
               ifelse(pub,"pub","decision"),"_",
               ifelse(exo,"exo","noexo"),"_",
               paste0("J",J),
               paste0("Q",Q),
               ifelse(sv,"sv",""),
               ifelse(svm,"-m",""),"_",
               ifelse(lev,"lev","iid"),"_",
               ifelse(svf,"facsv","fachom"),"_",
               paste0("h",formatC(0,width=2,flag="0")))
fnam <- paste0("irf_mcmc","/",spec,".rda")
if(file.exists(fnam)){
  load(fnam)
  tmp_ag <- data.frame(est$varexpl_ag)
  tmp_ciss <- data.frame(est$varexpl_ciss)
  
  colnames(tmp_ag) <- colnames(tmp_ciss) <- c("RF","AAF","ID")
  tmp_ag$id <- "All"
  tmp_ciss$id <- "CISS"
  tmp_ag$date <- tmp_ciss$date <- as.character(levels(F_melt$date))
  var_ruling_post <- bind_rows(tmp_ag,tmp_ciss) %>% 
    melt(id = c("id","date")) %>%
    select(date,variable,id,value)
}else{
  var_ruling_post <- NA
}

rm(sig2_store,cov_store,si_store,var_store)
gc()

# extract principal components from the factor
F_median <- F_post[3,,,"lvl"]
F_pca_R <- prcomp(F_median,center=TRUE,scale=TRUE)$x[,1]*(-1)
F_pca_h <- prcomp(t(F_median),center=TRUE,scale=TRUE)$x[,1]

F_pca_R <- (F_pca_R-mean(F_pca_R))/sd(F_pca_R)
F_pca_h <- (F_pca_h-mean(F_pca_h))/sd(F_pca_h)

F_pca_R <- data.frame("date"=names(F_pca_R),"value"=F_pca_R)
F_pca_h <- data.frame("horizon"=names(F_pca_h),"value"=F_pca_h)
rownames(F_pca_R) <- rownames(F_pca_h) <- NULL

F_point <- data.frame("date"=rownames(F_median),F_median)
rownames(F_point) <- NULL

F_ls <- list("point"=F_point,"pca"=F_pca_R)

# -----------------------------------------------------------------------------------------------------
# output all other objects of interest
save.loc <- "plots"
dir.create(save.loc,showWarnings = FALSE)
spec.loc <- spec0
dir.create(paste0(save.loc,"/",spec.loc),showWarnings = FALSE)
dir.create(paste0(save.loc,"/",spec.loc,"/irfs"),showWarnings = FALSE)
dir.create(paste0(save.loc,"/",spec.loc,"/other"),showWarnings = FALSE)

# subsetting
F_tmp <- F_full %>% subset(moment %in% "lvl") %>%
  subset(horizon %in% "h0")
F_tmp2 <- F_melt %>% subset(moment %in% "lvl") %>%
  subset(horizon %in% "h0")
H_tmp <- F_melt %>% subset(moment %in% "vola") %>%
  subset(horizon %in% "h0")
H_tmp_long <- H_tmp %>% melt(id=c("date","horizon","moment"))

pp_f <- F_tmp %>% ggplot() +
  geom_violin(aes(x=date,y=value)) + geom_hline(yintercept=0) + xlab("") + ylab(expression(f[t])) +
  geom_errorbar(aes(x=date,ymin=p5,ymax=p95),width=0.1,data=F_tmp2) +
  geom_errorbar(aes(x=date,ymin=p16,ymax=p84),width=0.3,data=F_tmp2) +
  geom_point(aes(x=date,y=p50),size=1.5,data=F_tmp2) +
  scale_x_discrete() +
  theme_cowplot() + 
  coord_cartesian(expand=TRUE,clip="off",ylim=c(-2,3)) +
  theme(text = element_text(family = font.type),
        axis.line.x=element_blank(),
        # axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        panel.grid.major.x = element_line(linewidth=0.7,color="grey90",linetype="longdash"))

pp_h <- H_tmp %>% ggplot(aes(x=date)) +
  geom_point(aes(y=p50),size=1.5) +
  geom_point(aes(y=p5),size=0.5) +
  geom_point(aes(y=p95),size=0.5) +
  geom_point(aes(y=p16),size=1) +
  geom_point(aes(y=p84),size=1) +
  geom_hline(yintercept=1) +
  geom_line(aes(x=date,y=value,group=variable,linewidth=variable),data=H_tmp_long) +
  scale_linewidth_manual(values=c(0.3,0.5,0.7,0.5,0.3)) +
  scale_x_discrete() +
  theme_cowplot() + ylab(expression(exp(h[t]))) + xlab("") +
  coord_cartesian(expand=TRUE,clip="off") +
  theme(text = element_text(family = font.type),
        axis.line.x=element_blank(),
        legend.position="none",
        # axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        panel.grid.major.x = element_line(linewidth=0.7,color="grey90",linetype="longdash"))

if(!any(is.na(var_ruling_post))){
  pp_var <- var_ruling_post %>% ggplot(aes(x=date)) +
    geom_bar(aes(y=100*value,group=id,fill=variable),width=0.5,stat="identity") +
    geom_hline(yintercept = 0) +
    facet_grid(id~.,switch="y") +
    theme_cowplot() + ylab("Share explained (%)") + xlab("Ruling") +
    coord_cartesian(expand=TRUE,clip="off") +
    scale_fill_manual(values = c(col_dark2,col_dark1,col_lght1),name="Component") +
    theme(text = element_text(family = font.type),
          axis.line.x=element_blank(),
          legend.position="bottom",
          axis.text.x = element_text(angle=90,vjust=0.5),
          panel.spacing.y = unit(0.5,"cm"),
          strip.placement = "outside",strip.background = element_blank())
  
  # including variance shares
  pdf(file = paste0(save.loc,"/",spec.loc,"/other/","factor.pdf"),width=10,height=8)
  print(plot_grid(pp_f,pp_h,pp_var,ncol=1,rel_heights = c(0.7,0.5,0.8),align="v",axis="bl"))
  dev.off()
}else{
  # excluding variance shares
  pdf(file = paste0(save.loc,"/",spec.loc,"/other/","factor.pdf"),width=10,height=8)
  print(plot_grid(pp_f,pp_h,ncol=1,rel_heights = c(0.7,0.5),align="hv",axis="trbl"))
  dev.off()
}

F_tmp <- F_melt %>% subset(moment %in% "lvl")
F_tmp$horizon <- as.numeric(gsub("h","",as.character(F_tmp$horizon)))

# sv learning process
sv_tmp <- sv_melt %>% 
  subset(horizon %in% "h0")

sv_phi <- sv_tmp %>% subset(parameter %in% "phi")
sv_sig2 <- sv_tmp %>% subset(parameter %in% "sig2")
sv_rho <- sv_tmp %>% subset(parameter %in% "rho")

phi_post <- quantile(sv_phi$value,c(0.16,0.5,0.84),na.rm=TRUE)
sig2_post <- quantile(sv_sig2$value,c(0.16,0.5,0.84),na.rm=TRUE)
rho_post <- quantile(sv_rho$value,c(0.16,0.5,0.84),na.rm=TRUE)

sig2_lim <- quantile(sv_sig2$value,c(0.9),na.rm=TRUE)

if(sv){
  dbetatfm <- function(x,shape1,shape2, ncp = 0, log=FALSE){
    return(dbeta((x + 1)/2, shape1 = shape1, shape2 = shape2, ncp = ncp, log = log))
  }
  
  para_low <- format(round(phi_post[1],digits=2),nsmall=2)
  para_high <- format(round(phi_post[3],digits=2),nsmall=2)
  pp_phi <- ggplot(sv_phi) + 
    stat_function(fun = dbetatfm, args = list(shape1 = 3, shape2 = 3),linetype="dashed") +
    geom_density(aes(x=value),fill="black",alpha=0.2,linewidth=0) + 
    coord_cartesian(expand=FALSE,clip="off") +
    scale_x_continuous(limits=c(-1,1),
                       breaks = c(-1,as.numeric(phi_post[1]),as.numeric(phi_post[2]),as.numeric(phi_post[3]),1),
                       labels = c(-1,"",format(round(as.numeric(phi_post[2]),digits=2),nsmall=2),"",1)) +
    # geom_vline(xintercept = phi_post,linewidth=c(0.1,0.3,0.7,0.3,0.1)) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = phi_post[2],linewidth=1) +
    ylab("Density") + xlab(bquote(phi[h] ~ " (" * .(para_low) ~","~ .(para_high)*")")) +
    theme_cowplot() + theme(text = element_text(family = font.type),axis.line.y = element_blank())
  
  para_low <- format(round(sig2_post[1],digits=2),nsmall=2)
  para_high <- format(round(sig2_post[3],digits=2),nsmall=2)
  pp_sig2 <- ggplot(sv_sig2) + 
    stat_function(fun = dinvgamma, args = list(shape = 1, rate = 0.1),linetype="dashed") +
    geom_density(aes(x=value),fill="black",alpha=0.2,linewidth=0) + 
    coord_cartesian(expand=FALSE,clip="off") +
    scale_x_continuous(limits=c(0,sig2_lim),
                       breaks = as.numeric(sig2_post),
                       labels = c("",format(round(as.numeric(sig2_post[2]),digits=2),nsmall=2),"")) +
    # geom_vline(xintercept = phi_post,linewidth=c(0.1,0.3,0.7,0.3,0.1)) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = sig2_post[2],linewidth=1) +
    ylab("") + xlab(bquote(sigma[h]^2 ~ " (" * .(para_low) ~","~ .(para_high)*")")) +
    theme_cowplot() + theme(text = element_text(family = font.type),axis.line.y = element_blank())
  
  para_low <- format(round(rho_post[1],digits=2),nsmall=2)
  para_high <- format(round(rho_post[3],digits=2),nsmall=2)
  pp_rho <- ggplot(sv_rho) + 
    stat_function(fun = dbetatfm, args = list(shape1 = 5, shape2 = 5),linetype="dashed") +
    geom_density(aes(x=value),fill="black",alpha=0.2,linewidth=0) + 
    coord_cartesian(expand=FALSE,clip="off") +
    scale_x_continuous(limits=c(-1,1),
                       breaks = c(-1,as.numeric(rho_post),1),
                       labels = c(-1,"",format(round(as.numeric(rho_post[2]),digits=2),nsmall=2),"",1)) +
    # geom_vline(xintercept = phi_post,linewidth=c(0.1,0.3,0.7,0.3,0.1)) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = rho_post[2],linewidth=1) +
    ylab("") + xlab(bquote(rho ~ " (" * .(para_low) ~","~ .(para_high)*")")) +
    # xlab(expression(rho)) +
    theme_cowplot() + 
    theme(text = element_text(family = font.type),axis.line.y = element_blank())
  
  pdf(file=paste0(save.loc,"/",spec.loc,"/other/","svpara.pdf"),width=10,height=2.5)
  print(plot_grid(pp_phi,pp_sig2,pp_rho,nrow=1))
  dev.off()
}

# -----------------------------------------------------------------------------------------------------
# output impulse response functions
ciss.sub <- varnames[grepl("CISS",varnames)]
ciss.mean <- irf_melt %>% subset(horizon %in% paste0("h",10)) %>%
  subset(moment %in% "beta") %>%
  subset(variable %in% ciss.sub) %>%
  subset(conf %in% "p50") %>%
  group_by(conf) %>%
  summarize(value=mean(value))
normalizer <- as.numeric(ciss.mean$value)

for(j in 1:nrow(scaleinfo)){
  sl.var <- scaleinfo$variable[j]
  sl.var.clean <- var.labs$clean[var.labs$code %in% sl.var]
  sl.var.scale <- scaleinfo$scale[scaleinfo$variable==sl.var]
  sl.var.lab <- scaleinfo$label[scaleinfo$variable==sl.var]
  
  sl.group.var <- scaleinfo[scaleinfo$group==scaleinfo$group[scaleinfo$variable==sl.var],"variable"]
  scale_tmp <- irf_melt %>% subset(variable %in% sl.group.var)
  scale_tmp_beta <- scale_tmp %>% subset(moment %in% "beta")
  scale_tmp_gamma <- scale_tmp %>% subset(moment %in% "gamma")
  
  ylimgroup_beta <- sort((-1)*sl.var.scale*c(min(scale_tmp_beta$value)/normalizer*ciss_sd,max(scale_tmp_beta$value)/normalizer*ciss_sd))
  ylimgroup_gamma <- sort((-1)*sl.var.scale*c(min(scale_tmp_gamma$value)/normalizer*ciss_sd,max(scale_tmp_gamma$value)/normalizer*ciss_sd))
  
  tmp.mix <- irf_melt %>% subset(variable %in% sl.var) %>%
    pivot_wider(names_from = conf,values_from = value)
  
  pp1 <- tmp.mix %>% ggplot() +
    geom_ribbon(aes(ymin=sl.var.scale*p5*ciss_sd/normalizer*(-1),
                    ymax=sl.var.scale*p95/normalizer*ciss_sd*(-1),
                    x=as.numeric(gsub("h","",horizon))),linewidth=0,
                fill=col_dark1,data=subset(tmp.mix,moment %in% "beta"),alpha=0.4) +
    geom_ribbon(aes(ymin=sl.var.scale*p16/normalizer*ciss_sd*(-1),
                    ymax=sl.var.scale*p84/normalizer*ciss_sd*(-1),
                    x=as.numeric(gsub("h","",horizon))),linewidth=0,
                fill=col_dark1,data=subset(tmp.mix,moment %in% "beta"),alpha=0.6) +
    geom_line(aes(y=sl.var.scale*p50/normalizer*ciss_sd*(-1),
                  x=as.numeric(gsub("h","",horizon)),color=moment,linetype=moment),
              data=subset(tmp.mix,moment %in% "beta")) +
    
    geom_hline(yintercept=0,linewidth=1) +
    xlab("Horizon (days)") + ylab(paste0(sl.var.clean," [",sl.var.lab,"]")) +
    scale_size_manual(values=c(0.3,0.5,0.7,0.5,0.3)) +
    scale_color_manual(values=c("black","black")) +
    coord_cartesian(expand=FALSE,clip="off",ylim=ylimgroup_beta) +
    theme_cowplot() +
    theme(text = element_text(family = font.type),
          legend.position="none",axis.line.x=element_blank(),
          panel.grid.major = element_line(linewidth=0.3,color="grey90"))
  
  pdf(file=paste0(save.loc,"/",spec.loc,"/irfs/",sl.var,".pdf"),width=3,height=2.5)
  print(pp1)
  dev.off()
  
  if(svm){
    pp2 <- tmp.mix %>% ggplot() +
      geom_ribbon(aes(ymin=sl.var.scale*p5,
                      ymax=sl.var.scale*p95,
                      x=as.numeric(gsub("h","",horizon))),linewidth=0,
                  fill=col_grey1,data=subset(tmp.mix,moment %in% "gamma"),alpha=0.4) +
      geom_ribbon(aes(ymin=sl.var.scale*p16,
                      ymax=sl.var.scale*p84,
                      x=as.numeric(gsub("h","",horizon))),linewidth=0,
                  fill=col_grey1,data=subset(tmp.mix,moment %in% "gamma"),alpha=0.6) +
      geom_line(aes(y=sl.var.scale*p50,
                    x=as.numeric(gsub("h","",horizon)),color=moment,linetype=moment),
                data=subset(tmp.mix,moment %in% "gamma")) +
      
      geom_hline(yintercept=0,linewidth=1) +
      xlab("Horizon (days)") + ylab(paste0(sl.var.clean," [",sl.var.lab,"]")) +
      scale_size_manual(values=c(0.3,0.5,0.7,0.5,0.3)) +
      scale_x_continuous(breaks = c(0,20,40,60)) +
      scale_color_manual(values=c("black","black")) +
      coord_cartesian(expand=FALSE,clip="off",ylim=ylimgroup_gamma) +
      theme_cowplot() +
      theme(text = element_text(family = font.type),
            legend.position="none",axis.line.x=element_blank(),
            axis.title = element_text(size=11),#axis.text = element_text(size=8),
            panel.grid.major = element_line(linewidth=0.3,color="grey90"))
    
    pdf(file=paste0(save.loc,"/",spec.loc,"/svm/",sl.var,".pdf"),width=2,height=1.9)
    print(pp2)
    dev.off()
    
    # both moments in one chart
    pp3 <- tmp.mix %>% ggplot() +
      geom_ribbon(aes(ymin=sl.var.scale*p5,
                      ymax=sl.var.scale*p95,
                      x=as.numeric(gsub("h","",horizon))),linewidth=0,
                  fill=col_grey1,data=subset(tmp.mix,moment %in% "gamma"),alpha=0.3) +
      geom_ribbon(aes(ymin=sl.var.scale*p16,
                      ymax=sl.var.scale*p84,
                      x=as.numeric(gsub("h","",horizon))),linewidth=0,
                  fill=col_grey1,data=subset(tmp.mix,moment %in% "gamma"),alpha=0.5) +
      
      geom_ribbon(aes(ymin=sl.var.scale*p5*ciss_sd/normalizer*(-1),
                      ymax=sl.var.scale*p95/normalizer*ciss_sd*(-1),
                      x=as.numeric(gsub("h","",horizon))),linewidth=0,
                  fill=col_dark1,data=subset(tmp.mix,moment %in% "beta"),alpha=0.3) +
      geom_ribbon(aes(ymin=sl.var.scale*p16/normalizer*ciss_sd*(-1),
                      ymax=sl.var.scale*p84/normalizer*ciss_sd*(-1),
                      x=as.numeric(gsub("h","",horizon))),linewidth=0,
                  fill=col_dark1,data=subset(tmp.mix,moment %in% "beta"),alpha=0.5) +
      
      geom_line(aes(y=sl.var.scale*p50/normalizer*ciss_sd*(-1),
                    x=as.numeric(gsub("h","",horizon))),linetype="solid",
                data=subset(tmp.mix,moment %in% c("beta"))) +
      geom_line(aes(y=sl.var.scale*p50,
                    x=as.numeric(gsub("h","",horizon))),linetype="dashed",
                data=subset(tmp.mix,moment %in% c("gamma"))) +
      
      geom_hline(yintercept=0,linewidth=1) +
      xlab("Horizon (days)") + ylab(paste0(sl.var.clean," [",sl.var.lab,"]")) +
      scale_size_manual(values=c(0.3,0.5,0.7,0.5,0.3)) +
      scale_color_manual(values=c("black","black")) +
      coord_cartesian(expand=FALSE,clip="off",ylim=ylimgroup_beta) +
      theme_cowplot() +
      theme(text = element_text(family = font.type),
            legend.position="none",axis.line.x=element_blank(),
            panel.grid.major = element_line(linewidth=0.3,color="grey90"))
    
    pdf(file=paste0(save.loc,"/",spec.loc,"/irfboth/",sl.var,".pdf"),width=2,height=1.9)
    print(pp3)
    dev.off()
  }
}

# ----------------------------------------------------------------------------------
# output regime classification
sig2_sub <- sig2_melt %>% 
  subset(conf %in% "p50")
si_sub <- si_melt
  
tmp <- F_melt %>% 
  subset(moment %in% "lvl") %>%
  subset(horizon %in% "h0")

sig2_sub <- left_join(sig2_sub,select(tmp,!moment))
sig2_sub[is.na(sig2_sub)] <- 0
sig2_sub$ma5 <- runner::runner(sig2_sub$value,f=function(x) mean(x),k=5)

pp0 <- ggplot(sig2_sub) +
  geom_hline(yintercept = 0) +
  geom_vline(aes(xintercept=as.Date(date),linetype=case),color="red",data=ruling.dates) +
  geom_bar(aes(x=as.Date(date),y=p50),stat="identity",width=15,fill="black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ylab(expression(f[t])) + xlab("") +
  coord_cartesian(expand=0) +
  theme_cowplot() + theme(text = element_text(family = font.type),
                          axis.line.x = element_blank(),legend.position="none")

pp <- ggplot() +
  geom_ribbon(aes(x=as.Date(date),ymin=-Inf,ymax=max(sig2_sub$value,na.rm=TRUE)*value),data=si_sub,fill="grey90") +
  # geom_hline(yintercept=1,color="grey50",size=1) +
  geom_vline(aes(xintercept=as.Date(date),linetype=case),color="red",data=ruling.dates) +
  geom_vline(aes(xintercept=as.Date(date)),color=middleblue2,data=add.dates) +
  geom_text(aes(x=as.Date(date),y=max(sig2_sub$value,na.rm=TRUE),label=event),
            color=middleblue2,size=3,hjust=1,vjust=-0.25,angle=90,data=add.dates,family=font.type) +
  geom_line(aes(x=as.Date(date),y=value),linewidth=0.3,data=sig2_sub,color="grey60") +
  geom_line(aes(x=as.Date(date),y=ma5),linewidth=0.5,data=sig2_sub) +
  ylab(expression(theta[t]^2)) + xlab("") +
  scale_linetype_discrete(name="Case") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(expand=0) +
  theme_cowplot() + theme(text = element_text(family = font.type),
                          legend.position="none",legend.key.height=unit(0.75,"cm"))
pp2 <- pp + coord_cartesian(expand=0,ylim=c(0,2)) + theme(legend.position="bottom",legend.key.height=unit(0.75,"cm"))

pdf(file=paste0(save.loc,"/",spec.loc,"/other/","minvariance.pdf"),width=10,height=6)
print(plot_grid(pp0,pp,pp2,ncol=1,align="hv",axis="l",rel_heights = c(0.25,0.4,0.35)))
dev.off()

# variance shares
var_melt$horizon <- as.numeric(gsub("h","",as.character(var_melt$horizon)))
var_melt$value <- var_melt$value*100
var_melt_sub <- var_melt %>% subset(variable %in% scaleinfo$variable)
var_melt_sub$variable <- factor(as.character(var_melt_sub$variable),levels=var.labs$code)
levels(var_melt_sub$variable) <- var.labs$clean
var_melt_sub$variable <- factor(as.character(var_melt_sub$variable),levels=rev(var.order))
var_melt <- var_melt_sub

pp_varexpl <- var_melt %>% 
  subset(conf %in% "p50") %>%
  ggplot(aes(x=horizon,y=variable,fill=value,color=value)) +
  geom_tile() + ylab("") + xlab("Horizon (days)") +
  scale_fill_gradient2(low="white",mid="#1f78b4",high=red,midpoint=50,limits=c(0,100),
                       name="Explained by factor\n(in percent)") +
  scale_color_gradient2(low="white",mid="#1f78b4",high=red,midpoint=50,limits=c(0,100),
                        name="Explained by factor\n(in percent)") +
  scale_x_continuous(breaks = seq(0,60,by=5)) +
  coord_cartesian(expand=FALSE) +
  theme_cowplot() + theme(text = element_text(family = font.type),
                          legend.position="bottom",legend.key.width=unit(2,"cm"),axis.text.y = element_text(size=10,hjust=0),
                          strip.background = element_blank(),strip.text = element_text(face="bold"),
                          panel.spacing = unit(0.5,"cm"))

pdf(file=paste0(save.loc,"/",spec.loc,"/other/","varshares.pdf"),width=10,height=6.5)
print(pp_varexpl)
dev.off()

# ----------------------------------------------------------------------------------
# covariance structure across R/NE/NR
cov_sub <- cov_melt %>% 
  subset(horizon %in% "h0") %>% 
  subset(variable_x %in% scaleinfo$variable) %>%
  subset(variable_y %in% scaleinfo$variable)
  
cov_sub$variable_x <- factor(as.character(cov_sub$variable_x),levels=var.labs$code)
cov_sub$variable_y <- factor(as.character(cov_sub$variable_y),levels=var.labs$code)
levels(cov_sub$variable_x) <- var.labs$clean
levels(cov_sub$variable_y) <- var.labs$clean
cov_sub$variable_x <- factor(cov_sub$variable_x,levels=(var.order))
cov_sub$variable_y <- factor(cov_sub$variable_y,levels=rev(var.order))

pp_ls <- list()
icount <- 0
for(ii in c("R_all","NR_full")){
  icount <- icount + 1
  if(icount == 1){
    gtitle <- "Ruling dates"
  }else if(icount == 5){
    gtitle <- "Ruling dates, full"
  }else if(icount == 2){
    gtitle <- "Non-ruling dates"
  }
  pp_ls[[ii]] <- cov_sub %>%
    subset(type %in% ii) %>%
    subset(!(variable_x %in% "Oil prices")) %>% subset(!(variable_y %in% "Oil prices")) %>%
    ggplot() +
    geom_tile(aes(x=variable_x,y=variable_y,fill=value)) +
    ggtitle(gtitle) +
    scale_fill_gradient2(mid="white",high=col_dark1,low=col_dark2,midpoint=0,na.value = "white",name="") +
    scale_color_gradient2(mid="white",high=col_dark1,low=col_dark2,midpoint=0,na.value = "white",name="") +
    coord_equal(expand=FALSE,clip="off") +
    xlab("") + ylab("") +
    theme_cowplot() + theme(text = element_text(family = font.type),
                            legend.position = "bottom",
                            panel.border = element_rect(linewidth=0.5,color="black"),
                            # axis.text = element_blank(),
                            axis.text.y = element_text(vjust=0.5,hjust=0,size=7.5),
                            axis.text.x = element_text(angle = 90,hjust=0,vjust=0.5,size=7.5),
                            axis.ticks = element_blank(),
                            legend.key.width = unit(1.5,"cm"),
                            legend.key.height = unit(0.1,"cm"))

}

pdf(file=paste0(save.loc,"/",spec.loc,"/other/","covariances.pdf"),width=10.5,height=6.5)
print(plot_grid(plotlist = pp_ls,nrow=1))
dev.off()

