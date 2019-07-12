rm(list=ls())
libraries_check <- c("data.table","xtable")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))

a01 <- 3; a1 <- 1.2; aU1 <- 2
a02 <- 2; a2 <- 1.6; aU2 <- -2; e12 <- 2
s1 <- 1; s2 <- 1; sY <- 1
b0 <- 1.6; bA <- 0.4; b1 <- 0.6; b2 <- 1.2

# no hidden confounding to check that NE estimates are unbiased
if (seed > 1000) aU1 <- aU2 <- 0

OneData <- function(n=400) {
  A <- rbinom(n, size = 1, prob = 0.5)
  U <- rnorm(n, mean = 1)
  M1 <- rnorm(n, mean = a01 + a1*A + aU1*U, sd = s1)
  M2 <- rnorm(n, mean = a02 + a2*A + aU2*U + e12*M1, sd = s2)
  Y <- rnorm(n, mean = b0 + bA*A + b1*M1 + b2*M2, sd = sY)
  return(data.frame(id = 1:n, A, M1, M2, Y))
}

true_effects <- list()
true_effects[["theta0"]] <- bA
true_effects[["theta1"]] <- b1*a1
true_effects[["theta2"]] <- b2*(a2+a1*e12)
true_effects[["gamma0"]] <- true_effects$theta0
true_effects[["gamma1"]] <- true_effects$theta1 + true_effects$theta2
true_effects[["ne0"]] <- true_effects$theta0
true_effects[["ne1"]] <- b1*a1 + b2*(a1*e12)
true_effects[["ne2"]] <- b2*a2

# helper function to create duplicated data
Dupdata <- function(t,A) {
  out <- matrix(A,nrow=t+1,ncol=t+1)
  diag(out) <- 1-A
  out <- data.frame(rbind(rep(A,times=(t+1)),out))
  colnames(out) <- paste0("a",0:t)
  out <- cbind("a.i"=1:(t+2),out)
  return(out)
}

OneMCestimator <- function(data,mc_draws=1e3) {
  # data=OneData();mc_draws=4
  res <- list()

  ## fit outcome and mediator models ------------------------------------------
  fitM1 <- glm(M1 ~ A, family = gaussian("identity"), data = data)
  fitM2 <- glm(M2 ~ A + M1, family = gaussian("identity"), data = data)
  fitM2_A <- glm(M2 ~ A, family = gaussian("identity"), data = data)
  fitY <- glm(Y ~ A + M1 + M2, family = gaussian("identity"), data = data)

  ## direct and joint indirect effects ----------------------------------------
  #### joint mediators as a single mediator
  #### both combinations of a0,a1
  dat <- data.table(data)
  setkey(dat)
  alevels <- dat[,as.data.table(cbind(
    "a.i"=1:2,
    "a0"=c(A,1-A),
    "a1"=c(A,A))),
    by=id]
  setkey(alevels)
  dat <- merge(dat,alevels)
  setkey(dat)
  rm(alevels)

  # observed mediator and outcome values
  dat[a1==A, paste0("M",1:2,".a1") := list(M1,M2)]
  dat[a0==A & a1==A, "Y.a" := Y] # observed value of Y
  dat[a0!=A & a1==A, "Y.a" := predict.glm( # predicted value of Y
    fitY, newdata=data.frame("A"=a0,"M1"=M1,"M2"=M2),
    type="response")]

  fit_mod1 <- glm(Y.a ~ a0 + a1, family = gaussian("identity"),
                  data = dat)

  gamma <- coef(fit_mod1)
  names(gamma) <- paste0("g",lapply(
    strsplit(names(gamma),"a"),paste0,collapse=""))
  res[["gamma"]] <- gamma
  rm(gamma,fit_mod1,dat)

  ## mediator-specific indirect effects ---------------------------------------
  t <- length(grep("M",colnames(data))) # number of distinct mediators
  #### all t+2 combinations of a(0),...,a(t)
  dat <- data.table(data)
  setkey(dat)
  alevels <- dat[,as.data.table(Dupdata(t,A)),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels)
  setkey(dat)
  rm(alevels)

  #### impute (average) potential outcomes using marginal mediator dist. only
  ## largest mediator index
  obsM <- apply(dat[,list(A,a1,a2)],1,function(x) {
    ifelse(any(x==x["A"]),max(which(x==x["A"]))-1,0)
  })
  ## set to observed value
  dat[a.i %in% which(obsM==1), "M1.a1" := M1] # observed values of M1
  dat[a.i %in% which(obsM==2), "M2.a2" := M2] # observed values of M2
  
  SampleMmarginal <- function(mydt) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)]
    setkey(mydt_mc)
    ## sample M1 given a1
    mydt_mc[is.na(M1.a1),"M1.a1" := predict.glm(
      fitM1, newdata=data.frame("A"=a1),type="response")+
        rnorm(n=mc_draws,sd=sqrt(summary(fitM1)$dispersion)),
      by=list(a.i,a1)]
    ## sample M2 given a2
    mydt_mc[is.na(M2.a2), "M2.a2" := predict.glm(
      fitM2_A, newdata=data.frame("A"=a2),type="response")+
        rnorm(n=mc_draws,sd=sqrt(summary(fitM2_A)$dispersion)),
      by=list(a.i,a2)]
    # predicted value of Y
    mydt_mc[, "Y.a" := predict.glm(
      fitY, newdata=data.frame(
        "A"=a0,"M1"=M1.a1,"M2"=M2.a2),
      type="response")]
    setkey(mydt_mc)
    
    ## average over MC draws
    mydt_mc_means <- mydt_mc[,lapply(.SD, mean),by=list(a.i,a0,a1,a2)]
    setkey(mydt_mc_means)
    # rows of duplicated data to impute
    mydt_a1 <- mydt[,list(a.i,a0,a1,a2)]
    setkey(mydt_a1)
    
    return(merge(mydt_a1,mydt_mc_means,sort=FALSE)[,list(
      M1.a1,M2.a2,Y.a)])
  }
  dat[, c("M1.a1","M2.a2","Y.a") := SampleMmarginal(.SD),
      by=id]
  setkey(dat)

  fit_mod2 <- glm(Y.a ~ a0 + a1 + a2, family = gaussian("identity"), 
                  data = dat)

  theta <- coef(fit_mod2)
  names(theta) <- paste0("t",lapply(
    strsplit(names(theta),"a"),paste0,collapse=""))
  res[["theta"]] <- theta
  rm(theta,fit_mod2)

  res[["ie12"]] <- as.numeric(res$gamma["g1"]-sum(res$theta[paste0("t",1:2)]))
  
  ## path analysis estimators -------------------------------------------------
  ## function of parameters in models fitted to observed data -----------------
  alpha1 <- coef(fitM1)
  alpha2 <- coef(fitM2)
  beta <- coef(fitY)
  a_1 <- alpha1["A"]
  a_2 <- alpha2["A"]
  e_12 <- alpha2["M1"]
  b_A <- beta["A"]
  b_1 <- beta["M1"]
  b_2 <- beta["M2"]
  
  gamma_lm <- c(b_A,b_1*a_1+b_2*(a_2+a_1*e_12))
  names(gamma_lm) <- names(res$gamma)[-1]
  
  theta_lm <- c(b_A,b_1*a_1,b_2*(a_2+a_1*e_12))
  names(theta_lm) <- names(res$theta)[-1]
  
  ne_lm <- c(b_1*a_1 + b_2*(a_1*e_12),b_2*a_2)
  names(ne_lm) <- paste0("ne",1:2)
  
  res[["gamma_lm"]] <- gamma_lm
  res[["theta_lm"]] <- theta_lm
  res[["ne_lm"]] <- ne_lm
  
  return(res)
}

source("NEmodel-2M-wrapper.R")
nsims <- 100
ptm=proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  res_list <- list()
  data=OneData()
  
  res <- OneMCestimator(data,mc_draws=1e3)
  res_list <- c(res_list,res)
  rm(res)
  
  res <- OneSteenEstimator(data)
  res_list <- c(res_list,res)
  rm(res)
  
  return(res_list)
})
proc.time()[3]-ptm
# 20s per sim

save(simres,file=paste0("interventional-sim5-MConly-2M-hiddenU-",
                        seed,".Rdata"))
q()

# results ---------------------------------------------------------------------
setwd("sim6-2M-hiddenU/")
# setwd("sim6-2M-hiddenU/noU/") # no hidden confounding to check NE estimates
myfiles <- list.files()
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
sim_results_list <- NULL
for (ll in myfiles) {
  load(ll)
  sim_results_list <- rbind(sim_results_list,
                            do.call(rbind,lapply(simres, unlist)))
  rm(simres)
}
res <- data.table(sim_results_list)
setkey(res); rm(sim_results_list)

# root filename for plots
plot_rootfilename <- strsplit(ll,"-")[[1]]
plot_rootfilename <- paste(plot_rootfilename[-length(plot_rootfilename)],
                           collapse="-")

# number of sims & MC error
nrow(res)
mc_err <- sqrt((1/nrow(res))*(1-1/nrow(res)))

library("xtable")
for (para in c(paste0("g",0:1),paste0("t",0:2))) {
  (meths <- names(res)[grep(para,names(res))])
  names(meths) <- c("MC","Path","NE model","NE model (o)")[1:length(meths)]
  if (strsplit(para,"")[[1]][1]=="g") {
    true_val <- true_effects[[paste0("gamma",strsplit(para,"")[[1]][2])]]  
  } else if (strsplit(para,"")[[1]][1]=="t") {
    true_val <- true_effects[[paste0("theta",strsplit(para,"")[[1]][2])]]
    if (any(grepl("NE",names(meths)))) {
      true_val <- c(
        rep(true_val,2),
        rep(true_effects[[paste0("ne",strsplit(para,"")[[1]][2])]],2))
    }
  }
  print(xtable(cbind(
    "method"=names(meths),
    "truth"=true_val, 
    "mean"=colMeans(res[,meths,with=FALSE]),
    "ese"=apply(res[,meths,with=FALSE], 2, sd))))
}


# function to plot histograms
OneHistogram <- function(est,true,main_,xlim_=NA) {
  xlab_ <- "estimate"
  if(any(is.na(xlim_))) {
    hist(est, probability=TRUE, main=main_, xlab=xlab_)
  } else {
    xlim_[1] <- ifelse(xlim_[1]>0,xlim_[1]/1.2,xlim_[1]*1.2)
    xlim_[2] <- ifelse(xlim_[2]>0,xlim_[2]*1.2,xlim_[1]/1.2)
    hist(est, probability=TRUE, xlim=xlim_, main=main_, xlab=xlab_)
  }
  axis(1, at=true, labels="", lwd=2, col.ticks=4)
  abline(v=mean(est),col=2,lwd=2)
}

# Model 1: direct and joint indirect effects
for (tt in 0:1) {
  (meths <- names(res)[grep(paste0("g",tt),names(res))])
  names(meths) <- c("MC","lm")[1:length(meths)]
  
  xlim_ <- range(res[,meths,with=FALSE])
  true_val <- true_effects[[paste0("gamma",tt)]]
  
  filename <- paste0(plot_rootfilename,"-Model1-",tt,".png")
  
  png(filename,width=800*length(meths)/2,height=600,pointsize=25)
  par(mfrow=c(1,length(meths)))
  for (nn in 1:length(meths)) {
    OneHistogram(unlist(res[,meths[nn],with=FALSE]),true=true_val,
                 main=names(meths[nn]),xlim_=xlim_)
  }
  dev.off()
}

# Model 2: direct and indirect effects
for (tt in 0:2) {
  (meths <- names(res)[grepl(paste0("t",tt),names(res))])
  names(meths) <- c("MC","lm","NE","NE.o")[1:length(meths)]
  
  filename <- paste0(plot_rootfilename,"-Model2-",tt,".png")
  
  png(filename,width=800*length(meths)/2,height=600,pointsize=25)
  par(mfrow=c(1,length(meths)))
  for (nn in 1:length(meths)) {
    if (grepl("NE",names(meths)[nn])) {
      xlim_ <- range(res[,meths[grep("NE",names(meths))],with=FALSE])
      true_val <- true_effects[[paste0("ne",tt)]]
    } else {
      xlim_ <- range(res[,meths[grep("NE",names(meths),invert=TRUE)],with=FALSE])
      true_val <- true_effects[[paste0("theta",tt)]]
    }
    OneHistogram(unlist(res[,meths[nn],with=FALSE]),true=true_val,
                 main=names(meths[nn]),xlim_=xlim_)
  }
  dev.off()
}

# Model 1+2: gamma_2-(theta_1+theta_2+theta_12)
(meths <- names(res)[grepl("ie12",names(res))])
names(meths) <- c("MC")

xlim_ <- range(res[,meths,with=FALSE])
true_val <- 0

filename <- paste0(plot_rootfilename,"-Model2-ie12.png")

png(filename,width=800*length(meths)/2,height=600,pointsize=25)
par(mfrow=c(1,length(meths)))
for (nn in 1:length(meths)) {
  OneHistogram(unlist(res[,meths[nn],with=FALSE]),true=true_val,
               main=names(meths[nn]),xlim_=xlim_)
}
dev.off()
