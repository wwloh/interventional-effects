rm(list=ls())
libraries_check <- c("data.table","mgcv","ranger")
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

a01 <- 1
a02 <- -2; a2 <- 1; e12 <- 1.5
s2 <- 1; sY <- 1
b0 <- 1.6; bA <- 0.4; b1 <- 1.2; b2 <- 0.6

OneData <- function(n=400) {
  A <- rbinom(n, size = 1, prob = 0.5)
  M1 <- rpois(n, lambda = a01)
  M2 <- rnorm(n, mean = a02 + a2*A + e12*M1, sd = s2)
  M2_probit <- pnorm(M2) # plot(M2,M2_probit,col=A+1)
  M2 <- rbinom(n,size=1,prob=M2_probit)
  Y <- rnorm(n, mean = b0 + bA*A + b1*sqrt(M1) + b2*(M1^3), sd = sY)
  return(data.frame(id = 1:n, A, M1, M2, Y))
}

true_effects <- list()
true_effects[["theta0"]] <- bA
true_effects[["theta1"]] <- 0
true_effects[["theta2"]] <- 0
true_effects[["gamma0"]] <- true_effects$theta0
true_effects[["gamma1"]] <- true_effects$theta1 + true_effects$theta2

# helper function to create duplicated data
Dupdata <- function(t,A) {
  out <- matrix(A,nrow=t+1,ncol=t+1)
  diag(out) <- 1-A
  out <- data.frame(rbind(rep(A,times=(t+1)),out))
  colnames(out) <- paste0("a",0:t)
  out <- cbind("a.i"=1:(t+2),out)
  return(out)
}

OneMCestimator <- function(data,mc_draws=1e3,y_fit) {
  # data=OneData();mc_draws=4;y_fit="ranger"
  res <- list()

  ## fit outcome and mediator models ------------------------------------------
  fitM1 <- glm(M1 ~ A, family = poisson(link="log"), data = data)
  fitM2 <- glm(M2 ~ A+M1, family = binomial(link="probit"), data = data)
  if (y_fit=="true") {
    fitY <- glm(Y ~ A + I(sqrt(M1)) + I(M1^3),
                family = gaussian("identity"), data = data)
  } else if (y_fit=="add") {
    fitY <- glm(Y ~ A + M1 + M2, family = gaussian("identity"), data = data)
  } else {
    data_fitY <- data
    data_fitY[,"M1sqrt"] <- sqrt(data$M1)
    if (y_fit=="gam") {
      fitY <- gam(Y ~ A + te(M1,M1sqrt) + M2,
                  family = gaussian("identity"), data = data_fitY)
      # gam.check(fitY)
    } else if (y_fit=="ranger") {
      data_fitY[,"M1pow3"] <- data$M1^3
      fitY <- ranger(Y ~ A+M1+M2+M1sqrt+M1pow3, data = data_fitY)
    }
  }
  # helper function to predict Y for different outcome models
  PredictY <- function(onedat) {
    # onedat: data.table with columns A,M1,M2
    onedat <- onedat[, list(A,M1,M2,sqrt(M1),M1^3)]
    setnames(onedat, c("A","M1","M2","M1sqrt","M1pow3"))
    if (y_fit=="true" || y_fit=="add") {
      Ya <- predict.glm(fitY, type="response", newdata=onedat)
    } else if (y_fit=="gam") {
      Ya <- predict.gam(fitY, newdata=onedat)
    } else if (y_fit=="ranger") {
      Ya <- predictions(predict(fitY,type="response",data=onedat))
    }
    return(Ya)
  }

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
  dat[a0!=A & a1==A, "Y.a" := PredictY(
    onedat=dat[a0!=A & a1==A, list("A"=a0,"M1"=M1.a1,"M2"=M2.a1)])]

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
    mydt_mc[is.na(M1.a1),"M1.a1" := rpois(
      n=mc_draws,
      lambda=predict.glm(fitM1, newdata=data.frame("A"=a1),type="response")),
      by=list(a.i,a1)] # group by a1 since not subsetting by a1
    ## M2 from marginal distribution
    #### first sample M1 given a2
    mydt_mc[is.na(M2.a2),"M1.a2" := rpois(
      n=mc_draws,
      lambda=predict.glm(fitM1, newdata=data.frame("A"=a2),type="response")),
      by=list(a.i,a2)]
    #### then sample M2 given a2, M1.a2
    mydt_mc[is.na(M2.a2),"M2.a2" := rbinom(
      n=mc_draws,size=1,
      prob=predict.glm(fitM2, newdata=data.frame("A"=a2,"M1"=M1.a2),
                       type="response")),
      by=list(a.i,a2)]
    mydt_mc[, "M1.a2" := NULL]
    
    mydt_mc[, "Y.a" := PredictY(
      onedat=mydt_mc[, list("A"=a0,"M1"=M1.a1,"M2"=M2.a2)])]
    setkey(mydt_mc)
    
    ## average over MC draws
    mydt_mc_means <- mydt_mc[,lapply(.SD, mean),by=list(a.i,a0,a1,a2)]
    setkey(mydt_mc_means)
    # rows of duplicated data to impute
    mydt_a1 <- mydt[,list(a.i,a0,a1,a2)]
    setkey(mydt_a1)
    
    return(merge(mydt_a1,mydt_mc_means,sort=FALSE)[,Y.a])
  }
  dat[, "Y.a" := SampleMmarginal(.SD), by=id]
  setkey(dat)

  fit_mod2 <- glm(Y.a ~ a0 + a1 + a2, family = gaussian("identity"), 
                  data = dat)

  theta <- coef(fit_mod2)
  names(theta) <- paste0("t",lapply(
    strsplit(names(theta),"a"),paste0,collapse=""))
  res[["theta"]] <- theta
  rm(theta,fit_mod2)

  res[["ie12"]] <- as.numeric(res$gamma["g1"]-sum(res$theta[paste0("t",1:2)]))
  
  return(res)
}

nsims <- 10
ptm=proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  res_list <- list()
  data=OneData()
  
  # true Y model
  res <- OneMCestimator(data,mc_draws=1e3,y_fit="true")
  
  names(res) <- paste0(names(res),".true")
  res_list <- c(res_list,res)
  rm(res)
  
  # fitted GAM for Y
  res <- OneMCestimator(data,mc_draws=1e3,y_fit="gam")
  
  names(res) <- paste0(names(res),".GAM")
  res_list <- c(res_list,res)
  rm(res)
  
  # fitted Random Forests for Y
  res <- OneMCestimator(data,mc_draws=1e3,y_fit="ranger")
  
  names(res) <- paste0(names(res),".ranger")
  res_list <- c(res_list,res)
  rm(res)
  
  # misspecified additive Y model
  res <- OneMCestimator(data,mc_draws=1e3,y_fit="add")
  
  names(res) <- paste0(names(res),".add")
  res_list <- c(res_list,res)
  rm(res)
  
  return(res_list)
})
proc.time()[3]-ptm
# 130s per sim

save(simres,file=paste0("interventional-sim5-MConly-2M-noIEs-",
                        seed,".Rdata"))
q()

# results ---------------------------------------------------------------------
setwd("sim6-2M-noIEs/")
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

# check which seeds failed to complete
saved <- sapply(myfiles, function(x) {
  xsplit = strsplit(x,"-")[[1]]
  as.integer(strsplit(xsplit[length(xsplit)],".Rdata")[[1]])
})
for (s in 1:max(saved)) if (!(s %in% sort(saved))) cat(s,"\n")

# number of sims & MC error
nrow(res)
mc_err <- sqrt((1/nrow(res))*(1-1/nrow(res)))

library("xtable")
for (para in c(paste0("g",0:1),paste0("t",0:2),"ie12")) {
  (meths <- names(res)[grep(para,names(res))])
  names(meths) <- c("True","GAM","Random Forests","Additive")[1:length(meths)]
  if (strsplit(para,"")[[1]][1]=="g") {
    true_val <- true_effects[[paste0("gamma",strsplit(para,"")[[1]][2])]]  
  } else if (strsplit(para,"")[[1]][1]=="t") {
    true_val <- true_effects[[paste0("theta",strsplit(para,"")[[1]][2])]]
  } else {
    true_val <- 0
  }
  print(xtable(cbind(
    "method"=names(meths),
    "truth"=true_val, 
    "mean"=round(colMeans(res[,meths,with=FALSE]),2),
    "ese"=round(apply(res[,meths,with=FALSE], 2, sd),2))))
}