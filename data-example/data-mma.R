rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# load dataset
load(file="data-mma.Rdata")

# specify mediator and outcome models -----------------------------------------
covariates <- as.character("age+car+gotosch+numpeople+race")
models <- list()
models[["M1"]] <- paste0("M1~A+",covariates)
for (s in 2:7) {
  # marginal models
  models[[paste0("M",s,"_A")]] <- gsub(
    pattern="M1", replacement=paste0("M",s), x=models$M1)
}
for (s in 6:7) {
  # conditional models
  models[[paste0("M",s)]] <- gsub(
    pattern="A",
    replacement=paste0("A+",paste(paste0("M",1:(s-1)),collapse="+")), 
    x=models[[paste0("M",s,"_A")]])
}
models[["Y"]] <- paste0("Y~A+",paste(paste0("M",1:7),collapse="+"),
                        "+",covariates)

models$M6_A <- NULL
models$M7_A <- NULL

# for creating interactions terms
## all subsets from 1:7
allcombi <- expand.grid(lapply(1:7, function(x) 0:1))

## M6 mean model
fitM6_int2_AM <- paste(paste0("A:M",1:5),collapse="+")
fitM6_int2_M <-
  unique(apply(allcombi[rowSums(allcombi[,1:5])==2,1:5],1,function(x) {
    paste(paste0("M",which(x==1)),collapse=":")
  }))
fitM6_int3_AM <- paste0("A:",fitM6_int2_M)
models$M6 <- paste(models$M6,
                   fitM6_int2_AM,
                   paste(fitM6_int2_M,collapse="+"),
                   paste(fitM6_int3_AM,collapse="+"),
                   sep="+")
rm(fitM6_int2_AM,fitM6_int2_M,fitM6_int3_AM)

## M7 mean model
fitM7_int2_AM <- paste(paste0("A:M",1:6),collapse="+")
fitM7_int2_M <-
  unique(apply(allcombi[rowSums(allcombi[,1:6])==2,1:6],1,function(x) {
    paste(paste0("M",which(x==1)),collapse=":")
  }))
fitM7_int3_AM <- paste0("A:",fitM7_int2_M)
models$M7 <- paste(models$M7,
                   fitM7_int2_AM,
                   paste(fitM7_int2_M,collapse="+"),
                   paste(fitM7_int3_AM,collapse="+"),
                   sep="+")
rm(fitM7_int2_AM,fitM7_int2_M,fitM7_int3_AM)

## outcome model 
fitY_int2_AM <- paste("A",paste0("M",1:7),sep=":",collapse="+")
fitY_int2 <- paste(apply(allcombi[rowSums(allcombi)==2,],1,function(x) {
  paste(paste0("M",which(x==1)),collapse=":")
}),collapse="+")
fitY_int3_AM <- paste0("A:",gsub(pattern="[+]M",replacement="+A:M",x=fitY_int2))
models$Y <- paste(models$Y,fitY_int2,fitY_int2_AM,fitY_int3_AM,sep="+")

(models <- lapply(models, as.formula))

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
  ## fit mediator and outcome models ------------------------------------------
  fits <- lapply(names(models), function(mod) {
    if (mod %in% c("M6","M7")) {
      glm(models[[mod]], family = binomial("logit"), data = data)
    } else {
      glm(models[[mod]], family = gaussian("identity"), data = data)
    }
  })
  names(fits) <- names(models)
  
  # rank-deficient fit or failure to converge?
  if(any(unlist(lapply(fits, function(x) any(is.na(summary(x)$coef[,1:2]))))) |
     !fits$M6$converged | !fits$M7$converged) {
    break()
  }

  # significant predictors
  print(lapply(fits, function(x) summary(x)$coef[summary(x)$coef[,4]<.1,]))
  
  res <- list()
  
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
  dat[a1==A, paste0("M",1:7,".a1") := list(M1,M2,M3,M4,M5,M6,M7)]
  dat[a0==A & a1==A, "Y.a" := Y] # observed value of Y
  setkey(dat)
  dat[a0!=A & a1==A, "Y.a" := predict.glm( # predicted value of Y
    fits$Y, newdata=data.frame("A"=a0,"M1"=M1,"M2"=M2,"M3"=M3,
                               "M4"=M4,"M5"=M5,"M6"=M6,"M7"=M7,
                               age,car,gotosch,numpeople,race),
    type="response")]
  setkey(dat)
  
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
  obsM <- apply(dat[,list(A,a1,a2,a3,a4,a5,a6,a7)],1,function(x) {
    ifelse(any(x==x["A"]),max(which(x==x["A"]))-1,0)
  })
  ## set to observed values
  for (s in 1:7) {
    dat[a.i %in% which(obsM==s), paste0("M",s,".a",s) := 
          dat[a.i %in% which(obsM==s), paste0("M",s), with=FALSE]]
  }
  
  SampleMmarginal <- function(mydt) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)]
    setkey(mydt_mc)
    # M1 - M5 are fitted with linear and additive mean models
    # => sample from marginals
    ## sample M1 given a1
    mydt_mc[is.na(M1.a1),M1.a1 := predict.glm(
      fits$M1,
      newdata=data.frame("A"=a1,age,car,gotosch,numpeople,race),
      type="response")+
        rnorm(n=mc_draws,
              sd=sqrt(summary(fits$M1)$dispersion)),
      by=list(a.i,a1)]
    ## sample M2 given a2
    mydt_mc[is.na(M2.a2),M2.a2 := predict.glm(
      fits$M2_A,
      newdata=data.frame("A"=a2,age,car,gotosch,numpeople,race),
      type="response")+
        rnorm(n=mc_draws,
              sd=sqrt(summary(fits$M2_A)$dispersion)),
      by=list(a.i,a2)]
    ## sample M3 given a3
    mydt_mc[is.na(M3.a3),M3.a3 := predict.glm(
      fits$M3_A,
      newdata=data.frame("A"=a3,age,car,gotosch,numpeople,race),
      type="response")+
        rnorm(n=mc_draws,
              sd=sqrt(summary(fits$M3_A)$dispersion)),
      by=list(a.i,a3)]
    ## sample M4 given a4
    mydt_mc[is.na(M4.a4),M4.a4 := predict.glm(
      fits$M4_A,
      newdata=data.frame("A"=a4,age,car,gotosch,numpeople,race),
      type="response")+
        rnorm(n=mc_draws,
              sd=sqrt(summary(fits$M4_A)$dispersion)),
      by=list(a.i,a4)]
    ## sample M5 given a5
    mydt_mc[is.na(M5.a5),M5.a5 := predict.glm(
      fits$M5_A,
      newdata=data.frame("A"=a5,age,car,gotosch,numpeople,race),
      type="response")+
        rnorm(n=mc_draws,
              sd=sqrt(summary(fits$M5_A)$dispersion)),
      by=list(a.i,a5)]
    
    # M6 - M7 are fitted with non-linear models
    ## M6 from marginal distribution
    #### first sample M1-M5 given a6
    mydt_mc[is.na(M6.a6),
            c("M1.a6","M2.a6","M3.a6","M4.a6","M5.a6") := list(
              predict.glm(
                fits$M1,
                newdata=data.frame("A"=a6,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M1)$dispersion)),
              predict.glm(
                fits$M2_A,
                newdata=data.frame("A"=a6,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M2_A)$dispersion)),
              predict.glm(
                fits$M3_A,
                newdata=data.frame("A"=a6,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M3_A)$dispersion)),
              predict.glm(
                fits$M4_A,
                newdata=data.frame("A"=a6,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M4_A)$dispersion)),
              predict.glm(
                fits$M5_A,
                newdata=data.frame("A"=a6,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M5_A)$dispersion))),
            by=list(a.i,a6)]
    #### then sample M6 given a6, M1.a6, ..., M5.a6
    mydt_mc[is.na(M6.a6), M6.a6 := rbinom(
      n=mc_draws,size=1,
      prob=predict.glm(
        fits$M6,
        newdata=data.frame("A"=a6,"M1"=M1.a6,"M2"=M2.a6,"M3"=M3.a6,"M4"=M4.a6,
                           "M5"=M5.a6,age,car,gotosch,numpeople,race),
        type="response")),
      by=list(a.i,a6)]
    mydt_mc[,c("M1.a6","M2.a6","M3.a6","M4.a6","M5.a6") := NULL]
    setkey(mydt_mc)
    
    
    ## M7 from marginal distribution
    #### first sample M1-M6 given a7
    mydt_mc[is.na(M7.a7),
            c("M1.a7","M2.a7","M3.a7","M4.a7","M5.a7") := list(
              predict.glm(
                fits$M1,
                newdata=data.frame("A"=a7,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M1)$dispersion)),
              predict.glm(
                fits$M2_A,
                newdata=data.frame("A"=a7,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M2_A)$dispersion)),
              predict.glm(
                fits$M3_A,
                newdata=data.frame("A"=a7,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M3_A)$dispersion)),
              predict.glm(
                fits$M4_A,
                newdata=data.frame("A"=a7,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M4_A)$dispersion)),
              predict.glm(
                fits$M5_A,
                newdata=data.frame("A"=a7,age,car,gotosch,numpeople,race),
                type="response")+
                rnorm(n=mc_draws,
                      sd=sqrt(summary(fits$M5_A)$dispersion))),
            by=list(a.i,a7)]
    #### then sample M6 given a7, M1.a7, ..., M5.a7
    mydt_mc[is.na(M7.a7), M6.a7 := rbinom(
      n=mc_draws,size=1,
      prob=predict.glm(
        fits$M6,
        newdata=data.frame("A"=a7,"M1"=M1.a7,"M2"=M2.a7,"M3"=M3.a7,"M4"=M4.a7,
                           "M5"=M5.a7,age,car,gotosch,numpeople,race),
        type="response")),
      by=list(a.i,a7)]
    #### then sample M7 given a7, M1.a7, ..., M6.a7
    mydt_mc[is.na(M7.a7),"M7.a7" := rbinom(
      n=mc_draws,size=1,
      prob=predict.glm(
        fits$M7,
        newdata=data.frame("A"=a7,"M1"=M1.a7,"M2"=M2.a7,"M3"=M3.a7,"M4"=M4.a7,
                           "M5"=M5.a7,"M6"=M6.a7,
                           age,car,gotosch,numpeople,race),
        type="response")),
      by=list(a.i,a7)]
    mydt_mc[,c("M1.a7","M2.a7","M3.a7","M4.a7","M5.a7","M6.a7") := NULL]
    setkey(mydt_mc)
    
    # predicted value of Y
    mydt_mc[, "Y.a" := predict.glm(
      fits$Y, 
      newdata=data.frame("A"=a0,"M1"=M1.a1,"M2"=M2.a2,"M3"=M3.a3,
                         "M4"=M4.a4,"M5"=M5.a5,"M6"=M6.a6,"M7"=M7.a7,
                         age,car,gotosch,numpeople,race),
      type="response")]
    setkey(mydt_mc)
    
    ## average over MC draws
    mydt_mc_means <- mydt_mc[,mean(Y.a),by=list(a.i,a0,a1,a2,a3,a4,a5,a6,a7)]
    setnames(mydt_mc_means,"V1","Y.a")
    setkey(mydt_mc_means)
    # rows of duplicated data to impute
    mydt_a1 <- mydt[,list(a.i,a0,a1,a2,a3,a4,a5,a6,a7)]
    setkey(mydt_a1)
    
    return(merge(mydt_a1,mydt_mc_means,sort=FALSE)[,Y.a])
  }
  
  dat[, "Y.a" := SampleMmarginal(.SD), by=id]
  setkey(dat)
  
  fit_mod2 <- glm(Y.a ~ a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7, 
                  family = gaussian("identity"), 
                  data = dat)
  
  theta <- coef(fit_mod2)
  names(theta) <- paste0("t",lapply(
    strsplit(names(theta),"a"),paste0,collapse=""))
  res[["theta"]] <- theta
  rm(theta,fit_mod2)
  
  res[["ie12"]] <- as.numeric(res$gamma["g1"]-sum(res$theta[paste0("t",1:7)]))
  
  return(res)
}

ptm=proc.time()[3]
OneMCestimator(data=data, mc_draws=1e4) # 1e4 draws: 12 mins
proc.time()[3]-ptm

# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))

# IDs from sampling with replacement
boot_id <- seed
perm_ids <- sample(x=nrow(data),replace=TRUE)
# extract resampled dataset for each resample 
onedat_boot <- data[perm_ids]
onedat_boot[, "id" := 1:nrow(onedat_boot)] # ensure unique IDs
setkey(onedat_boot)
res <- OneMCestimator(data=onedat_boot, mc_draws=1e4)
res <- unlist(res)
res <- c(res,"boot_id"=boot_id)

save(res,file=paste0("data-mma-boots-",seed,".Rdata"))
q()

# results ---------------------------------------------------------------------
setwd("boots/")
myfiles <- list.files()
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
boot_results <- NULL
for (ll in myfiles) {
  load(ll)
  boot_results[[ll]] <- res
  rm(res)
}
res <- data.table(do.call(rbind,boot_results))
setkey(res); rm(boot_results)

# check how many bootstrap samples failed
saved <- sapply(myfiles, function(x) {
  xsplit = strsplit(x,"-")[[1]]
  as.integer(strsplit(xsplit[length(xsplit)],".Rdata")[[1]])
})
(max(saved)-length(saved))/max(saved)

res_summary <- t(rbind(apply(res, 2, sd),
                       apply(res, 2, quantile, probs=c(.025,.975))))
apply(res_summary, 1, function(x) {
  x_ <- round(x,2)
  paste0(x_[1]," & (",x_[2],",",x_[3],")")
})
