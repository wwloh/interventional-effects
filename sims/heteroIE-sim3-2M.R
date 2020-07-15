rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# simulation settings =========================================================
simsettings <- expand.grid(
  "contY"=c(FALSE), # continuous or binary outcome
  "a2"=c(0,0.8), # A --> M2
  "e21"=c(0,0.8), # M1 --> M2
  "b1"=c(0), # M1 --> Y
  "n"=c(5e2,5e4) # population or sample
)
simsettings <- simsettings[simsettings$a2==0 | simsettings$e21==0,]
row.names(simsettings) <- NULL
simsettings

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[
    c(rep(which(simsettings$n<5e4),each=20),which(simsettings$n>5e2)),]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

# specific setting
ss <- seed
(contY <- simsettings[ss,"contY"])
(a2 <- as.numeric(simsettings[ss,"a2"]))
(e21 <- as.numeric(simsettings[ss,"e21"]))
(b1 <- as.numeric(simsettings[ss,"b1"]))
(n <- as.numeric(simsettings[ss,"n"]))
(nsims <- ifelse(n>500,1,50))

source("funs-MC_Wt.R")

OneData <- function() {
  
  # generate population data ==================================================
  L1 <- rnorm(n)
  L2 <- rbinom(n,size=1,prob=0.5)
  ps <- 0.7*L2 # propensity score
  A <- rbinom(n,size=1,prob=exp(ps)/(1+exp(ps)))
  rm(ps)
  a1c <- 1
  M1 <- a1c*A*L2 - L2 + L1 + rnorm(n) # A --> M1
  M2 <- a2*A*M1 + e21*M1 + L1 + L2 + rnorm(n) # M1 --> M2
  b2 <- 1
  Yast <- b1*M1 + b2*M2 - 0.4*M1*M2 + L1 + L2 # M1, M2 --> Y
  if (contY) {
    Y <- Yast + rnorm(n)  
  } else {
    Y <- rbinom(n,size=1,prob=exp(Yast)/(1+exp(Yast)))
  }
  popdata <- data.frame("id"=1:n,L1,L2,A,M1,M2,Y)
  rm(L1,L2,A,M1,M2,Y)
  
  # fit mediator and outcome models to population data ========================
  (m_names <- grep("M",colnames(popdata),value=TRUE))
  p <- length(m_names) # number of mediators
  l_names <- paste0("L",1:2) # confounders
  fit_MY <- list()
  for (varname in c(paste0("M",1:p),"Y")) {
    # exposure group-specific fitted models
    for (aa in 0:1) {
      # conditional on confounders L
      if (varname == "Y") {
        # outcome
        modelformula <- paste0(varname,"~",paste(
          c(paste0("M",1:p),l_names,paste0("M",1:p,":L2"),"M1:M2","M1:M2:L2"),
          collapse="+"))
      } else {
        ## marginal for each mediator unconditional on other mediators
        modelformula <- paste0(varname,"~",paste(l_names,collapse="+"))
      }
      
      # subset is evaluated in the same way as variables in formula
      if (contY || varname!="Y") {
        fit_MY[[paste0(varname,"_A",aa)]] <- 
          glm(formula=as.formula(modelformula),
              family = gaussian("identity"),
              data=popdata,subset=A==aa)
      } else {
        fit_MY[[paste0(varname,"_A",aa)]] <- 
          glm(formula=as.formula(modelformula),
              family = binomial("logit"),
              data=popdata,subset=A==aa)
      }
      rm(modelformula)
        
      # conditional on other mediators (for joint mediator distribution)
      if (varname %in% paste0("M",2:p)) {
        m_idx <- as.integer(strsplit(varname,"M")[[1]][2])
        modelformula <- paste0(varname,"~",paste(c(
          paste(c(paste0("M",1:(m_idx-1)),"L2"),collapse="*"),
          l_names),collapse="+"))
        fit_MY[[paste0(varname,"condMs_A",aa)]] <- 
          glm(formula=as.formula(modelformula),
              family = gaussian("identity"),
              data=popdata,subset=A==aa)
        rm(modelformula)
      }
    }
  }
  # PS model
  fit_A <- glm(paste0("A~",paste(l_names,collapse="+")), 
               family = binomial("logit"), data = popdata)
  
  # misspecified outcome model
  fit_MY_Ymis <- fit_MY
  fit_MY_Ymis[grep("Y",names(fit_MY_Ymis))] <- NULL
  for (aa in 0:1) {
    varname <- "Y"
    modelformula <- paste0(varname,"~",paste(
      c(paste0("M",1:p),l_names,paste0("M",1:p,":L2")),collapse="+"))
    # subset is evaluated in the same way as variables in formula
    if (contY) {
      fit_MY_Ymis[[paste0(varname,"_A",aa)]] <- 
        glm(formula=as.formula(modelformula),
            family = gaussian("identity"),
            data=popdata,subset=A==aa)
    } else {
      fit_MY_Ymis[[paste0(varname,"_A",aa)]] <- 
        glm(formula=as.formula(modelformula),
            family = binomial("logit"),
            data=popdata,subset=A==aa)
    }
    rm(modelformula,varname)
  }
  # lapply(fit_MY, "[[", "formula")
  # lapply(fit_MY_Ymis, "[[", "formula")
  
  # effect model posited ======================================================
  eff_mod <- as.formula(Y.a ~ a1:I(1-J) + a1:L2:I(1-J) + 
                          a2:I(1-J) + a2:L2:I(1-J) + L1 + L2 + 
                          J + a0:J + a0:J:L2 + a1:a2:J + a1:a2:J:L2)
  Effects_from_coefs <- function(regcoef) {
    ie1 <- as.numeric(cumsum(regcoef[c("a1:I(1 - J)","a1:I(1 - J):L2")]))
    ie2 <- as.numeric(cumsum(regcoef[c("I(1 - J):a2","I(1 - J):L2:a2")]))
    ie_jt <- as.numeric(cumsum(regcoef[c("a1:a2:J","a1:L2:a2:J")]))
    ie12 <- ie_jt-ie1-ie2
    ie.out <- c(ie1,ie2,ie12,ie_jt,as.numeric(regcoef["J:a0"]))
    names(ie.out) <- c("ie1","ieC1","ie2","ieC2","ie_mu","ie_Cmu",
                       "ie_jt","ie_Cjt","de")
    return(ie.out)
  }
  
  # observed sample ===========================================================
  onedata <- popdata
  
  res_ie <- list() # for saving results
  
  # weighted estimator
  res <- OneWTestimator(Data=onedata,
                        fit_M=fit_MY[grep("M",names(fit_MY),value=TRUE)],
                        fit_A=fit_A,Lnames=l_names,eff_mod)
  res_ie[["wt"]] <- Effects_from_coefs(res$eff_coef)
  res_weights <- res$weights
  rm(res)
  
  # MC estimator
  res <- OneMCestimator(Data=onedata,mc_draws=100,
                        fit_MY=fit_MY,Lnames=l_names,eff_mod)
  res_ie[["mc"]] <- Effects_from_coefs(res)
  rm(res)
  
  if (n<5e4) {
    res <- OneMCestimator(Data=onedata,mc_draws=100,
                          fit_MY=fit_MY_Ymis,Lnames=l_names,eff_mod)
    res_ie[["mc_ymis"]] <- Effects_from_coefs(res)
  } else {
    res <- rep(NA,length(res_ie[["wt"]]))
    names(res) <- names(res_ie[["wt"]])
    res_ie[["mc_ymis"]] <- res
  }
  rm(res)
  
  return( list(res_ie, res_weights) )
}

ptm=proc.time()[3]
simres <- replicate(n=nsims,OneData(),simplify=FALSE)
proc.time()[3]-ptm
# 30 sec per sim

simres_weights <- lapply(simres, "[[", 2)
simres_weights <- do.call(rbind,simres_weights)
simres_weights <- cbind(simsettings[ss,],simres_weights)
simres <- lapply(simres, "[[", 1)
simres <- do.call(rbind,lapply(simres, unlist))
simres <- cbind(simsettings[ss,],simres)
filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[ss,]))),collapse="_")
save(simres,simres_weights,
     file=paste0("heteroIE-sim3-2M-",filename,"-",seed,".Rdata"))
q()

# results =====================================================================
subfolder <- "sim3/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
resweights_list <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  res_list <- c(res_list,list(simres))
  resweights_list <- c(resweights_list, list(simres_weights))
  rm(simres,simres_weights)
  cat(ll,"\n")
}
simres <- data.table(do.call(rbind,res_list))
setkey(simres)

effects <- c("ie1","ieC1","ie2","ieC2","ie_mu","ie_Cmu","ie_jt","ie_Cjt","de")
names(effects) <- c("IE1","IE1_L2","IE2","IE2_L2","Mutual IE","Mutual IE_L2",
                    "Joint IE","Joint IE_L2","DE")
est_types <- c("Wt","MC","MC (Y mis)")
resout <- list()
for (i in 1:length(effects)) {
  eff <- effects[i]
  # true effects
  true_eff <- simres[n==5e4, .SD, by=list(contY,b1,a2,e21), 
                     .SDcols=grep(eff,names(simres),value=TRUE)[2]]
  setnames(true_eff,ncol(true_eff),"true")
  setkey(true_eff)
  # sample estimates
  resout.eff <- simres[n<5e4, lapply(.SD,function(x) c(mean(x),sd(x))), 
                       by=list(contY,b1,a2,e21), 
                       .SDcols=grep(eff,names(simres),value=TRUE)]
  setnames(resout.eff, 5:ncol(resout.eff), est_types)
  resout.eff[, est := rep(c("pt","se"),times=nrow(resout.eff)/2)]
  resout.eff[, eff := names(eff)]
  setkey(resout.eff)
  resout.eff <- dcast.data.table(resout.eff, 
                                 formula=contY+b1+e21+a2+eff ~ est,
                                 value.var=est_types)
  resout.eff <- merge(resout.eff,true_eff,all.x=TRUE)
  setcolorder(resout.eff,c(1:5,ncol(resout.eff),6:7,10:11,8:9))
  resout[[eff]] <- resout.eff
}
resout <- rbindlist(resout)
library("xtable")
print(xtable(resout[,-c(1:2)]),include.rownames=FALSE)
