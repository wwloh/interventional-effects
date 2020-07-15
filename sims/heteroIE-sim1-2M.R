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
  "a2"=c(0,0.4,0.8), # A --> M2
  "e21"=c(0,0.4,0.8), # M1 --> M2
  "b1"=c(0,0.8), # M1 --> Y
  "n"=c(5e2,5e4) # population or sample
)
simsettings <- simsettings[simsettings$a2==0 | simsettings$e21==0,]
rownames(simsettings) <- NULL
simsettings

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[
    c(rep(which(simsettings$n<5e4),each=10),which(simsettings$n>5e2)),]
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
(nsims <- ifelse(n>500,1,100))

source("funs-MC_Wt.R")
source("funs-Wt_M2M1.R") # overwrite with joint density assuming M2 --> M1

OneData <- function() {
  
  # generate population data ==================================================
  L <- rnorm(n)
  ps <- 0.7*L # propensity score
  A <- rbinom(n,size=1,prob=exp(ps)/(1+exp(ps)))
  rm(ps)
  M1 <- A - 2*L + rnorm(n) # A --> M1
  M2 <- a2*A + e21*M1 + L + rnorm(n) # M1 --> M2
  Yast <- b1*M1 + M2 + L # M1, M2 --> Y
  if (contY) {
    Y <- Yast + rnorm(n)  
  } else {
    Y <- rbinom(n,size=1,prob=exp(Yast)/(1+exp(Yast)))
  }
  popdata <- data.frame("id"=1:n,L,A,M1,M2,Y)
  rm(L,A,M1,M2,Y)
  
  # fit mediator and outcome models to population data ========================
  (m_names <- grep("M",colnames(popdata),value=TRUE))
  p <- length(m_names) # number of mediators
  l_names <- "L" # confounders
  fit_MY <- list()
  for (varname in c(paste0("M",1:p),"Y")) {
    # exposure group-specific fitted models
    for (aa in 0:1) {
      # conditional on confounders L
      if (varname == "Y") {
        # outcome
        modelformula <- paste0(varname,"~",paste(c(paste0("M",1:p),l_names),
                                                 collapse="+"))
      } else {
        ## marginal for each mediator unconditional on other mediators
        modelformula <- paste0(varname,"~",paste(l_names,collapse="*"))
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
      if (varname=="M1") {
        fit_MY[[paste0(varname,"condMs_A",aa)]] <- 
          glm(formula=M1~M2+L,
              family = gaussian("identity"),
              data=popdata,subset=A==aa)
      }
    }
  }
  # lapply(fit_MY, "[[", "formula")
  # lapply(fit_MY, "[[", "coefficients")
  # PS model
  fit_A <- glm(paste0("A~",paste(l_names,collapse="*")), 
               family = binomial("logit"), data = popdata)
  
  # effect model posited ======================================================
  eff_mod <- as.formula(Y.a ~ a1:I(1-J) + a2:I(1-J) + J + a0:J + a1:a2:J + L)
  Effects_from_coefs <- function(regcoef) {
    c("ie1"=as.numeric(regcoef["a1:I(1 - J)"]),
      "ie2"=as.numeric(regcoef["I(1 - J):a2"]),
      "ie_jt"=as.numeric(regcoef["a1:a2:J"]),
      "de"=as.numeric(regcoef["J:a0"]))
  }
  
  # observed sample ===========================================================
  onedata <- popdata
  
  res_ie <- list() # for saving results
  
  # weighted estimator
  res <- OneWTestimator(Data=onedata,
                        fit_M=fit_MY[grep("M",names(fit_MY),value=TRUE)],
                        fit_A=fit_A,Lnames="L",eff_mod)
  res_ie[["wt"]] <- Effects_from_coefs(res$eff_coef)
  res_weights <- res$weights
  rm(res)
  
  # MC estimator
  res <- OneMCestimator(Data=onedata,mc_draws=100,
                        fit_MY=fit_MY,Lnames="L",eff_mod)
  res_ie[["mc"]] <- Effects_from_coefs(res)
  rm(res)
  
  return( list(res_ie, res_weights) )
}

ptm=proc.time()[3]
simres <- replicate(n=nsims,OneData(),simplify=FALSE)
proc.time()[3]-ptm
# 12 sec per sim

simres_weights <- lapply(simres, "[[", 2)
simres_weights <- do.call(rbind,simres_weights)
simres_weights <- cbind(simsettings[ss,],simres_weights)
simres <- lapply(simres, "[[", 1)
simres <- do.call(rbind,lapply(simres, unlist))
simres <- cbind(simsettings[ss,],simres)
filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[ss,]))),collapse="_")
save(simres,simres_weights,
     file=paste0("heteroIE-sim1-2M-",filename,"-",seed,".Rdata"))
q()

# results =====================================================================
subfolder <- "sim1/"
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

effects <- c("ie1","ie2","ie_jt","de")
names(effects) <- c("IE1","IE2","Joint IE","DE")
est_types <- c("Wt","MC")
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
  setcolorder(resout.eff,c(1:5,ncol(resout.eff),6:(ncol(resout.eff)-1)))
  resout[[eff]] <- resout.eff
}
resout <- rbindlist(resout)
library("xtable")
print(xtable(resout[,-c(1)]),include.rownames=FALSE)

simwts <- data.table(do.call(rbind,resweights_list))
setkey(simwts)
print(xtable(dcast.data.table(
  simwts[W>0 & n==500, lapply(.SD,sd), 
         by=list(contY,a2,e21,b1,a.i), .SDcols="W"],
  formula=b1+e21+a2~a.i,
  value.var="W")),include.rownames=FALSE)
