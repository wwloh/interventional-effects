rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# simulation settings =========================================================
simsettings <- expand.grid(
  "contY"=c(TRUE), # continuous or binary outcome
  "a2"=c(0,1.6), # A --> M2
  "e21"=c(0,1.6), # M1 --> M2
  "b1"=c(0,1.6) # M1 --> Y
)
simsettings <- simsettings[simsettings$a2==0 | simsettings$e21==0,]
row.names(simsettings) <- NULL
nsims <- 10
simsettings

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[rep(1:nrow(simsettings),each=1000/nsims),]
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
(n <- 100)

source("funs-MC_Wt.R")

OneDataset <- function() {
  
  # generate population data ==================================================
  L1 <- rnorm(n)
  L2 <- rbinom(n,size=1,prob=0.5)
  ps <- 0.7*L2 # propensity score
  A <- rbinom(n,size=1,prob=exp(ps)/(1+exp(ps)))
  rm(ps)
  a1c <- 1
  M1 <- a1c*A*L2 - L2 + L1 + rnorm(n) # A --> M1
  M2 <- a2*A + e21*M1 + L1 + L2 + rnorm(n) # M1 --> M2
  b2 <- 1
  Yast <- b1*M1 + b2*M2 + L1 + L2 # M1, M2 --> Y
  if (contY) {
    Y <- Yast + rnorm(n)  
  } else {
    Y <- rbinom(n,size=1,prob=exp(Yast)/(1+exp(Yast)))
  }
  popdata <- data.frame("id"=1:n,L1,L2,A,M1,M2,Y)
  return(popdata)
}

OnePtEst <- function(popdata) {
  
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
          c(paste0("M",1:p),l_names,paste0("M",1:p,":L2")),collapse="+"))
      } else {
        ## marginal for each mediator unconditional on other mediators
        modelformula <- paste0(varname,"~",paste(l_names,collapse="+"))
      }
      fit_MY[[paste0(varname,"_A",aa)]] <- 
        glm(formula=as.formula(modelformula),
            family = gaussian("identity"),
            data=popdata,subset=A==aa)
      # subset is evaluated in the same way as variables in formula
      rm(modelformula)
        
      # conditional on other mediators (for joint mediator distribution)
      if (varname %in% paste0("M",2:p)) {
        modelformula <- paste0(varname,"~",paste(
          c(paste0("M",1:(p-1)),l_names,paste0("M",1:(p-1),":L2")),collapse="+"))
        fit_MY[[paste0(varname,"condMs_A",aa)]] <- 
          glm(formula=as.formula(modelformula),
              family = gaussian("identity"),
              data=popdata,subset=A==aa)
        rm(modelformula)  
      }
    }
  }
  # lapply(fit_MY, "[[", "formula")
  # lapply(fit_MY, "[[", "coefficients")
  # PS model
  fit_A <- glm(paste0("A~",paste(l_names,collapse="+")), 
               family = binomial("logit"), data = popdata)
  
  # effect model posited ======================================================
  eff_mod <- as.formula(Y.a ~ a1 + a1:L2 + a2 + a2:L2 + L1 + L2 + J + 
                          a0:J + a1:a2:J + a1:a2:J:L2 + J:L1 + J:L2 + a0:J)
  Effects_from_coefs <- function(regcoef) {
    ie1 <- cumsum(regcoef[c("a1","a1:L2")])
    ie2 <- cumsum(regcoef[c("a2","L2:a2")])
    ie12 <- cumsum(regcoef[c("a1:a2:J","a1:L2:a2:J")])
    ie_jt <- ie1+ie2+ie12
    ie.out <- c(ie1,ie2,ie_jt,as.numeric(regcoef["J:a0"]))
    names(ie.out) <- c("ie1","ieC1","ie2","ieC2","ie_jt","ie_Cjt","de")
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
  res <- OneMCestimator(Data=onedata,mc_draws=5,
                        fit_MY=fit_MY,Lnames=l_names,eff_mod)
  res_ie[["mc"]] <- Effects_from_coefs(res)
  rm(res)
  
  # # product-of-coef
  # ie.poc.betas <- coef(lm(Y~A+M1+M2+L1+L2+M1:L2++M2:L2,data=onedata))[
  #   c("M1","M1:L2","M2","M2:L2","A")]
  # ie.poc.alphas <- c(coef(lm(M1~A*L2+L1,data=onedata))[c("A","A:L2")],
  #                    coef(lm(M2~A*L2+L2,data=onedata))[c("A","A:L2")])
  # ie.poc <- c(
  #   ie.poc.betas["M1"]*ie.poc.alphas["A"],
  #   sum(ie.poc.betas[c("M1","M1:L2")])*sum(ie.poc.alphas[1:2]),
  #   ie.poc.betas["M2"]*ie.poc.alphas["A"],
  #   sum(ie.poc.betas[c("M2","M2:L2")])*sum(ie.poc.alphas[3:4]))
  # ie.poc <- c(ie.poc, sum(ie.poc[c(1,3)]), sum(ie.poc[c(2,4)]),
  #             ie.poc.betas["A"])
  # names(ie.poc) <- names(res_ie[["wt"]])
  # res_ie[["poc"]] <-  ie.poc
  
  return( list(res_ie, res_weights) )
}

OneData <- function() {
  # generate observed data ====================================================
  Data <- OneDataset()
  # point estimates ===========================================================
  res <- OnePtEst(popdata=Data)
  # bootstrap CIs =============================================================
  res_boot <- sapply(1:n, function(idx) {
    bootdat <- Data[sample(n,n,TRUE),]
    bootdat$id <- 1:nrow(bootdat)
    unlist(OnePtEst(popdata=bootdat)[[1]])
  })
  return( c(res,list(res_boot)) )
}

ptm=proc.time()[3]
simres <- replicate(n=nsims,OneData(),simplify=FALSE)
proc.time()[3]-ptm
# 5 mins per sim

warnings() # print warnings

simres_weights <- lapply(simres, "[[", 2)
simres_weights <- do.call(rbind,simres_weights)
simres_weights <- cbind(simsettings[ss,],simres_weights)
simres_boot <- lapply(simres, "[[", 3)
simres_boot <- lapply(simres_boot, function(cis)
  cbind(simsettings[ss,],"eff"=rownames(cis),data.frame(cis)))
simres <- lapply(simres, "[[", 1)
simres <- do.call(rbind,lapply(simres, unlist))
simres <- cbind(simsettings[ss,],simres)

filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[ss,]))),collapse="_")
save(simres,simres_weights,simres_boot,
     file=paste0("heteroIE-sim2-2M-",filename,"-",seed,".Rdata"))
q()

# results =====================================================================
subfolder <- "sim2/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
resweights_list <- list()
resboot_list <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  res_list <- c(res_list,list(simres))
  resweights_list <- c(resweights_list, list(simres_weights))
  resboot_list <- c(resboot_list, simres_boot)
  rm(simres,simres_weights,simres_boot)
  cat(ll,"\n")
}
simres <- data.table(do.call(rbind,res_list))
setkey(simres)
simres_boot <- data.table(do.call(rbind,resboot_list))
setkey(simres_boot)
setnames(simres_boot,"eff","est.type")

effects <- c("ie1","ieC1","ie2","ieC2","ie_jt","ie_Cjt","de")
names(effects) <- c("IE1","IE1_L2","IE2","IE2_L2","Joint IE","Joint IE_L2","DE")
est_types <- c("Wt","MC")
resout <- list()
library("coxed")
for (i in 1:length(effects)) {
  eff <- effects[i]
  # point estimates ===========================================================
  resout.eff <- simres[, lapply(.SD,function(x) c(mean(x),sd(x))), 
                       by=list(contY,b1,a2,e21), 
                       .SDcols=grep(eff,names(simres),value=TRUE)]
  setnames(resout.eff, 5:ncol(resout.eff), est_types)
  resout.eff[, est := rep(c("pt","se"),times=nrow(resout.eff)/2)]
  resout.eff[, eff := names(eff)]
  setkey(resout.eff)
  resout.eff <- dcast.data.table(resout.eff, 
                                 formula=contY+b1+e21+a2+eff ~ est,
                                 value.var=est_types)
  resout.eff[,"true" := 0]
  if (eff=="ieC1") {
    resout.eff[,"true" := b1]
  } else if (eff=="ie2") {
    resout.eff[,"true" := a2]
  } else if (eff=="ieC2") {
    resout.eff[,"true" := a2+e21]
  } else if (eff=="ie_jt") {
    resout.eff[,"true" := a2]
  } else if (eff=="ie_Cjt") {
    resout.eff[,"true" := b1+a2+e21]
  }
  setcolorder(resout.eff,c(1:5,ncol(resout.eff),6:(ncol(resout.eff)-1)))
  
  # 95% CIs ===================================================================
  resout.ci <- simres_boot[est.type %in% grep(eff,names(simres),value=TRUE)]
  # estimator: IW or MC
  resout.ci$est.type <- unlist(lapply(
    strsplit(as.character(resout.ci$est.type),"[.]"),"[",1))
  # calculate 95% bca CIs
  resout.ci <- cbind(
    resout.ci[, list(contY,b1,a2,e21,est.type)],
    t(apply(resout.ci[,grep("X",names(resout.ci),value=TRUE),with=FALSE], 1, 
            function(boots) quantile(boots,probs=c(.025,.975),na.rm=TRUE)))
  )
  setnames(resout.ci,c("2.5%","97.5%"),c("l","u"))
  # insert true values
  resout.ci <- merge(resout.ci,resout.eff[,1:6],by=c("contY","b1","a2","e21"))
  setkey(resout.ci)
  # calculate coverage
  resout.ci[, cover := (l <= true) & (u >= true)]
  resout.ci <- resout.ci[,mean(cover), by=list(contY,b1,a2,e21,est.type)]
  resout.ci <- dcast.data.table(resout.ci, 
                                 formula=contY+b1+e21+a2 ~ est.type,
                                 value.var="V1")
  setcolorder(resout.ci,c(1:(ncol(resout.ci)-2),ncol(resout.ci)-(0:1)))
  resout[[eff]] <- merge(resout.eff,resout.ci)
  rm(eff)
}
resout <- rbindlist(resout)
library("xtable")
print(xtable(resout[,-c(1)]),include.rownames=FALSE)
