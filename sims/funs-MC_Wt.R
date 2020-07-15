# helper function to create duplicated data for one individual
Dupdata <- function(t, Ai, dup_type="MC") {
  # for indirect effects via each mediator
  out <- do.call(rbind,lapply(1:t, function(s) {
    a <- rep(0L,t)
    a[1:s] <- 1L
    return(a)
  }))
  out <- rbind(0,out)  # first row all 0s
  out <- cbind(0,out) # first column all zeroes for a0
  # for direct, joint, and mutual effects
  if (dup_type=="MC") {
    out <- rbind(out, c(1-Ai, rep(Ai,t)), rep(Ai,t+1))
  } else {
    out <- rbind(out, c(Ai, rep(1-Ai,t)), rep(Ai,t+1))
  }
  colnames(out) <- paste0("a",0:t)
  out <- cbind("a.i"=1:nrow(out),out,"J"=c(rep(0,t+1),rep(1,2)))
  return(out)
}

OneMCestimator <- function(Data,mc_draws=1e2,fit_MY,Lnames,eff_mod) {
  # helper function to predict Y for different outcome models
  PredictY <- function(onedat) {
    setnames(onedat,"a0","A")
    Ya <- rep(NA, nrow(onedat))
    # exposure group-specific fitted models
    for (aa in 0:1) {
      fitY_Aa <- fit_MY[[paste0("Y","_A",aa)]]
      # relevant rows in duplicated data
      indexY_Aa <- as.vector(onedat$A==aa)
      Ya[indexY_Aa] <- predict.glm(
        object=fitY_Aa,newdata=data.frame(onedat[indexY_Aa,]),type="response")
    }
    return(Ya)
  }
  
  # mediator column names
  Mnames <- grep("M",colnames(Data),value=TRUE)
  t <- length(Mnames) # number of distinct mediators
  
  # duplicated data for each individual =======================================
  dat <- data.table(Data)
  setkey(dat)
  alevels <- dat[,as.data.table(Dupdata(t=t,Ai=A,dup_type="MC")),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  SampleMs <- function(mydt,av_mc) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    ## only duplicate rows where J==0
    mydt_mc <- mydt[c(
      rep(which(mydt$J==0),each=mc_draws),which(mydt$J==1))]
    mydt_mc <- cbind(mydt_mc,"mc"=c(
      rep(1:mc_draws,time=sum(mydt$J==0)),rep(1,sum(mydt$J==1))))
    setkey(mydt_mc)
    
    # initialize counterfactual mediator draws
    ## observed values for J==1
    mydt_mc[J==1, paste0(Mnames,".a") := mydt_mc[J==1, Mnames, with=FALSE]]
    setkey(mydt_mc)
    
    # newdata for imputation: same observed L for all duplicated rows
    newdataMs <- data.frame(mydt_mc[1,Lnames,with=FALSE])
    
    # marginal distributions of mediators
    for (s in 1:t) {
      # exposure group-specific fitted models
      for (aa in 0:1) {
        fitMs_Aa <- fit_MY[[paste0("M",s,"_A",aa)]]
        # relevant rows in duplicated data
        s_aa <- mydt_mc[, J==0] & 
          (as.vector(unlist(mydt_mc[, paste0("a",s), with=FALSE]))==aa)
        sampMs_Aa <- sum(s_aa)
        # predicted counterfactual mediator
        meanMs_Aa <- predict.glm(
          object=fitMs_Aa,newdata=newdataMs,type="response")
        if (fitMs_Aa$family$family=="gaussian") {
          # normal distribution
          drawnMs_Aa <- rnorm(n=sampMs_Aa,mean=meanMs_Aa,
                              sd=sqrt(summary(fitMs_Aa)$dispersion))
        } else if (fitMs_Aa$family$family=="binomial") {
          # binomial distribution
          drawnMs_Aa <- rbinom(n=sampMs_Aa,size=1,prob=meanMs_Aa)
        } else if (fitMs_Aa$family$family=="poisson") {
          # Poisson distribution
          drawnMs_Aa <- rpois(n=sampMs_Aa,lambda=meanMs_Aa)
        }
        mydt_mc[s_aa, paste0("M",s,".a") := drawnMs_Aa]
        rm(drawnMs_Aa)
      }
    }
    setkey(mydt_mc)
    
    ## predict potential outcomes using sampled mediator values
    mydt_mc[, (Mnames) := NULL]
    setnames(mydt_mc,paste0(Mnames,".a"),Mnames)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(
      onedat=mydt_mc[, c("a0",Lnames,Mnames), with=FALSE])]
    setkey(mydt_mc)
    mydt_mc[J==1 & a0==A & a1==A, Y.a := Y] # set to observed outcome
    
    if (av_mc==TRUE) {
      ## average over MC draws
      mydt_mc_means <- mydt_mc[,lapply(.SD, mean),
                               by=c("a.i",paste0("a",0:t),"J"),.SDcols="Y.a"]
      setkey(mydt_mc_means)
      return(mydt_mc_means[,Y.a])
    } else {
      return(mydt_mc)  
    }
  }
  
  # average potential outcomes for each individual only
  dat[, "Y.a" := SampleMs(.SD,av_mc=TRUE), by=id]
  setkey(dat)
  
  if (fit_MY$Y_A0$family$family=="gaussian") {
    fit.eff <- lm(eff_mod, data = dat)
  } else if (fit_MY$Y_A0$family$family=="binomial") {
    fit.eff <- glm(eff_mod, family = binomial("logit"), data = dat)
  }
  
  return( coef(fit.eff) )
}


OneWTestimator <- function(Data,fit_M,fit_A,Lnames,eff_mod) {
  # mediator column names
  Mnames <- grep("M",colnames(Data),value=TRUE)
  t <- length(Mnames) # number of distinct mediators
  
  # duplicated data for each individual =======================================
  dat <- data.table(Data)
  setkey(dat)
  # IPW weights for observed treatment
  dat[,"pA1" := predict.glm(fit_A,newdata=dat,type="response")]
  dat[,"w.A" := A/pA1 + (1-A)/(1-pA1)]
  dat[,pA1 := NULL]
  setkey(dat)
  
  alevels <- dat[,as.data.table(Dupdata(t=t,Ai=A,dup_type="Wt")),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  # initialize weights
  dat[, paste0("w.",Mnames) := NA*1.0]
  setkey(dat)
  for (s in 1:t) {
    # exposure group-specific fitted models
    for (aa in 0:1) {
      # different mediator densities 
      for (mtype in c("marginal","joint","observed")) {
        if (mtype=="marginal") {
          # marginal mediator densities =======================================
          fitMs_Aa <- fit_M[[paste0("M",s,"_A",aa)]] # conditional on L
          # relevant rows in duplicated data
          s_aa <- dat[, J==0] & (dat[, paste0("a",s), with=FALSE]==aa)[,1]
        } else {
          # joint mediator densities ==========================================
          if (s==1) {
            fitMs_Aa <- fit_M[[paste0("M1_A",aa)]]
          } else {
            fitMs_Aa <- fit_M[[paste0("M",s,"condMs_A",aa)]]
          }
          if (mtype=="joint") {
            # relevant rows in duplicated data
            s_aa <- dat[, J==1] & (dat[, paste0("a",s), with=FALSE]==aa)[,1]  
          } else if (mtype=="observed") {
            # observed mediators
            s_aa <- (dat[, A]==aa)
          }
        }
        if (fitMs_Aa$family$family=="gaussian") {
          # normal distribution
          drawnMs_Aa <- dnorm(
            x=unlist(dat[s_aa,paste0("M",s),with=FALSE]),
            mean=predict.glm(
              object=fitMs_Aa,
              newdata=data.frame(dat[s_aa,c(Mnames,Lnames),with=FALSE]),
              type="response"),
            sd=sqrt(summary(fitMs_Aa)$dispersion))
        }
        
        if (mtype=="marginal" || mtype=="joint") {
          dat[s_aa, paste0("w.M",s) := drawnMs_Aa]
        } else {
          dat[s_aa, paste0("w.M",s,".obs") := drawnMs_Aa]
        }
        rm(s_aa,fitMs_Aa,drawnMs_Aa)
      }
    }
  }
  setkey(dat)
  
  dat[, W := (a0==A)*w.A*(w.M1*w.M2)/(w.M1.obs*w.M2.obs)]
  setkey(dat)
  setnames(dat,old="Y",new="Y.a")
  setkey(dat)

  binaryY <- all((0 <= dat$Y.a) & (dat$Y.a <= 1))
  if (!binaryY) {
    # weights evaluated in the same way as variables in formula: first in data
    fit.eff <- lm(eff_mod, data=dat, weights=W)
  } else if (binaryY) {
    fit.eff <- glm(eff_mod, family = binomial("logit"), data = dat, weights=W)
  }
  return( list("eff_coef"=coef(fit.eff),
               "weights"=dat[,c(grep("a[.]i",names(dat))+(0:(t+2)),ncol(dat)),
                             with=FALSE]) )
}
