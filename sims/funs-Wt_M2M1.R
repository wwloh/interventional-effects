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
          if (s==2) {
            fitMs_Aa <- fit_M[[paste0("M2_A",aa)]]
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
               "weights"=dat[,list(a.i,a0,a1,a2,J,W)]) )
}
