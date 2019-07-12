rm(list=ls())
libraries_check <- c("data.table", "mvtnorm", "xtable")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

a0 <- 3:0
a <- .4*(3:6)
aC <- 1
sigma <- matrix(.85, nrow=length(a), ncol=length(a)); diag(sigma) <- 1
sY <- 1
b <- c(1.6,0.4,0.6*c(1:2,1:0),-2)

OneData <- function(n=400) {
  C <- rnorm(n)
  A <- rbinom(n, size = 1, prob = pnorm(C))
  M <- sapply(1:n, function(i) {
    rmvnorm(1,mean=a0+a*A[i]+aC*C[i],sigma=sigma)
  })
  M <- t(M)
  colnames(M) <- paste0("M",1:ncol(M))
  Y <- rnorm(n, mean = cbind(1,A,M,C) %*% b, sd = sY)
  return(data.frame(id = 1:n, C, A, M, Y))
}

true_effects <- list()
true_effects[["theta0"]] <- b[2]
for (tt in 1:length(a)) {
  true_effects[[paste0("theta",tt)]] <- b[tt+2]*a[tt]
}

OnePAestimator <- function(data) {
  # data=OneData()
  res <- list()

  ## fit outcome and mediator models ------------------------------------------
  fitMs <- lapply(1:4, function(s) {
    glm(formula=as.formula(paste0("M",s, "~A+C")), 
        family = gaussian("identity"), data = data)
  })
  fitY <- glm(Y ~ A + M1 + M2 + M3 + M4 + C, family = gaussian("identity"), 
              data = data)
  
  ## path analysis estimators -------------------------------------------------
  ## function of parameters in models fitted to observed data -----------------
  deltas <- lapply(fitMs, function(lm.Ms) coef(lm.Ms)["A"])
  beta <- coef(fitY)
  theta_lm <- beta[paste0("M",1:4)]*unlist(deltas)
  theta_lm <- c(beta["A"],theta_lm)
  names(theta_lm) <- paste0("t",0:(length(theta_lm)-1))
  return(theta_lm)
}

nsims <- 20000
ptm=proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  res_list <- list()
  data=OneData()
  
  res <- OnePAestimator(data)
  res_list <- c(res_list,res)
  rm(res)
  
  return(res_list)
})
proc.time()[3]-ptm

save(simres,file="interventional-sim5-lm-4M.Rdata")
q()

# results ---------------------------------------------------------------------
load("interventional-sim5-lm-4M.Rdata")
res <- rbindlist(simres)
setkey(res)
nrow(res)

xtable(cbind(
  "truth"=true_effects, 
  "mean"=colMeans(res),
  "ese"=apply(res, 2, sd)), digits=3)
