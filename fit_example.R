rm(list=ls())

library(tidyverse)
library(lubridate)
library(cvTools)
library(mice)
library(tcltk)
source("functions_simul_mvpbma.R")

#### Data import and specifying variables type
dat <- read_csv("data/example_data.csv")
dat$AGEG <- as.factor(dat$AGEG)

covar_name <- c("FEVER", "MYAL", "HACH",
                "JPAIN", "CHILL", "NSEA", "VOMIT", "DIRHEA", "COUGH",
                "DIZY", "SORE", "CONUL", "JAUN", "SKINR", "RUNNO",
                "MOUSO", "BOIL", "EYEDIS", "SWOTON", "ABPAIN", "CONSTI",
                "NECKSF", "LOSSCO", "LETHA", "ODEMA", "BODYIC", "DYSPH", "RAIN")
covar_cat <- c("PLACE", "AGEG", "GENDA", "BNET")
covar_cont <- c("TEMP", "WBC", "RBC", "HGB", "PLT", "RDW",
                "MPV", "PDW", "LYM", "MON", "GRA")
resp_name <- c("MIC", "URI")

#### Covariates
options(na.action = 'na.pass')
X <- model.matrix(~ ., dat[,c(covar_name, covar_cat, covar_cont)])
X <- X[,-1]
options(na.action = 'na.omit')

id <- which(!colnames(X) %in% covar_cont) # Index of discrete columns
ic <- which(colnames(X) %in% covar_cont) # Index of continuous columns

#### Response
Y <- dat[,resp_name]

#### Model fitting
## Constants and data
params <- list(
  id = id,
  ic = ic,
  Y = as.matrix(Y),
  X = as.matrix(X),
  n = nrow(X),
  p = ncol(X),
  r = ncol(Y)
)

## Start RJMCMC with intercept only model
params$H <- list()
for (i in resp_name) params$H[[i]] <- list(In = 1, Out = c(id, ic) + 1)

## Initialize parameter at some sensible values
params <- within(params, {
  Beta <- rnorm(r*(p+1)) %>%
    matrix(p+1, r)
  Beta[H[[1]]$Out,] <- 0
  nu <- rnorm(p)
  Sigma <- diag(p)
  mu <- nu
  Omega <- Sigma
  Yst <- ifelse(Y, 1, -1)
})

params$Z <- params$X
params$Z[is.na(params$Z)] <- 0
W_i <- ifelse(params$Z, 1, -1)
W_i[,params$ic] <- params$Z[,params$ic]
params$W <- W_i

## MCMC behaviour parameters
skip <- 20
iter <- 10000
burnin <- iter/2
use <- burnin + 1
burnin <- 1:burnin
use <- use:iter
use <- use[use %% skip == 0]

## Storage for MCMC output
storage <- list(mu = matrix(NA, iter, params$p),
                Omega = matrix(NA, iter, params$p^2),
                Beta = matrix(NA, iter, (params$p + 1) * params$r))

## Progress bar
pb <- tkProgressBar(title = "Gibbs Progress", label = "", 0, iter, 0)

## MCMC loops
for (i in 1:iter) {
  setTkProgressBar(pb, i)
  
  params$nu <- update_nu(params)
  params$Sigma <- update_Sigma(params)
  
  tmp <- decomp_Sig(params$Sigma)
  params$V.5 <- (1/diag(tmp$D.5)) %>% diag
  storage$mu[i,] <- params$mu <- params$V.5 %*% params$nu
  storage$Omega[i,] <- params$Omega <- tmp$Omega
  
  tmp <- update_w(params)
  params$W <- tmp$W
  params$Z <- tmp$Z
  params$H <- update_H(params)
  storage$Beta[i,] <- params$Beta <- update_Beta(params)
  params$Yst <- update_yst(params)
}


#### Predictions (Demonstrate using training data as "new" data here)
## Setting up "new" data (X) with intial values for Z
XNew <- X
MNew <- is.na(XNew)
ZNew <- ifelse(XNew, 1, -1)
ZNew[,params$ic] <- XNew[,params$ic]
ZNew[MNew] <- 0

## Creating storage for predictions
preds <- list()
for (l in 1:params$r) preds[[resp_name[l]]] <- matrix(NA, nrow(ZNew), length(use))

## Looping through the posterior samples of parameters
for (i in 1:length(use)) {
  mu <- storage$mu[use[i],]
  Omega <- storage$Omega[use[i],] %>% matrix(length(mu), length(mu))
  lo <- ifelse(XNew, 0, -Inf)
  lo[is.na(lo)] <- -Inf
  hi <- ifelse(XNew, Inf, 0)
  hi[is.na(hi)] <- Inf
  
  ## Cycling through each covariates (impute missing covariates)
  for (p in 1:ncol(ZNew)) {
    Om12 <- Omega[p,-p]
    Om22_Inv <- solve(Omega[-p, -p])
    diffMat <- (ZNew[,-p] - rep(1, nrow(XNew)) %*% t(mu[-p]))
    mu_cond <- mu[p] + as.vector(t(Om12) %*% Om22_Inv %*% t(diffMat))
    Sig_cond <- as.vector(Omega[p, p] - t(Om12) %*% Om22_Inv %*% Om12)
    
    if (p %in% id) {
      ZNew[,p] <- rtnorm(nrow(XNew), mu_cond, sqrt(Sig_cond), lo[,p], hi[,p])
    } else {
      ZNew[,p] <- XNew[,p]
      cond <- is.na(ZNew[,p])
      ZNew[cond,p] <- rnorm(sum(cond), mu_cond[cond], sqrt(Sig_cond))
    }
  }
  
  ZtNew <- ZNew
  ZtNew[,params$id] <- ZtNew[,params$id] > 0
  B <- storage$Beta[use[i],] %>%
    matrix(params$p + 1, params$r)
  pred <- (cbind(1, ZtNew) %*% B) %>% pnorm
  
  for (l in 1:params$r) preds[[l]][,i] <- pred[,l]
}

preds[["MIC"]] # Posterior prediction (probabilities) of "MIC", row = observation, column = replicate
preds[["TACPLAS"]] # Posterior prediction (probabilities) of "TACPLAS"