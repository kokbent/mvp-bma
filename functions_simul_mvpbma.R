### DEV MODE
library(truncnorm)
library(Rfast)

## Gibbs Sampler functions
rtnorm <- function(n, m, s, lo, hi) {
  u_lo <- pnorm(lo, m, s)
  u_hi <- pnorm(hi, m, s)
  
  u <- runif(n, u_lo, u_hi)
  tn <- qnorm(u, m, s)
  
  if (any(tn == Inf | tn == -Inf)) {
    tn <- ifelse(tn == Inf, 10, tn)
    tn <- ifelse(tn == -Inf, -10, tn)
  }
  
  return(tn)
}

decomp_Sig <- function (Sigma) {
  D.5 <- diag(sqrt(diag(Sigma)))
  Omega <- solve(D.5) %*% Sigma %*% solve(D.5)
  
  return(list(D.5 = D.5, Omega = Omega))
}

update_w <- function (params) {
  W_upd <- params$W
  Z_upd <- params$Z
  
  for (j in 1:params$p) {
    Z_upd[,j] <- params$X[,j]
    M <- is.na(Z_upd[,j])
    
    Om22 <- params$Omega[-j, -j]
    Om22Inv <- solve(Om22)
    Om12 <- params$Omega[j, -j, drop = F]
    res <- t(params$Z[,-j]) - params$mu[-j]
    eta <- (params$mu[j] + Om12 %*% Om22Inv %*% res) %>%
      as.vector
    psi <- (params$Omega[j, j] - Om12 %*% Om22Inv %*% t(Om12)) %>%
      as.vector
    
    if (j %in% params$id) {
      if (sum(M) > 0) {
        Z_M0 <- Z_M1 <- Z_upd[M,,drop=F]
        Z_M0[,j] <- 0
        Z_M1[,j] <- 1
        
        pn <- 1 - pnorm(0, eta[M], psi)
        prob1 <- with(params, exp(-0.5 * (Yst[M,] - cbind(1, Z_M1) %*% Beta)^2)) %>%
          apply(1, prod, na.rm = T) * pn
        prob0 <- with(params, exp(-0.5 * (Yst[M,] - cbind(1, Z_M0) %*% Beta)^2)) %>%
          apply(1, prod, na.rm = T) * (1 - pn)
        prob <- prob1/(prob0+prob1)
        prob[is.na(prob)] <- 0.5
        Z_upd[M,j] <- (runif(sum(M)) < prob) %>% as.numeric
      }
      
      lo <- ifelse(Z_upd[,j], 0, -Inf)
      hi <- ifelse(Z_upd[,j], Inf, 0)
      
      W_upd[,j] <- with(params, rtnorm(n, eta, psi, lo, hi))
    } else {
      if (sum(M) > 0) {
        A <- params$Beta[j+1,]
        B <- with(params, Yst[M,] - cbind(1, Z_upd[M,-j, drop = F]) %*% Beta[-(j+1),])
        C <- colSums(t(B) * A, na.rm = T)
        D <- (t(A)^2 %*% t(!is.na(B))) %>% as.vector
        
        Mu <- (eta[M] + psi * C) / (1 + psi * D)
        Var <- psi / (1 + psi * D)
        Z_upd[M, j] <- rnorm(sum(M), Mu, sqrt(Var))
      }
      
      W_upd[,j] <- Z_upd[,j]
    }
    # print(paste(max(Mu), Var))
  }
  
  return(list(W = W_upd, Z = Z_upd))
}

update_nu <- function (params) {
  SigmaInv <- params$Sigma %>% solve()
  Sig <- solve(diag(params$p) + params$n * SigmaInv)
  W_bar <- with(params, colMeans(W))
  m <- with(params, n * Sig %*% SigmaInv %*% W_bar)
  
  return(rmvnorm(1, m, Sig) %>% as.vector)
}

update_Sigma <- function (params) {
  resid_W <- t(params$W) - params$nu
  resid_mat <- apply(resid_W, 2, function (x) x %*% t(x)) %>%
    rowSums %>%
    matrix(params$p, params$p)
  S <- solve(diag(params$p) + resid_mat)
  new_Sigma <- rWishart(1, params$n + params$p + 2, S)[,,1] %>% solve()
  return(new_Sigma)
}

update_yst <- function (params) {
  Yst_new <- params$Yst
  lo <- with(params, ifelse(Y, 0, -Inf))
  hi <- with(params, ifelse(Y, Inf, 0))
  
  Mu <- with(params, cbind(1, Z) %*% Beta)
  Sig <- 1
  
  for (j in 1:ncol(Mu)) {
    cond <- !is.na(params$Yst[,j])
    Yst_new[cond, j] <- rtnorm(sum(cond), Mu[cond, j], Sig, 
                          lo[cond, j], hi[cond, j])
  }
  
  return(Yst_new)
}

update_H <- function (params, prob = c(1/3, 1/3, 1/3), upperp = params$p) {
  H_new <- params$H
  
  for (j in 1:length(params$H)) {
    Hin_old <- params$H[[j]]$In
    pin <- length(Hin_old)
    rand1 <- runif(1)
    p0 <- 1
    if (pin == 1) {
      Hin_new <- with(params, birth(H[[j]]$In, H[[j]]$Out))
      move <- "B"
      p0 <- 1/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1. 
    }
    
    if (pin == upperp) {
      Hin_new <- death(params$H[[j]]$In)
      move <- "D"
      p0 <- 1/3
    }
    
    if (1 < pin & pin < upperp) {
      if (rand1 < 1/3) {
        Hin_new <- with(params, birth(H[[j]]$In, H[[j]]$Out))
        move <- "B"
      }
      if (1/3 <= rand1 & rand1 < 2/3) {
        Hin_new <- death(params$H[[j]]$In)
        move="D"
        if (pin==2) p0=3 #birth prob from 1 -> 2 is 1 and death prob from 2 -> 1 is 1/3
      }
      if (2/3 <= rand1) {
        Hin_new <- swap(params$H[[j]]$In, params$H[[j]]$Out)
        move="S"
      }
    }
    
    pold <- with(params, log_marg_l(cbind(1, Z)[,Hin_old, drop=F], Yst[,j]))
    pnew <- with(params, log_marg_l(cbind(1, Z)[,Hin_new, drop=F], Yst[,j]))+log(p0)
    prob <- exp(pnew-pold)
    rand2 <- runif(1)
    
    if (rand2<prob) {
      Allind <- with(params, c(H[[j]]$In, H[[j]]$Out) %>% sort)
      Hin <- Hin_new
      Hout <- Allind[!Allind %in% Hin]
      H_new[[j]] <- list(In = Hin, Out = Hout)
    } 
  }
  
  return(H_new)
}

log_marg_l <- function(cov, yst, lambda = 1){
  cond <- !is.na(yst)
  yst <- yst[cond]
  cov <- cov[cond,,drop = F]
  
  pp <- ncol(cov)
  
  Tinv <- diag(x=1/lambda, pp)
  Tinv[1,1] <- 1/10
  prec <- t(cov) %*% cov + Tinv
  var1 <- solve(prec)
  mu <- var1 %*% t(cov) %*% yst
  
  diag1 <- c(10, rep(lambda, pp-1))
  
  (1/2)*(-sum(log(diag1))-t(yst)%*%yst+t(mu)%*%prec%*%mu+determinant(var1)$modulus[1])
}

death <- function(indin){
  if (length(indin) == 2) return(indin[-2])
  
  k <- sample(2:length(indin),size=1) 
  if (indin[k] == 1) stop("Something wrong: Proposing death of intercept")
  return(indin[-k])
}

swap <- function(indin, indout){
  if (length(indin) == 2) k_out <- 2
  if (length(indin) > 2) k_out <- sample(2:length(indin), size=1)  
  tmp <- indin[-k_out]
  if (length(indout) == 1) k_in <- indout[1]
  if (length(indout) > 1) k_in <- sample(indout, size=1)
  # print(paste("In:", paste(indin, collapse = " ")))
  # print(paste("Out:", paste(indout, collapse = " ")))
  # if ((k_in %in% indin)) {
  #   print(paste("in",indin))
  #   print(paste("out",indout))
  #   stop("Having trouble in swap")
  # }
  if (indin[k_out] == 1) stop("Something wrong: Proposing swap of intercept")
  sort(c(tmp, k_in))
}

birth <- function(indin, indout){
  k <- sample(indout, size=1)
  sort(c(indin,k))
}

update_Beta <- function (params) {
  Beta_new <- with(params, matrix(NA, p+1, r))
  for (j in 1:ncol(params$Beta)) {
    Z1 <- cbind(1, params$Z)[,params$H[[j]]$In, drop = F]
    Tau <- diag(1, length(params$H[[j]]$In))
    Tau[1, 1] <- 1/10
    
    yst <- params$Yst[,j]
    cond <- !is.na(yst)
    yst <- yst[cond]
    Z1 <- Z1[cond,]
    
    Var <- solve(t(Z1) %*% Z1 + Tau)
    Mu <- Var %*% t(Z1) %*% yst
    beta_new <- rep(0, params$p + 1)
    beta_new[params$H[[j]]$In] <- rmvnorm(1, Mu, Var) %>% as.vector
    Beta_new[,j] <- beta_new
  }
  
  return(Beta_new)
}
