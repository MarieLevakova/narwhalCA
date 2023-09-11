# Wrapper functions for estimating VECM model
# Output here is more user friendly and usable.


# Wrapper function for the Johansen cpp function to fit a VECM model
johansen <- function(Y = NULL, Z = NULL, Z2 = NULL, r = 1, dt = 0.1){
  N <- nrow(Y)
  p <- ncol(Y)

  intercept <- TRUE
  if(is.null(Z2)){
    intercept <- FALSE
    Z2 <- matrix(1, nrow = N, ncol = 1) # "Fake" Z2 to match argument type, not used anywhere later
  }

  out <- johansenCpp(Y, Z, Z2, r, dt, intercept)

  p2 <- dim(Z2)[2]

  if(r > 0){
    a.hat = matrix(out[1:r,], nr = p, byrow = TRUE)
    b.hat = matrix(out[(r+1):(2*r),], nr = p, byrow = TRUE)
    rownames(a.hat) = paste0("x",1:p)
    colnames(a.hat) = paste0("r",1:r)
    rownames(b.hat) = paste0("x",1:p)
    colnames(b.hat) = paste0("r",1:r)
  }
  P.hat = out[(2*r+1):(2*r+p2),]
  O.hat = out[(2*r+p2+1):(2*r+p2+p),]
  test  = out[(2*r+p2+p+1),]
  eigs  = out[(2*r+p2+p+2),]
  res   = out[(2*r+p2+p+3):(2*r+p2+p+N+2),]

  res0  = Y - matrix(apply(Y, 2, mean), ncol = p, nrow = N, byrow = T)
  R2    = 1 - sum(res^2)/sum(res0^2)

  names(P.hat)    = paste0("x",1:p)
  rownames(O.hat) = paste0("x",1:p)
  colnames(O.hat) = paste0("x",1:p)
  names(test)     = paste0("r=",0:(p-1))
  names(eigs)     = paste0("l",1:p)
  colnames(res)   = paste0("x",1:p)

  if(r==0){
    a.hat = b.hat = rep(0, p)
  }
  return(list(N = N, p = p, r = r, alpha = a.hat, beta = b.hat,
              Psi = P.hat, Omega = O.hat,
              test = test, lambda = eigs, dt = dt, r2 = R2))
}

# A function that fits the VECM model and inserts the missing and reference channel into
# the fitted model
johansen.full <- function(Y = NULL, Z = NULL, Z2 = NULL, r = 1, dt = 0.1, intercept = TRUE,
                          reference = 64, excluded = NULL){
  p <- ncol(Y)
  N <- nrow(Y)

  channel.sel <- (1:p)[-c(excluded, reference)] # channels used to fit the model
  channel.sel.Z2 <- c()
  excluded.Z2 <- c()
  if(!is.null(Z2)){
    k <- ncol(Z2) %/% p
    for(ii in 0:(k-1)){
      channel.sel.Z2 <- c(channel.sel.Z2, 1+ii*p+channel.sel)
      excluded.Z2 <- c(excluded.Z2, 1+ii*p+excluded)
    }
  }

  # Calculation of centered variables
  Ystd <- Y
  meanY <- apply(Y, 2, mean)
  Ystd <- Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))
  meanZ <- apply(Z, 2, mean)
  Zstd <- Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))

  # model <- johansen(Y[,channel.sel], Z[,channel.sel], r = r,
  #                   dt = dt, intercept = intercept)

  out <- johansenCpp(Y[,channel.sel], Z[,channel.sel], Z2 =  Z2[,channel.sel.Z2],
                     r = r, dt = dt,
                     intercept = T)

  # Extract beta and add the row for the reference channel and create a new set of stationary predictors
  BETA.p <- t(out[(r+1):(2*r),])
  BETA.p.mean <- length(channel.sel)/(length(channel.sel)+1)*apply(BETA.p, 2, mean)
  BETA.full <- rbind(BETA.p, -BETA.p.mean) -
    matrix(c(rep(1, length(channel.sel)), 0), ncol = 1) %*% t(BETA.p.mean)
  if(length(excluded>=1)){
    Z.r <- matrix(Zstd[, -excluded] %*% BETA.full, nrow = N, ncol = r)
    Z2.prep <- Z2[,-excluded.Z2]
  } else {
    Z.r <- matrix(Zstd %*% BETA.full, nrow = N, ncol = r)
    Z2.prep <- Z2
  }

  all.param  <- t(coef(lm(Ystd ~ Z.r + Z2.prep - 1)))
  ALPHA <- all.param[,1:r]
  MU <- matrix(all.param[,r+1], ncol = 1)
  GAMMAS <- all.param[,(r+2):ncol(all.param)]

  if(length(excluded)>=1){
    ALPHA[excluded,] <- 0
  }

  BETA.final <- matrix(0, nrow = p, ncol = r)

  if(length(excluded)>=1){
    BETA.final[-excluded,] <- BETA.full
  }

  PI <- ALPHA %*% t(BETA.final)

  # MU <- meanY - PI %*% as.matrix(meanZ, ncol = 1)

  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(MU) - Z %*% t(PI) - Z2[,-1] %*% t(GAMMAS)
  if(length(excluded)>=1){
    res[,excluded] <- NA
  }
  OMEGA <- (t(res) %*% res)/N

  res0  = Y - matrix(meanY, nrow = N, ncol = p, byrow = T)
  if(length(excluded)>=1){
    res0[,excluded] <- NA
  }
  R2    = 1 - sum(res^2, na.rm = T)/sum(res0^2, na.rm = T)

  if(length(excluded)>=1){
    MU[excluded,] <- NA
  }

  test  = out[(2*r+2+length(channel.sel)),]
  eigs  = out[2*r+2+length(channel.sel+1),]

  return(list(N = N, p = p, r = r,
         alpha = ALPHA/dt, beta = BETA.final, Pi = PI/dt, Omega = OMEGA/dt,
         Psi = MU/dt, test = test, lambda = eigs, dt = dt, r2 = R2))
}

# Johansen procedure with penalization of alpha
pen.alpha <- function(Y, Z, Z2 = NULL, r, dt = 1, equal.penalty = F, n.penalty = 100,
                      n.cv = 10, lasso.prop = 1, excluded = NULL, reference = ncol(Y)){

  ## INPUT
  # Y: differenced time series
  # Z: lagged time series
  # Z2: lagged differences
  # r: cointegration rank
  # dt: timestep
  # equal.penalty: logical value -  should the penalty parameter be the same for all rows of alpha?
  # n.penalty: number of penalty values to use
  # n.cv: number of repetitions of the crossvalidation procedure
  # lasso.prop: the proportion of the lasso penalty vs. the ridge penalty
  #             (the default is 1, meaning the penalty is purely of lasso type)
  # excluded: numbers of excluded channels in case of corrupted signals
  # reference: the reference channel

  ## OUTPUT
  # N: sample size
  # p: dimension of the multivariate system
  # r: cointegration rank
  # alpha: estimate of adjustment coefficients
  # beta: estimate of cointegrating vectors
  # Pi: estimate of Pi
  # Omega: estimate of covariance matrix
  # Psi: estimate of the intercept
  # test: trace test statistics
  # lambda: eigenvalues
  # r2: coefficient of determination
  # penalty: the value of the penalty parameter used in the final model
  # lasso.prop: the proportion of the lasso penalty in the final fit

  ## START CODE
  N <- dim(Y)[1]
  p <- dim(Y)[2]

  channel.sel.Z2 <- c()
  excluded.Z2 <- c()
  if(!is.null(Z2)){
    k <- ncol(Z2) %/% p
    for(ii in 0:(k-1)){
      channel.sel.Z2 <- c(channel.sel.Z2, ii*p+channel.sel)
      excluded.Z2 <- c(excluded.Z2, ii*p+excluded)
    }
  }

  # Centering variables, to remove the effect of intercept
  Ystd <- Y
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))

  meanZ2 <- apply(Z2[,-1], 2, mean)
  sdZ2 <- apply(Z2[,-1], 2, sd)
  Z2std <- Z[,-1] - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ2))

  # Remove excluded channels + one extra (reference) channel because of the linear dependence
  Ystd.p <- Y[,-c(excluded, reference)]
  Zstd.p <- Z[,-c(excluded, reference)]
  Z2std.p <- Z2[, channel.sel.Z2]

  # Fit by the standard Johansen procedure
  fit0 <- johansen(Y = Ystd.p, Z = Zstd.p, Z2 = Z2std.p,
                   r = r, dt = dt, intercept = F)

  # Extract beta, add the row for the reference channel and create a new set of stationary predictors
  BETA.p <- fit0$beta
  BETA.p.mean <- (p-length(excluded)-1)/(p-length(excluded))*apply(BETA.p, 2, mean)
  BETA.full <- rbind(BETA.p, -BETA.p.mean) -
    matrix(c(rep(1, p-length(excluded)-1), 0), ncol = 1) %*% t(BETA.p.mean)
  Z.r <- matrix(Zstd[, -excluded] %*% BETA.full, nrow = N, ncol = r)

  # Estimate alpha with LASSO penalty
  ALPHA.Sparse <- matrix(0, nrow = p, ncol = r)
  if((r == 1) & is.null(Z2)){
    ALPHA.Sparse  <- matrix(coef(lm(Ystd ~ Z.r - 1)), ncol = 1)
    ALPHA.Sparse[excluded,] <- NA
  } else {

    if(equal.penalty){
      penalty.opt <- cv.opt <- rep(NA, length(lasso.prop))
      for(i in 1:length(lasso.prop)){
        lasso.prop.i <- lasso.prop[i]
        # Determine the sequence of tuning parameters
        penalty.max <- 0
        penalty.min <- NA
        for(i.r in (1:p)[-excluded]){
          determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = cbind(Z.r, Z2std.p),
                                             intercept = F,
                                             family = "gaussian",
                                             alpha = lasso.prop.i)
          penalty.max <- max(c(penalty.max, determine_lambdasequence$lambda))
          penalty.min <- min(c(penalty.min, determine_lambdasequence$lambda), na.rm = T)
        }

        penalty.seq <- exp(seq(log(penalty.max), log(penalty.min), length = n.penalty))
        cv <- rep(0, n.penalty)
        for(i.cv in 1:n.cv){
          for(i.r in (1:p)[-excluded]){
            determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = cbind(Z.r, Z2std.p),
                                          intercept = F,
                                          family = "gaussian", lambda = penalty.seq,
                                          alpha = lasso.prop.i)
            cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
          }
        }
        penalty.opt[i] <- penalty.seq[which.min(cv)]
        cv.opt[i] <- min(cv)
      }

      # Fit the final model
      lasso.prop.opt <- lasso.prop[which.min(cv.opt)]
      penalty.final <- penalty.opt[which.min(cv.opt)]

      for(i.r in (1:p)[-excluded]){
        LASSOfinal <- glmnet(y = Ystd[,i.r], x = cbind(Z.r, Zstd.p),
                             intercept = F,
                             lambda = penalty.seq, family = "gaussian",
                             alpha = lasso.prop.opt)
        ALPHA.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = penalty.final), nrow = 1)[-1]
      }
    } else {

      lasso.prop.opt <- penalty.final <- rep(NA, p)
      for(i.r in (1:p)[-excluded]){
        for(i in 1:length(lasso.prop)){
          lasso.prop.i <- lasso.prop[i]
          penalty.opt <- cv.opt <- rep(NA, length(lasso.prop.opt))

          determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = cbind(Z.r, Zstd.p),
                                             intercept = F,
                                             family = "gaussian",
                                             alpha = lasso.prop.i)
          penalty.max <- max(determine_lambdasequence$lambda)
          penalty.min <- min(determine_lambdasequence$lambda)

          penalty.seq <- exp(seq(log(penalty.max), log(penalty.min), length = n.penalty))

          cv <- rep(0, n.penalty)
          for(i.cv in 1:n.cv){
            determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = cbind(Z.r, Zstd.p),
                                          intercept = F,
                                          family = "gaussian", lambda = penalty.seq,
                                          alpha = lasso.prop.i)
            cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
          }
          penalty.opt[i] <- determine_lambda$lambda[which.min(cv)]
          cv.opt[i] <- min(cv)
        }
        lasso.prop.opt[i.r] <- lasso.prop[which.min(cv.opt, na.rm = T)]
        penalty.final[i.r] <- penalty.opt[which.min(cv.opt, na.rm = T)]

        # Fit the final model
        LASSOfinal <- glmnet(y = Ystd[,i.r], x = cbind(Z.r, Zstd.p), intercept = F,
                             family = "gaussian",
                             alpha = lasso.prop.opt[i.r])
        ALPHA.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = penalty.final[i.r]), nrow = 1)[-1]
      }
    }
  }

  # Calculation of the remaining parameters and adding NAs for excluded channels
  BETA.final <- matrix(0, nrow = p, ncol = r)
  BETA.final[-excluded,] <- BETA.full

  PI <- matrix(ALPHA.Sparse[,1:r], ncol = r) %*% t(BETA.final)
  if(!is.null(Z2)){
    GAMMAS <- ALPHA.Sparse[,(r+1):ncol(ALPHA.Sparse)]
  }

  MU <- meanY - PI %*% as.matrix(meanZ, ncol = 1) -
    GAMMAS %*% as.matrix(meanZ2, ncol = 1)

  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(MU) - Z %*% t(PI) - Z2 %*% t(GAMMAS)
  OMEGA <- (t(res) %*% res)/N

  res0  = Y - matrix(meanY, nrow = N, ncol = p, byrow = T)
  res0[,excluded] <- NA
  MU[excluded, ] <- NA
  GAMMAS[excluded, ] <- NA
  R2    = 1 - sum(res^2, na.rm = T)/sum(res0^2, na.rm = T)

  return(list(N = N, p = p, r = r,
              alpha = ALPHA.Sparse/dt, beta = BETA.final, Pi = PI/dt, Omega = OMEGA/dt,
              Psi = MU/dt, test = fit0$test, lambda = fit0$lambda, dt = dt, r2 = R2,
              penalty = penalty.final, lasso.prop = lasso.prop.opt))
}

# Wrapper function for fitting a VECM model with structural breaks
johansen.sb <- function(Y = NULL, Z = NULL, Z2 = NULL, H = NULL, h = NULL,
                        A.0 = NULL, B.0 = NULL, C.0 = NULL, Omega.0 = NULL,
                        break.pts, r = 1, n.iter = 20){
  N <- nrow(Y)
  p <- ncol(Y)
  p2 <- ncol(Z2)
  m <- length(break.pts)+1
  M <- ncol(B.0)
  p1 <- nrow(B.0)/M
  K <- ncol(H)/(p-1)

  out <- glsConstCpp(Y, Z, Z2, H, h, A.0, B.0, C.0, Omega.0, break.pts, r, n.iter)

  a.hat <- t(out[1:M,1:p])
  b.hat <- t(out[M+(1:M),1:p1])
  c.hat <- t(out[M+M+(1:p2),1:p])
  O.hat <- t(out[M+M+p2+(1:p),1:p])
  resid <- out[M+M+p2+p+(1:N),1:p]
  Lmax <- out[M+M+p2+p+N+(1:n.iter),1]

  # Put the parameters in the form that captures individual segments

  # alpha.hat <- vector("list", m)
  # beta.hat <- vector("list", m)
  # Omega.hat <- vector("list", m)
  # for(mm in 1:m){
  #   alpha.hat[[mm]] <- a.hat[,(mm-1)*r+(1:r)]
  #   beta.hat[[mm]] <- b.hat[,(mm-1)*r+(1:r)]
  #   Omega.hat[[mm]] <- O.hat[,(mm-1)*p+(1:p)]
  # }

  # Get the block diagonal form of beta
  # B <- bdiag(beta.hat)

  # return(out)

  return(list(A = a.hat, B = b.hat, C = c.hat, Omega = O.hat,
              resid = resid, Lmax = Lmax))
}
