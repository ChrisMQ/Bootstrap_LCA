##-------------------------------------------------------------------------------------#
## What Can Bootstrap Do in Latent Class Analysis? A Whole-Process Bootstrap Framework
## Author: Meng Qiu
## Last update: February 2024
##-------------------------------------------------------------------------------------#
## R functions for conducting the empirical example
##-------------------------------------------------------------------------------------#

#=======================================
# Functions for the Empirical Example
#=======================================
# Create a frequency table that allows zeros
table.pk <- function(vec, K) {
  v <- sum(vec==1)/length(vec)
  for(k in 2:K) {
    vv <- sum(vec==k)/length(vec)
    v <- c(v, vv)
  }
  return(v)
}

# Compute Bernulli density 
lik_j <- function(x, y) {
  z <- (x^y)*((1-x)^(1-y))
  return(z)
}

# Compute asymptotic CI
ml.ci <- function(theta.est, se.est, conf.level=.95) {
  alpha <- 1 - conf.level
  cv <- qnorm(p=alpha/2, mean=0, sd=1, lower.tail=F)
  ci <- c(theta.est - cv*se.est, theta.est + cv*se.est)
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  labs <- c(alpha/2, 1-alpha/2)*100
  labs <- paste0(labs, "%")
  names(ci) <- labs
  return(ci)
}

# Compute Student's t CI
st.ci <- function(theta.est, se.bt, N, conf.level=.95) {
  alpha <- 1 - conf.level
  tn <- qt(p=alpha/2, df=N-1, lower.tail=F)
  tp <- tn*se.bt
  low <- theta.est - tp
  high <- theta.est + tp
  ci <- c(low, high)
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  labs <- c(alpha/2, 1-alpha/2)*100
  labs <- paste0(labs, "%")
  names(ci) <- labs
  return(ci)
}

# Compute percentile CI
pt.ci <- function(theta.bt, conf.level=.95, bound=T) {
  alpha <- 1 - conf.level
  low <- alpha/2
  high <- 1 - low
  ci <- quantile(theta.bt, probs = c(low, high))
  if(bound == T) {
    if(ci[1] < 0) ci[1] <- 0
    if(ci[2] > 1) ci[2] <- 1
    return(ci)
  } else {
    return(ci)
  }
}

# Compute bias-corrected CI 
bc.ci <- function(theta.est, theta.bt, conf.level=.95){
  low <- (1 - conf.level)/2
  high <- 1 - low
  sims <- length(theta.bt)
  z.inv <- length(theta.bt[theta.bt < theta.est])/sims
  z <- qnorm(z.inv)
  lower.inv <-  pnorm(z + z + qnorm(low))
  lower <- quantile(theta.bt, lower.inv, names=FALSE)
  upper.inv <-  pnorm(z + z + qnorm(high))
  upper <- quantile(theta.bt, upper.inv, names=FALSE)
  ci <- c(lower, upper)
  labs <- c(low, high)*100
  labs <- paste0(labs, "%")
  names(ci) <- labs
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  return(ci)
}

# Compute three types of entropy
calc.entropy <- function(fit) {
  N <- fit$Nobs
  S <- length(fit$P)
  dat <- fit$y
  names(dat) <- paste0("Y", 1:ncol(dat))
  post.p <- fit$posterior
  if(S > 1) {
    ei <- round(1-apply(post.p,1,max),digits=3) # individual uncertainty: closer to 0 is better
    indi.ei <- cbind(dat, ei)
    pat.ei <- (indi.ei) %>% group_by_all %>% count
    pat.ei <- as.data.frame(pat.ei)
    func <- function(x)-x*log(x)
    entropy <- sum(apply(post.p, 2, func)) # entropy: lower the better
    rel.entropy <- 1 - entropy/(N*log(S)) # relative entropy: closer to 1 is better
  } else {
    pat.ei <- entropy <- rel.entropy <- NA
    message("Warning: Entropy is undefined for K=1.")
  }
  rtn <- list(pat.ei=pat.ei, entropy=entropy, rel.entropy=rel.entropy)
  return(rtn)
}

# Part 1: model selection and classification uncertainty
LCA.bt.NB <- function(data, nrep, nrep.bt, maxiter, K, B, conf.level=.95, verbose=T) {
  N <- nrow(data)
  J <- ncol(data)
  colnames(data) <- paste0("Y", 1:J)
  
  fit.bt.list <- list()
  count <- matrix(0, nrow=B, ncol=K)
  aics <- aiccs <- caics <- bics <- ssbics <- dbics <- c()
  aics.bt <- aiccs.bt <- caics.bt <- bics.bt <- ssbics.bt <- dbics.bt <- matrix(0, nrow=B, ncol=K)
  aic.pk <- aicc.pk <- caic.pk <- bic.pk <- ssbic.pk <- dbic.pk <- c()
  select.r <- indx.m <- indx.sd <- indx.min <- indx.max <- matrix(0, nrow=6, ncol=K,
                                                                  dimnames=list(c("aic","aicc","caic","bic","ssbic","dbic"),paste0("K=",1:K)))
  f <- as.formula( paste( "cbind(", paste( paste0("Y", 1:J), collapse="," ), ")", "~1" ) )
  fit.org <- list()
  for(k in 1:K) {
    fit <- poLCA(f, data, nclass=k, maxiter=maxiter, nrep=nrep, verbose=F)
    npar <- fit$npar
    aics[k] <- fit$aic
    aiccs[k] <- fit$aic + (2*npar)*(npar+1)/(N-npar-1)
    caics[k] <- -2*(fit$llik) + npar*(log(N)+1)
    bics[k] <- fit$bic
    ssbics[k] <- -2*(fit$llik) + npar*log((N+2)/24)
    dbics[k] <- -2*(fit$llik) + npar*(log(N)-log(2*pi))
    fit.org[[k]] <- fit
  }
  
  # container for fit information across bootstrap
  fit.bt.save <- list()
  fit.bt.names <- paste0("K",1:K,".fits")
  for(k in 1:K) {
    fit.bt.save[[k]] <- assign(fit.bt.names[k], list())
  }
  names(fit.bt.save) <- fit.bt.names
  
  dat.bt.list <- list()
  if(verbose)
    pb <- txtProgressBar(min=0, max=B, style=3, width=40, char="=")
  for(b in 1:B) {
    idx <- sample(1:N, size=N, replace=T)
    dat.bt <- data[idx, ]
    dat.bt.list[[b]] <- dat.bt
    for(k in 1:K) {
      post0 <- fit.org[[k]]$posterior
      inits <- post0[idx, ]
      fit.bt <- poLCA.useprior(f, dat.bt, nclass=k, nrep=nrep.bt, maxiter=maxiter,
                               pp.start=inits, verbose=F, graphs=F)
      fit.bt.save[[k]][[b]] <- fit.bt
      if(fit.bt$numiter <= fit.bt$maxiter) {
        count[b, k] <- 0
        npar <- fit.bt$npar
        aics.bt[b,k] <- fit.bt$aic
        aiccs.bt[b,k] <- fit.bt$aic + (2*npar)*(npar+1)/(N-npar-1)
        caics.bt[b,k] <- -2*(fit.bt$llik) + npar*(log(N)+1)
        bics.bt[b,k] <- fit.bt$bic
        ssbics.bt[b,k] <- -2*(fit.bt$llik) + npar*log((N+2)/24)
        dbics.bt[b,k] <- -2*(fit.bt$llik) + npar*(log(N)-log(2*pi))
      } else {
        count[b, k] <- 1
      }
    }
    aic.pk[b] <- which.min(aics.bt[b,])
    aicc.pk[b] <- which.min(aiccs.bt[b,])
    caic.pk[b] <- which.min(caics.bt[b,])
    bic.pk[b] <- which.min(bics.bt[b,])
    ssbic.pk[b] <- which.min(ssbics.bt[b,])
    dbic.pk[b] <- which.min(dbics.bt[b,])
    if(verbose) setTxtProgressBar(pb, b)
  }
  if(verbose) close(pb)
  
  # compute EIC's (i.e., AIC variants)
  eic1 <- eic2 <- eic3 <- eic4 <- eic5 <- c()
  for(k in 1:K) {
    fit.use <- fit.bt.save[[k]]
    llk <- fit.org[[k]]$llik
    llk1 <- llk2 <- llk3 <- c()
    for(b in 1:B) {
      # compute two types of log-likelihood
      dat.org <- data - 1
      dat.bt <- dat.bt.list[[b]] - 1
      P.org <- fit.org[[k]]$P
      P.bt <- fit.use[[b]]$P
      KK <- length(P.bt)
      JJ <- ncol(data)
      NN <- nrow(data)
      ip.org <- fit.org[[k]]$probs
      ip.bt <- fit.use[[b]]$probs
      pi.org <- sapply(ip.org, FUN=function(S){S[,2]}) # ans=1
      qi.org <- sapply(ip.org, FUN=function(S){S[,1]}) # ans=0
      pi.bt <- sapply(ip.bt, FUN=function(S){S[,2]}) # ans=1
      qi.bt <- sapply(ip.bt, FUN=function(S){S[,1]}) # ans=0
      if(length(P.org)==1) {
        f_i_1 <- f_i_2 <- f_i_3 <- c()
        f_ij_1 <- f_ij_2 <- f_ij_3 <- matrix(0, nrow=NN, ncol=JJ)
        for (i in 1:NN){
          for (j in 1:JJ){
            f_ij_1[i,j] <- (pi.bt[j]^dat.org[i,j])*(qi.bt[j]^(1-dat.org[i,j]))
            f_ij_2[i,j] <- (pi.bt[j]^dat.bt[i,j])*(qi.bt[j]^(1-dat.bt[i,j]))
            f_ij_3[i,j] <- (pi.org[j]^dat.bt[i,j])*(qi.org[j]^(1-dat.bt[i,j]))
          }
        }
        f_i_1 <- apply(f_ij_1, 1, prod)
        f_i_2 <- apply(f_ij_2, 1, prod)
        f_i_3 <- apply(f_ij_3, 1, prod)
        llk1[b] <- sum(log(f_i_1))
        llk2[b] <- sum(log(f_i_2))
        llk3[b] <- sum(log(f_i_3))
      } else {
        f_ik_1 <- f_ik_2 <- f_ik_3 <- matrix(0, nrow=NN, ncol=KK)
        f_ijk_1 <- f_ijk_2 <- f_ijk_3 <- array(0, dim=c(NN,JJ,KK))
        for (i in 1:NN){
          for (c in 1:KK){
            for (j in 1:JJ){
              f_ijk_1[i,j,c] <- (pi.bt[c,j]^dat.org[i,j])*(qi.bt[c,j]^(1-dat.org[i,j]))
              f_ijk_2[i,j,c] <- (pi.bt[c,j]^dat.bt[i,j])*(qi.bt[c,j]^(1-dat.bt[i,j]))
              f_ijk_3[i,j,c] <- (pi.org[c,j]^dat.bt[i,j])*(qi.org[c,j]^(1-dat.bt[i,j]))
            }
          }
        }
        for (c in 1:KK) {
          f_ik_1[,c] <- P.bt[c]*apply(f_ijk_1[,,c], 1, prod)
          f_ik_2[,c] <- P.bt[c]*apply(f_ijk_2[,,c], 1, prod)
          f_ik_3[,c] <- P.org[c]*apply(f_ijk_3[,,c], 1, prod)
        }
        llk1[b] <- sum(log(rowSums(f_ik_1)))
        llk2[b] <- sum(log(rowSums(f_ik_2)))
        llk3[b] <- sum(log(rowSums(f_ik_3)))
      }
    }
    llk.mat <- cbind(llk1, llk2, llk3)
    mm <- llk.mat[is.finite(rowSums(llk.mat)), ]
    n.inf <- B - nrow(mm) # number of cases with infinity llk
    llk1 <- mm[,1]; llk2 <- mm[,2]; llk3 <- mm[,3]
    eic1[k] <- -2*llk + 2*mean(llk2 - llk1)
    eic2[k] <- -2*llk + 4*mean(llk - llk1)
    eic3[k] <- -2*llk + 4*mean(llk2 - llk3)
    eic4[k] <- -2*llk + 4*mean(llk3 - llk1)
    eic5[k] <- -2*llk + 4*mean(llk2 - llk)
  }
  indx.list <- list(aics, aiccs, caics, bics, ssbics, dbics,
                    eic1, eic2, eic3, eic4, eic5)
  indx.mat <- rbind(aics, aiccs, caics, bics, ssbics, dbics,
                    eic1, eic2, eic3, eic4, eic5)
  indx.s <- unlist(lapply(indx.list, FUN=which.min))
  dif.best2 <- unlist(lapply(indx.list, function(x) sort(x)[2]-sort(x)[1]))
  mod.best2 <- do.call(rbind, lapply(indx.list, function(x) c(order(x)[1], order(x)[2])))
  indx.names <- c("aic","bic","aicc","caic","ssbic","dbic",
                  "eic1","eic2","eic3","eic4","eic5")
  names(indx.s) <- names(dif.best2) <- indx.names
  rownames(indx.mat) <- indx.names
  colnames(indx.mat) <- paste0("K=", 1:K)
  rownames(mod.best2) <- indx.names
  colnames(mod.best2) <- c("1st.best","2nd.best")
  
  indx.bt.list <- list(aics.bt, aiccs.bt, caics.bt, bics.bt, ssbics.bt, dbics.bt)
  for(d in 1:6) {
    indx.m[d,] <- colMeans(indx.bt.list[[d]])
    indx.sd[d,] <- apply(indx.bt.list[[d]], 2, sd)
    indx.min[d,] <- apply(indx.bt.list[[d]], 2, min)
    indx.max[d,] <- apply(indx.bt.list[[d]], 2, max)
  }
  select.r[1,] <- table.pk(aic.pk, K)
  select.r[2,] <- table.pk(aicc.pk, K)
  select.r[3,] <- table.pk(caic.pk, K)
  select.r[4,] <- table.pk(bic.pk, K)
  select.r[5,] <- table.pk(ssbic.pk, K)
  select.r[6,] <- table.pk(dbic.pk, K)
  select.r <- round(select.r, digits=4)
  
  non.conv <- colSums(count)/nrow(count)
  names(non.conv) <- paste0("K=", 1:K)
  res <- list(select.single=indx.s, dif.best2=dif.best2, mod.best2=mod.best2,
              indx.value=indx.mat, select.rate=select.r,
              indx.mean=indx.m, indx.sd=indx.sd,
              indx.min=indx.min, indx.max=indx.max,
              fit.save=fit.org, fit.bt.save=fit.bt.save,
              indx.bt.list=indx.bt.list, non.conv=non.conv)
  return(res)
}

# Part 2: parameter estimation & classification uncertainty
LCA.bt.EST <- function(MS.out, nclass, conf.level) {
  
  if (nclass == 1) stop("Error: The number of classes must be greater than 1.")
  
  # subset fitted models for further usage
  fit.org.use <- MS.out$fit.save[[nclass]]
  fit.bt.use <- MS.out$fit.bt.save[[nclass]]
  
  # organize bootstrap parameter estimates and standard errors
  K <- nclass
  N <- fit.org.use$Nobs
  J <- length(fit.org.use$probs)
  B <- length(fit.bt.use)
  bt.paras <- matrix(0, nrow=B, ncol=(K+J*K)) # a B by (K+K*J) container for all bt parameters
  P.names <- paste0("P",1:K)
  K.ind <- rep(1:K, each=J)
  J.ind <- rep(1:J, K)
  ip.names <- c()
  for(i in 1:(J*K)) {
    tt <- paste0("ip.", J.ind[i], K.ind[i])
    ip.names <- c(ip.names, tt)
  }
  para.names <- c(P.names, ip.names)
  colnames(bt.paras) <- para.names 
  for(b in 1:B) {
    P <- fit.bt.use[[b]]$P
    ip <- sapply(fit.bt.use[[b]]$probs, FUN=function(S){S[,2]})
    ip <- c(t(ip))
    bt.paras[b, ] <- c(P, ip)
  }
  P.est.ML <- fit.org.use$P
  ip.est.ML <- sapply(fit.org.use$probs, FUN=function(S){S[,2]})
  P.se.ML <- fit.org.use$P.se
  ip.se.ML <- sapply(fit.org.use$probs.se, FUN=function(S){S[,2]})
  ML.est <- c(P.est.ML, c(t(ip.est.ML)))
  ML.se <- c(P.se.ML, c(t(ip.se.ML)))
  names(ML.est) <- names(ML.se) <- para.names
  bt.est <- colMeans(bt.paras)
  bt.se <- apply(bt.paras, 2, sd)
  bc.bt.paras <- sweep(bt.paras, 2, bt.est) # bias-correct all bootstrap est.
  bc.se <- apply(bc.bt.paras, 2, sd) # standard error for bias-corrected bootstrap est.
  bc.est <- 2*ML.est - bt.est # bias-corrected bootstrap estimates
  
  # compute all kinds of CI's
  ci.ML <- t(mapply(ml.ci, ML.est, ML.se, conf.level))
  ci.st <- t(mapply(st.ci, bt.est, bt.se, N, conf.level))
  ci.pt <- t(apply(bt.paras, 2, pt.ci, conf.level))
  ci.bc <- matrix(0, nrow=(K+J*K), ncol=2)
  rownames(ci.bc) <- para.names
  alpha <- 1-conf.level
  labs <- c(alpha/2, 1-alpha/2)*100
  labs <- paste0(labs, "%")
  colnames(ci.bc) <- labs
  for(d in 1:(K+J*K)) {
    ci.bc[d, ] <- bc.ci(ML.est[d], bt.paras[,d], conf.level=conf.level)
  }
  
  # organize parameter estimation information
  est.MLnBT <- rbind(ML.est, bt.est, bc.est) # put ML est. and bias-corrected bootstrap est. together
  est.MLnBT[est.MLnBT > 1] <- 1
  est.MLnBT[est.MLnBT < 0] <- 0
  se.MLnBT <- rbind(ML.se, bt.se) # put ML se. and bootstrap se. together
  Estimates <- list(est.MLnBT=est.MLnBT, se.MLnBT=se.MLnBT)
  Estimates <- lapply(Estimates, round, 4)
  CI <- list(ML=ci.ML, Student_t=ci.st, Percentile=ci.pt, Bias_corrected=ci.bc) # put all CI's together
  CI <- lapply(CI, round, 4)
  
  # compute entropy information
  # 1. ML entropy estimates
  en.ML <- calc.entropy(fit.org.use)
  n.pat <- nrow(en.ML$pat.ei)        # number of response patterns
  # 2. CI for entropy estimates
  ei <- list()
  en <- rel.en <- c()
  for(b in 1:B) {
    temp <- calc.entropy(fit.bt.use[[b]])
    ei[[b]] <- temp$pat.ei$ei
    en[b] <- temp$entropy
    rel.en[b] <- temp$rel.entropy
  }
  indx.use <- which(is.nan(en) == FALSE & unlist(lapply(ei, function(x) length(x) == n.pat)) == TRUE)
  ei <- ei[indx.use]
  ei.mat <- matrix(unlist(ei), nrow=n.pat)
  en <- en[indx.use]
  rel.en <- rel.en[indx.use]
  # BT mean
  ei.m.bt <- apply(ei.mat, 1, mean)
  en.m.bt <- mean(en)
  rel.en.m.bt <- mean(rel.en)
  en.Mean <- list(ei.m=ei.m.bt, en.m=en.m.bt, rel.en.m=rel.en.m.bt)
  # BT se
  ei.se.bt <- apply(ei.mat, 1, sd)
  en.se.bt <- sd(en)
  rel.en.se.bt <- sd(rel.en)
  en.SE <- list(ei.se=ei.se.bt, en.se=en.se.bt, rel.en.se=rel.en.se.bt)
  # 95% percentile CI
  ei.ci.pt <- apply(ei.mat, 1, pt.ci)
  en.ci.pt <- pt.ci(en, bound=F)
  rel.en.ci.pt <- pt.ci(rel.en)
  en.PT <- list(ei.ci=ei.ci.pt, en.ci=en.ci.pt, rel.en.ci=rel.en.ci.pt)
  en.PT <- lapply(en.PT, round, 4)
  
  Entropy <- list(en.ML=en.ML, en.Mean=en.Mean, en.SE=en.SE, en.PT=en.PT)
  res <- list(Estimates=Estimates, CI=CI, Entropy=Entropy)
  return(res)
}

#======================================================
# Functions are for obtaining p-value for the BVR
#======================================================
# Creating unique patterns for values in two vectors
unique.pat <- function(v1, v2) {
  mat <- expand.grid(1:v1, 1:v2)
  dup <- apply(mat, 1, function(row) length(unique(row)) == 1)
  df <- mat[-which(dup==TRUE), ]
  pat <- as.data.frame((df <- t(apply(df, 1, sort)))[!duplicated(df), ])
  return(pat)
}
# Computing BVR statistic for a pair of items
compute.bvr <- function(dat, c1, c2) {
  N <- nrow(dat)
  ctab <- crosstab(dat[, c1], dat[, c2], plot=F)
  of <- ctab$tab
  ef <- chisq.test(ctab$tab)$expected
  out <- sum(mapply(function(v1,v2) (v1-v2)^2/v2, c(of), c(ef)))
  return(out)
}
# Performing parametric bootstrap to obtain p-values for each pair of items
bvr.pVal <- function(data, nclass, B=100, nrep=20, alpha=0.05, seed=123) {
  N <- nrow(data)
  nc <- ncol(data)
  names <- paste0("V",1:nc)
  colnames(data) <- names
  pair <- unique.pat(nc, nc)
  npair <- nrow(pair)
  f <- as.formula(paste("cbind(",paste(names,collapse=","),")","~1"))
  # compute observed BVR values
  set.seed(seed)
  Mod <- poLCA(f, data, nclass, nrep, verbose = F) 
  prop <- Mod$P
  ip <- Mod$probs
  bvr.obs <- c()
  for(n in 1:npair) {
    bvr.obs[n] <- compute.bvr(data, pair[n,1], pair[n,2])
  }
  # compute p-values for the BVR statistic 
  bvr.pb <- matrix(0, nrow = B, ncol = npair)
  set.seed(seed)
  for(b in 1:B) {
    dat.sim <- poLCA.simdata(N = N, probs = ip, P = prop)$dat
    for(n in 1:npair) {
      bvr.pb[b,n] <- compute.bvr(dat.sim, pair[n,1], pair[n,2])
    }
  }
  p.val <- c()
  for(i in 1:npair) {
    p.val[i] <- sum(bvr.pb[,i] >= bvr.obs[i])/B
  }
  select.pair <- which(p.val < alpha)
  names.pair <- apply(pair, 1, function(x) paste("(",x[1],",",x[2],")", sep=""))
  if(length(select.pair) == 0 ) {
    identified.pair <- "NA"
  } else {
    identified.pair <- data.frame(names.pair[select.pair], bvr.obs[select.pair], p.val[select.pair])
    colnames(identified.pair) <- c("Pair", "Observed BVR", "p-value")
  }
  all.pair <- data.frame(names.pair, bvr.obs, p.val)
  colnames(all.pair) <- c("Pair", "Observed BVR", "p-value")
  out <- list("Identified Pairs"=identified.pair, "All Pairs"=all.pair)
  return(out)
}

#===========================================================
# The following function is adopted from the poLCA package, 
# with the modification of incorporating a user-defined
# initialization matrix. 
#===========================================================
## poLCA using user-defined initial classification probabilities
poLCA.useprior <- 
  function(formula,data,nclass=2,maxiter=1000,graphs=FALSE,tol=1e-10,
           na.rm=TRUE,probs.start=NULL,nrep=1,verbose=TRUE,calc.se=TRUE,pp.start) {
    starttime <- Sys.time()
    mframe <- model.frame(formula,data,na.action=NULL)
    mf <- model.response(mframe)
    if (any(mf<1,na.rm=TRUE) | any(round(mf) != mf,na.rm=TRUE)) {
      cat("\n ALERT: some manifest variables contain values that are not
    positive integers. For poLCA to run, please recode categorical
    outcome variables to increment from 1 to the maximum number of
    outcome categories for each variable. \n\n")
      ret <- NULL
    } else {
      data <- data[rowSums(is.na(model.matrix(formula,mframe)))==0,]
      if (na.rm) {
        mframe <- model.frame(formula,data)
        y <- model.response(mframe)
      } else {
        mframe <- model.frame(formula,data,na.action=NULL)
        y <- model.response(mframe)
        y[is.na(y)] <- 0
      }
      if (any(sapply(lapply(as.data.frame(y),table),length)==1)) {
        y <- y[,!(sapply(apply(y,2,table),length)==1)]
        cat("\n ALERT: at least one manifest variable contained only one
    outcome category, and has been removed from the analysis. \n\n")
      }
      x <- model.matrix(formula,mframe)
      N <- nrow(y)
      J <- ncol(y)
      K.j <- t(matrix(apply(y,2,max)))
      R <- nclass
      S <- ncol(x)
      if (S>1) { calc.se <- TRUE }
      eflag <- FALSE
      probs.start.ok <- TRUE
      ret <- list()
      if (R==1) {
        ret$probs <- list()
        for (j in 1:J) {
          ret$probs[[j]] <- matrix(NA,nrow=1,ncol=K.j[j])
          for (k in 1:K.j[j]) { ret$probs[[j]][k] <- sum(y[,j]==k)/sum(y[,j]>0) }
        }
        ret$probs.start <- ret$probs
        ret$P <- 1
        ret$posterior <- ret$predclass <- prior <- matrix(1,nrow=N,ncol=1)
        ret$llik <- sum(log(poLCA.ylik.C(poLCA.vectorize(ret$probs),y)) - log(.Machine$double.xmax))
        if (calc.se) {
          se <- poLCA.se(y,x,ret$probs,prior,ret$posterior)
          ret$probs.se <- se$probs           # standard errors of class-conditional response probabilities
          ret$P.se <- se$P                   # standard errors of class population shares
        } else {
          ret$probs.se <- NA
          ret$P.se <- NA
        }
        ret$numiter <- 1
        ret$probs.start.ok <- TRUE
        ret$coeff <- NA
        ret$coeff.se <- NA
        ret$coeff.V <- NA
        ret$eflag <- FALSE
        if (S>1) {
          cat("\n ALERT: covariates not allowed when nclass=1; will be ignored. \n \n")
          S <- 1
        }
      } else {
        if (!is.null(probs.start)) { # error checking on user-inputted probs.start
          if ((length(probs.start) != J) | (!is.list(probs.start))) {
            probs.start.ok <- FALSE
          } else {
            if (sum(sapply(probs.start,dim)[1,]==R) != J) probs.start.ok <- FALSE
            if (sum(sapply(probs.start,dim)[2,]==K.j) != J) probs.start.ok <- FALSE
            if (sum(round(sapply(probs.start,rowSums),4)==1) != (R*J)) probs.start.ok <- FALSE
          }
        }
        ret$llik <- -Inf
        ret$attempts <- NULL
        for (repl in 1:nrep) { # automatically reestimate the model multiple times to locate the global max llik
          error <- TRUE; firstrun <- TRUE
          probs <- probs.init <- probs.start
          while (error) { # error trap
            error <- FALSE
            b <- rep(0,S*(R-1))
            # prior <- poLCA.updatePrior(b,x,R)
            prior <- pp.start
            if ((!probs.start.ok) | (is.null(probs.start)) | (!firstrun) | (repl>1)) { # only use the specified probs.start in the first nrep
              probs <- list()
              for (j in 1:J) { 
                probs[[j]] <- matrix(runif(R*K.j[j]),nrow=R,ncol=K.j[j])
                probs[[j]] <- probs[[j]]/rowSums(probs[[j]]) 
              }
              probs.init <- probs
            }
            vp <- poLCA.vectorize(probs)
            iter <- 1
            llik <- matrix(NA,nrow=maxiter,ncol=1)
            llik[iter] <- -Inf
            dll <- Inf
            while ((iter <= maxiter) & (dll > tol) & (!error)) {
              iter <- iter+1
              rgivy <- poLCA.postClass.C(prior,vp,y)      # calculate posterior
              vp$vecprobs <- poLCA.probHat.C(rgivy,y,vp)  # update probs
              if (S>1) {
                dd <- poLCA.dLL2dBeta.C(rgivy,prior,x)
                b <- b + ginv(-dd$hess) %*% dd$grad     # update betas
                prior <- poLCA.updatePrior(b,x,R)       # update prior    
              } else {
                prior <- matrix(colMeans(rgivy),nrow=N,ncol=R,byrow=TRUE)
              }
              llik[iter] <- sum(log(rowSums(prior*poLCA.ylik.C(vp,y))) - log(.Machine$double.xmax))
              dll <- llik[iter]-llik[iter-1]
              if (is.na(dll)) {
                error <- TRUE
              } else if ((S>1) & (dll < -1e-7)) {
                error <- TRUE
              }
            }
            if (!error) { 
              if (calc.se) {
                se <- poLCA.se(y,x,poLCA.unvectorize(vp),prior,rgivy)
              } else {
                se <- list(probs=NA,P=NA,b=NA,var.b=NA)
              }
            } else {
              eflag <- TRUE
            }
            firstrun <- FALSE
          } # finish estimating model without triggering error
          ret$attempts <- c(ret$attempts,llik[iter])
          if (llik[iter] > ret$llik) {
            ret$llik <- llik[iter]             # maximum value of the log-likelihood
            ret$probs.start <- probs.init      # starting values of class-conditional response probabilities
            ret$probs <- poLCA.unvectorize(vp) # estimated class-conditional response probabilities
            ret$probs.se <- se$probs           # standard errors of class-conditional response probabilities
            ret$P.se <- se$P                   # standard errors of class population shares
            ret$posterior <- rgivy             # NxR matrix of posterior class membership probabilities
            ret$predclass <- apply(ret$posterior,1,which.max)   # Nx1 vector of predicted class memberships, by modal assignment
            ret$P <- colMeans(ret$posterior)   # estimated class population shares
            ret$numiter <- iter-1              # number of iterations until reaching convergence
            ret$probs.start.ok <- probs.start.ok # if starting probs specified, logical indicating proper entry format
            if (S>1) {
              b <- matrix(b,nrow=S)
              rownames(b) <- colnames(x)
              rownames(se$b) <- colnames(x)
              ret$coeff <- b                 # coefficient estimates (when estimated)
              ret$coeff.se <- se$b           # standard errors of coefficient estimates (when estimated)
              ret$coeff.V <- se$var.b        # covariance matrix of coefficient estimates (when estimated)
            } else {
              ret$coeff <- NA
              ret$coeff.se <- NA
              ret$coeff.V <- NA
            }
            ret$eflag <- eflag                 # error flag, true if estimation algorithm ever needed to restart with new initial values
          }
          if (nrep>1 & verbose) { cat("Model ",repl,": llik = ",llik[iter]," ... best llik = ",ret$llik,"\n",sep=""); flush.console() }
        } # end replication loop
      }
      names(ret$probs) <- colnames(y)
      if (calc.se) { names(ret$probs.se) <- colnames(y) }
      ret$npar <- (R*sum(K.j-1)) + (R-1)                  # number of degrees of freedom used by the model (number of estimated parameters)
      if (S>1) { ret$npar <- ret$npar + (S*(R-1)) - (R-1) }
      ret$aic <- (-2 * ret$llik) + (2 * ret$npar)         # Akaike Information Criterion
      ret$bic <- (-2 * ret$llik) + (log(N) * ret$npar)    # Schwarz-Bayesian Information Criterion
      ret$Nobs <- sum(rowSums(y==0)==0)                   # number of fully observed cases (if na.rm=F)
      if (all(rowSums(y==0)>0)) { # if no rows are fully observed
        ret$Chisq <- NA
        ret$Gsq <- NA
        ret$predcell <- NA
      } else {
        compy <- poLCA.compress(y[(rowSums(y==0)==0),])
        datacell <- compy$datamat
        rownames(datacell) <- NULL
        freq <- compy$freq
        ylik <- poLCA.ylik.C(poLCA.vectorize(ret$probs),datacell) 
        if (!na.rm) {
          fit <- matrix(ret$Nobs/.Machine$double.xmax * (ylik  %*% ret$P))
          ret$Chisq <- sum((freq-fit)^2/fit) + (ret$Nobs-sum(fit)) # Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables
        } else {
          fit <- matrix(N/.Machine$double.xmax * (ylik %*% ret$P))
          ret$Chisq <- sum((freq-fit)^2/fit) + (N-sum(fit))
        }
        ret$predcell <- data.frame(datacell,observed=freq,expected=round(fit,3)) # Table that gives observed vs. predicted cell counts
        ret$Gsq <- 2 * sum(freq*log(freq/fit))  # Likelihood ratio/deviance statistic
      }
      y[y==0] <- NA
      ret$y <- data.frame(y)             # outcome variables
      ret$x <- data.frame(x)             # covariates, if specified
      for (j in 1:J) {
        rownames(ret$probs[[j]]) <- paste("class ",1:R,": ",sep="")
        if (is.factor(data[,match(colnames(y),colnames(data))[j]])) {
          lev <- levels(data[,match(colnames(y),colnames(data))[j]])
          colnames(ret$probs[[j]]) <- lev
          ret$y[,j] <- factor(ret$y[,j],labels=lev)
        } else {
          colnames(ret$probs[[j]]) <- paste("Pr(",1:ncol(ret$probs[[j]]),")",sep="")
        }
      }
      ret$N <- N                         # number of observations
      ret$maxiter <- maxiter             # maximum number of iterations specified by user
      ret$resid.df <- min(ret$N,(prod(K.j)-1))-ret$npar # number of residual degrees of freedom
      class(ret) <- "poLCA"
      if (graphs) plot.poLCA(ret)
      if (verbose) print.poLCA(ret)
      ret$time <- Sys.time()-starttime   # how long it took to run the model
    }
    ret$call <- match.call()
    return(ret)
  }

#=============================================================
# The following functions are adopted from the poLCA package
#=============================================================

poLCA.dLL2dBeta.C <-
  function(rgivy,prior,x) {
    classes <- dim(prior)[2]
    numx <- dim(x)[2]
    ret <-  .C("d2lldbeta2",
               as.double(t(rgivy)),
               as.double(t(prior)),
               as.double(t(x)),
               as.integer(dim(x)[1]),
               as.integer(classes),
               as.integer(numx),
               grad = double((classes-1)*numx),
               hess = double(((classes-1)*numx)^2)                
    )
    return(list(grad=ret$grad,hess=-matrix(ret$hess,ncol=((classes-1)*numx),byrow=TRUE)))
  }

print.poLCA <-
  function(x, ...) {
    R <- length(x$P)
    S <- ifelse(is.na(x$coeff[1]),1,nrow(x$coeff))
    cat("Conditional item response (column) probabilities,\n by outcome variable, for each class (row) \n \n")
    print(lapply(x$probs,round,4))
    cat("Estimated class population shares \n", round(x$P,4), "\n \n")
    cat("Predicted class memberships (by modal posterior prob.) \n",round(table(x$predclass)/x$N,4), "\n \n")
    cat("========================================================= \n")
    cat("Fit for", R, "latent classes: \n")
    cat("========================================================= \n")
    if (S>1) {
      for (r in 2:R) {
        cat(r,"/ 1 \n")
        disp <- data.frame(coeff=round(x$coeff[,(r-1)],5),
                           se=round(x$coeff.se[,(r-1)],5),
                           tval=round(x$coeff[,(r-1)]/x$coeff.se[,(r-1)],3),
                           pr=round(1-(2*abs(pt(x$coeff[,(r-1)]/x$coeff.se[,(r-1)],x$resid.df)-0.5)),3))
        colnames(disp) <- c("Coefficient"," Std. error"," t value"," Pr(>|t|)")
        print(disp)
        cat("========================================================= \n")
      }
    }
    cat("number of observations:", x$N, "\n")
    if(x$N != x$Nobs) cat("number of fully observed cases:", x$Nobs, "\n")
    cat("number of estimated parameters:", x$npar, "\n")
    cat("residual degrees of freedom:", x$resid.df, "\n")
    cat("maximum log-likelihood:", x$llik, "\n \n")
    cat("AIC(",R,"): ",x$aic,"\n",sep="")
    cat("BIC(",R,"): ",x$bic,"\n",sep="")
    if (S==1) cat("G^2(",R,"): ",x$Gsq," (Likelihood ratio/deviance statistic) \n",sep="")
    cat("X^2(",R,"): ",x$Chisq," (Chi-square goodness of fit) \n \n",sep="")
    if (x$numiter==x$maxiter) cat("ALERT: iterations finished, MAXIMUM LIKELIHOOD NOT FOUND \n \n")
    if (!x$probs.start.ok) cat("ALERT: error in user-specified starting values; new start values generated \n \n")
    if (x$npar>x$N) cat("ALERT: number of parameters estimated (",x$npar,") exceeds number of observations (",x$N,") \n \n")
    if (x$resid.df<0) cat("ALERT: negative degrees of freedom; respecify model \n \n")
    if (x$eflag) cat("ALERT: estimation algorithm automatically restarted with new initial values \n \n")
    flush.console()
    invisible(x)
  }

poLCA.updatePrior <-
  function(b,x,R) {
    b <- matrix(b,ncol=(R-1))
    exb <- exp(x %*% b)
    p <- cbind(1,exb)/(rowSums(exb)+1)
    return(p)
  }

poLCA.vectorize <-
  function(probs) {
    classes <- nrow(probs[[1]])
    vecprobs <- unlist(lapply(probs,t))
    numChoices <- sapply(probs,ncol)
    return(list(vecprobs=vecprobs,numChoices=numChoices,classes=classes))
  }

poLCA.postClass.C <-
  function(prior,vp,y) {
    ret <-  .C("postclass",
               as.double(t(prior)),
               as.double(vp$vecprobs),
               as.integer(t(y)),
               as.integer(length(vp$numChoices)),
               as.integer(dim(y)[1]),
               as.integer(vp$numChoices),
               as.integer(vp$classes),
               posterior = double(dim(y)[1]*vp$classes)
    )
    ret$posterior <- matrix(ret$posterior,ncol=vp$classes,byrow=TRUE)
    return(ret$posterior)
  }

poLCA.ylik.C <-
  function(vp,y) {
    ret <-  .C("ylik",
               as.double(vp$vecprobs),
               as.integer(t(y)),
               as.integer(dim(y)[1]),
               as.integer(length(vp$numChoices)),
               as.integer(vp$numChoices),
               as.integer(vp$classes),
               lik = double(dim(y)[1]*vp$classes)
    )
    ret$lik <- matrix(ret$lik,ncol=vp$classes,byrow=TRUE)
    return(ret$lik)
  }

poLCA.se <-
  function(y,x,probs,prior,rgivy) {
    J <- ncol(y)
    R <- ncol(prior)
    K.j <- sapply(probs,ncol)
    N <- nrow(y)
    ymat <- y
    y <- list()
    for (j in 1:J) {
      y[[j]] <- matrix(0,nrow=N,ncol=K.j[j])
      y[[j]][cbind(c(1:N),ymat[,j])] <- 1
      y[[j]][ymat[,j]==0,] <- NA      # To handle missing values
    }
    s <- NULL
    # score matrix contains \sum_j [R(K.j-1)] columns correpsonding to log-odds response probs...
    for (r in 1:R) {
      for (j in 1:J) {
        s <- cbind(s,rgivy[,r] * t(t(y[[j]][,2:K.j[j]]) - probs[[j]][r,2:K.j[j]]))
      }
    }
    # ...and (R-1)*ncol(x) columns corresponding to coefficients (betas)
    ppdiff <- rgivy-prior
    if (R>1) for (r in 2:R) { s <- cbind(s,x*ppdiff[,r]) }
    
    s[is.na(s)] <- 0      # To handle missing values
    info <- t(s) %*% s    # Information matrix
    VCE <- ginv(info)     # VCE matrix of log-odds response probs and covariate coefficients
    
    # Variance of class conditional response probs using delta fn. transformation with 
    # Jacobian a block diagonal matrix (across r) of block diagonal matrices (across j)
    VCE.lo <- VCE[1:sum(R*(K.j-1)),1:sum(R*(K.j-1))]
    Jac <- matrix(0,nrow=nrow(VCE.lo)+(J*R),ncol=ncol(VCE.lo))
    rpos <- cpos <- 1
    for (r in 1:R) {
      for (j in 1:J) {
        Jsub <- -(probs[[j]][r,] %*% t(probs[[j]][r,]))
        diag(Jsub) <- probs[[j]][r,]*(1-probs[[j]][r,])
        Jsub <- Jsub[,-1]
        Jac[rpos:(rpos+K.j[j]-1),cpos:(cpos+K.j[j]-2)] <- Jsub
        rpos <- rpos+K.j[j]
        cpos <- cpos+K.j[j]-1
      }
    }
    VCE.probs <- Jac %*% VCE.lo %*% t(Jac)
    
    maindiag <- diag(VCE.probs)
    maindiag[maindiag<0] <- 0 # error trap
    se.probs.vec <- sqrt(maindiag)
    se.probs <- list()
    for (j in 1:J) { se.probs[[j]] <- matrix(0,0,K.j[j]) }
    pos <- 1
    for (r in 1:R) {
      for (j in 1:J) { 
        se.probs[[j]] <- rbind(se.probs[[j]],se.probs.vec[pos:(pos+K.j[j]-1)])
        pos <- pos+K.j[j]
      }
    }
    
    # Variance of mixing proportions (priors) and coefficients (betas)
    if (R>1) {
      VCE.beta <- VCE[(1+sum(R*(K.j-1))):dim(VCE)[1],(1+sum(R*(K.j-1))):dim(VCE)[2]]
      se.beta <- matrix(sqrt(diag(VCE.beta)),nrow=ncol(x),ncol=(R-1))
      
      ptp <- array(NA,dim=c(R,R,N))
      for (n in 1:N) {
        ptp[,,n] <- -(prior[n,] %*% t(prior[n,]))
        diag(ptp[,,n]) <- prior[n,] * (1-prior[n,])
      }
      Jac.mix <- NULL
      for (r in 2:R) {
        for (l in 1:ncol(x)) {
          Jac.mix <- cbind(Jac.mix,colMeans(t(ptp[,r,]) * x[,l]))
        }
      }
      VCE.mix <- Jac.mix %*% VCE.beta %*% t(Jac.mix)
      se.mix <- sqrt(diag(VCE.mix))
    } else {
      VCE.beta <- se.beta <- se.mix <- NULL
    }
    return( list(probs=se.probs,P=se.mix,b=se.beta,var.b=VCE.beta) )
  }

poLCA.unvectorize <-
  function(vp) {
    probs <- list()
    idx <- c(0,cumsum(vp$numChoices*vp$classes))
    for (i in 1:length(vp$numChoices)){
      probs[[i]] <- matrix(vp$vecprobs[(idx[i]+1):idx[i+1]],nrow=vp$classes,byrow=TRUE)
    }
    return(probs)
  }

poLCA.compress <-
  function(y) {
    ym.sorted <- y[do.call(order,data.frame(y)),]
    vars <- ncol(ym.sorted)
    datamat <- ym.sorted[1,]
    freq <- 1
    curpos <- 1
    for (i in 2:nrow(ym.sorted)) {
      if (sum(ym.sorted[i,] == ym.sorted[i-1,])==vars) {
        freq[curpos] <- freq[curpos]+1
      } else {
        datamat <- rbind(datamat,ym.sorted[i,])
        freq <- c(freq,1)
        curpos <- curpos+1
      }
    }
    rownames(datamat) <- c(1:length(freq))
    ret <- list(datamat=datamat,freq=freq)
    return(ret)
  }

poLCA.probHat.C <-
  function(rgivy,y,vp) {
    ret <-  .C("probhat",
               as.integer(t(y)),
               as.double(t(rgivy)),
               as.integer(length(vp$numChoices)),
               as.integer(dim(y)[1]),
               as.integer(vp$numChoices),
               as.integer(vp$classes),
               ph = double(sum(vp$numChoices)*vp$classes)
    )
    return(ret$ph)
  }
