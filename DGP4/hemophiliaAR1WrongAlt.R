rm(list=ls())
library(Rcpp)
library(rlang)
library(permute)
library(hemophiliaAR1Hetero)

args <- commandArgs(trailingOnly = TRUE)
nSeed <- as.integer(args[1])

#####
set.seed(nSeed)
# D1
#####
Ts1 <- 120
alpha1 <- matrix(60,nrow = 1,ncol = 1)
n1 <- 10
theta1 <- c(20,1,10,-0.1,8)

# wrong D
Tsw <- 120
alphaw <- matrix(60,nrow = 1,ncol = 1)
nw <- 10
thetaw <- c(10,4,20,-1,12)

tList1 <- list()
nList1 <- list()
XList1 <- list()
YList1 <- list()

# First external data
for(i in 1:n1){
  ni <- rpois(1,lambda = 20)
  tList1[[i]] <- matrix(sort(runif(ni,min = 0,max = Ts1)),ncol = 1)
  nList1[[i]] <- ni
}

XList1 <- getXList(tList1,alpha1[1,1],Ts1)
thetaMat1 <- matrix(0,nrow = n1 + nw,ncol = 5)
SigMat1 <- diag(c(abs(c(25,0.1,9,0.1,4))))

for(i in 1:n1){
  ni <- length(tList1[[i]])
  thetai <- tmvtnorm::rtmvnorm(1,mean = theta1,sigma = SigMat1)
  thetaMat1[i,] <- thetai
  epsMat <- arMat(rho = 0.5,tListi = tList1[[i]])
  eps <- mvtnorm::rmvnorm(1,mean = rep(0,ni),sigma = 4*epsMat)
  
  YList1[[i]] <- matrix(c(XList1[[i]]%*%matrix(thetai,ncol = 1)) + c(eps), ncol = 1) 
  
}
# Do z-transform
phaseExt1 <- matrix(0,nrow = sum(unlist(nList1)),ncol = 4)

n1count <- 0
for(i in 1:n1){
  ni <- nList1[[i]]
  phaseExt1[(n1count + 1):(n1count + ni),1] <- as.numeric(rep(i,ni))
  phaseExt1[(n1count + 1):(n1count + ni),2] <- as.numeric(tList1[[i]])
  phaseExt1[(n1count + 1):(n1count + ni),3] <- as.numeric(YList1[[i]])
  phaseExt1[(n1count + 1):(n1count + ni),4] <- rep("External",ni)
  n1count <- n1count + ni
}

phaseExt1 <- as.data.frame(phaseExt1)
colnames(phaseExt1) <- c("psuedo_id","TimeW","Value","Phase")
phaseExt1$psuedo_id <- as.numeric(phaseExt1$psuedo_id)
phaseExt1$TimeW <- as.numeric(phaseExt1$TimeW)
phaseExt1$Value <- as.numeric(phaseExt1$Value)

phaseExt1$Value <- (phaseExt1$Value - mean(phaseExt1$Value[which(phaseExt1$TimeW <= 10)]))/
  sd(phaseExt1$Value[which(phaseExt1$TimeW <= 10)])

# Second external data
tList1 <- list()
nList1 <- list()
XList1 <- list()
YList1 <- list()

for(i in 1:nw){
  ni <- rpois(1,lambda = 20)
  tList1[[i]] <- matrix(sort(runif(ni,min = 0,max = Tsw)),ncol = 1)
  nList1[[i]] <- ni
}

XList1 <- getXList(tList1,alphaw[1,1],Tsw)
SigMatw <- diag(c(abs(c(25,0.1,9,0.1,4))))

for(i in 1:nw){
  ni <- length(tList1[[i]])
  thetai <- tmvtnorm::rtmvnorm(1,mean = thetaw,sigma = SigMatw)
  thetaMat1[i,] <- thetai
  
  epsMat <- arMat(rho = 0.5,tListi = tList1[[i]])
  eps <- mvtnorm::rmvnorm(1,mean = rep(0,ni),sigma = 4*epsMat)
  
  YList1[[i]] <- matrix(c(XList1[[i]]%*%matrix(thetai,ncol = 1)) + c(eps), ncol = 1) 
  
}

# Do z-transform
phaseExt2 <- matrix(0,nrow = sum(unlist(nList1)),ncol = 4)

n1count <- 0
for(i in 1:nw){
  ni <- nList1[[i]]
  phaseExt2[(n1count + 1):(n1count + ni),1] <- as.numeric(rep(i + n1,ni))
  phaseExt2[(n1count + 1):(n1count + ni),2] <- as.numeric(tList1[[i]])
  phaseExt2[(n1count + 1):(n1count + ni),3] <- as.numeric(YList1[[i]])
  phaseExt2[(n1count + 1):(n1count + ni),4] <- rep("External",ni)
  n1count <- n1count + ni
}

phaseExt2 <- as.data.frame(phaseExt2)
colnames(phaseExt2) <- c("psuedo_id","TimeW","Value","Phase")
phaseExt2$psuedo_id <- as.numeric(phaseExt2$psuedo_id)
phaseExt2$TimeW <- as.numeric(phaseExt2$TimeW)
phaseExt2$Value <- as.numeric(phaseExt2$Value)

phaseExt2$Value <- (phaseExt2$Value - mean(phaseExt2$Value[which(phaseExt2$TimeW <= 10)]))/
  sd(phaseExt2$Value[which(phaseExt2$TimeW <= 10)])

phaseExt <- rbind(phaseExt1,phaseExt2)

#####
# D2
Ts2 <- 120
alpha2 <- matrix(60,nrow = 1,ncol = 1)
n2 <- 40
theta2 <- c(20,1,10,-0.1,8)

tList2 <- list()
nList2 <- list()
XList2 <- list()
YList2 <- list()

for(i in 1:n2){
  ni <- rpois(1,lambda = 20)
  # truncated before 40
  tList2[[i]] <- matrix(sort(runif(ni,min = 0,max = 60)),ncol = 1)
  nList2[[i]] <- ni
}

XList2 <- getXList(tList2,alpha2[1,1],Ts2)
thetaMat2 <- matrix(0,nrow = n2,ncol = 5)
SigMat2 <- diag(c(abs(c(25,0.1,9,0.1,4))))

for(i in 1:n2){
  ni <- length(tList2[[i]])
  thetai <- tmvtnorm::rtmvnorm(1,mean = theta2,sigma = SigMat2)
  thetaMat2[i,] <- thetai
  
  epsMat <- arMat(rho = 0.5,tListi = tList2[[i]])
  eps <- mvtnorm::rmvnorm(1,mean = rep(0,ni),sigma = 4*epsMat)
  
  YList2[[i]] <- matrix(c(XList2[[i]]%*%matrix(thetai,ncol = 1)) + c(eps), ncol = 1) 
  
}

phaseInt <- matrix(0,nrow = sum(unlist(nList2)),ncol = 4)

n2count <- 0
for(i in 1:n2){
  ni <- nList2[[i]]
  phaseInt[(n2count + 1):(n2count + ni),1] <- as.numeric(rep(i,ni))
  phaseInt[(n2count + 1):(n2count + ni),2] <- as.numeric(tList2[[i]])
  phaseInt[(n2count + 1):(n2count + ni),3] <- as.numeric(YList2[[i]])
  phaseInt[(n2count + 1):(n2count + ni),4] <- rep("Internal",ni)
  n2count <- n2count + ni
}

phaseInt <- as.data.frame(phaseInt)
colnames(phaseInt) <- c("psuedo_id","TimeW","Value","Phase")
phaseInt$psuedo_id <- as.numeric(phaseInt$psuedo_id)
phaseInt$TimeW <- as.numeric(phaseInt$TimeW)
phaseInt$Value <- as.numeric(phaseInt$Value)

# conduct z-transform
mVal <- mean(phaseInt$Value[which(phaseInt$TimeW <= 10)])
sdVal <- sd(phaseInt$Value[which(phaseInt$TimeW <= 10)])

phaseInt$Value <- (phaseInt$Value - mean(phaseInt$Value[which(phaseInt$TimeW <= 10)]))/
  sd(phaseInt$Value[which(phaseInt$TimeW <= 10)])

tList1Sub <- list()
YList1Sub <- list()
tList2Sub <- list()
YList2Sub <- list()

TsThres <- 60
nExt <- n1 + nw
count <- 1
for(i in 1:nExt){
  tList1[[i]] <- matrix(phaseExt$TimeW[which(phaseExt$psuedo_id == i)],ncol = 1)
  YList1[[i]] <- matrix(phaseExt$Value[which(phaseExt$psuedo_id == i)],ncol = 1)
  
  idx <- which(tList1[[i]] <= TsThres)
  if(length(idx) > 0){
    tList1Sub[[count]] <- tList1[[i]][idx,,drop = FALSE]
    YList1Sub[[count]] <- YList1[[i]][idx,,drop = FALSE]
    count <- count + 1
  }
}

nInt <- n2
count <- 1
for(i in 1:nInt){
  tList2[[i]] <- matrix(phaseInt$TimeW[which(phaseInt$psuedo_id == i)],ncol = 1)
  YList2[[i]] <- matrix(phaseInt$Value[which(phaseInt$psuedo_id == i)],ncol = 1)
  
  idx <- which(tList2[[i]] <= TsThres)
  if(length(idx) > 0){
    tList2Sub[[count]] <- tList2[[i]][idx,,drop = FALSE]
    YList2Sub[[count]] <- YList2[[i]][idx,,drop = FALSE]
    count <- count + 1
  }
}

set.seed(1)
niter <- 2e3
nburn <- 1e3
niterOut <- 1e3
nSamp <- 1e3

Ts <- 120
beta0 <- matrix(0,nrow = 5,ncol = 1)
nu0 <-  1e-2
Psi0 <- diag(5) * 1e2
a0 <- 1e-2
b0 <- 1e-2

# train both subset
betaStarSub <- matrix(0,nrow = 5,ncol = 1)
SigmaStarSub <- diag(5)
sig201Sub <- matrix(1,nrow = 1,ncol = 1)
sig202Sub <- matrix(1,nrow = 1,ncol = 1)
alphaSub <- matrix(60,nrow = 1,ncol = 1)
rhoSub <- matrix(0.2,nrow = 1,ncol = 1)

alphaCollectSub <- list()
betaStarCollectSub <- list()
SigmaStarCollectSub <- list()
sig201CollectSub <- list()
sig202CollectSub <- list()
rhoCollectSub <- list()

# and complete data
betaStar <- matrix(0,nrow = 5,ncol = 1)
SigmaStar <- diag(5)

sig201 <- matrix(1,nrow = 1,ncol = 1)
sig202 <- matrix(1,nrow = 1,ncol = 1)

alpha <- matrix(60,nrow = 1,ncol = 1)
rho <- matrix(0.2,nrow = 1,ncol = 1)

alphaCollect <- list()
betaStarCollect <- list()
SigmaStarCollect <- list()
sig201Collect <- list()
sig202Collect <- list()
rhoCollect <- list()

SigRhoiList1 <- list()
SigRhoiList1Sub <- list()

for(i in 1:nExt){
  SigRhoiList1[[i]] <- arMat(rho = rho[1,1],tListi = tList1[[i]])
  SigRhoiList1Sub[[i]] <- arMat(rho = rhoSub[1,1],tListi = tList1Sub[[i]])
}

SigRhoiList2 <- list()
SigRhoiList2Sub <- list()

for(i in 1:nInt){
  SigRhoiList2[[i]] <- arMat(rho = rho[1,1],tListi = tList2[[i]])
  SigRhoiList2Sub[[i]] <- arMat(rho = rhoSub[1,1],tListi = tList2Sub[[i]])
}

####
# initialization
for(i in 1:niter){
  
  XList1 <- getXList(tList = tList1, alpha = alpha[1,1], Ts = Ts)
  
  XList2 <- getXList(tList = tList2, alpha = alpha[1,1], Ts = Ts)
  
  XList1Sub <- getXList(tList = tList1Sub, alpha = alphaSub[1,1], Ts = Ts)
  
  XList2Sub <- getXList(tList = tList2Sub, alpha = alphaSub[1,1], Ts = Ts)
  
  # update everything except for alpha
  getSampMixEffectHetero(YListS = YList1, XListS = XList1,
                         YListT = YList2, XListT = XList2,
                         beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                         a0 = a0, b0 = b0, 
                         SigRhoiListS = SigRhoiList1, SigRhoiListT = SigRhoiList2, 
                         betaStar = betaStar, SigmaStar = SigmaStar, 
                         sig20S = sig201, sig20T = sig202)
  
  getSampMixEffectHetero(YListS = YList1Sub,XListS = XList1Sub,
                         YListT = YList2Sub,XListT = XList2Sub,
                         beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                         a0 = a0, b0 = b0, 
                         SigRhoiListS = SigRhoiList1Sub, SigRhoiListT = SigRhoiList2Sub, 
                         betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                         sig20S = sig201Sub, sig20T = sig202Sub)
  
  # update Alpha
  getAlphaMixEffectHetero(YListS = YList1, XListS = XList1, tListS = tList1, 
                          YListT = YList2, XListT = XList2, tListT = tList2, 
                          Ts = Ts, alpha = alpha, 
                          SigRhoiListS = SigRhoiList1, SigRhoiListT = SigRhoiList2, 
                          betaStar = betaStar, SigmaStar = SigmaStar, 
                          sig20S = sig201, sig20T = sig202, deltaAlpha = Ts/60)
  
  getAlphaMixEffectHetero(YListS = YList1Sub, XListS = XList1Sub, tListS = tList1Sub, 
                          YListT = YList2Sub, XListT = XList2Sub, tListT = tList2Sub, 
                          Ts = Ts, alpha = alphaSub, 
                          SigRhoiListS = SigRhoiList1Sub, SigRhoiListT = SigRhoiList2Sub, 
                          betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                          sig20S = sig201Sub, sig20T = sig202Sub, deltaAlpha = Ts/60)
  
  # update rho
  SigRhoiList1 <- getRhoMixEffect(YList = YList1, tList = tList1, Ts = Ts, 
                                  rho = rho, SigRhoiList = SigRhoiList1, 
                                  alpha = alpha, betaStar = betaStar, SigmaStar = SigmaStar, 
                                  sig20 = sig201, deltaRho = 0.03)
  
  SigRhoiList2 <- getRhoMixEffect(YList = YList2, tList = tList2, Ts = Ts, 
                                  rho = rho, SigRhoiList = SigRhoiList2, 
                                  alpha = alpha, betaStar = betaStar, SigmaStar = SigmaStar, 
                                  sig20 = sig202, deltaRho = 0.03)
  # Sub
  SigRhoiList1Sub <- getRhoMixEffect(YList = YList1Sub, tList = tList1Sub, Ts = Ts,
                                     rho = rhoSub, SigRhoiList = SigRhoiList1Sub,
                                     alpha = alphaSub, betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                                     sig20 = sig201Sub, deltaRho = 0.03)
  
  SigRhoiList2Sub <- getRhoMixEffect(YList = YList2Sub, tList = tList2Sub, Ts = Ts,
                                     rho = rhoSub, SigRhoiList = SigRhoiList2Sub,
                                     alpha = alphaSub, betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                                     sig20 = sig202Sub, deltaRho = 0.03)
  
  betaStarCollect[[i]] <- duplicate(betaStar,shallow = TRUE)
  SigmaStarCollect[[i]] <- duplicate(SigmaStar,shallow = TRUE)
  sig201Collect[[i]] <- duplicate(sig201,shallow = TRUE)
  sig202Collect[[i]] <- duplicate(sig202,shallow = TRUE)
  alphaCollect[[i]] <- duplicate(alpha,shallow = TRUE)
  rhoCollect[[i]] <- duplicate(rho,shallow = TRUE)
  # Sub
  betaStarCollectSub[[i]] <- duplicate(betaStarSub,shallow = TRUE)
  SigmaStarCollectSub[[i]] <- duplicate(SigmaStarSub,shallow = TRUE)
  sig201CollectSub[[i]] <- duplicate(sig201Sub,shallow = TRUE)
  sig202CollectSub[[i]] <- duplicate(sig202Sub,shallow = TRUE)
  alphaCollectSub[[i]] <- duplicate(alphaSub,shallow = TRUE)
  rhoCollectSub[[i]] <- duplicate(rhoSub,shallow = TRUE)
}

####
# data selection
set.seed(1)
p1L <- c()
Z1L <- list()

alphaCollectSub <- list()
betaStarCollectSub <- list()
SigmaStarCollectSub <- list()
sig201CollectSub <- list()
sig202CollectSub <- list()
rhoCollectSub <- list()

alphaCollect <- list()
betaStarCollect <- list()
SigmaStarCollect <- list()
sig201Collect <- list()
sig202Collect <- list()
rhoCollect <- list()

# Z1 <- sample(c(0,1),size = nExt,replace = TRUE)
Z1 <- c(rep(1,nExt/2),rep(0,nExt/2))
cpst <- 1

pb = txtProgressBar(min = 0, max = niterOut, style = 3)

for(i in 1:niterOut){
  ## when phase 1 is external
  Z0Trans <- Z1
  Z1Trans <- Z1
  
  ratio <- runif(1,0,0.5)
  idx0 <- which(Z0Trans == 0)
  r0 <- ceiling(length(idx0) * ratio)
  if(length(idx0) > 1){
    idx0New <- sample(idx0,r0,replace = FALSE)
  }else{
    idx0New <- idx0
  }
  Z0Trans[idx0New] <- 1
  idxZ0Trans <- which(Z0Trans == 1)
  
  idx1 <- which(Z1Trans == 1)
  r1 <- ceiling(length(idx1) * ratio)
  if(length(idx1) > 1){
    idx1New <- sample(idx1,r1,replace = FALSE)
    Z1Trans[idx1New] <- 0
  }else{
    # if too few samples remain, try increasing the selected external data
    idx0 <- which(Z1Trans == 0)
    r0 <- ceiling(length(idx0) * ratio)
    idx0New <- sample(idx0,r0,replace = FALSE)
    Z1Trans[idx0New] <- 1
  }
  idxZ1Trans <- which(Z1Trans == 1)
  
  idxOrigin <- which(Z1 == 1)
  
  ###
  idxZ0Trans <- which(Z0Trans == 1)
  idxZ1Trans <- which(Z1Trans == 1)
  idxOrigin <- which(Z1 == 1)
  
  set.seed(1)
  tList1SubNew <- tList1Sub[idxOrigin]
  YList1SubNew <- YList1Sub[idxOrigin]
  
  mlOrigin <- getMLCondParamsHetero(YListS = YList1SubNew, tListS = tList1SubNew,
                                    YListT = YList2Sub, tListT = tList2Sub,
                                    Ts = Ts, nburn = nburn, niter = niter, deltaAlpha = Ts/60,
                                    beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                                    a0 = a0, b0 = b0, rho = rhoSub, 
                                    alpha = alphaSub, betaStar = betaStarSub, 
                                    SigmaStar = SigmaStarSub, 
                                    sig20S = sig201Sub, sig20T = sig202Sub, deltaRho = 0.03)
  
  set.seed(1)
  tList1SubNew <- tList1Sub[idxZ0Trans]
  YList1SubNew <- YList1Sub[idxZ0Trans]
  
  ml0 <- getMLCondParamsHetero(YListS = YList1SubNew, tListS = tList1SubNew,
                               YListT = YList2Sub, tListT = tList2Sub,
                               Ts = Ts, nburn = nburn, niter = niter, deltaAlpha = Ts/60,
                               beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                               a0 = a0, b0 = b0, rho = rhoSub, 
                               alpha = alphaSub, betaStar = betaStarSub, 
                               SigmaStar = SigmaStarSub, 
                               sig20S = sig201Sub, sig20T = sig202Sub, deltaRho = 0.03)
  
  set.seed(1)
  tList1SubNew <- tList1Sub[idxZ1Trans]
  YList1SubNew <- YList1Sub[idxZ1Trans]
  
  ml1 <- getMLCondParamsHetero(YListS = YList1SubNew, tListS = tList1SubNew,
                               YListT = YList2Sub, tListT = tList2Sub, 
                               Ts = Ts, nburn = nburn, niter = niter, deltaAlpha = Ts/60,
                               beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                               a0 = a0, b0 = b0, rho = rhoSub, 
                               alpha = alphaSub, betaStar = betaStarSub, 
                               SigmaStar = SigmaStarSub, 
                               sig20S = sig201Sub, sig20T = sig202Sub, deltaRho = 0.03)
  
  logprob <- c(cpst * (ml0 - mlOrigin),
               cpst * (ml1 - mlOrigin),
               0)
  
  prob <- exp(logprob - max(logprob))
  
  label <- sample(1:3,1,prob = prob)
  
  if(label == 1){
    Z1 <- Z0Trans
  }else if(label == 2){
    Z1 <- Z1Trans
  }else{
    # pass
  }
  
  tList1New <- tList1[which(Z1 == 1)]
  YList1New <- YList1[which(Z1 == 1)]
  
  if(length(tList1New) > 0){
    SigRhoiList1 <- list()
    for(iAr in 1:length(tList1New)){
      SigRhoiList1[[iAr]] <- arMat(rho = rho[1,1],tListi = tList1New[[iAr]])
    }
    
    tList1NewSub <- tList1Sub[which(Z1 == 1)]
    YList1NewSub <- YList1Sub[which(Z1 == 1)]
    
    SigRhoiList1Sub <- list()
    for(iAr in 1:length(tList1NewSub)){
      SigRhoiList1Sub[[iAr]] <- arMat(rho = rhoSub[1,1],tListi = tList1NewSub[[iAr]])
    }    
  }
  
  
  for(j in 1:nSamp){
    
    XList1New <- getXList(tList = tList1New, alpha = alpha[1,1], Ts = Ts)
    
    XList2 <- getXList(tList = tList2, alpha = alpha[1,1], Ts = Ts)
    
    XList1NewSub <- getXList(tList = tList1NewSub, alpha = alphaSub[1,1], Ts = Ts)
    
    XList2Sub <- getXList(tList = tList2Sub, alpha = alphaSub[1,1], Ts = Ts)
    
    # update everything except for alpha
    getSampMixEffectHetero(YListS = YList1New, XListS = XList1New,
                           YListT = YList2, XListT = XList2,
                           beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                           a0 = a0, b0 = b0, 
                           SigRhoiListS = SigRhoiList1, SigRhoiListT = SigRhoiList2, 
                           betaStar = betaStar, SigmaStar = SigmaStar, 
                           sig20S = sig201, sig20T = sig202)
    
    getSampMixEffectHetero(YListS = YList1NewSub,XListS = XList1NewSub,
                           YListT = YList2Sub,XListT = XList2Sub,
                           beta0 = beta0, nu0 = nu0, Psi0 = Psi0,
                           a0 = a0, b0 = b0, 
                           SigRhoiListS = SigRhoiList1Sub, SigRhoiListT = SigRhoiList2Sub, 
                           betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                           sig20S = sig201Sub, sig20T = sig202Sub)
    
    # update Alpha
    getAlphaMixEffectHetero(YListS = YList1New, XListS = XList1New, tListS = tList1New, 
                            YListT = YList2, XListT = XList2, tListT = tList2, 
                            Ts = Ts, alpha = alpha, 
                            SigRhoiListS = SigRhoiList1, SigRhoiListT = SigRhoiList2, 
                            betaStar = betaStar, SigmaStar = SigmaStar, 
                            sig20S = sig201, sig20T = sig202, deltaAlpha = Ts/60)
    
    getAlphaMixEffectHetero(YListS = YList1NewSub, XListS = XList1NewSub, tListS = tList1NewSub, 
                            YListT = YList2Sub, XListT = XList2Sub, tListT = tList2Sub, 
                            Ts = Ts, alpha = alphaSub, 
                            SigRhoiListS = SigRhoiList1Sub, SigRhoiListT = SigRhoiList2Sub, 
                            betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                            sig20S = sig202Sub, sig20T = sig201Sub, deltaAlpha = Ts/60)
    
    # update rho
    SigRhoiList1 <- getRhoMixEffect(YList = YList1New, tList = tList1New, Ts = Ts, 
                                    rho = rho, SigRhoiList = SigRhoiList1, 
                                    alpha = alpha, betaStar = betaStar, SigmaStar = SigmaStar, 
                                    sig20 = sig201, deltaRho = 0.03)
    
    SigRhoiList2 <- getRhoMixEffect(YList = YList2, tList = tList2, Ts = Ts, 
                                    rho = rho, SigRhoiList = SigRhoiList2, 
                                    alpha = alpha, betaStar = betaStar, SigmaStar = SigmaStar, 
                                    sig20 = sig202, deltaRho = 0.03)
    # Sub
    SigRhoiList1Sub <- getRhoMixEffect(YList = YList1NewSub, tList = tList1NewSub, Ts = Ts,
                                       rho = rhoSub, SigRhoiList = SigRhoiList1Sub,
                                       alpha = alphaSub, betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                                       sig20 = sig201Sub, deltaRho = 0.03)
    
    SigRhoiList2Sub <- getRhoMixEffect(YList = YList2Sub, tList = tList2Sub, Ts = Ts,
                                       rho = rhoSub, SigRhoiList = SigRhoiList2Sub,
                                       alpha = alphaSub, betaStar = betaStarSub, SigmaStar = SigmaStarSub, 
                                       sig20 = sig202Sub, deltaRho = 0.03)
    
  }
  
  alphaCollect[[i]] <- duplicate(alpha,shallow = TRUE)
  betaStarCollect[[i]] <- duplicate(betaStar,shallow = TRUE)
  SigmaStarCollect[[i]] <- duplicate(SigmaStar,shallow = TRUE)
  sig201Collect[[i]] <- duplicate(sig201,shallow = TRUE)
  sig202Collect[[i]] <- duplicate(sig202,shallow = TRUE)
  rhoCollect[[i]] <- duplicate(rho,shallow = TRUE)
  
  alphaCollectSub[[i]] <- duplicate(alphaSub,shallow = TRUE)
  betaStarCollectSub[[i]] <- duplicate(betaStarSub,shallow = TRUE)
  SigmaStarCollectSub[[i]] <- duplicate(SigmaStarSub,shallow = TRUE)
  sig201CollectSub[[i]] <- duplicate(sig201Sub,shallow = TRUE)
  sig202CollectSub[[i]] <- duplicate(sig202Sub,shallow = TRUE)
  rhoCollectSub[[i]] <- duplicate(rhoSub,shallow = TRUE)
  
  Z1L[[i]] <- duplicate(Z1,shallow = TRUE)
  p1L[i] <- sum(Z1L[[i]])/length(Z1L[[i]])

  setTxtProgressBar(pb,i)
}

close(pb)

result <- list()
result$alphaList <- alphaCollect
result$betaStarList <- betaStarCollect
result$SigmaStarList <- SigmaStarCollect
result$sig201List <- sig201Collect
result$sig202List <- sig202Collect
result$rhoList <- rhoCollect
result$p1L <- p1L
result$Z1L <- Z1L

fileName <- paste("hemophiliaWrongPlateauRho05n10Seed",nSeed,sep = "")
save(result,file = fileName)
