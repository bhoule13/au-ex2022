
################
# install.packages("MortCast")
 library(MortCast)


### first all the functions
# for the a0 in the first year of life 
Calculate_a0 <- function(m0,sex) {
  #Andreev-Kingkade formulas for computing a0 given m0
  # HMD full protocol Table 1 pg37
  #Males
  if(sex=="m"){
    if (m0<0.02300) {a0<-0.14929-(1.99545*m0)}
    if ((0.0230<= m0)&(m0<0.08307)) {a0<-0.02832+(3.26021*m0)}
    if (0.08307<= m0) {a0<-0.29915}
  }
  if (sex=="f"){
    #Females
    if (m0<0.01724) {a0<-0.14903-(2.05527*m0)}
    if ((0.01724 <= m0)&(m0< 0.06891)) {a0<-0.04667+(3.88089*m0)}
    if (0.06891<= m0) {a0<-0.31411}
  }
  return(a0) }


lifetable.mx<-function(mx,sex){
  
  N<-length(mx)
  AgeI<-rep(1,N)
  a0<-Calculate_a0(mx[1],sex)
  ax<-c(a0,rep(0.5,(N-1)))
  if(mx[N]>0){ax[N]<-1/mx[N]}
  qx<-mx/(1+(1-ax)*mx)
  qx[N]<-1             
  
  px<-1-qx
  
  lx<-1
  
  for(y in 1:(N-1)){          
    lx[y+1]<-lx[y]*px[y]
  }
  
  dx<-lx*qx
  dx[N]<-lx[N]
  
  Lx<-lx+(ax-AgeI)*dx
  Lx[N]<-lx[N]*ax[N]                 
  
  Tx<-c()
  for(y in 1:N){
    Tx[y]<-sum(Lx[y:N])
  }
  
  ex<-Tx/lx
  Age<-0:(N-1) 
  AgeI<-rep(1,N)
  ALL<-cbind(Age,AgeI,ax,mx,qx,lx,dx,Lx,Tx,ex)
  return(ALL)
}


Kannisto<-function(mx){
  names(mx) <- 0:100
  mx1<-kannisto(mx, est.ages = seq(80, 100), proj.ages = seq(85, 110))
  return(mx1)}


RxiMatrix<-function(PP,Cum){
  #PP<-B
  NumC<-dim(PP)[2]
  NumR<-dim(PP)[1]
  
  G<-PP
  FD<-colSums(PP)
  
  if (Cum==1){
    G<-t(apply(G,1,cumsum))
    FD<-cumsum(FD)
  } 
  
  FRx3<-t(matrix(rep(FD,(NumR)),NumC))
  FRx<-G/FRx3
  FRx[is.infinite(FRx)]<-0
  return(FRx)
}


AgeDecomp<-function(LT1,LT2){
  N<-dim(LT1)[1]
  
  lx1<-LT1[,6]
  lx2<-LT2[,6]
  Lx1<-LT1[,8]
  Lx2<-LT2[,8] 
  Tx1<-LT1[,9]
  Tx2<-LT2[,9] 

  AgeD<-lx1[-N]*((Lx2[-N]/lx2[-N])-(Lx1[-N]/lx1[-N]))+Tx2[-1]*((lx1[-N]/lx2[-N])-(lx1[-1]/lx2[-1]))
  AgeD<-c(AgeD,lx1[N]*((Tx2[N]/lx2[N])-(Tx1[N]/lx1[N])))
  return(AgeD)
  }


CauseDecomp<-function(LT1,LT2,ICD1,ICD2,SEX){
  a<-AgeDecomp(LT1,LT2)
  AgeGD<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  
  a<-LT1[,7]
  dx1<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  a<-LT1[,8]
  Lx1<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  mx1<-dx1/Lx1
  
  a<-LT2[,7]
  dx2<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  a<-LT2[,8]
  Lx2<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  mx2<-dx2/Lx2
  
  ICD1f<-ICD1[ICD1$Males==SEX,-1]
  R1<-ICD1f[-c(1,13),-c(1,2)]/matrix(rep(rowSums(ICD1f[-c(1,13),-c(1,2)]),6),11)
  ICD2f<-ICD2[ICD2$Males==SEX,-1]
  R2<-ICD2f[-c(1,13),-c(1,2)]/matrix(rep(rowSums(ICD2f[-c(1,13),-c(1,2)]),6),11)
  
  CauseD<-AgeGD*(R2*mx2-R1*mx1)/(mx2-mx1)
  return(colSums(CauseD))}

CauseDecomp2<-function(LT1,LT2,ICD1,ICD2){
  a<-AgeDecomp(LT1,LT2)
  AgeGD<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  
  a<-LT1[,7]
  dx1<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  a<-LT1[,8]
  Lx1<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  mx1<-dx1/Lx1
  
  a<-LT2[,7]
  dx2<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  a<-LT2[,8]
  Lx2<-c(a[1],sum(a[2:15]),sum(a[16:25]),sum(a[26:35]),sum(a[36:45]),sum(a[46:55]),sum(a[56:65]),sum(a[66:75]),sum(a[76:85]),sum(a[86:95]),sum(a[96:111]))
  mx2<-dx2/Lx2
  
  
  R1<-ICD1
  R2<-ICD2
  
  CauseD<-AgeGD*(R2*mx2-R1*mx1)/(mx2-mx1)
  return(colSums(CauseD))}

CauseDecomp3<-function(LT1,LT2,ICD1,ICD2){
  a<-AgeDecomp(LT1,LT2)
  AgeGD<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
           sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))
  
  a<-LT1[,7]
  dx1<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  a<-LT1[,8]
  Lx1<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  mx1<-dx1/Lx1
  
  a<-LT2[,7]
  dx2<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  a<-LT2[,8]
  Lx2<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  mx2<-dx2/Lx2
  
  
  R1<-ICD1
  R2<-ICD2
  
  CD<-AgeGD*(R2*mx2-R1*mx1)/(mx2-mx1)
  
  CauseD<-rbind(colSums(CD[1:6,]),colSums(CD[7:8,]),colSums(CD[9:10,]))
  
  return(CauseD)}


CauseDecomp4<-function(LT1,LT2,ICD1,ICD2){
  a<-AgeDecomp(LT1,LT2)
  AgeGD<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
           sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))
  
  a<-LT1[,7]
  dx1<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  a<-LT1[,8]
  Lx1<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  mx1<-dx1/Lx1
  
  a<-LT2[,7]
  dx2<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  a<-LT2[,8]
  Lx2<-c(sum(a[1:10]),sum(a[11:20]),sum(a[21:30]),sum(a[31:40]),sum(a[41:50]),
         sum(a[51:60]),sum(a[61:70]),sum(a[71:80]),sum(a[81:90]),sum(a[91:111]))  
  mx2<-dx2/Lx2
  
  
  R1<-ICD1
  R2<-ICD2
  
  CD<-AgeGD*(R2*mx2-R1*mx1)/(mx2-mx1)
  
  CauseD<-rbind(colSums(CD[2:6,]),colSums(CD[7:8,]),colSums(CD[9:10,]))
  
  return(CauseD)}

# functions to ungroup 
pclm <- function(y, C, X, lambda = 1, deg = 2, show = F){
  
  # Fit a PCLM (estimate b in ) E(y) = C %*% exp(X %*% b)
  
  # y = the vector of observed counts of length i
  
  # C = the composition matrix of dimension IxJ
  
  # X = the identity matrix of dimension JxJ; or B-spline basis
  
  # lambda = smoothing parameter
  
  # deg = order of differences of the components of b
  
  # show = indicates whether iteration details should be shown
  
  # Fit the penalized composite link model
  
  # Some preparations
  
  nx <- dim(X)[2]
  
  D <- diff(diag(nx), diff=deg)
  
  la2 <- sqrt(lambda)
  
  it <- 0
  
  bstart <- log(sum(y) / nx);
  
  b <- rep(bstart, nx);
  
  # Perform the iterations
  
  for (it in 1:50) {
    
    b0 <- b
    
    eta <- X %*% b
    
    gam <- exp(eta)
    
    mu <- C %*% gam
    
    w <- c(1 / mu, rep(la2, nx - deg))
    
    Gam <- gam %*% rep(1, nx)
    
    Q <- C %*% (Gam * X)
    
    z <- c(y - mu + Q %*% b, rep(0, nx - deg))
    
    Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F)
    
    b <- Fit$coef
    
    db <- max(abs(b - b0))
    
    if (show) cat(it, " ", db, "\n")
    
    if (db < 1e-6) break
    
  }
  
  cat(it, " ", db, "\n")
  
  # Regression diagnostic
  
  R <- t(Q) %*% diag(c(1 / mu)) %*% Q
  
  H <- solve(R + lambda * t(D) %*% D) %*% R
  
  fit <- list()
  
  fit$trace <- sum(diag(H))
  
  ok <- y > 0 & mu > 0
  
  fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
  
  fit$gamma <- gam
  
  fit$aic<- fit$dev + 2 * fit$trace
  
  fit$mu <- mu
  
  fit
  
}

UnGroup<-function(y){
  #### now ungrouping deaths age 100+ 
  ##  y<-D2020m
  
  x<-c(0:110)
  
  # Make C matrix and (trivial) basis B
  n<-length(y)
  m<-length(x)
  C1<-diag(1,100,100)
  C <- matrix(0, n, m)
  
  C[1:100, 1:100]<-C1
  C[101, 101:111] <- 1
  
  B <- diag(m)
  lambda <- 10^7
  mod <- pclm(y, C, B,lambda = lambda, deg = 2)
  
  DX<-as.numeric(as.character(c(y[1:100],y[101]*(mod$gamma[101:111]/sum(mod$gamma[101:111])))))
  return(DX)}

UnGroup22<-function(y){
  #### now ungrouping deaths age 100+ 
   #y<-D22[1:5,2]
   
   x<-c(0:100)
   
 # x<-c(0:44,45:64,65:74,75:84,85:110)
  
  # Make C matrix and (trivial) basis B
  n<-length(y)
  m<-length(x)
  C <- matrix(0, n, m)
  
  C[1, 1:45]<-1
  C[2, 46:65]<-1
  C[3, 66:75]<-1
  C[4, 76:85]<-1
  C[5, 86:101] <- 1
  
  B <- diag(m)
  lambda <- 10^7
  mod <- pclm(y, C, B,lambda = lambda, deg = 2)
  
  Dis<-function(h,g){
  Per<-h/sum(h)
  return(Per*y[g])
  }
  
  DX<-c(Dis(mod$gamma[1:45],1),Dis(mod$gamma[46:65],2),
        Dis(mod$gamma[66:75],3),Dis(mod$gamma[76:85],4),
        Dis(mod$gamma[86:101],5))
  
    return(DX)}


UnGroup23<-function(y){
  #### now ungrouping deaths age 100+ 
  #y<-D22[1:5,2]
  
  x<-c(0:100)
  
  # x<-c(0:44,45:64,65:74,75:84,85:110)
  
  # Make C matrix and (trivial) basis B
  n<-length(y)
  m<-length(x)
  C <- matrix(0, n, m)
  
  C[1, 1:5]<-1
  C[2, 6:10]<-1
  C[3, 11:15]<-1
  C[4, 16:20]<-1
  C[5, 21:25]<-1
  C[6, 26:30]<-1
  C[7, 31:35]<-1
  C[8, 36:40]<-1
  C[9, 41:45]<-1
  C[10, 46:50]<-1
  C[11, 51:55]<-1
  C[12, 56:60]<-1
  C[13, 61:65]<-1
  C[14, 66:70]<-1
  C[15, 71:75]<-1
  C[16, 76:80]<-1
  C[17, 81:85]<-1
  C[18, 86:90]<-1
  C[19, 91:95]<-1
  C[20, 96:100]<-1
  C[21, 101] <- 1
  
  B <- diag(m)
  lambda <- 10^7
  mod <- pclm(y, C, B,lambda = lambda, deg = 2)
  
  Dis<-function(h,g){
    Per<-h/sum(h)
    return(Per*y[g])
  }
  
  DX<-c(Dis(mod$gamma[1:5],1),Dis(mod$gamma[6:10],2),
        Dis(mod$gamma[11:15],3),Dis(mod$gamma[16:20],4),
        Dis(mod$gamma[21:25],5),Dis(mod$gamma[26:30],6),
        Dis(mod$gamma[31:35],7),Dis(mod$gamma[36:40],8),
        Dis(mod$gamma[41:45],9),Dis(mod$gamma[46:50],10),
        Dis(mod$gamma[51:55],11),Dis(mod$gamma[56:60],12),
        Dis(mod$gamma[61:65],13),Dis(mod$gamma[66:70],14),
        Dis(mod$gamma[71:75],15),Dis(mod$gamma[76:80],16),
        Dis(mod$gamma[81:85],17),Dis(mod$gamma[86:90],18),
        Dis(mod$gamma[91:95],19),Dis(mod$gamma[96:100],20),
        Dis(mod$gamma[101],21))
  
  return(DX)}

# CI for e0  
LTci <- function(Dx,mx,sex) {
  # Dx<-UnGroup(D2020m)[1:111]
  # mx<-R2020m
  # sex<-"m"
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  fun.ex <- function(mx) {
    return(lifetable.mx(mx,sex)[1,10])
  }
 
  exsim.ex <- apply(MX, 2, fun.ex)
 
  ## confidence interval
  CI.ex <- quantile(exsim.ex,
                    probs = c((1-0.95)/2,0.5,
                              1 - (1-0.95)/2))
 # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CI.ex
     # exsim=exsim
  )
  return(out)
}

# CI for differences in e0
LTci2 <- function(Dx,mx,Dx2,mx2,sex) {
  # Dx<-UnGroup(D2020m)[1:111]
  # mx<-R2020m
  # Dx2<-UnGroup(D2019m)[1:111]
  # mx2<-R2019m
  # sex<-"m"
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  
  Ntil2 <- round(Dx2/mx2)
  Y2 <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil2,
                                      mx2),
                               m, 1000))
  MX2 <- Y2/Ntil2
  
  fun.ex <- function(mx) {
    return(lifetable.mx(mx,sex)[1,10])
  }
  
  exsim.ex <- apply(MX, 2, fun.ex)
  exsim.ex2 <- apply(MX2, 2, fun.ex)
  
  ## confidence interval
  CI.ex <- quantile((exsim.ex-exsim.ex2),
                    probs = c((1-0.95)/2,0.5,
                              1 - (1-0.95)/2))
  # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CI.ex
    # exsim=exsim
  )
  return(out)
}

# CI for e0  
LTcie60 <- function(Dx,mx,sex) {
  # Dx<-UnGroup(D2020m)[1:111]
  # mx<-R2020m
  # sex<-"m"
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  fun.ex <- function(mx) {
    return(lifetable.mx(mx,sex)[61,10])
  }
  
  exsim.ex <- apply(MX, 2, fun.ex)
  
  ## confidence interval
  CI.ex <- quantile(exsim.ex,
                    probs = c((1-0.95)/2,0.5,
                              1 - (1-0.95)/2))
  # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CI.ex
    # exsim=exsim
  )
  return(out)
}

# CI for differences in e60
LTci2b <- function(Dx,mx,Dx2,mx2,sex) {
  # Dx<-UnGroup(D2020m)[1:111]
  # mx<-R2020m
  # Dx2<-UnGroup(D2019m)[1:111]
  # mx2<-R2019m
  # sex<-"m"
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  
  Ntil2 <- round(Dx2/mx2)
  Y2 <- suppressWarnings(matrix(rbinom(m * 1000,
                                       Ntil2,
                                       mx2),
                                m, 1000))
  MX2 <- Y2/Ntil2
  
  fun.ex <- function(mx) {
    return(lifetable.mx(mx,sex)[61,10])
  }
  
  exsim.ex <- apply(MX, 2, fun.ex)
  exsim.ex2 <- apply(MX2, 2, fun.ex)
  
  ## confidence interval
  CI.ex <- quantile((exsim.ex-exsim.ex2),
                    probs = c((1-0.95)/2,0.5,
                              1 - (1-0.95)/2))
  # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CI.ex
    # exsim=exsim
  )
  return(out)
}

# CI for Age-decomp in e0
LTci2c <- function(Dx,mx,Dx2,mx2,sex) {
  # Dx<-UnGroup(D2019m)[1:111]
  # mx<-R2019m
  # Dx2<-UnGroup(D2020m)[1:111]
  # mx2<-R2020m
  # sex<-"m"
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  
  Ntil2 <- round(Dx2/mx2)
  Y2 <- suppressWarnings(matrix(rbinom(m * 1000,
                                       Ntil2,
                                       mx2),
                                m, 1000))
  MX2 <- Y2/Ntil2
  
  fun.AD <- function(mxA,mxB) {
    A<-AgeDecomp(lifetable.mx(mxA,sex),lifetable.mx(mxB,sex))
    return(c(sum(A[1:55]),sum(A[56:85]),sum(A[86:111])))
  }
  
  MXA<-lapply(seq_len(ncol(MX)), function(x) MX[,x])
  
  MXB<-lapply(seq_len(ncol(MX2)), function(x) MX2[,x])
  
  exsim.AD <-mapply(fun.AD, mxA=MXA, mxB=MXB)
  
  Mat <- matrix(unlist(exsim.AD), 3)
  
  
  ## confidence interval
 CIlyl1<-c()
  for (i in 1:3){
    CI.lyli3 <- quantile(Mat[i,],
                         probs = c((1-0.95)/2,.5,
                                   1 - (1-0.95)/2),na.rm = TRUE)
    CIlyl1<-rbind(CIlyl1,CI.lyli3)} 
  # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CIlyl1
    # exsim=exsim
  )
  return(out)
}

# CI for Cause-decomp in e0
LTci2d <- function(Dx,mx,B,Dx2,mx2,B2,sex) {
  # Dx<-UnGroup(DD2019m)[1:111]
  # mx<-R2019m
  # Dx2<-UnGroup(D2020m)[1:111]
  # mx2<-R2020m
  # sex<-"m"
  # B<-
  # B2<-
 
  SEX<-c()
  SEX[sex=="m"]<-"Males"
  SEX[sex=="f"]<-"Females"
  
  b<-B[B$Males==SEX,-1]
  mxi<-b[-c(1,13),-c(1,2)]/matrix(rep(rowSums(b[-c(1,13),-c(1,2)]),6),11)
  b2<-B2[B2$Males==SEX,-1]
  mxi2<-b2[-c(1,13),-c(1,2)]/matrix(rep(rowSums(b2[-c(1,13),-c(1,2)]),6),11)
  
  mxi[is.na(mxi)]<-0
  mxi2[is.na(mxi2)]<-0
  Nmxi<-dim(mxi)[2]  # causes of death
  Ncol<-dim(mxi)[1]
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 10000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  
  Ntil2 <- round(Dx2/mx2)
  Y2 <- suppressWarnings(matrix(rbinom(m * 10000,
                                       Ntil2,
                                       mx2),
                                m, 1000))
  MX2 <- Y2/Ntil2
  
  #rmultinom
  BM<-list() 
  BMa<-BM
  BM2<-BM
  BMa2<-BM
  n=1000
  
  for (t in 1:Ncol){
    xx<-rmultinom(n, size =10000, prob = mxi[t,])/10000
    BM[[t]]<-lapply(seq_len(ncol(xx)), function(i) xx[,i])
    xx2<-rmultinom(n, size =10000, prob = mxi2[t,])/10000
    BM2[[t]]<-lapply(seq_len(ncol(xx2)), function(i) xx2[,i])
  }
  BMa<- Map(rbind,BM[[1]])
  BMa2<- Map(rbind,BM2[[1]])
  for (tt in 2:Ncol){
    BMa<- Map(rbind,BMa,BM[[tt]])   
    BMa2<- Map(rbind,BMa2,BM2[[tt]])
  }
  
  fun.CD <- function(mxA,mxB,BB1,BB2) {
    return(CauseDecomp2(lifetable.mx(mxA,sex),lifetable.mx(mxB,sex),BB1,BB2))
    }
  
  MXA<-lapply(seq_len(ncol(MX)), function(x) MX[,x])
  
  MXB<-lapply(seq_len(ncol(MX2)), function(x) MX2[,x])
  
  exsim.CD <-mapply(fun.CD, mxA=MXA, mxB=MXB,BB1=BMa,BB2=BMa2)
  
  Mat <- matrix(unlist(exsim.CD), 6)
  
  
  ## confidence interval
  CIlyl1<-c()
  for (i in 1:6){
    CI.lyli3 <- quantile(Mat[i,],
                         probs = c((1-0.95)/2,.5,
                                   1 - (1-0.95)/2),na.rm = TRUE)
    CIlyl1<-rbind(CIlyl1,CI.lyli3)} 
  # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CIlyl1
    # exsim=exsim
  )
  return(out)
}

# CI for Cause-decomp in e0
LTci2e <- function(Dx,mx,B,Dx2,mx2,B2,sex) {
  # Dx<-UnGroup(DD2019m)[1:111]
  # mx<-R2019m
  # Dx2<-UnGroup(D2020m)[1:111]
  # mx2<-R2020m
  # sex<-"m"
  # B<-IICD19m
  # B2<-CoDABS(CoD20m)
  
  SEX<-c()
  SEX[sex=="m"]<-"Males"
  SEX[sex=="f"]<-"Females"
  
  mxi<-B
  mxi2<-B2
  
  mxi[is.na(mxi)]<-0
  mxi2[is.na(mxi2)]<-0
  Nmxi<-dim(mxi)[2]  # causes of death
  Ncol<-dim(mxi)[1]
  
  m <- length(mx)
  Ntil <- round(Dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 10000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  
  Ntil2 <- round(Dx2/mx2)
  Y2 <- suppressWarnings(matrix(rbinom(m * 10000,
                                       Ntil2,
                                       mx2),
                                m, 1000))
  MX2 <- Y2/Ntil2
  
  #rmultinom
  BM<-list() 
  BMa<-BM
  BM2<-BM
  BMa2<-BM
  n=1000
  
  for (t in 1:Ncol){
    xx<-rmultinom(n, size =10000, prob = mxi[t,])/10000
    BM[[t]]<-lapply(seq_len(ncol(xx)), function(i) xx[,i])
    xx2<-rmultinom(n, size =10000, prob = mxi2[t,])/10000
    BM2[[t]]<-lapply(seq_len(ncol(xx2)), function(i) xx2[,i])
  }
  BMa<- Map(rbind,BM[[1]])
  BMa2<- Map(rbind,BM2[[1]])
  for (tt in 2:Ncol){
    BMa<- Map(rbind,BMa,BM[[tt]])   
    BMa2<- Map(rbind,BMa2,BM2[[tt]])
  }
  
  fun.CD <- function(mxA,mxB,BB1,BB2) {
    return(CauseDecomp2(lifetable.mx(mxA,sex),lifetable.mx(mxB,sex),BB1,BB2))
  }
  
  MXA<-lapply(seq_len(ncol(MX)), function(x) MX[,x])
  
  MXB<-lapply(seq_len(ncol(MX2)), function(x) MX2[,x])
  
  exsim.CD <-mapply(fun.CD, mxA=MXA, mxB=MXB,BB1=BMa,BB2=BMa2)
  
  Mat <- matrix(unlist(exsim.CD), 6)
  
  
  ## confidence interval
  CIlyl1<-c()
  for (i in 1:6){
    CI.lyli3 <- quantile(Mat[i,],
                         probs = c((1-0.95)/2,.5,
                                   1 - (1-0.95)/2),na.rm = TRUE)
    CIlyl1<-rbind(CIlyl1,CI.lyli3)} 
  # output
  out <- data.frame(
    #meanex=mean(exsim),
    CIex=CIlyl1
    # exsim=exsim
  )
  return(out)
}

##############################

#############################


####################### average 2017-2019


### DATA
#ANU
#setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/Australia Covid19/Re-do2022/Data")
#home
# setwd("C:/Users/u1019088/Desktop/Papers/Australia Covid19/Re-do22/Data")

P<-read.table("../data/PopO.csv",header=TRUE,fill=TRUE,sep=",")
PQ<-read.table("../data/PopQ.csv",header=TRUE,fill=TRUE,sep=",")
D<-read.table("../data/DeathsO.csv",header=TRUE,fill=TRUE,sep=",")
COD2<-read.table("../data/CausesNational3.csv",header=TRUE,fill=TRUE,sep=",", skip = 1)
DS<-read.table("../data/DeathsStates2.csv",header=TRUE,fill=TRUE,skip=2,sep=",")
#ANU
#setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/Australia Covid19/Re-do2022/Results")
#home
# setwd("C:/Users/u1019088/Desktop/Papers/Australia Covid19/Re-do22/Results")

# now i select the data by year and type (deaths D, Population P and rates r)
# this first just to make the file smaller (too many repetitions in columns)
POP<-P[,c(3:5,7,8)]
POPQ<-PQ[,c(3:5,7,8)]

A<-unlist(strsplit(as.character(POP$AGE..Age),":"))[seq(2,2*length(POP$AGE..Age),by=2)]
A[A==" 100+ years"]<-100
#Ages with A99 as the 100+ and TT as the total
POP$Age<-as.numeric(A)
POPb<-POP[(POP$REGION..Region=="AUS: Australia"),]

POPb$Time<-POPb$TIME_PERIOD..Time.Period

POPb$SEX_ABS<-as.character(POPb$SEX..Sex)
POPb$SEX_ABS[POPb$SEX_ABS=="2: Females"]<-"2"
POPb$SEX_ABS[POPb$SEX_ABS=="1: Males"]<-"1"
POPb$SEX_ABS<-as.numeric(POPb$SEX_ABS)
###

A<-unlist(strsplit(as.character(POPQ$AGE..Age),":"))[seq(2,2*length(POPQ$AGE..Age),by=2)]
A[A==" All ages"]<-999
#All ages is the total
POPQ$Age<-as.numeric(A)
POPbQ<-POPQ[(POPQ$REGION..Region=="AUS: Australia")&
            (POPQ$Age<101),]

POPbQ$Time<-unlist(strsplit(as.character(POPbQ$TIME_PERIOD..Time.Period),"-"))[seq(1,2*length(POPbQ$AGE..Age),by=2)]

POPbQ$SEX_ABS<-as.character(POPbQ$SEX..Sex)
POPbQ$SEX_ABS[POPbQ$SEX_ABS=="2: Females"]<-"2"
POPbQ$SEX_ABS[POPbQ$SEX_ABS=="1: Males"]<-"1"
POPbQ$SEX_ABS<-as.numeric(POPbQ$SEX_ABS)


###now with deaths
DEA<-D[,c(3:5,7,8)]

B<-unlist(strsplit(as.character(DEA$AGE..Age),":"))[seq(2,2*length(DEA$AGE..Age),by=2)]
B[B==" 100+ years"]<-100
B[B==" 105+ years"]<-105
#Ages with A99 as the 100+ and TT as the total
DEA$Age<-as.numeric(B)
DEAb<-DEA[(DEA$REGION..Region=="AUS: Australia"),]
DEAb<-DEAb[DEAb$Age<101,]
DEAb<-DEAb[DEAb$AGE..Age!="100: 100",]

DEAb$Time<-DEAb$TIME_PERIOD..Time.Period

DEAb$SEX_ABS<-as.character(DEAb$SEX..Sex)
DEAb$SEX_ABS[DEAb$SEX_ABS=="2: Females"]<-"2"
DEAb$SEX_ABS[DEAb$SEX_ABS=="1: Males"]<-"1"

# P2020m = population-P year

P2022m<-POPbQ[(POPbQ$Time==2022)&(POPbQ$SEX_ABS==1),]
P2022f<-POPbQ[(POPbQ$Time==2022)&(POPbQ$SEX_ABS==2),]

P2021m<-POPb[(POPb$Time==2021)&(POPb$SEX_ABS==1),]
P2021f<-POPb[(POPb$Time==2021)&(POPb$SEX_ABS==2),]
D2021m<-DEAb[(DEAb$Time==2021)&(DEAb$SEX_ABS==1),]
D2021f<-DEAb[(DEAb$Time==2021)&(DEAb$SEX_ABS==2),]

P2020m<-POPb[(POPb$Time==2020)&(POPb$SEX_ABS==1),]
P2020f<-POPb[(POPb$Time==2020)&(POPb$SEX_ABS==2),]
D2020m<-DEAb[(DEAb$Time==2020)&(DEAb$SEX_ABS==1),]
D2020f<-DEAb[(DEAb$Time==2020)&(DEAb$SEX_ABS==2),]

P2019m<-POPb[(POPb$Time==2019)&(POPb$SEX_ABS==1),]
P2019f<-POPb[(POPb$Time==2019)&(POPb$SEX_ABS==2),]
D2019m<-DEAb[(DEAb$Time==2019)&(DEAb$SEX_ABS==1),]
D2019f<-DEAb[(DEAb$Time==2019)&(DEAb$SEX_ABS==2),]

P2018m<-POPb[(POPb$Time==2018)&(POPb$SEX_ABS==1),]
P2018f<-POPb[(POPb$Time==2018)&(POPb$SEX_ABS==2),]
D2018m<-DEAb[(DEAb$Time==2018)&(DEAb$SEX_ABS==1),]
D2018f<-DEAb[(DEAb$Time==2018)&(DEAb$SEX_ABS==2),]

P2017m<-POPb[(POPb$Time==2017)&(POPb$SEX_ABS==1),]
P2017f<-POPb[(POPb$Time==2017)&(POPb$SEX_ABS==2),]
D2017m<-DEAb[(DEAb$Time==2017)&(DEAb$SEX_ABS==1),]
D2017f<-DEAb[(DEAb$Time==2017)&(DEAb$SEX_ABS==2),]


############## Ordering them

P2022m<-P2022m[order(P2022m$Age,decreasing=FALSE),5]
P2022f<-P2022f[order(P2022f$Age,decreasing=FALSE),5]


  D2021mb<-DS[(DS$Year==2021)&(DS$Sex=="Males"),12][1:21]
  D2021fb<-DS[(DS$Year==2021)&(DS$Sex=="Females"),12][1:21]
  
  D2022mb<-DS[(DS$Year==2022)&(DS$Sex=="Males"),12][1:21]
  D2022fb<-DS[(DS$Year==2022)&(DS$Sex=="Females"),12][1:21]
  
  D2022m<-UnGroup23(as.numeric(as.character(D2022mb)))
  D2022f<-UnGroup23(as.numeric(as.character(D2022fb)))
  D2021m<-UnGroup23(as.numeric(as.character(D2021mb)))
  D2021f<-UnGroup23(as.numeric(as.character(D2021fb)))



P2021m<-P2021m[order(P2021m$Age,decreasing=FALSE),5]
P2020m<-P2020m[order(P2020m$Age,decreasing=FALSE),5]
P2019m<-P2019m[order(P2019m$Age,decreasing=FALSE),5]
P2018m<-P2018m[order(P2018m$Age,decreasing=FALSE),5]
P2017m<-P2017m[order(P2017m$Age,decreasing=FALSE),5]

#D2021m<-D2021m[order(D2021m$Age,decreasing=FALSE),5]
D2020m<-D2020m[order(D2020m$Age,decreasing=FALSE),5]
D2019m<-D2019m[order(D2019m$Age,decreasing=FALSE),5]
D2018m<-D2018m[order(D2018m$Age,decreasing=FALSE),5]
D2017m<-D2017m[order(D2017m$Age,decreasing=FALSE),5]


P2021f<-P2021f[order(P2021f$Age,decreasing=FALSE),5]
P2020f<-P2020f[order(P2020f$Age,decreasing=FALSE),5]
P2019f<-P2019f[order(P2019f$Age,decreasing=FALSE),5]
P2018f<-P2018f[order(P2018f$Age,decreasing=FALSE),5]
P2017f<-P2017f[order(P2017f$Age,decreasing=FALSE),5]

#D2021f<-D2021f[order(D2021f$Age,decreasing=FALSE),5]
D2020f<-D2020f[order(D2020f$Age,decreasing=FALSE),5]
D2019f<-D2019f[order(D2019f$Age,decreasing=FALSE),5]
D2018f<-D2018f[order(D2018f$Age,decreasing=FALSE),5]
D2017f<-D2017f[order(D2017f$Age,decreasing=FALSE),5]



#sum(D2021m)+sum(D2021f)

#sum(P2021m)+sum(P2021f)
##  Kannisto smooth the death rates

## I combine the deaths and populations for the 3 years
##  I did not change the name of D2019 and P2019 in the rest of the 
## code but it refers now to the addition of the 3 years
dD2019m<-D2019m
dD2019f<-D2019f
dD2018m<-D2018m
dD2018f<-D2018f
dD2017m<-D2017m
dD2017f<-D2017f

D2019m<-D2019m+D2018m+D2017m
P2019m<-P2019m+P2018m+P2017m
D2019f<-D2019f+D2018f+D2017f
P2019f<-P2019f+P2018f+P2017f


D202012m<-D2020m+D2021m+D2022m
P202012m<-P2020m+P2021m+(P2022m)
D202012f<-D2020f+D2021f+D2022f
P202012f<-P2020f+P2021f+(P2022f)

R202012m<-Kannisto(D202012m/P202012m) 
R202012f<-Kannisto(D202012f/P202012f) 
  
R2022m<-Kannisto(D2022m/(P2022m)) 
R2022f<-Kannisto(D2022f/(P2022f))
R2021m<-Kannisto(D2021m/P2021m)
R2021f<-Kannisto(D2021f/P2021f)
R2020m<-Kannisto(D2020m/P2020m)
R2020f<-Kannisto(D2020f/P2020f)
R2019m<-Kannisto(D2019m/P2019m)
R2019f<-Kannisto(D2019f/P2019f)
R2018m<-Kannisto(D2018m/P2018m)
R2018f<-Kannisto(D2018f/P2018f)
R2017m<-Kannisto(D2017m/P2017m)
R2017f<-Kannisto(D2017f/P2017f)


#####################


#####################


#####################


## life tables including CI
## LT20m = life table-LT year-2020 males-m


LT2012m<-lifetable.mx(R202012m,"m")
CI2012m<-unlist(LTci(UnGroup(D202012m)[1:111],R202012m,"m"))

LT2012f<-lifetable.mx(R202012f,"f")
CI2012f<-unlist(LTci(UnGroup(D202012f)[1:111],R202012f,"f"))


LT22m<-lifetable.mx(R2022m,"m")
CI22m<-unlist(LTci(UnGroup(D2022m)[1:111],R2022m,"m"))


LT22f<-lifetable.mx(R2022f,"f")
CI22f<-unlist(LTci(UnGroup(D2022f)[1:111],R2022f,"f"))

LT21m<-lifetable.mx(R2021m,"m")
CI21m<-unlist(LTci(UnGroup(D2021m)[1:111],R2021m,"m"))


LT21f<-lifetable.mx(R2021f,"f")
CI21f<-unlist(LTci(UnGroup(D2021f)[1:111],R2021f,"f"))


LT20m<-lifetable.mx(R2020m,"m")
CI20m<-unlist(LTci(UnGroup(D2020m)[1:111],R2020m,"m"))


LT20f<-lifetable.mx(R2020f,"f")
CI20f<-unlist(LTci(UnGroup(D2020f)[1:111],R2020f,"f"))

LT19m<-lifetable.mx(R2019m,"m")
CI19m<-unlist(LTci(UnGroup(D2019m)[1:111],R2019m,"m"))

LT19f<-lifetable.mx(R2019f,"f")
CI19f<-unlist(LTci(UnGroup(D2019f)[1:111],R2019f,"f"))

# the Tables area arranged by value, followed by the ci-L low and then ci-H high
# first for females and then for males
Table1<-round(   cbind(
  rbind(c(LT19f[1,10],CI19f[c(1,3)]),
        c(LT20f[1,10],CI20f[c(1,3)]),
        c(LT21f[1,10],CI21f[c(1,3)]),
        c(LT22f[1,10],CI22f[c(1,3)]),
        c(LT2012f[1,10],CI2012f[c(1,3)])),
  rbind(c(LT19m[1,10],CI19m[c(1,3)]),
        c(LT20m[1,10],CI20m[c(1,3)]),
        c(LT21m[1,10],CI21m[c(1,3)]),
        c(LT22m[1,10],CI22m[c(1,3)]),
        c(LT2012m[1,10],CI2012m[c(1,3)])))
  ,2)





dimnames(Table1)[[1]]<-c("ave2017-19","2020", "2021","2022","ave2020-22")

dimnames(Table1)[[2]]<-c("mean e0 fem","CI_l fem","CI_h fem","mean e0 male","CI_l male","CI_h male")




########################
Dif2012m<-LT2012m[1,10]-LT19m[1,10]
CID2012m<-unlist(LTci2(UnGroup(D202012m)[1:111],R202012m,UnGroup(D2019m)[1:111],R2019m,"m"))
Dif2012f<- LT2012f[1,10]-LT19f[1,10] 
CID2012f<-unlist(LTci2(UnGroup(D202012f)[1:111],R202012f,UnGroup(D2019f)[1:111],R2019f,"f"))


Dif22m<-LT22m[1,10]-LT19m[1,10]
CID22m<-unlist(LTci2(UnGroup(D2022m)[1:111],R2022m,UnGroup(D2019m)[1:111],R2019m,"m"))
Dif22f<- LT22f[1,10]-LT19f[1,10] 
CID22f<-unlist(LTci2(UnGroup(D2022f)[1:111],R2022f,UnGroup(D2019f)[1:111],R2019f,"f"))

Dif21m<-LT21m[1,10]-LT19m[1,10]
CID21m<-unlist(LTci2(UnGroup(D2021m)[1:111],R2021m,UnGroup(D2019m)[1:111],R2019m,"m"))
Dif21f<- LT21f[1,10]-LT19f[1,10] 
CID21f<-unlist(LTci2(UnGroup(D2021f)[1:111],R2021f,UnGroup(D2019f)[1:111],R2019f,"f"))

Dif20m<-LT20m[1,10]-LT19m[1,10]
CID20m<-unlist(LTci2(UnGroup(D2020m)[1:111],R2020m,UnGroup(D2019m)[1:111],R2019m,"m"))
Dif20f<-LT20f[1,10]-LT19f[1,10]
CID20f<-unlist(LTci2(UnGroup(D2020f)[1:111],R2020f,UnGroup(D2019f)[1:111],R2019f,"m"))


Table2<-round(   
 rbind( 
   c(Dif20f[1],CID20f[c(1,3)],Dif20m[1],CID20m[c(1,3)]),
   c(Dif21f[1],CID21f[c(1,3)],Dif21m[1],CID21m[c(1,3)]),
   c(Dif22f[1],CID22f[c(1,3)],Dif22m[1],CID22m[c(1,3)]),
   c(Dif2012f[1],CID2012f[c(1,3)],Dif2012m[1],CID2012m[c(1,3)]))
  ,2)




dimnames(Table2)[[1]]<-c("2020vs2017-19","2021vs2017-19","2022vs2017-19","2020-22vs2017-19")



dimnames(Table2)[[2]]<-c("mean dif e0 fem","CI_l fem","CI_h fem","mean dif e0 male","CI_l male","CI_h male")





############################


########################


Dif22mb<-LT22m[1,10]-LT21m[1,10]
CID22mb<-unlist(LTci2(UnGroup(D2022m)[1:111],R2022m,UnGroup(D2021m)[1:111],R2021m,"m"))
Dif22fb<- LT22f[1,10]-LT21f[1,10] 
CID22fb<-unlist(LTci2(UnGroup(D2022f)[1:111],R2022f,UnGroup(D2021f)[1:111],R2021f,"f"))


Dif21mb<-LT21m[1,10]-LT20m[1,10]
CID21mb<-unlist(LTci2(UnGroup(D2021m)[1:111],R2021m,UnGroup(D2020m)[1:111],R2020m,"m"))
Dif21fb<- LT21f[1,10]-LT20f[1,10] 
CID21fb<-unlist(LTci2(UnGroup(D2021f)[1:111],R2021f,UnGroup(D2020f)[1:111],R2020f,"f"))



Table2b<-round(   
  rbind( 
    c(Dif20f[1],CID20f[c(1,3)],Dif20m[1],CID20m[c(1,3)]),
    c(Dif21fb[1],CID21fb[c(1,3)],Dif21mb[1],CID21mb[c(1,3)]),
    c(Dif22fb[1],CID22fb[c(1,3)],Dif22mb[1],CID22mb[c(1,3)]),
    c(Dif2012f[1],CID2012f[c(1,3)],Dif2012m[1],CID2012m[c(1,3)]))
  ,2)




dimnames(Table2b)[[1]]<-c("2020vs2017-19","2021vs2020","2022vs2021","2020-22vs2017-19")



dimnames(Table2b)[[2]]<-c("mean dif e0 fem","CI_l fem","CI_h fem","mean dif e0 male","CI_l male","CI_h male")





############# Age contribution


## age decomposition including CI
# males
cbind(LT21m[1,10]-LT19m[1,10],(LT20m[1,10]-LT19m[1,10]))

# 2019 to 2020
A<-AgeDecomp(LT19m,LT20m)
A2<-AgeDecomp(LT19m,LT21m)
A3<-AgeDecomp(LT19m,LT22m)
A4<-AgeDecomp(LT19m,LT2012m)

# age-groups = 0-59, 60-79, 80+
Agm<-c(sum(A[1:60]),sum(A[61:80]),sum(A[81:111]))
Agm2<-c(sum(A2[1:60]),sum(A2[61:80]),sum(A2[81:111]))
Agm3<-c(sum(A3[1:60]),sum(A3[61:80]),sum(A3[81:111]))
Agm4<-c(sum(A4[1:60]),sum(A4[61:80]),sum(A4[81:111]))

# and their CI
CIAgm<-t(matrix(unlist(LTci2c(UnGroup(D2019m)[1:111],R2019m,UnGroup(D2020m)[1:111],R2020m,"m")),3))
CIAgm2<-t(matrix(unlist(LTci2c(UnGroup(D2019m)[1:111],R2019m,UnGroup(D2021m)[1:111],R2021m,"m")),3))
CIAgm3<-t(matrix(unlist(LTci2c(UnGroup(D2019m)[1:111],R2019m,UnGroup(D2022m)[1:111],R2022m,"m")),3))

# now for females
cbind(LT21f[1,10]-LT19f[1,10],LT20f[1,10]-LT19f[1,10])
A<-AgeDecomp(LT19f,LT20f)
A2<-AgeDecomp(LT19f,LT21f)
A3<-AgeDecomp(LT19f,LT22f)
A4<-AgeDecomp(LT19f,LT2012f)


# age-groups = 0-59, 60-79, 80+
Agf<-c(sum(A[1:60]),sum(A[61:80]),sum(A[81:111]))
Agf2<-c(sum(A2[1:60]),sum(A2[61:80]),sum(A2[81:111]))
Agf3<-c(sum(A3[1:60]),sum(A3[61:80]),sum(A3[81:111]))
Agf4<-c(sum(A4[1:60]),sum(A4[61:80]),sum(A4[81:111]))


CIAgf<-t(matrix(unlist(LTci2c(UnGroup(D2019f)[1:111],R2019f,UnGroup(D2020f)[1:111],R2020f,"f")),3))
CIAgf2<-t(matrix(unlist(LTci2c(UnGroup(D2019f)[1:111],R2019f,UnGroup(D2021f)[1:111],R2021f,"f")),3))
CIAgf3<-t(matrix(unlist(LTci2c(UnGroup(D2019f)[1:111],R2019f,UnGroup(D2022f)[1:111],R2022f,"f")),3))


Table3<-round(   
  rbind(cbind(c(Agf,sum(Agf)),c(Agf2,sum(Agf2)),c(Agf3,sum(Agf3)),c(Agf4,sum(Agf4))),
        c(rep(0,4)), 
        cbind(c(Agm,sum(Agm)),c(Agm2,sum(Agm2)),c(Agm3,sum(Agm3)),c(Agm4,sum(Agm4))))
  ,2)



dimnames(Table3)[[1]]<-c("fem 0-59","fem 60-79","fem 80+","Total",0, "male 0-59","male 60-79","male 80+","Total")


dimnames(Table3)[[2]]<-c("2020vs2017-19", "2021vs2017-19", "2022vs2017-19", "2020-22vs2017-19")


############################


# males
# 2019 to 2020, 20 to 21 and 21 to 22
A2<-AgeDecomp(LT20m,LT21m)
A3<-AgeDecomp(LT21m,LT22m)

# age-groups = 0-59, 60-79, 80+
Agm2<-c(sum(A2[1:60]),sum(A2[61:80]),sum(A2[81:111]))
Agm3<-c(sum(A3[1:60]),sum(A3[61:80]),sum(A3[81:111]))

# now for females
A2<-AgeDecomp(LT20f,LT21f)
A3<-AgeDecomp(LT21f,LT22f)

# age-groups = 0-59, 60-79, 80+
Agf2<-c(sum(A2[1:60]),sum(A2[61:80]),sum(A2[81:111]))
Agf3<-c(sum(A3[1:60]),sum(A3[61:80]),sum(A3[81:111]))



Table3b<-round(   
  rbind(cbind(c(Agf,sum(Agf)),c(Agf2,sum(Agf2)),c(Agf3,sum(Agf3)),c(Agf4,sum(Agf4))),
        c(rep(0,4)), 
        cbind(c(Agm,sum(Agm)),c(Agm2,sum(Agm2)),c(Agm3,sum(Agm3)),c(Agm4,sum(Agm4))))
  ,2)



dimnames(Table3b)[[1]]<-c("fem 0-59","fem 60-79","fem 80+","Total",0, "male 0-59","male 60-79","male 80+","Total")


dimnames(Table3b)[[2]]<-c("2020vs2017-19", "2021vs2020", "2022vs2021", "2020-22vs2017-19")

# checking that they match nicely the total
#sum(AgeDecomp(LT19m,LT20m))
#sum(AgeDecomp(LT15m,LT19m))/4

###############################################################
############ causes of death
#################################################################
## D2022m D2022f D2021m D2021f D2020m D2020f D2019m D2019f D2018m D2018f D2017m D2017f
#unique(COD2$Age)
# [1] "0-4"   "5-9"   "10-14" "15-19" "20-24" "25-29" "30-34" "35-39" "40-44" "45-49" "50-54"
# [12] "55-59" "60-64" "65-69" "70-74" "75-79" "80-84" "85-89" "90-94" "95+" 

COD2[is.na(COD2)]<-0 

Causes<-function(D20,S,Y){
  # D20<-D2022f
  # S<-"Females"
  # Y<-2022
  
D<-c(sum(D20[1:5]),sum(D20[6:10]),sum(D20[11:15]),sum(D20[16:20]),sum(D20[21:25]),sum(D20[26:30]),
       sum(D20[31:35]),sum(D20[36:40]),sum(D20[41:45]),sum(D20[46:50]),sum(D20[51:55]),sum(D20[56:60]),
       sum(D20[61:65]),sum(D20[66:70]),sum(D20[71:75]),sum(D20[76:80]),sum(D20[81:85]),sum(D20[86:90]),
       sum(D20[91:95]),sum(D20[96:101]))
COD<-COD2[(COD2$Year==Y)&(COD2$Sex==S),][,-(1:3)]*D
#COD<-cbind(COD[,5],rowSums(COD[,6:7]),rowSums(COD[,2:4]),COD[,10],COD[,1],D-rowSums(COD[,c(1:7,10)]))
Others<-D-rowSums(COD[,c(1:10)])
COD<-cbind(COD[,1:10],Others)
CODn<-rbind(colSums(COD[1:2,]),colSums(COD[3:4,]),colSums(COD[5:6,]),colSums(COD[7:8,]),
            colSums(COD[9:10,]),colSums(COD[11:12,]),colSums(COD[13:14,]),colSums(COD[15:16,]),
            colSums(COD[17:18,]),colSums(COD[19:20,]))
return(CODn)
}




CoD17m<-Causes(dD2017m,"Males",2017)
CoD18m<-Causes(dD2018m,"Males",2018)
CoD19m<-Causes(dD2019m,"Males",2019)
CoD20m<-Causes(D2020m,"Males",2020)
CoD21m<-Causes(D2021m,"Males",2021)
CoD22m<-Causes(D2022m,"Males",2022)
CoD17f<-Causes(dD2017f,"Females",2017)
CoD18f<-Causes(dD2018f,"Females",2018)
CoD19f<-Causes(dD2019f,"Females",2019)
CoD20f<-Causes(D2020f,"Females",2020)
CoD21f<-Causes(D2021f,"Females",2021)
CoD22f<-Causes(D2022f,"Females",2022)

cbind(colSums(CoD19f),colSums(CoD20f),colSums(CoD21f),colSums(CoD22f))
cbind(colSums(CoD19m),colSums(CoD20m),colSums(CoD21m),colSums(CoD22m))


CoDABS<-function(Mm){
mm<-Mm
n<-dim(mm)[2]
mm3<-mm/matrix(rep(rowSums(mm),n),10)
return(mm3)}


CoDABS3<-function(MM1,MM2,MM3){
    # MM1<-CoD20m
    # MM2<-CoD21m
    # MM3<-CoD22m 
  mm<-MM1+MM2+MM3
  n<-dim(mm)[2]
  mm3<-mm/matrix(rep(rowSums(mm),n),10)
  return(mm3)}


IICD19m<-CoDABS3(CoD17m,CoD18m,CoD19m)
IICD19f<-CoDABS3(CoD17f,CoD18f,CoD19f)

#2020

LT1<-lifetable.mx(R2019m,"m")
LT2<-lifetable.mx(R2020m,"m")
ICD1<-IICD19m
ICD2<-CoDABS(CoD20m)

T2020m<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



LT1<-lifetable.mx(R2019f,"f")
LT2<-lifetable.mx(R2020f,"f")
ICD1<-IICD19f
ICD2<-CoDABS(CoD20f)

T2020f<-CauseDecomp3(LT1,LT2,ICD1,ICD2)

### 2021


LT1<-lifetable.mx(R2019m,"m")
LT2<-lifetable.mx(R2021m,"m")
ICD1<-IICD19m
ICD2<-CoDABS(CoD21m)

T2021m<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



LT1<-lifetable.mx(R2019f,"f")
LT2<-lifetable.mx(R2021f,"f")
ICD1<-IICD19f
ICD2<-CoDABS(CoD21f)

T2021f<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



### 2022

LT1<-lifetable.mx(R2019m,"m")
LT2<-lifetable.mx(R2022m,"m")
ICD1<-IICD19m
ICD2<-CoDABS(CoD22m)

T2022m<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



LT1<-lifetable.mx(R2019f,"f")
LT2<-lifetable.mx(R2022f,"f")
ICD1<-IICD19f
ICD2<-CoDABS(CoD22f)

T2022f<-CauseDecomp3(LT1,LT2,ICD1,ICD2)


########## Average 2020-22 

IICD2012m<-CoDABS3(CoD20m,CoD21m,CoD22m)
IICD2012f<-CoDABS3(CoD20f,CoD21f,CoD22f)


LT1<-lifetable.mx(R2019m,"m")
LT2<-lifetable.mx(R202012m,"m")
ICD1<-IICD19m
ICD2<-IICD2012m

T202012m<-CauseDecomp4(LT1,LT2,ICD1,ICD2)


LT1<-lifetable.mx(R2019f,"f")
LT2<-lifetable.mx(R202012f,"f")
ICD1<-IICD19f
ICD2<-IICD2012f

T202012f<-CauseDecomp4(LT1,LT2,ICD1,ICD2)



#######################

Table4f<-rbind(rep(0,11),T2020f,rep(0,11),T2021f,rep(0,11),T2022f,rep(0,11),T202012f)
              


dimnames(Table4f)[[1]]<-c("2020vs2017-19","0-59","60-79","80+",
                         "2021vs2017-19","0-59","60-79","80+",
                         "2022vs2017-19","0-59","60-79","80+",
                         "2020-22vs2017-19","0-59","60-79","80+")

Table4m<-rbind(rep(0,11),T2020m,rep(0,11),T2021m,rep(0,11),T2022m,rep(0,11),T202012m)



dimnames(Table4m)[[1]]<-c("2020vs2017-19","0-59","60-79","80+",
                          "2021vs2017-19","0-59","60-79","80+",
                          "2022vs2017-19","0-59","60-79","80+",
                          "2020-22vs2017-19","0-59","60-79","80+")


Table5<-round(rbind(c(colSums(T2020f),sum(T2020f)),c(colSums(T2021f),sum(T2021f)),c(colSums(T2022f),sum(T2022f)),c(colSums(T202012f),sum(T202012f))
              ,rep(0,12),
              c(colSums(T2020m),sum(T2020m)),c(colSums(T2021m),sum(T2021m)),c(colSums(T2022m),sum(T2022m)),c(colSums(T202012m),sum(T202012m))),3)


dimnames(Table5)[[1]]<-c("F-2020vs2017-19","F-2021vs2017-19","F-2022vs2017-19","F-2020-22vs2017-19"," ",
                         "M-2020vs2017-19","M-2021vs2017-19","M-2022vs2017-19","M-2020-22vs2017-19")


###################



### 2021


LT1<-lifetable.mx(R2020m,"m")
LT2<-lifetable.mx(R2021m,"m")
ICD1<-CoDABS(CoD20m)
ICD2<-CoDABS(CoD21m)

T2021mb<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



LT1<-lifetable.mx(R2020f,"f")
LT2<-lifetable.mx(R2021f,"f")
ICD1<-CoDABS(CoD20f)
ICD2<-CoDABS(CoD21f)

T2021fb<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



### 2022

LT1<-lifetable.mx(R2021m,"m")
LT2<-lifetable.mx(R2022m,"m")
ICD1<-CoDABS(CoD21m)
ICD2<-CoDABS(CoD22m)

T2022mb<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



LT1<-lifetable.mx(R2021f,"f")
LT2<-lifetable.mx(R2022f,"f")
ICD1<-CoDABS(CoD21f)
ICD2<-CoDABS(CoD22f)

T2022fb<-CauseDecomp3(LT1,LT2,ICD1,ICD2)



Table4bf<-rbind(rep(0,11),T2020f,rep(0,11),
              T2021fb,rep(0,11),T2022fb,rep(0,11),T202012f)


dimnames(Table4bf)[[1]]<-c("2020vs2017-19","0-59","60-79","80+",
                         "2021vs2020","0-59","60-79","80+",
                         "2022vs2021","0-59","60-79","80+",
                         "2020-22vs2017-19","0-59","60-79","80+")


Table4bm<-rbind(rep(0,11),T2020m,rep(0,11),
                T2021mb,rep(0,11),T2022mb,rep(0,11),T202012m)


dimnames(Table4bm)[[1]]<-c("2020vs2017-19","0-59","60-79","80+",
                           "2021vs2020","0-59","60-79","80+",
                           "2022vs2021","0-59","60-79","80+",
                           "2020-22vs2017-19","0-59","60-79","80+")

Table5b<-round(rbind(c(colSums(T2020f),sum(T2020f)),
                     c(colSums(T2021fb),sum(T2021fb)),
                     c(colSums(T2022fb),sum(T2022fb)),
                     c(colSums(T202012f),sum(T202012f))
                    ,rep(0,12),
                    c(colSums(T2020m),sum(T2020m)),c(colSums(T2021mb),sum(T2021mb)),c(colSums(T2022mb),sum(T2022mb)),c(colSums(T202012m),sum(T202012m))),3)

Covid<-Table5b[,1]
Respiratory<-rowSums(Table5b[,2:4])
Cancers<-Table5b[,5]
IHD<-rowSums(Table5b[,6:7])
Others<-rowSums(Table5b[,8:11])
Table5c<-cbind(Covid,Respiratory,Cancers,IHD,Others)


dimnames(Table5b)[[1]]<-c("F-2020vs2017-19","F-2021vs2020","F-2022vs2021","F-2020-22vs2017-19"," ",
                         "M-2020vs2017-19","M-2021vs2020","M-2022vs2021","M-2020-22vs2017-19")


dimnames(Table5c)[[1]]<-c("F-2020vs2017-19","F-2021vs2020","F-2022vs2021","F-2020-22vs2017-19"," ",
                          "M-2020vs2017-19","M-2021vs2020","M-2022vs2021","M-2020-22vs2017-19")


#############################

####### Figure 1 comes from data in Table1
####### Figure 2 is yours
####### Figure 3 comes from data in Table3b
####### Figure 4 comes from data in Table5c


