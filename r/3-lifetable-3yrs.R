
## all from HMD

######## functions
library(MortCast)

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
  names(mx) <- 0:110
  mx1<-kannisto(mx, est.ages = seq(80, 110), proj.ages = seq(85, 110))
  return(mx1)}

Kannisto2<-function(mx){
  names(mx) <- 0:110
  mx[is.nan(mx)]<-0
  mx1<-kannisto(mx, est.ages = seq(80, 100), proj.ages = seq(85, 110))
  return(mx1)}
##############



#Country<-c("USA","GBRCENW","FRATNP","NOR") 
Country <- c(
  "BGR"
  , "CZE"
  , "DNK"
  , "FIN"
  , "JPN"
  , "LUX"
  , "NZL_NP"
  , "NOR"
  , "PRT"
  , "SWE"
  , "CHE"
  , "USA"
)



E0f20<-c()
E0m20<-c()


E0f21<-c()
E0m21<-c()

### only denmark,norway,and sweden have data for 2022 
E0f22<-rep(0,length(Country))
E0m22<-rep(0,length(Country))





for(i in 1:length(Country)){
data.path <- "../data/E0per/"

  Ex<-read.table(paste(data.path, Country[i], ".E0per.txt", sep = ""),header=TRUE,fill=TRUE,skip=1)

  E0f20<-c(E0f20,Ex$Female[Ex$Year==2020])
  E0m20<-c(E0m20,Ex$Male[Ex$Year==2020])
  
  E0f21<-c(E0f21,Ex$Female[Ex$Year==2021])
  E0m21<-c(E0m21,Ex$Male[Ex$Year==2021])
 
   if((i==3)|(i==8)|(i==10)){
  E0f22[i]<-Ex$Female[Ex$Year==2022]
  E0m22[i]<-Ex$Male[Ex$Year==2022]
  }
}


E0f17to19<-c()
E0m17to19<-c()  
E0f17to19b<-c()
E0m17to19b<-c()  

for(i in 1:length(Country)){
data.path <- "../data/Deaths_1x1/"

  Dx<-read.table(paste(data.path, Country[i], ".Deaths_1x1.txt", sep = ""),header=TRUE,fill=TRUE,skip=1)

data.path <- "../data/Exposures_1x1/"

  Px<-read.table(paste(data.path, Country[i], ".Exposures_1x1.txt", sep = ""),header=TRUE,fill=TRUE,skip=1)
  
  
  Dx1<-Dx[Dx$Year%in%c(2017,2018,2019),]
  Px1<-Px[Px$Year%in%c(2017,2018,2019),]
  
  if((i!=3)&(i!=4)&(i!=6)&(i!=9)&(i!=11)){
  LTF<- lifetable.mx(Kannisto(rowSums(matrix(Dx1$Female,111))/rowSums(matrix(Px1$Female,111))),"f")
  LTM<- lifetable.mx(Kannisto(rowSums(matrix(Dx1$Male,111))/rowSums(matrix(Px1$Male,111))),"m")
  }
  if((i==3)|(i==4)|(i==6)|(i==9)|(i==11)){
  LTF<- lifetable.mx(Kannisto2(rowSums(matrix(Dx1$Female,111))/rowSums(matrix(Px1$Female,111))),"f")
  LTM<- lifetable.mx(Kannisto2(rowSums(matrix(Dx1$Male,111))/rowSums(matrix(Px1$Male,111))),"m")
  }
  E0f17to19<- c(E0f17to19,LTF[1,10])
  E0m17to19<- c(E0m17to19,LTM[1,10])
 
  ### now without Kannisto
  mx.f<-rowSums(matrix(Dx1$Female,111))/rowSums(matrix(Px1$Female,111))
  mx.f[is.nan(mx.f)]<-0
  
  mx.m<-rowSums(matrix(Dx1$Male,111))/rowSums(matrix(Px1$Male,111))
  mx.m[is.nan(mx.m)]<-0
  
  LTF2<- lifetable.mx( mx.f,"f")
  LTM2<- lifetable.mx(mx.m,"m")
  
  E0f17to19b<- c(E0f17to19b,LTF2[1,10])
  E0m17to19b<- c(E0m17to19b,LTM2[1,10])
}



E0f<-rbind(E0f22,E0f21,E0f20,E0f17to19)
E0m<-rbind(E0m22,E0m21,E0m20,E0m17to19)

dimnames(E0f)[[2]]<-Country
dimnames(E0m)[[2]]<-Country


### ...
diff.f<-E0f[-4,]-E0f[-1,]
diff.m<-E0m[-4,]-E0m[-1,]
diff.f[1,c(-3,-8,-10)]<-0
diff.m[1,c(-3,-8,-10)]<-0


dimnames(diff.f)[[1]]<-c("22-21","21-20","20-17to19")
dimnames(diff.m)[[1]]<-c("22-21","21-20","20-17to19")
