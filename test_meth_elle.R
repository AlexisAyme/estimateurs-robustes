source("methode_ellipsoide.R")
source("spectralMethod.R")
source("contamination.R")

#X= oneParameterContamination(n=1000,s=100,p=3,mean1=cbind(c(1,1,1)),mean2=cbind(c(3,3,3)))

#X<- oneParameterContamination(n=10000,s=1000,p=32,mean1=cbind(rep(1,32)),mean2=cbind(rep(3,32)))
#print (diakonikolasInequality(X,100,0.05))
#
#print(system.time(muSp <- spectralMethod(X,C=0.7)))
#
#muEll =ellipsoideMethode(X,tau=0.2,s=100)
#
#
#
#muEmp= cbind(rowMeans(X))
#
#
#print(list( diffEll= norm(muEll-cbind(rep(1,3)),type="F"),
#           diffSp= norm(muSp-cbind(rep(1,3)),type="F"),
#           diffEmp= norm(muEmp-cbind(rep(1,3)),type="F")))
#



traceResult <- function (X,s,mu,a,b,pt){
  C= a *rep(1,pt) + ((b-a)/pt)*c(1:pt)
  res= rep(0,pt)
  muEmp= cbind(rowMeans(X))
  for (i in 1:pt){
    diffSp <-  norm(spectralMethod(X,C[i])-mu,type="F")
    res[i] <- diffSp
  }
  
  plot(C,res,type="l",xlim=c(a,b),ylim=c(min(res),max(res)), xlab = "C", ylab = "distance to mu",main = "teststst")
  abline(b = laiInequality(X,s), a = 0, col = 2)
  text(0.5,0.5*laiInequality(X,s), "Lai inequality", col = 2, adj = c(-.1, -.1))
  
  abline(b = 0, a = norm(muEmp-mu,type="F"), col = 3)
  text(1, norm(muEmp-mu,type="F"), "distance with mean emp", col = 3, adj = c(-.1, -.1))
  print(list(distanceWithMuEmp= norm(muEmp-mu,type="F")))
  
}

findConstante <- function(X,mu,a,b,s,eps){
  
  repeat{
    c <- (a+b)/2 
    if (b-a < eps){
      return (c)
    }
    
    if (norm(spectralMethod(X,c)-mu,type="F")> c*laiInequality(X,s)){
      a <- c 
    }else {
      b <- c
    }
  }
  
}
#
#X<- oneParameterContamination(n=10000,s=1000,p=32,mean1=cbind(rep(1,32)),mean2=cbind(rep(3,32)))
#traceResult(X,1000,cbind(rep(1,32)),0.3,2,100)
#
#X<- oneParameterContamination(n=9000,s=500,p=32,mean1=cbind(rep(1,32)),mean2=cbind(rep(3,32)))
#traceResult(X,500,cbind(rep(1,32)),0.3,2,100)
#
matMeans <- matrix (data= rep(1,32*3), ncol=3,nrow=32)
mean1<- cbind(rep(1,32))
mean1[1:10,1]<- rep(2,10)
mean2 <- cbind(rep(1,32)+ rnorm(32,0,5))


mean3 <- cbind(rep(3,32))

matMeans[,1]<- mean1
matMeans[,2]<- mean2
matMeans[,3]<- mean3

X <- uniParameterContamination(n=10000,s=500,p=32,mean1=cbind(rep(1,32)),matMeans = matMeans)
print(findConstante (X,mu=cbind(rep(1,32)),a=0.3,b=2,s=500,eps=0.05))
#traceResult(X,500,cbind(rep(1,32)),0.3,2,100)

X<- oneParameterContamination2(n=10000,s=1000,p=32,mean1=cbind(rep(1,32)),mean2=cbind(rep(1.5,32)),sd=0.3)
traceResult(X,1000,cbind(rep(1,32)),0.3,2,100)
#X<- oneParameterContamination(n=10000,s=1000,p=32,mean1=cbind(rep(1,32)),mean2=cbind(rep(1.5,32)))
#traceResult(X,1000,cbind(rep(1,32)),0.3,2,100)
