oneParameterContamination <- function (n,s,p,mean1,mean2){
  ## on creer un echantillon de donnée dont les vecteur suivent N(0,I_p) 
  X<- matrix(data= rnorm(n*p), nrow=p,ncol=n)
  ## on somme aves means1 pour obtenir un echnatillon de N(means1,I_p)
  
  X <- X + mean1%*% rep(1,n)
  
  ## on choisi s observation:
  
  S<- sample (1:n,s)
  ## on les contamine :
  
  for (i in S){
    X[,i] <- cbind(X[,i]) +mean2-mean1
  }
  
  return (X)
  
}

uniParameterContamination <- function (n,s,p,mean1,matMeans){
  
  k <- ncol(matMeans) # le nombre de parametres différents
  
  ## on creer un echantillon de donnée dont les vecteur suivent N(0,I_p) 
  X<- matrix(data= rnorm(n*p), nrow=p,ncol=n)
  ## on somme aves means1 pour obtenir un echnatillon de N(means1,I_p)
  
  X <- X + mean1%*% rep(1,n)
  
  ## on choisi s observation:
  
  S<- sample (1:n,s)
  ## on les contamine :
  
  for (i in S){
    a <- as.integer(runif(1,1,k))[1]
    X[,i] <- cbind(X[,i]) +cbind(matMeans[,a])-mean1
  }
  
  return (X)
  
}

oneParameterContamination2 <- function (n,s,p,mean1,mean2,sd){
  ## on creer un echantillon de donnée dont les vecteur suivent N(0,I_p) 
  X<- matrix(data= rnorm(n*p), nrow=p,ncol=n)
  ## on somme aves means1 pour obtenir un echnatillon de N(means1,I_p)
  
  X <- X + mean1%*% rep(1,n)
  
  ## on choisi s observation:
  
  S<- sample (1:n,s)
  ## on les contamine :
  
  for (i in S){
    X[,i] <- (cbind(X[,i]) -mean1)*sd +mean2
  }
  
  return (X)
  
}


cont1 <- function(n,p,eps){
  return(list(X= oneParameterContamination(n,s=as.integer(n*eps),p=p,mean1=cbind(rep(1,p)),mean2=cbind(rep(3,p))),
              mu = cbind(rep(1,p))))
}
