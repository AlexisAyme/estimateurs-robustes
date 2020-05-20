projector <- function(X,V){
  ### project the database X on subspace V which is defined by columns of V 
  ### we implemente a sub space V such as a matrix V (columns are vectors of bases of V, 
  ### so dim(V) is the number of columns of matrix V)
  return (t(V)%*% X)
}


outlierRemoval <- function (X,C){
  ### remove observations (columns of matrix X) too far from componentwise median 
  
  n= ncol(X) # number of observations 
  p=nrow(X) # dimension of observations
  
  ## componentwise median,
  
  compMed <-  cbind(apply(X,1,median))
  
  ## remove observation too far,
  
  Y <- X - compMed %*% rep(1,n)
  # vector of distance^2 to componentwise median 
  d2 <- rep(1,n)
  for (i in 1:n){
    d2 [i] <- norm (cbind(Y[,i]),type="F")^2
    
  }
  
  
  return (X[,d2< (C^2)*p*log(n)])
  
}

spectralMethod <- function(X,C){
  ## use spectral methode to estimate mean of adversaly noisy gausian 
  ## we implemente a sub space V such as a matrix V (columns are vectors of bases of V, 
  ## so dim(V) is the number of columns of matrix V)
  
  n <- ncol(X) # number of observations 
  p <- nrow(X) # dimension of observations
  
  ### data-splitting,
  
  k <- as.integer(log(p,base=2))
  j= as.integer(n/k) ##we will split data in sub set of j observations 
  
  ### initialise subspace (as R^p) W and mean 
  
  dimW <- p # dimension of W 
  W<-  matrix ( data= rep(0,p*p), nrow=p, ncol=p)
  diag(W) <- rep(1,p)
  
  muSp = cbind(rep(0,p))
  
  ### main loop  
  
  for (i in 1:k){

    ## In sub space W
    Dp <- projector(X[,((i-1)*j+1):(i*j)], W) # projection on W of D_i
    Dp <- outlierRemoval(Dp, C)
    Y <- Dp - cbind(rowMeans(Dp))%*%rep(1,ncol(Dp)) # centered data 
    
    matCov <-  Y%*%t(Y)
    vp <- eigen (matCov)
  
    dimWp <- as.integer(dimW/2)
    Wp <- cbind(vp$vectors[,1:dimWp]) # we select dimWp top singular vectors 
    Vp <- cbind(vp$vectors[,(dimWp+1):dimW])
    
    meansVp <- rowMeans(t(Vp)%*%Dp ) # means of observations projeted on Vp, 
                                    ## this is the projection of muSp on Vp
    
    ## In R^p 
    
    muSp = muSp + W%*%Vp%*%meansVp # Vp%*%meansVp in W and W%*%Vp%*%meansVp in R^p
                                  ## part in R^p of the means of meansVp
    
    ## define new sub space W and its dimension 
    
    W <- W%*%Wp 
    dimW <- ncol(W)
  }
  
  ### compute the  median on D_k 
  Dp <- projector(X[,((k-1)*j+1):(k*j)],W) 
  medV_k <- median(Dp) #in W 
  muSp <- muSp + W%*% medV_k
  
  return (muSp)
}

laiInequality <- function (X,s){
  ## deteminate the variable factor in Lai inequality,
  ## we have norm(muSp - mu)< C*laiInequality(X,s) 
  
  n <- ncol(X) # number of observations 
  p <- nrow(X) # dimension of observations
  
  return (p*log(p)*log(p)*log(n)/n + ((s/n)^2)*log(p))
}