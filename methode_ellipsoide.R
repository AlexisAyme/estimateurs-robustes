estDansOmega <- function (omega, n, s ) {
  ### fonction qui indique si un vecteur poid omega est dans l'ensemble Omega 
  ### dans le cas contraire la fonction retourne la direction de l'ensemble oméga 
  ### n le nombre d'observation, s le nombre de contamination 
  
  
  ## on verifie d'abord si les n-1 composantes sont dans [0,1/(n-s)] :
  for(i in 1:(n-1)){
    if (omega[i,1]<0 ){
      a <- cbind(rep(0,n-1))
      a[i,1]<- -1  # la contrainte n'est pas respecté, on retourne la direction de Omega 
      return (list(hyperVect=a,bool = FALSE) )
    }
    else{ 
      if (omega[i,1]>1/(n-2*s) ){
        a <- cbind(rep(0,n-1))
        a[i,1]<- 1  # meme chose mais de l'autre signe  
        return (list(hyperVect=a,bool = FALSE))
      }
    }
    
  } 
  
  ## a ce stade les n-1 premiere composante sont dans le bon intervalle, il reste a verifier 
  ## que la derniere forme affine des autres est bien dans [0,1/(n-s)]
  
  if (omega[n,1] < 0) { # somme supérieur a 1
    return  (list(hyperVect=cbind(rep(1,n-1)),bool = FALSE))
  }else {
    if (omega[n,1]> 1/(n-2*s)){
      return  (list(hyperVect=-cbind(rep(1,n-1)),bool = FALSE))
    }
    else {## toutes les contraintes sont verifiées, on est bien dans Oméga 
      return (list(bool=TRUE))
    }
  }
}

#omega <- cbind(c(1/8,1/4,1/4,1/8,1/8,1/8))
#print(estDansOmega(omega,6,1))
#omega <- cbind(c(1/8,1/8,1/8,1/8,1/8,3/8))
#print(estDansOmega(omega,6,1))



separationOracle <- function (omega, X, tau,s) {
  ### Fonction qui indique si le point omega est dans l'ensemble Omega* et
  ### dans le cas inverse donne la forme affine qui situe cette espace dans
  ### le demi espace négatif pour cette forme affine 
  
  p <- nrow(X) # la dimension des observations
  n <- ncol(X) # la dimension des observations
  
  ### Etape 1 : on verifie si on est bien dans Oméga 
  
  dansOmega <- estDansOmega(omega,n,s)
  if (dansOmega$bool == FALSE ){
    return (list(vect =dansOmega$hyperVect, bool = FALSE))
  }
  
  ### Etape 2: on situe Omega* 
  
  hat_mu <-  X %*% omega  # La moyenne estimée avec les poids omega
  Y <- X - hat_mu%*% rep(1,n) # Les observations centrées avec hat_mu
  
  ## definition de la matrice de covariance empirique estimé avec omega :
  
  mat <- matrix(data= rep(0,p*p),nrow=p,ncol=p)
  diag(mat)<- rep(-1,p)
  for (i in 1:n) {
    mat <- mat + omega[i,1] * (Y[,i]%o%Y[,i])
  }
  ## on cherche la valeur propre de module le plus grand (norme spectrale)
  ## et son vecteur propre :
  vp <- eigen(mat)
  i <-  which.max(abs(vp$values))
  lambda <- vp$values[i] # la valeur propre de plus grand module
  v <- cbind(vp$vectors[,i])

  if (abs(lambda) < tau){ # dans ce cas omega appartient à Omega*
    return (list(bool=TRUE))
  }
  else if (lambda> tau) {
    # on calcule les coordonnées de la forme lineaire dans sa base duale (en dimension n-1) qui
    # defini l'hyperplan (coordonnées sous forme de vecteur ligne) :
    
    a <- ((t(Y[,1:(n-1)]) %*% v )^2 - cbind(rep(1,n-1))) #on calcule sur le n-1 premiere composante 
    
    ## la derniere composante est vu comme une forme affine des autres, donc elle intervient dans les autres 
    ## composantes : 
    
    
    a <- a - (Y[,n]%*% v )[1,1]^2 * cbind(rep(1,n-1)) 
    
    return (list(vect=a,bool=FALSE))
    
  }
  else {
    # la meme chose mais avec un moins 
    
    a <- ((t(Y[,1:(n-1)]) %*% v )^2 - cbind(rep(1,n-1)))
    a <- -a + (Y[,n]%*% v)[1,1]^2 * cbind(rep(1,n-1))
    
    return (list(vect=a,bool=FALSE))
  }
}

##omega = cbind(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
##X=matrix(data= rep(1,30),nrow=3,ncol=10)
##separationOracle(omega,X,0.2,s=1)




newEllipsoide <- function (centre,matEll,vectHyperplan){
  ###Prend un ellipsoide de matrice matEll et de centre centre, un hyperplan
  ###affine de vecteur vectHyperplan et d'origine centre.
  ###retourne l'ellipsoide sous forme de couple (centre, matrice) qui contient
  ###la demi ellipsoide et qui soit de volume minimum
  
  n <- nrow(centre)
  b <- (matEll %*% vectHyperplan)
  b <- b / sqrt((t(vectHyperplan) %*% matEll %*% vectHyperplan)[1,1])
  newCentre <- centre - (1/(n+1))* b
  newMatEll <- (n^2/(n^2-1))*(matEll- (2/(n+1))*b%*%t(b))

  
  return (list(newCentre=newCentre, newMatEll=newMatEll,qte= (t(vectHyperplan) %*% matEll %*% vectHyperplan)[1,1] ))
}

ellipsoideMethode <- function (X, tau,s){
  ### utilise la methode des ellipsoides pour envoyer une estimation de la moyenne 
  ### on delimite l'ensemble cible convexe Omega* avec tau :
  ### l'ensemble cible est l'ensemble des poids omega qui permet d'avoir 
  ### une variance empirique distant de Ip d'au plus sqrt(tau) en norme spectrale
  ### pour que Oméga est un volume on réduit d'une dimmension, la derniere composante 
  ###d'un vecteur omega sera toujours 1-somme des autres composantes 
  
  p <- nrow(X) # la dimension des observations
  n <- ncol(X) # la dimension des observations
  
  ## IINITIALISATION : On prend comme ellipsoide initial la boule de centre 0 et de rayon sqrt(n) :
  
  matEll <- matrix(data= rep(0,(n-1)*(n-1)), ncol=n-1,nrow=n-1) 
  diag(matEll) <- sqrt(1/(n-s)) *rep(1,n-1)
  omega <- cbind(rep(0,n))
  omega [n,1]<- 1 # obligatoire car on prend les autres composantes nulles 
  
  ## Iteration de la méthode des ellipsoides : 
  it <-0
  
  
  T3 <- system.time(f<-1)
  
  repeat {
    sep <- separationOracle(omega,X,tau,s)
    if (sep$bool){
      break
    }
    else {
      ## on calcule l'ellipsoide pour le n-1 premieres composantes
      newEll <- newEllipsoide(cbind(omega[1:(n-1),1]),matEll,sep$vect)
      omega[1:(n-1), 1]<- newEll$newCentre
      omega[n,1] <- 1- sum(omega[1:(n-1),1])
      matEll <- newEll$newMatEll
      
      it<- it+1
      #if (it%%100 == 0){
      #  print(list(iteration=it,Omega=(estDansOmega(omega,n,s)$bool),volume=determinant.matrix(matEll)$modulus))
      #}
  
      
      
      
    }
  }
  
  ## on retrourne la moyenne empirique estimé avec le vecteur poids omega 
  return (X %*% omega)
  
}

diakonikolasInequality <- function (X,s,delta){
  ## return the variable multiplier 
  ## norm(muEll - mu)^2 < C* DiakonikolasInequality(X,s,delta)
  p <- nrow(X)
  n <- ncol(X)
  
  eps <- s/n 
  return (p/n + (eps^2)*log(1/eps)+ log(1/delta)/n )
}

ellipsoideMethode2 <- function(X,s,delta){
  return (ellipsoideMethode(X,10*diakonikolasInequality(X,s,delta),s))
}

