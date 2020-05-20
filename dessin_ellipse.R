source("methode_ellipsoide.R")

dessinerEllipse <- function (centre, matrice,n){
  ## donne les x,y pour représenté l'ellipse representé par sa matrice def positive
  
  theta <- 2*pi /n * c(1:n)
  x <- rep(0,n)  # on simule les point de la sphere unité 
  y <- rep (0,n)
  for (i in 1:n){
    v = chol(matrice)%*% cbind(c(cos(theta[i]),sin(theta[i]))) + centre # la decomposition de choleski nous donne 
    ## une representation possible de l'ellipse sous forme d'apllication affine inversible 
    x[i] <- v[1,1]
    y[i]<- v[2,1]
  }
  return (list(x=x,y=y))
}

f<- function(x){return (x)}

mat = matrix( data = c(1,0,0,1), nrow=2, ncol=2)
centre =cbind(c(0,0))
l=dessinerEllipse(centre, mat,1000)
curve(f)
plot(l$x,l$y,type="l",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
pranel.first=grid()

new_ell =newEllipsoide(centre,mat,cbind(c(2,0)))

lp=dessinerEllipse(new_ell$newCentre, new_ell$newMatEll,1000)
points(lp$x,lp$y,type="l",col="red")



new_ell2 =newEllipsoide(new_ell$newCentre, new_ell$newMatEll,cbind(c(1,0)))

lp2=dessinerEllipse(new_ell2$newCentre, new_ell2$newMatEll,1000)
points(lp2$x,lp2$y,type="l",col="blue")

