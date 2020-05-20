
source("contamination.R")
source("methode_ellipsoide.R")
source("spectralMethod.R")


graphForObservation<- function(eps,p,contamination,n_iter){
  n = 90*rep(1,20) + 100*seq(1,20)
  diff = rep(0,20)
  diff2 = rep(0,20)
  time = rep(0,20)
  time2 = rep(0,20)
  diffEmp = rep(0,20)
  for (i in 1:20){
    
    print(i)
    diffp = rep(0,n_iter)
    diff2p = rep(0,n_iter)
    timep = rep(0,n_iter)
    time2p = rep(0,n_iter)
    diffEmpp = rep(0,n_iter)
    for (j in 1:n_iter){
      C <- contamination(n[i],p,eps)
      X <- C$X
      mu <- C$mu 
      timep [j] <- system.time(diffp[j] <- norm(spectralMethod(X,0.5)-mu,type="F"))[3]
      time2p [j] <- system.time(diff2p[j] <- norm(ellipsoideMethode2(X,as.integer(n[i]*eps),0.05)-mu,type="F"))
      diffEmpp[j] = norm(cbind(rowMeans(X))-mu,type="F")
      
    }
    diff[i] <- mean(diffp)
    diff2[i] <- mean(diff2p)
    diffEmp[i] <- mean(diffEmpp)
    time[i]<- mean(timep)
    time2[i]<- mean(time2p)
    
    title = paste(paste(paste("Efficacité des différents estimateurs pour p=",deparse(p)),"et epsilon ="),deparse(eps))
    
    ## graph of difference 
    
    plot(n[1:i],diffEmp[1:i],
         type="l",
         ylim=c(0,max(c(diff,diff2,diffEmp))),
         ylab = "distance ",
         xlab= "n",
         col=3,
         main =title)
    text(200,diffEmp[1], "Moyenne empirique", col = 3, adj = c(-.1, -.1))
    points(n[1:i],diff[1:i],type="l",col=2)
    text(200,diff[1], "Moyenne spectrale", col = 2, adj = c(-.1, -.1))
    points(n[1:i],diff2[1:i],type="l",col=1)
    text(200,diff2[1], "Moyenne ellipsoide", col = 1, adj = c(-.1, -.1))
    
    ## graph of empiric complexity
    
    title = paste(paste(paste("Compléxité de la méthode spectrale p=",deparse(p)),"et epsilon ="),deparse(eps))
    
    plot(n[1:i],time[1:i],
         type="l",
         ylab = "Temps ",
         xlab= "n",
         col=3,
         main =title)
    
    title = paste(paste(paste("Compléxité de la méthode des ellipsoides p=",deparse(p)),"et epsilon ="),deparse(eps))
    
    plot(n[1:i],sqrt(time2[1:i]),
         type="l",
         ylab = "racine carrée du temps ",
         xlab= "n",
         col=3,
         main =title)
    
  }

  
  
  
}

graphForVariable <- function(eps,n,contamination,n_iter){
  p = c(32,100,200,400,600)
  diff = rep(0,5)
  diff2 = rep(0,5)
  time = rep(0,5)
  time2 = rep(0,5)
  diffEmp = rep(0,5)
  for (i in 1:5){
    
    print(i)
    diffp = rep(0,n_iter)
    diff2p = rep(0,n_iter)
    timep = rep(0,n_iter)
    time2p = rep(0,n_iter)
    diffEmpp = rep(0,n_iter)
    for (j in 1:n_iter){

      C <- contamination(n,p[i],eps)
      X <- C$X
      mu <- C$mu 
      timep [j] <- system.time(diffp[j] <- norm(spectralMethod(X,0.5)-mu,type="F"))[3]
      time2p [j] <- system.time(diff2p[j] <- norm(ellipsoideMethode2(X,as.integer(n*eps),0.05)-mu,type="F"))[3]
      diffEmpp[j] = norm(cbind(rowMeans(X))-mu,type="F")
      
    }
    diff[i] <- mean(diffp)
    diff2[i] <- mean(diff2p)
    diffEmp[i] <- mean(diffEmpp)
    time[i]<- mean(timep)
    time2[i]<- mean(time2p)
    
    title = paste(paste(paste("Efficacité des différents estimateurs pour n=",deparse(n)),"et epsilon ="),deparse(eps))
    
    ## graph of difference 
    
    plot(p[1:i],diffEmp[1:i],
         type="l",
         ylim=c(0,max(c(diff,diff2,diffEmp))),
         ylab = "distance ",
         xlab= "p",
         col=3,
         main =title)
    text(50,diffEmp[1], "Moyenne empirique", col = 3, adj = c(-.1, -.1))
    points(p[1:i],diff[1:i],type="l",col=2)
    text(50,diff[1], "Moyenne spectrale", col = 2, adj = c(-.1, -.1))
    points(p[1:i],diff2[1:i],type="l",col=1)
    text(50,diff2[1], "Moyenne ellipsoide", col = 1, adj = c(-.1, -.1))
    
    ## graph of empiric complexity
    
    title = paste(paste(paste("Compléxité de la méthode spectrale n=",deparse(n)),"et epsilon ="),deparse(eps))
    
    plot(p[1:i],time[1:i],
         type="l",
         ylab = "Temps ",
         xlab= "p",
         col=3,
         main =title)
    
    title = paste(paste(paste("Compléxité de la méthode des ellipsoides n=",deparse(n)),"et epsilon ="),deparse(eps))
    
    plot(p[1:i],sqrt(time2[1:i]),
         type="l",
         ylab = "Racine carrée du temps ",
         xlab= "p",
         col=3,
         main =title)
    
  }
  
}

graphForEpsilon <- function(n,p,contamination,n_iter){
  eps = 0.05 *seq(1,5)
  diff = rep(0,5)
  diff2 = rep(0,5)
  diffEmp = rep(0,5)
  for (i in 1:5){
    
    print(i)
    diffp = rep(0,n_iter)
    diff2p = rep(0,n_iter)
    diffEmpp = rep(0,n_iter)
    for (j in 1:n_iter){
      C <- contamination(n,p,eps[i])
      X <- C$X
      mu <- C$mu 
      diffp[j] <- norm(spectralMethod(X,0.5)-mu,type="F")
      diff2p[j] <- norm(ellipsoideMethode2(X,as.integer(n*eps[i]),0.05)-mu,type="F")
      diffEmpp[j] = norm(cbind(rowMeans(X))-mu,type="F")
      
    }
    diff[i] <- mean(diffp)
    diff2[i] <- mean(diff2p)
    diffEmp[i] <- mean(diffEmpp)
    
    title = paste(paste(paste("Efficacité des différents estimateurs pour n=",deparse(n)),"et p ="),deparse(p))
    
    ## graph of difference 
    
    plot(eps[1:i],diffEmp[1:i],
         type="l",
         ylim=c(0,max(c(diff,diff2,diffEmp))),
         ylab = "distance ",
         xlab= "eps",
         col=3,
         main =title)
    text(0.05,diffEmp[1], "Moyenne empirique", col = 3, adj = c(-.1, -.1))
    points(eps[1:i],diff[1:i],type="l",col=2)
    text(0.05,diff[1], "Moyenne spectrale", col = 2, adj = c(-.1, -.1))
    points(eps[1:i],diff2[1:i],type="l",col=1)
    text(0.05,diff2[1], "Moyenne ellipsoide", col = 1, adj = c(-.1, -.1))
    
    
  }
  
  
  
}

#graphForObservation(0.1,32,cont1,10)
graphForVariable(0.1,600,cont1,1)
#graphForEpsilon(500,32,cont1,4)

graphForSpectralMethod <- function(eps,p,contamination,n_iter){
  n = 1000*seq(1,20)
  diff = rep(0,20)
  diffEmp = rep(0,20)
  np = c()
  diffp=c()
  sigma= rep(0,20)
  for (i in 1:20){
    print(i)
    diffpp = rep(0,n_iter)
    diffEmpp = rep(0,n_iter)
    for (j in 1:n_iter){
      C <- contamination(n[i],p,eps)
      X <- C$X
      mu <- C$mu 
      diffpp[j] <- norm(spectralMethod(X,0.5)-mu,type="F")
      
      diffEmpp[j] = norm(cbind(rowMeans(X))-mu,type="F")
      
    }
    diff[i] <- mean(diffpp)
    diffEmp[i] <- mean(diffEmpp)
    sigma [i] <-sqrt(var(diffpp))
    diffp = c(diffp,diffpp)
    np=c(np,rep(n[i],n_iter))
    
  }
  
  title = paste(paste(paste("Efficacité de la méthode spectrale pour p=",deparse(p)),"et epsilon ="),deparse(eps))
  
  plot(np,diffp,ylab = "distance ",
       xlab= "n",
       ylim=c(0,max(diffEmp)),
       col=1,
       main =title)
  #points(n,diff,type="l",col=2)
  #text(200,diff[1], "Moyenne spectrale", col = 2, adj = c(-.1, -.1))
  #points(n,diff+sigma,type="l",col=4)
  #points(n,diff-sigma,type="l",col=4)
  #
  points(n,diffEmp,type="l",col=3)
  text(200,diffEmp[1], "Moyenne empirique", col = 3, adj = c(-.1, -.1))
  
}

#graphForSpectralMethod(0.1,32,cont1,50)

