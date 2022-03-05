U1.Clayton<-function(x.obs,y.obs,dx,dy,lower=0.001,upper=50,U.plot=TRUE){
    n=length(x.obs)
    
    U.func=function(theta){
      u=0
      for(i in 1:(n-1)){
        xi=x.obs[i];yi=y.obs[i];dxi=dx[i];dyi=dy[i]
        for(j in (i+1):n){
          xj=x.obs[j];yj=y.obs[j];dxj=dx[j];dyj=dy[j]
          Ox=((xi<=xj)&(dxi))|((xj<=xi)&(dxj))
          Oy=((yi<=yj)&(dyi))|((yj<=yi)&(dyj))
          if(Ox&Oy){
            xij=min(xi,xj);yij=min(yi,yj)
            Rij=sum( (x.obs>=xij)&(y.obs>=yij) )
            Wij=(theta+2)/(theta+1)/(Rij+theta)
            u=u+Wij*(  ((xi-xj)*(yi-yj)>0)-(theta+1)/(theta+2)  )
          }
        }
      }
      u
    }
    
    if((U.func(lower)<0)&(U.func(upper)<0)){ theta=lower }else{  
      theta=uniroot(U.func,lower=lower,upper=upper)$root
    }
    tau=theta/(theta+2)
    
    if(U.plot==TRUE){  
      curve(U.func,lower,upper,xlab="theta",ylab="U1(theta)")
      points(theta,U.func(theta),col="red",cex=1)
      abline(h=0,lty="dotted",col="blue")
    }
    
    return(c(theta=theta,tau=tau))
}
  
  