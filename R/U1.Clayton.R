U1.Clayton<-function(x.obs,y.obs,dx,dy,lower=0.001,upper=50){
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
    
    theta=uniroot(U.func,lower=lower,upper=upper)$root
    tau=theta/(theta+2)
    return(c(theta=theta,tau=tau))
}
  
  