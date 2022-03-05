U2.Gumbel<-function(x.obs,y.obs,dx,dy,lower=0.01,upper=50,U.plot=TRUE){
    
  n=length(x.obs)
  x.ox=x.obs[order(x.obs)]
  y.oy=y.obs[order(y.obs)]
  d.ox=dx[order(x.obs)]
  d.oy=dy[order(y.obs)]
  Sx.vec=cumprod(1-d.ox/n:1)
  Sy.vec=cumprod(1-d.oy/n:1)
  Sx.vec=c(1,Sx.vec[-n])
  Sy.vec=c(1,Sy.vec[-n])
    
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
          Sx=Sx.vec[sum(x.ox<=xij)];Sy=Sy.vec[sum(y.oy<=yij)]
          Aij=(-log(Sx))^(theta+1)+(-log(Sy))^(theta+1)
          Sij=exp(-Aij^( 1/(theta+1) ))
          Eij=(log(Sij)-theta)/(2*log(Sij)-theta)
          u=u+(((xi-xj)*(yi-yj)>0)-Eij)
        }
      }
    }
    u
  }
  
  if((U.func(lower)<0)&(U.func(upper)<0)){ theta=lower }else{  
    theta=uniroot(U.func,lower=lower,upper=upper)$root
  }
  tau=theta/(theta+1)
  
  if(U.plot==TRUE){  
    curve(U.func,lower,upper,xlab="theta",ylab="U2(theta)")
    points(theta,U.func(theta),col="red",cex=1)
    abline(h=0,lty="dotted",col="blue")
  }
  
  return(c(theta=theta,tau=tau))
}
