U1.Clayton<-function(x.obs,y.obs,dx,dy,lower=0.001,upper=50,U.plot=TRUE){
    n=length(x.obs)

    U.func=function(alpha){
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
            Wij=(alpha+2)/(alpha+1)/(Rij+alpha)
            u=u+Wij*(  ((xi-xj)*(yi-yj)>0)-(alpha+1)/(alpha+2)  )
          }
        }
      }
      u
    }

    if((U.func(lower)<0)&(U.func(upper)<0)){ alpha=lower }else{
      alpha=uniroot(U.func,lower=lower,upper=upper)$root
    }
    tau=alpha/(alpha+2)

    if(U.plot==TRUE){
      curve(U.func,lower,upper,xlab="alpha",ylab="U1(alpha)")
      points(alpha,U.func(alpha),col="red",cex=1)
      abline(h=0,lty="dotted",col="blue")
    }

    return(c(alpha=alpha,tau=tau))
}

