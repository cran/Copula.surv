U2.Clayton<-function(x.obs,y.obs,dx,dy){
  n=length(x.obs)
  m=delta=0
  for(i in 1:(n-1)){
    xi=x.obs[i];yi=y.obs[i];dxi=dx[i];dyi=dy[i]
    for(j in (i+1):n){
      xj=x.obs[j];yj=y.obs[j];dxj=dx[j];dyj=dy[j]
      Ox=((xi<=xj)&(dxi))|((xj<=xi)&(dxj))
      Oy=((yi<=yj)&(dyi))|((yj<=yi)&(dyj))
      if(Ox&Oy){
        m=m+1
        delta=((xi-xj)*(yi-yj)>0)+delta
      }
    }
  }
  
  delta=delta/m
  theta=(2*delta-1)/(1-delta)
  tau=theta/(theta+2)
 
  return(c(theta=theta,tau=tau))
}
  