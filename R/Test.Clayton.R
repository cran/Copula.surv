Test.Clayton<-function(x.obs,y.obs,dx,dy,lower=0.001,upper=50){
  theta1=U1.Clayton(x.obs,y.obs,dx,dy)["theta"]
  theta2=U2.Clayton(x.obs,y.obs,dx,dy)["theta"]
  stat=log(theta1)-log(theta2)
  
  Stat<-function(x.obs,y.obs,dx,dy){
    theta1=U1.Clayton(x.obs,y.obs,dx,dy,lower,upper)["theta"]
    theta2=U2.Clayton(x.obs,y.obs,dx,dy)["theta"]
    log(theta1)-log(theta2)
  }
  n=length(x.obs)
  T.del=numeric(n)
  for(i in 1:n){
    T.del[i]=Stat(x.obs[-i],y.obs[-i],dx[-i],dy[-i])
  }
  V=(n-1)^2/n*var(T.del)
  Z=stat/sqrt(V)
  P=1-pchisq(Z^2,df=1)
  
  Res=c(theta1=unname(theta1),theta2=unname(theta2),Stat=stat,Z=Z,P=P)
  return(Res)
}
