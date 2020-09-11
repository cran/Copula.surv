Test.Gumbel<-function(x.obs,y.obs,dx,dy,lower=0.1,upper=50,U.plot=TRUE){
  alpha1=U1.Gumbel(x.obs,y.obs,dx,dy,lower,upper,U.plot)["alpha"]
  alpha2=U2.Gumbel(x.obs,y.obs,dx,dy,lower,upper,U.plot)["alpha"]
  stat=log(alpha1)-log(alpha2)
  stat=unname(stat)

  Stat<-function(x.del,y.del,dx.del,dy.del){
    alpha1=U1.Gumbel(x.del,y.del,dx.del,dy.del,lower,upper,U.plot=FALSE)["alpha"]
    alpha2=U2.Gumbel(x.del,y.del,dx.del,dy.del,lower,upper,U.plot=FALSE)["alpha"]
    stat=log(alpha1)-log(alpha2)
    unname(stat)
  }
  n=length(x.obs)
  T.del=numeric(n)
  for(i in 1:n){
    T.del[i]=Stat(x.obs[-i],y.obs[-i],dx[-i],dy[-i])
  }
  V=(n-1)^2/n*var(T.del)
  Z=stat/sqrt(V)
  P=1-pchisq(Z^2,df=1)

  Res=c(alpha1=unname(alpha1),alpha2=unname(alpha2),Stat=stat,Z=Z,P=P)
  return(Res)
}
