
simu.BB1reg <-
  function(n,alpha,delta=0,scale1=1,scale2=1,shape1=1,shape2=1,
           beta1=0,beta2=0,beta12=0,Z.dist=runif,...){

  X.vec=Y.vec=Z.vec=U.vec=V.vec=numeric(n)
  for (i in 1:n){
    Z.vec[i]=Z.dist(1,...)
    r1=scale1*exp(beta1*Z.vec[i])
    r2=scale2*exp(beta2*Z.vec[i])
    a=alpha*exp(beta12*Z.vec[i])
    dat=simu.BB1(n=1,alpha=a,delta=delta,scale1=r1,scale2=r2,
                 shape1=shape1,shape2=shape2)
    U.vec[i]=dat[,"U"]
    V.vec[i]=dat[,"V"]
    X.vec[i]=dat[,"X"]
    Y.vec[i]=dat[,"Y"]
  }
  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec,Z=Z.vec)

}

