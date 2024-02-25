simu.Frank=function(n,alpha,scale1=1,scale2=1,shape1=1,shape2=1,
                    Print=FALSE){
  if((0<=alpha)&(alpha<0.0001)){alpha=0.0001}
  if((0>alpha)&(-0.0001<alpha)){alpha=-0.0001}

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  for (i in 1:n){
    U=runif(1, 0, 1)
    W=runif(1, 0, 1)
    A=exp(-alpha*U)-W*exp(-alpha*U)+W*exp(-alpha)
    B=exp(-alpha*U)-W*exp(-alpha*U)+W
    V=-log(A/B)/alpha
    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }

  if(Print==TRUE){
    meanX=scale1^(-1/shape1)*gamma(1+1/shape1)
    meanY=scale2^(-1/shape2)*gamma(1+1/shape2)
    func1=function(x){x/(exp(x)-1)}
    Tau=1-4/alpha*(1-integrate(func1,0,alpha)$value/alpha)
    print(c(true_Kendall_tau=Tau,meanX=meanX,meanY=meanY))
  }
  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
