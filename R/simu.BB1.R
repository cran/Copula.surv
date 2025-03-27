simu.BB1=function(n,alpha,delta=0,scale1=1,scale2=1,shape1=1,shape2=1,
                  Print=FALSE){

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  a=alpha
  d=delta
  for (i in 1:n){
    U=runif(1, 0, 1)
    W=runif(1, 0, 1)
    U1=U^(-a)-1
    func=function(v){
      v1=v^(-a)-1
      A=(U1^(d+1)+v1^(d+1))^(1/(d+1))
      U^(-a-1)*U1^d*A^(-d)*(1+A)^(-(1+a)/a)-W
    }
    V=uniroot(func,lower=0.00000000001,upper=0.9999999999)$root

    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }

  if(Print==TRUE){
    meanX=scale1^(-1/shape1)*gamma(1+1/shape1)
    meanY=scale2^(-1/shape2)*gamma(1+1/shape2)
    print(
      c(true_Kendall_tau=1-2/(d+1)/(a+2),meanX=meanX,meanY=meanY)
    )
  }

  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
