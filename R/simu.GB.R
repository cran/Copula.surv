simu.GB=function(n,alpha,scale1=1,scale2=1,shape1=1,shape2=1,
                     Print=FALSE){

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  for (i in 1:n){
    U=runif(1, 0, 1)
    W=runif(1, 0, 1)
    GB = function(v){
      (1-alpha*log(v))*U^(-alpha*log(v))*v-W
    }
    V=uniroot(GB,lower=0.00000000001,upper=0.9999999999)$root

    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }

  if(Print==TRUE){
    meanX=scale1^(-1/shape1)*gamma(1+1/shape1)
    meanY=scale2^(-1/shape2)*gamma(1+1/shape2)
    func1=function(x){x*(1-alpha*log(x))*log(1-alpha*log(x))}
    if(alpha==0){Tau=0}else{Tau=1-4/alpha*integrate(func1,0,1)$value}
    print(c(true_Kendall_tau=Tau,meanX=meanX,meanY=meanY))
  }

  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
