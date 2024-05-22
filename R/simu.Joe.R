simu.Joe=function(n,alpha,scale1=1,scale2=1,shape1=1,shape2=1,
                      Print=FALSE){

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  for (i in 1:n){
    U=runif(1, 0, 1)
    W=runif(1, 0, 1)
    Joe=function(v){
      A=(1-U)^alpha+(1-v)^alpha-((1-U)^alpha)*((1-v)^alpha)
      W-A^(1/alpha-1)*(1-(1-v)^alpha)*(1-U)^(alpha-1)
    }
    V=uniroot(Joe,lower=0.00000000001,upper=0.9999999999)$root

    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }

  if(Print==TRUE){
    meanX=scale1^(-1/shape1)*gamma(1+1/shape1)
    meanY=scale2^(-1/shape2)*gamma(1+1/shape2)
    func1=function(x){x*(1-exp(-x))^(2/alpha-2)*exp(-2*x)}
    Tau=1-4/alpha^2*integrate(func1,0,Inf)$value
    print(
      c(true_Kendall_tau=Tau,meanX=meanX,meanY=meanY)
    )
  }

  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
