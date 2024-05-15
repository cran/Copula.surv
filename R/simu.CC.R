simu.CC=function(n,alpha,scale1=1,scale2=1,shape1=1,shape2=1,
                     Print=FALSE){

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  for (i in 1:n){
    U=runif(1, 0, 1)
    W=runif(1, 0, 1)
    CC = function(v){
      (v-alpha*U*v*(1-v))*exp(alpha*(1-U)*(1-v))-W
    }
    V=uniroot(CC,lower=0.00000000001,upper=0.9999999999)$root

    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }

  if(Print==TRUE){
    meanX=scale1^(-1/shape1)*gamma(1+1/shape1)
    meanY=scale2^(-1/shape2)*gamma(1+1/shape2)

    func1=function(u){
      func2=function(v){
        C=u*v*exp(alpha*(1-u)*(1-v))
        D=(1+(3*u*v-u-v)*alpha+u*v*(1-u)*(1-v)*alpha^2)
        C*D*exp(alpha*(1-u)*(1-v))
      }
      integrate(func2,0,1)$value
    }
    Tau=4*integrate(Vectorize(func1),0,1)$value-1
    print(c(true_Kendall_tau=Tau,meanX=meanX,meanY=meanY))
  }

  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
