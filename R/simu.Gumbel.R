simu.Gumbel=function(n,alpha,scale1=1,scale2=1,shape1=1,shape2=1){

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  for (i in 1:n){
    U=runif(1, 0, 1)
    W=runif(1, 0, 1)
    Gumbul = function(v){
      A=(-log(U))^alpha/U*((-log(U))^(alpha+1)+(-log(v))^(alpha+1))^(-alpha/(alpha+1))
      A*exp(-((-log(U))^(alpha+1)+(-log(v))^(alpha+1))^(1/(alpha+1)))-W
    }
    V=uniroot(Gumbul,lower=0.00000000001,upper=0.9999999999)$root

    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }
  print(c(true_Kendall_tau=alpha/(alpha+1)))
  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
