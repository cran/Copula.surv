simu.t=function(n,alpha,df=1,scale1=1,scale2=1,shape1=1,shape2=1,
                     Print=FALSE){

  U.vec=V.vec=X.vec=Y.vec=numeric(n)
  for (i in 1:n){
    dat=mvrnorm(mu=c(0,0),Sigma=matrix(c(1,alpha,alpha,1),ncol=2))
    C2=rchisq(1,df=df)
    X=dat[1]/sqrt(C2/df)
    Y=dat[2]/sqrt(C2/df)
    U=pt(X,df=df)
    V=pt(Y,df=df)
    U.vec[i]=U
    V.vec[i]=V
    X.vec[i]=(-log(U)/scale1)^(1/shape1)
    Y.vec[i]=(-log(V)/scale2)^(1/shape2)
  }

  if(Print==TRUE){
    meanX=scale1^(-1/shape1)*gamma(1+1/shape1)
    meanY=scale2^(-1/shape2)*gamma(1+1/shape2)
    Tau=2/pi*asin(alpha)
    print(c(true_Kendall_tau=Tau,meanX=meanX,meanY=meanY))
  }

  cbind(U=U.vec,V=V.vec,X=X.vec,Y=Y.vec)
}
