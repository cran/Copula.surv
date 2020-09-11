simu.Clayton=function(n,alpha,scale1=1,scale2=1,shape1=1,shape2=1){

  V = runif(n, 0, 1)
  W = runif(n, 0, 1)
  U = (V^(-alpha)*W^(-alpha/(alpha+1))-V^(-alpha)+1)^(-1/alpha)

  X=(-log(U)/scale1)^(1/shape1)
  Y=(-log(V)/scale2)^(1/shape2)

  print(c(true_Kendall_tau=alpha/(alpha+2)))
  cbind(U,V,X,Y)
}
