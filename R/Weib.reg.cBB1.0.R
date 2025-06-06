Weib.reg.cBB1.0=function(x.obs,y.obs,dx,dy,zx,zy,zxy,delta=0,
                          convergence.par=FALSE){

  x.obs[x.obs<=0]=min(1,min(x.obs[x.obs>0]))
  y.obs[y.obs<=0]=min(1,min(y.obs[y.obs>0]))

  tx = x.obs
  ty = y.obs
  zx = as.matrix(zx)
  zy = as.matrix(zy)
  zxy = as.matrix(zxy)
  px = ncol(zx)
  py = ncol(zy)
  pxy = ncol(zxy)
  d=delta
  N = length(x.obs)

  mx=sum(dx)
  my=sum(dy)
  mxy=sum(dx*dy)

  count=c(N=N,Event.X=mx,Event.Y=my)
  CEN=c(Censor.X=(N-mx)/N,Censor.Y=(N-my)/N)

  ## Likelihood function ##
  l.func=function(phi){
    gx=exp(pmax(pmin(phi[1:2], 500), -500))
    gy=exp(pmax(pmin(phi[3:4], 500), -500))
    a=min(exp(phi[5]), exp(3))
    bx=phi[(5+1):(5+px)]
    by=phi[(5+px+1):(5+px+py)]
    bxy=phi[(5+px+py+1):(5+px+py+pxy)]
    bzx=as.vector(zx%*%bx)
    bzy=as.vector(zy%*%by)
    bzxy=as.vector(zxy%*%bxy)
    a=pmin(a*exp(bzxy),exp(3))
    Sx=pmin(exp(a*exp(bzx)*gx[1]*tx^gx[2]),exp(500))
    Sy=pmin(exp(a*exp(bzy)*gy[1]*ty^gy[2]),exp(500))
    A=1+((Sx-1)^(d+1)+(Sy-1)^(d+1))^(1/(d+1))
    Ex=Sx/A*((Sx-1)/(A-1))^d
    Ey=Sy/A*((Sy-1)/(A-1))^d
    logR=log(a+1+a*d*A/(A-1))
    l=mx*(log(gx[1])+log(gx[2]))+(gx[2]-1)*sum(dx*log(tx))+sum(dx*bzx)
    l=l+my*(log(gy[1])+log(gy[2]))+(gy[2]-1)*sum(dy*log(ty))+sum(dy*bzy)
    l=l+sum(dx*dy*logR)+sum(dx*log(Ex))+sum(dy*log(Ey))-sum(1/a*log(A))
    -l
  }

  p0=rep(0,5+px+py+pxy)
  p0[c(1,3)]=p0[c(1,3)]-log(c(mean(tx),mean(ty)))
  res=nlm(l.func,p=p0,hessian=TRUE)
  ML=-res$minimum
  HL=-res$hessian
  V=solve(-HL,tol=10^(-50))
  DF=5+px+py+pxy
  AIC=2*DF-2*(-l.func(res$estimate))
  BIC=log(N)*DF-2*(-l.func(res$estimate))
  convergence_res=c(ML=ML,DF=DF,AIC=AIC,BIC=BIC,code=res$code,
                    No.of.iterations=res$iterations)
  est=c(exp(res$est[1:5]),res$est[(5+1):(5+px+py+pxy)])
  D=diag(c(est[1:5],rep(1,px+py+pxy)))
  est_var=D%*%V%*%D

  bx_est=res$est[(5+1):(5+px)]
  by_est=res$est[(5+px+1):(5+px+py)]
  bxy_est=res$est[(5+px+py+1):(5+px+py+pxy)]
  g_est=exp(res$est[1:2])
  h_est=exp(res$est[3:4])
  a_est=exp(res$est[5])
  tau_est=1-2/(a_est+2)/(d+1)
  bx_se=sqrt(diag(V)[(5+1):(5+px)])
  by_se=sqrt(diag(V)[(5+px+1):(5+px+py)])
  bxy_se=sqrt(diag(V)[(5+px+py+1):(5+px+py+pxy)])
  a_se=a_est*sqrt(diag(V)[5])
  tau_se=2/(a_est+2)^2/(d+1)*a_se
  g_var=diag(g_est)%*%V[1:2,1:2]%*%diag(g_est)
  h_var=diag(h_est)%*%V[3:4,3:4]%*%diag(h_est)
  g_se=sqrt(diag(g_var))
  h_se=sqrt(diag(h_var))

  bx_res=c(estimate=bx_est,SE=bx_se,
              Lower=bx_est-1.96*bx_se,Upper=bx_est+1.96*bx_se)
  by_res=c(estimate=by_est,SE=by_se,
              Lower=by_est-1.96*by_se,Upper=by_est+1.96*by_se)
  bxy_res=c(estimate=bxy_est,SE=bxy_se,
           Lower=bxy_est-1.96*bxy_se,Upper=bxy_est+1.96*bxy_se)
  a_Lower=a_est*exp(-1.96*sqrt(diag(V)[5]))
  a_Upper=a_est*exp(1.96*sqrt(diag(V)[5]))
  a_res=c(Estimate=a_est,SE=a_se,Lower=a_Lower,Upper=a_Upper)
  tau_res=c(Estimate=tau_est,SE=tau_se,
            Lower=tau_est-1.96*tau_se,Upper=tau_est+1.96*tau_se)
  rx=c(Estimate=g_est[1],SE=g_se[1],
       Lower=g_est[1]*exp(-1.96*sqrt(diag(V)[1:1])),
       Upper=g_est[1]*exp(+1.96*sqrt(diag(V)[1:1])))
  vx=c(Estimate=g_est[2],SE=g_se[2],
       Lower=g_est[2]*exp(-1.96*sqrt(diag(V)[2:2])),
       Upper=g_est[2]*exp(+1.96*sqrt(diag(V)[2:2])))
  ry=c(Estimate=h_est[1],SE=h_se[1],
       Lower=h_est[1]*exp(-1.96*sqrt(diag(V)[3:3])),
       Upper=h_est[1]*exp(+1.96*sqrt(diag(V)[3:3])))
  vy=c(Estimate=h_est[2],SE=h_se[2],
       Lower=h_est[2]*exp(-1.96*sqrt(diag(V)[4:4])),
       Upper=h_est[2]*exp(+1.96*sqrt(diag(V)[4:4])))

  if (convergence.par==FALSE){convergence.parameters = NULL}else{
    convergence.parameters=list(log_estimate=res$est,gradient=-res$gradient,log_var=V)
  }
  list(Number=count,Proportion=CEN,
       beta_x=bx_res,beta_y=by_res,beta_xy=bxy_res,
       alpha=a_res,tau=tau_res,
       scale_x=rx,shape_x=vx,scale_y=ry,shape_y=vy,
       convergence=convergence_res,
       convergence.parameters=convergence.parameters)
}
