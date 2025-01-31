Weib.reg.Frank=function(x.obs,y.obs,dx,dy,zx,zy,
                          convergence.par=FALSE){

  x.obs[x.obs<=0]=min(1,min(x.obs[x.obs>0]))
  y.obs[y.obs<=0]=min(1,min(y.obs[y.obs>0]))

  tx = x.obs
  ty = y.obs
  zx = as.matrix(zx)
  zy = as.matrix(zy)
  px = ncol(zx)
  py = ncol(zy)
  N = length(x.obs)

  mx=sum(dx)
  my=sum(dy)
  mxy=sum(dx*dy)

  count=c(N=N,Event.X=mx,Event.Y=my)
  CEN=c(Censor.X=(N-mx)/N,Censor.Y=(N-my)/N)

  ## Log-likelihood function ##
  l.func=function(phi){
    gx=exp(pmax(pmin(phi[1:2], 500), -500))
    gy=exp(pmax(pmin(phi[3:4], 500), -500))
    #a=phi[5]
    a=max(min(phi[5],100),-100)
    bx=phi[(5+1):(5+px)]
    by=phi[(5+px+1):(5+px+py)]
    bzx=as.vector(zx%*%bx)
    bzy=as.vector(zy%*%by)
    Sx=pmax(exp(-exp(bzx)*gx[1]*tx^gx[2]),0.00001)
    Sy=pmax(exp(-exp(bzy)*gy[1]*ty^gy[2]),0.00001)
    ex=pmin(exp(-a*Sx),exp(500))
    ey=pmin(exp(-a*Sy),exp(500))
    ea=exp(-a)
    S=-1/a*log(1+(ex-1)*(ey-1)/(ea-1))
    #D=Ea-1+(Ex-1)*(Ey-1)
    #Dx=Ex*(Ey-1)/D
    #Dy=Ey*(Ex-1)/D
    #Dxy=-a*(Ea-1)*Ex*Ey/(D^2)
    Ex=Sx*ex*(ey-1)/(ea-1)/exp(-a*S)/S
    Ey=Sy*ey*(ex-1)/(ea-1)/exp(-a*S)/S
    R=a*S/(1-exp(-a*S))
    #l=mx*(log(gx[1])+log(gx[2]))+(gx[2]-1)*sum(dx*log(tx))+sum(dx*bzx)
    #l=l+my*(log(gy[1])+log(gy[2]))+(gy[2]-1)*sum(dy*log(ty))+sum(dy*bzy)
    #l=l+sum(dx*dy*log(Dxy))+sum(dx*(1-dy)*log(Dx))+sum(dy*(1-dx)*log(Dy))
    #l=1+sum((1-dx)*(1-dy)*log(S))
    l=mx*(log(gx[1])+log(gx[2]))+(gx[2]-1)*sum(dx*log(tx))+sum(dx*bzx)
    l=l+my*(log(gy[1])+log(gy[2]))+(gy[2]-1)*sum(dy*log(ty))+sum(dy*bzy)
    l=l+sum(dx*dy*log(R))+sum(dx*log(Ex))+sum(dy*log(Ey))+sum(log(S))
    -l
  }

  p0=rep(0,5+px+py)
  p0[5]=1 ### copula parameter cannot be zero
  p0[c(1,3)]=p0[c(1,3)]-log(c(mean(tx),mean(ty)))
  res=nlm(l.func,p=p0,hessian=TRUE)
  ML=-res$minimum
  HL=-res$hessian
  V=solve(-HL,tol=10^(-50))
  DF=5+px+py
  AIC=2*DF-2*(-l.func(res$estimate))
  BIC=log(N)*DF-2*(-l.func(res$estimate))
  convergence_res=c(ML=ML,DF=DF,AIC=AIC,BIC=BIC,code=res$code,
                    No.of.iterations=res$iterations)
  est=c(exp(res$est[1:5]),res$est[(5+1):(5+px+py)])
  D=diag(c(est[1:5],rep(1,px+py)))
  est_var=D%*%V%*%D

  bx_est=res$est[(5+1):(5+px)]
  by_est=res$est[(5+px+1):(5+px+py)]
  g_est=exp(res$est[1:2])
  h_est=exp(res$est[3:4])
  a_est=res$est[5]
  func1=function(x){x/(exp(x)-1)}
  tau_est=1-4/a_est*(1-integrate(func1,0,a_est)$value/a_est)

  bx_se=sqrt(diag(V)[(5+1):(5+px)])
  by_se=sqrt(diag(V)[(5+px+1):(5+px+py)])
  a_se=sqrt(diag(V)[5])
  dtau=1/a_est*( 2*(1-tau_est)+4*(1/(exp(a_est)-1)-1/a_est) )
  tau_se=abs(dtau)*a_se
  g_var=diag(g_est)%*%V[1:2,1:2]%*%diag(g_est)
  h_var=diag(h_est)%*%V[3:4,3:4]%*%diag(h_est)
  g_se=sqrt(diag(g_var))
  h_se=sqrt(diag(h_var))

  bx_res=c(estimate=bx_est,SE=bx_se,
              Lower=bx_est-1.96*bx_se,Upper=bx_est+1.96*bx_se)
  by_res=c(estimate=by_est,SE=by_se,
              Lower=by_est-1.96*by_se,Upper=by_est+1.96*by_se)
  a_Lower=a_est-1.96*sqrt(diag(V)[5])
  a_Upper=a_est+1.96*sqrt(diag(V)[5])
  a_res=c(Estimate=a_est,SE=a_se,Lower=a_Lower,Upper=a_Upper)
  tau_res=c(Estimate=tau_est,SE=tau_se,
            Lower=tau_est-1.96*tau_se,Upper=tau_est+1.96*tau_se)
  rx=c(Estimate=g_est[1],SE=g_se[1],Lower=g_est[1]*exp(-1.96*sqrt(diag(V)[1:1])),
       Upper=g_est[1]*exp(+1.96*sqrt(diag(V)[1:1])))
  vx=c(Estimate=g_est[2],SE=g_se[2],Lower=g_est[2]*exp(-1.96*sqrt(diag(V)[2:2])),
       Upper=g_est[2]*exp(+1.96*sqrt(diag(V)[2:2])))
  ry=c(Estimate=h_est[1],SE=h_se[1],Lower=h_est[1]*exp(-1.96*sqrt(diag(V)[3:3])),
       Upper=h_est[1]*exp(+1.96*sqrt(diag(V)[3:3])))
  vy=c(Estimate=h_est[2],SE=h_se[2],Lower=h_est[2]*exp(-1.96*sqrt(diag(V)[4:4])),
       Upper=h_est[2]*exp(+1.96*sqrt(diag(V)[4:4])))

  if (convergence.par==FALSE){convergence.parameters = NULL}else{
    convergence.parameters=list(log_estimate=res$est,gradient=-res$gradient,log_var=V)
  }
  list(Number=count,Proportion=CEN,
       beta_x=bx_res,beta_y=by_res,alpha=a_res,tau=tau_res,
       scale_x=rx,shape_y=ry,scale_y=ry,shape_y=ry,
       convergence=convergence_res,convergence.parameters=convergence.parameters)
}
