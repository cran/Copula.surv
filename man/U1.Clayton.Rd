\name{U1.Clayton}
\alias{U1.Clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of an association parameter via the pseudo-likelihood}
\description{
Estimate the association parameter (alpha) of the Clayton copula using bivariate survival data.
 The estimator was derived by Clayton (1978) and reformulated by Emura, Lin and Wang (2010).
}
\usage{
U1.Clayton(x.obs,y.obs,dx,dy,lower=0.001,upper=50,U.plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.obs}{censored times for X}
  \item{y.obs}{censored times for Y}
  \item{dx}{censoring indicators for X}
  \item{dy}{censoring indicators for Y}
  \item{lower}{lower bound for the association parameter}
  \item{upper}{upper bound for the association parameter}
  \item{U.plot}{if TRUE, draw the plot of U_1(alpha)}
}

\details{
Details are seen from the references.
}

\value{
  \item{alpha}{association parameter}
  \item{tau}{Kendall's tau (=alpha/(alpha+2))}
}
\references{
Clayton DG (1978). A model for association in bivariate life tables and its
 application to epidemiological studies of familial tendency in chronic disease
 incidence. Biometrika 65: 141-51.

Emura T, Lin CW, Wang W (2010) A goodness-of-fit test for Archimedean copula models
 in the presence of right censoring, Compt Stat Data Anal 54: 3033-43
}
\author{Takeshi Emura}

\examples{
n=100
Dat=simu.Clayton(n=n,scale1=1,scale2=2,shape1=0.5,shape2=2,alpha=1)
x.obs=Dat[,"X"]
y.obs=Dat[,"Y"]
dx=dy=rep(1,n) ## uncensored data
U1.Clayton(x.obs,y.obs,dx,dy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Estimation }% __ONLY ONE__ keyword per line
