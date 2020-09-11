\name{U2.Clayton}
\alias{U2.Clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of an association parameter (alpha) via the unweighted estimator}
\description{
Estimate the association parameter of the Clayton copula using bivariate survival data.
 The estimator was defined as the unweighted estimator in Emura, Lin and Wang (2010).
}
\usage{
U2.Clayton(x.obs,y.obs,dx,dy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.obs}{censored times for X}
  \item{y.obs}{censored times for Y}
  \item{dx}{censoring indicators for X}
  \item{dy}{censoring indicators for Y}
}

\details{
Details are seen from the references.
}

\value{
  \item{alpha}{association parameter}
  \item{tau}{Kendall's tau (=alpha/(alpha+2))}
}
\references{
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
U2.Clayton(x.obs,y.obs,dx,dy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Estimation }% __ONLY ONE__ keyword per line
