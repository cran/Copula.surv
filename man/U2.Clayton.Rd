\name{U2.Clayton}
\alias{U2.Clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of an association parameter via the unweighted estimator}
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
  \item{theta}{association parameter}
  \item{tau}{Kendall's tau (=theta/(theta+2))}
}
\references{
Emura T, Lin CW, Wang W (2010) A goodness-of-fit test for Archimedean copula models
 in the presence of right censoring, Compt Stat Data Anal 54: 3033-43
}
\author{Takeshi Emura}

\examples{
n=200
theta_true=2 ## association parameter ##
r1_true=1 ## hazard for X
r2_true=1 ## hazard for Y

set.seed(1)
V1=runif(n)
V2=runif(n)
X=-1/r1_true*log(1-V1)
W=(1-V1)^(-theta_true)
Y=1/theta_true/r2_true*log(  1-W+W*(1-V2)^(-theta_true/(theta_true+1))  )
C=runif(n,min=0,max=5)

x.obs=pmin(X,C)
y.obs=pmin(Y,C)
dx=X<=C
dy=Y<=C

U2.Clayton(x.obs,y.obs,dx,dy)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Goodness-of-fit test }% __ONLY ONE__ keyword per line
