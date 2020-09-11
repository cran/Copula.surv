\name{Test.Gumbel}
\alias{Test.Gumbel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A goodness-of-fit test for the Gumbel copula}
\description{
Perform a goodness-of-fit test for the Gumbel copula based on Emura, Lin and Wang (2010).
}
\usage{
Test.Gumbel(x.obs,y.obs,dx,dy,lower=0.1,upper=50,U.plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.obs}{censored times for X}
  \item{y.obs}{censored times for Y}
  \item{dx}{censoring indicators for X}
  \item{dy}{censoring indicators for Y}
  \item{lower}{lower bound for the association parameter}
  \item{upper}{upper bound for the association parameter}
  \item{U.plot}{if TRUE, draw the plot of U_1(theta) and U_2(theta)}
}

\details{
See the references.
}

\value{
  \item{alpha1}{association parameter by the pseudo-likelihood estimator}
  \item{alpha2}{association parameter by the unweighted estimator}
  \item{Stat}{log(alpha1)-log(alpha2)}
  \item{Z}{Z-value of the goodness-of-fit for the Gumbel copula}
  \item{P}{P-value of the goodness-of-fit for the Gumbel copula}
}
\references{
Emura T, Lin CW, Wang W (2010) A goodness-of-fit test for Archimedean copula models
 in the presence of right censoring, Compt Stat Data Anal 54: 3033-43
}
\author{Takeshi Emura}

\examples{
x.obs=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
y.obs=c(2,1,4,5,6,8,3,7,10,9,11,12,13,14,15)
dx=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
dy=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
Test.Gumbel(x.obs,y.obs,dx,dy)

n=20
set.seed(1)
Dat=simu.Clayton(n=n,alpha=2)
C=runif(n,min=0,max=5)
x.obs=pmin(Dat[,"X"],C)
y.obs=pmin(Dat[,"Y"],C)
dx=Dat[,"X"]<=C
dy=Dat[,"Y"]<=C
Test.Gumbel(x.obs,y.obs,dx,dy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Gumbel copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Goodness-of-fit test }% __ONLY ONE__ keyword per line
