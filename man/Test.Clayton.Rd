\name{Test.Clayton}
\alias{Test.Clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A goodness-of-fit test for the Clayton copula}
\description{
Perform a goodness-of-fit test for the Clayton copula based on Emura, Lin and Wang (2010).
 The test is asymptotically equivalent to the test of Shih (1998).
}
\usage{
Test.Clayton(x.obs,y.obs,dx,dy,lower=0.001,upper=50,U.plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.obs}{censored times for X}
  \item{y.obs}{censored times for Y}
  \item{dx}{censoring indicators for X}
  \item{dy}{censoring indicators for Y}
  \item{lower}{lower bound for the association parameter}
  \item{upper}{upper bound for the association parameter}
  \item{U.plot}{if TRUE, draw the plot of U_1(theta)}
}

\details{
See the references.
}

\value{
  \item{theta1}{association parameter by the pseudo-likelihood estimator}
  \item{theta2}{association parameter by the unweighted estimator}
  \item{Stat}{log(theta1)-log(theta2)}
  \item{Z}{Z-value of the goodness-of-fit for the Clayton copula}
  \item{P}{P-value of the goodness-of-fit for the Clayton copula}
}
\references{
Emura T, Lin CW, Wang W (2010) A goodness-of-fit test for Archimedean copula models
 in the presence of right censoring, Compt Stat Data Anal 54: 3033-43

Shih JH (1998) A goodness-of-fit test for association in a bivariate survival model.
 Biometrika 85: 189-200
}
\author{Takeshi Emura}

\examples{
n=20
theta_true=2 ## association parameter ##
r1_true=2 ## hazard for X
r2_true=2 ## hazard for Y

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

Test.Clayton(x.obs,y.obs,dx,dy)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Goodness-of-fit test }% __ONLY ONE__ keyword per line
