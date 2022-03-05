\name{U1.Gumbel}
\alias{U1.Gumbel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of an association parameter via the unweighted estimator}
\description{
Estimate the association parameter of the Gumbel copula using bivariate survival data.
 The estimator was derived by Emura, Lin and Wang (2010).
}
\usage{
U1.Gumbel(x.obs,y.obs,dx,dy,lower=0.01,upper=50,U.plot=TRUE)
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
x.obs=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
y.obs=c(2,1,4,5,6,8,3,7,10,9,11,12,13,14,15)
dx=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
dy=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
U1.Gumbel(x.obs,y.obs,dx,dy)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Gumbel copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Goodness-of-fit test }% __ONLY ONE__ keyword per line
