\name{simu.BB1reg}
\alias{simu.BB1reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating data from the BB1 copula regression model}
\description{
n pairs of (U,V) are generated from the BB1 copula.
n paris of (X,Y) are generated from the corresponding bivariate survival model with the Weibull marginal distributions. The default parameters (scale1=scale2=shape1=shape2=1) give the unit exponential distributions.
}
\usage{
simu.BB1reg(n,alpha,delta=0,scale1=1,scale2=1,shape1=1,shape2=1,
 beta1=0,beta2=0,beta12=0,Z.dist=runif,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{alpha}{association (copula) parameter}
  \item{delta}{BB1 copula's departure parameter from the Clayton (0 is the default)}
  \item{scale1}{scale parameter for X}
  \item{scale2}{scale parameter for Y}
  \item{shape1}{shape parameter for X}
  \item{shape2}{shape parameter for Y}
  \item{beta1}{regression coefficient for X}
  \item{beta2}{regression coefficient for Y}
  \item{beta12}{regression coefficient for copula}
  \item{Z.dist}{distribution for covariates}
  \item{...}{parameters for Z.dist}
}

\details{
See Section 2.6 of Emura et al.(2019) for copulas and bivariate survival times.
}

\value{
  \item{U}{uniformly distributed on (0,1)}
  \item{V}{uniformly distributed on (0,1)}
  \item{X}{Weibull distributed (scale1, shape1)}
  \item{Y}{Weibull distributed (scale2, shape2)}
  \item{Z}{Covariates}
}
\references{
Emura T, Lin CW, Wang W (2010) A goodness-of-fit test for Archimedean copula models in the presence of right censoring, Compt Stat Data Anal 54: 3033-43

Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoints, Joint Frailty-Copula Models, JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura}

\examples{
n=10
simu.BB1reg(n=n,alpha=1,delta=2,scale1=1,scale2=2,shape1=0.5,beta1=1,beta2=-1,beta12=2,shape2=2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BB1 copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Simulation }% __ONLY ONE__ keyword per line
\keyword{ Weibull distribution }% __ONLY ONE__ keyword per line
