\name{simu.BB1}
\alias{simu.BB1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating data from the BB1 copula}
\description{
n pairs of (U,V) are generated from the BB1 copula.
n paris of (X,Y) are generated from the corresponding bivariate survival model with the Weibull marginal distributions. The default parameters (scale1=scale2=shape1=shape2=1) give the unit exponential distributions.
}
\usage{
simu.BB1(n,alpha,d=0,scale1=1,scale2=1,shape1=1,shape2=1,Print=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{alpha}{association (copula) parameter}
  \item{d}{BB1 copula's departure parameter from the Clayton (d=0 is the default)}
  \item{scale1}{scale parameter for X}
  \item{scale2}{scale parameter for Y}
  \item{shape1}{shape parameter for X}
  \item{shape2}{shape parameter for Y}
  \item{Print}{print Kendall's tau and means of X and Y if "TRUE"}
}

\details{
See Section 2.6 of Emura et al.(2019) for copulas and bivariate survival times.
}

\value{
  \item{U}{uniformly distributed on (0,1)}
  \item{V}{uniformly distributed on (0,1)}
  \item{X}{Weibull distributed (scale1, shape1)}
  \item{Y}{Weibull distributed (scale2, shape2)}
}
\references{
Emura T, Lin CW, Wang W (2010) A goodness-of-fit test for Archimedean copula models in the presence of right censoring, Compt Stat Data Anal 54: 3033-43

Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoints, Joint Frailty-Copula Models, JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura}

\examples{
n=100
Dat=simu.BB1(n=n,alpha=1,d=2,scale1=1,scale2=2,shape1=0.5,shape2=2)
plot(Dat[,"U"],Dat[,"V"])
cor(Dat[,"U"],Dat[,"V"],method="kendall")
plot(Dat[,"X"],Dat[,"Y"])
cor(Dat[,"X"],Dat[,"Y"],method="kendall")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BB1 copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Simulation }% __ONLY ONE__ keyword per line
\keyword{ Weibull distribution }% __ONLY ONE__ keyword per line
