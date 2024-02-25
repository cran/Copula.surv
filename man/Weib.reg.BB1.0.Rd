\name{Weib.reg.BB1.0}
\alias{Weib.reg.BB1.0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weibull regression under the BB1 copula with known "delta"}
\description{
See Section 2.6.3 of Emura et al. (2019).
}
\usage{
Weib.reg.BB1.0(x.obs,y.obs,dx,dy,zx,zy,delta=0,convergence.par=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.obs}{censored times for X}
  \item{y.obs}{censored times for Y}
  \item{dx}{censoring indicators for X}
  \item{dy}{censoring indicators for Y}
  \item{zx}{matrix of covariates for X}
  \item{zy}{matrix of covariates for Y}
  \item{delta}{known copula parameter (d>=0)}
  \item{convergence.par}{if TRUE, show the details}
}

\details{
Details are seen from the references.
}

\value{
  \item{beta_x}{regression coefficients for X}
  \item{beta_y}{regression coefficients for Y}
}
\references{
Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoints, Joint Frailty-Copula Models, JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura}

\examples{#TBA}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BB1 copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Regression }% __ONLY ONE__ keyword per line
