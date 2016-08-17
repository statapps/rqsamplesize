\name{rqfun}
\alias{rqfun}
\alias{rqfun}


\title{
  rqfun object.
}
\description{
  Create a rqfun object for power.rq.test.
  }
\usage{
  x = rqfun(mu=0, sd=1, dist='norm', term=c('1'), pos=2, method='exact', a =NA,  b=NA)
}
\arguments{
  \item{mu}{Mean of independent variable under univariate case.}
  \item{sd}{Standard deviation of independent variable under univariate case..}
  \item{dist}{The distribution which takes "norm" as normal and "unif" as uniform.}
  \item{term}{It takes a vector of strings of the form of integers or "log", "sqrt" and "exp" eg. term=c("1","-1","log", "sqrt") indicates there are 4 variables: x,1/x,log(x) and sqrt(x)}
  \item{pos}{The position of the regression coefficient(s) of interest (position of intercept is 1, but it is not important, so 1 is not recommended). Default is 2.}
  \item{method}{"exact" or "sim" which are the methods to compute variance of regression coefficients, "exact" method is recommended when variables are chosen from x,x^2 and x^3.}
  \item{a}{If dist is "unif", this is the lower limit of uniform distribution. It is optional if "mu" and "sd" are already specified.}
  \item{b}{If dist is "unif", this is the upper limit of uniform distribution. It is optional if "mu" and "sd" are already specified.}
  }

\examples{
    #Construct x
    x = rqfun(mu = 5, sd = 1.5, term = c('1','log','1/sqrt'))
    }
    