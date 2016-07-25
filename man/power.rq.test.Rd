\name{power.rq.test}
\alias{power.rq.test}
\alias{power.rq.test}


 \title{
  Power and sample size for quantile regression.
   }
 \description{

    Compute power and sample size of test under alternative hypothesis to obtain target power (same as power.anova.test)   }
 \usage{
     power.rq.test(x, n = NULL, sig.level = 0.05, power = NULL, 
     tau = 0.5, beta = 1, sd = 1, dist = "Norm", kernel.smooth = NULL, 
     bw = NULL, alternative = c("two.sided", "one.sided"))

 }
    \arguments{
     \item{x}{A \code{\link{rqfun}} object.See "Details" and "Examples".}
     \item{n}{Given sample size to compute power.}
     \item{sig.level}{Significance level of the test(Type I error probability). The default is 0.05.}
     \item{power}{Power of the test between 0 and 1(One minus Type II error probability). It is required to compute sample size.}
     \item{tau}{The desired regression quantile between 0 and 1. The default is 0.5.}
     \item{beta}{The desired quantile regression coefficient. It can be a number or a vector.}
     \item{sd}{Standard deviation or scale of error distribution.}
     \item{dist}{The error distribution for the purpose of power/sample size calculation. It takes "Norm","Cauchy","Gamma" or a vector of residuals when distribution is unknown and to be estimated. }
     \item{kernel.smooth}{It takes values 1,2,3,4 indicating the number of kernels used. 1 is for a uniform variable bounded by -1/2 and 1/2. 2 is for triangle density bounded by -1 and 1. 3 is for three-piece density bounded by -3/2 and 3/2. 4 is for four-piece density bounded by -2 and 2. Standard normal kernel is used by default if no values are specified.}
     \item{bw}{The bandwidth used in kernel smoothing. Small bandwidth results in under-smoothing whereas large bandwidth causes over-smoothing. Default is NULL and bandwidth is estimated by normal distribution approximation or Silverman's (1986) rule of thumb by default.}
     \item{alternative}{Choose either "one.sided" for one sided test or "two.sided" for two sided test.}
     }
\value{
Object of class "power.htest", a list of the arguments (including the computed one) augmented with method and note elements.
}
 \details{
A \code{\link{rqfun}} object needs to be defined in univariate regression situation. Here "univariate" means the independent variables are transformation of a single variable. It is defined as following: x = rqfun(mu=0, sd=1, dist='norm', term=c('1'), pos=2, method='exact', a =NA,  b=NA). It contains the information of the independent variable.

It is assumed that error distribution is iid, so the calculation of variance of regression coefficients are based on Section 3.2.2 of Quantile Regression(2005).
 }
 \note{
 Either n or power needs to be defined. They can not be NULL at the same time. If estimated sample size is desired, power should be given; if estimated power is desired, actual sample size n should be given.

uniroot is used to solve power equation for unknowns, so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
 }
\author{Zhenxian Gong}
\seealso{rq, rqfun, power.anova.test, uniroot}
 \references{
Roger Koenker. (2005)Quantile Regression. Cambridge University Press.

Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis. London: Chapman & Hall/CRC.
 }
 \examples{
 #Construct x
 x = rqfun(mu = 5, sd = 1.5, dist = "norm", pos = 2, term = c('1'),a=NA,b=NA, method = 'exact')
 #Compute estimated sample size when power is 0.8
 power.rq.test(x=x, power = 0.8, tau = 0.2, beta = c(1,0.5),
                sd = 10, dist = "Norm", alternative = 'two.sided')
 }
