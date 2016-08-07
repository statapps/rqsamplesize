############################################################################################
#
# Author:   Zhenxian Gong
# Crated:   May 3, 2016
# Revised:  Jun 14, 2016
#
# Versions history:  V0.1: sample size for univariate quantile regression
#                    v0.2: add sample size for multivariate quantile regression
#                          through a non-central chisquare test.
#
# Variance of the quanitle regression coefficient.
############################################################################################

source("./library.R")

############################################################################################
## Sample size for quantile regression

power.rq.test = function (x, n = NULL, sig.level = 0.05, power = NULL, 
                                  tau = 0.5, delta = 1, sd = 1, dist = "Norm", kernel.smooth = "norm", 
                                  bw = NULL, alternative = c("two.sided", "one.sided"))
{
  #
  #x        is data in vector or matrix form
  #sig.level    is Type I error and power is Type II error
  #tau      is the desired quantile i.e 0.1,0.25,0.5,0.75,0.9
  #beta       is a vector of all estimated regression coefficients in quantile regression
  #sigm     is variance?
  #
  # dist    is the distribution of error term: it can be Normal ('norm'), Cauchy ('cauchy')
  #         or mixed normal ('mixed')
  #
  #Note that since we test the coefficient is 0 under null hypothesis, beta[pos] is also the
  #difference, delta in Section 1.3.3 of Sample Size Calculations in Clinical Research 2008
  if (sum(sapply(list(n, power), is.null)) != 1)
  stop("exactly one of 'n', 'power' must be NULL")
  if (sig.level<0||sig.level>1)
    stop("Type I error is between 0 and 1.")
  if(!is.null(power)) {
    if (power<0||power>1)
      stop("Type II error is between 0 and 1.")
  }
  beta = c(1, delta)
  alternative = match.arg(alternative)
  tside = switch(alternative, one.sided = 1, two.sided = 2)

  NOTE = paste('Sample size based on ', x$method, 'method')
  METHOD = 'Sample size for quantile regression'
  pos = x$pos
  b   = beta[pos]
  Vq  = qrV(x, sd, tau, dist,kernel.smooth,bw)[pos, pos]
  p   = length(beta)
  df  = length(pos)

  if (df == 1) {
    if(is.null(n)) {
      p.body = quote({
        1-pt(qt(1-sig.level/tside,n-2), ncp = b/(sqrt(Vq/n)), df=n-2)
      })
      n = uniroot(function(n) eval(p.body) - power, c(3,1e+06))$root
    # n=Vq*(qnorm(1-sig.level/tside)+qnorm(power))^2/(delta^2)
    }
    if (is.null(power)) {
      if (n-2 <= 0)
        stop("Sample size has to be greater than 2.")
      power = 1-pt(qt(1-sig.level/tside,n-2), ncp = b/(sqrt(Vq/n)), df=n-2)
    }
  }
  if(df > 1) {
    lambda = t(b)%*%solve(Vq)%*%b
    p.body = quote({
      c=p*(n-1)/(n-p)
      #cat(c)
      1-c*pf(c*qf(1-sig.level,df1 = df,df2 = n-p), ncp = n*lambda, df1=df,df2 = n-p)
    })
    if(is.null(n)) {
      n = uniroot(function(n) eval(p.body) - power, c(p+1,1e+06))$root
    } else if(is.null(power)) {
      power = eval(p.body)
    }
  }
  n = ceiling(n)

  structure(list(n = n, delta = delta, pos = pos, sd = sd, sig.level = sig.level,
                 power = power, alternative = alternative, note = NOTE,
                 method = METHOD), class = "power.htest")
}



########################################Examples##################################################
#
# Below are examples of using sample size for quantile regression, please do not delete
#
##############################################################################################
#delta = c(.23, -.015)
#x = rqfun(dist = "unif", pos = c(2, 3), term = c('2','3'),a= 3,b= 8, method = 'exact')
#pw = power.rq.test(x=x, power = 0.8, delta=delta, sd = 2)
#print(pw)
#plot(x, beta = beta)

# test multivariate Q:
#
# x is a class that describe the distribution and structure of X
# for example if I add a component "term" to x
# say x$term = c('1', '-1', 'sqrt', '2', 'log')
# and let x$dim = 2, x$pos = 5
# that means X = (1, x, x^-1, sqrt(x), x^2, log(x))
#

#tau = c(0.1,0.3,0.5,0.6,0.8)
#x = qrfun(mu = 5, sd = 3, dist = "norm", pos = 2, term = c('1'),a=NA,b=NA, method = 'exact')


##########Examples of using the sample size function ###############################

#n = samplesize(x, sig.level, power, tau[1], beta, sd, dist = "Mix", w)
#cat('tau = ', tau[16], 'n = ', n, '\n')

#x$sd = 0.8
#x$dist = 'norm'
#sd = 2
#checkX(x)

#plot(tau, n.tau, type = 'l')

#x$mu = 10
#x$sd = 0.8
#x$dist = 'norm'
#x$method = 'sim'
#x$term = c('1', '-1')
#checkX(x)
#Q = getQs(x)
#print(Q)
#sd=c(2,3)
#n = samplesize(x, sig.level, power, tau = 0.4, beta=c(3, 5), sd, dist = "Mix", w)
