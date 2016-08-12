###################Define a Class of rqfun for qunatile regression function ###############
# example of usage:
#
#    x = rqfun(mu = 3, sd = 5)
#
#source("./kerneldens.R")
rqfun = function(x, ...) UseMethod("rqfun")

########################Initiate regression function#####################################
rqfun.default = function(mu=0, sd=1, dist='norm', term=c('1'), pos=2, method='exact',
                         a =NA,  b=NA) {
  x = list(mu=mu, sd=sd, dist=dist, term=term, pos=pos, method=method)
  if (x$dist == 'bin') {
    if (x$mu<0||x$mu>1)
      stop("Probability of binary distribution is between 0 and 1.")
    x$sd = sqrt(x$mu*(1-x$mu))
  }

  if (x$dist == 'unif') {
    if(is.na(a)|is.na(b)) {
      a = mu - sqrt(3)*sd
      b = mu + sqrt(3)*sd
    } else {
      x$mu = (a+b)/2
      x$sd = sqrt(1/12)*(b-a)
    }
    x$a = a
    x$b = b
  }
  class(x) = 'rqfun'
  return(x)
}

##################print summary of quantile regression function################################
# example of usage:
#
#    x = rqfun(a = 1, b = 10, dist = 'unif')
#    print(x)
#
print.rqfun = function(x, ...) {
  x.dim = length(x$term)
  pos = x$pos
  cat('The covriate X has the following', x.dim+1, 'components:\n')
  cat('    Intercept \n')
  for (i in 1:x.dim) {
    if (x$term[i] == '1')
      cat('    x     ')
    if (x$term[i] == '-1')
      cat('    1/x')
    if (x$term[i] == 'sqrt')
      cat('    sqrt(x)')
    if (x$term[i] == '2')
      cat('    x^2   ')
    if (x$term[i] == '-2')
      cat('    1/x^2   ')
    if (x$term[i] == '3')
      cat('    x^3   ')
    if (x$term[i] == '-3')
      cat('    1/x^3   ')
    if (x$term[i] == 'log')
      cat('    log(x) ')
    if (x$term[i] == 'exp')
      cat('    exp(x) ')
    for(j in 1:length(pos)) {
      if ((1+1) == pos[j]) cat('    *')
    }
    cat('\n')
  }
  if (x$dist == 'unif') {
    cat('where x has a uniform distribution with a = ', x$a, 'and b = ', x$b, '\n')
  } else {
    cat('where x has a ', x$dist, 'distribution with mean = ', x$mu, 'and variance = ', (x$sd)^2, '\n')
  }
  if (x$method == 'exact')
    cat('Exact ')
  else
    cat('Simulation ')
  cat('method is used to calculate matrix Q = E(t(X)%*%(X)).\n')
  cat('Term(s) with * will be used in sample size calculation')
}

#############Plot qunatile regression function curve #####################
# example of usage:
#
#    x = rqfun(mu = 3, sd = 5, term = c('sqrt', '1', '2'))
#    plot(x)
#
plot.rqfun = function(x, beta=NULL, B = 100, ...) {
  x.dim = length(x$term)
  if(is.null(beta))
    beta = rep(1, (x.dim+1))

  Qs = getQs(x, B = B, X.return = TRUE)
  X = Qs$X
  plot(X[, 1], X[, -1]%*%beta, type = 'l', xlab = 'x', ylab = 'g(x)')
  bc = paste(beta[1])
  for(i in 1:x.dim)
    bc = paste(bc, ', ', beta[i+1])
  title(paste('Quantile regression function for \nbeta = (', bc, ').'))
}

###########Find the matrix Q by simulation method#######################
getQs = function(x, B = 100000, X.return=FALSE) {
  # I removed set.seed because it was used in the main program.
  # We do not use set.seed multiple time in the simulation
  #set.seed(seed)
  if(x$dist == 'norm')
    xi = rnorm(B, x$mu, x$sd)
  if(x$dist == 'bin')
    xi = rbinom(B, 1, x$mu)
  if(x$dist == 'unif')
    xi = runif(B, x$a, x$b)

  x.dim    = length(x$term)
  X        = matrix(0, B, x.dim+1)
  x.lab    = rep(' ', x.dim+1)
  X[, 1]   = rep(1, B)
  x.lab[1] = 'Intercept'
  for (i in 1:x.dim) {
    j = i + 1
    
    if (x$term[i] == '-2') {
      # how to handle xi = 0?
      if (min(xi)==0)
        stop("Error: 0 is in denominator.")
      X[, j] = 1/xi^2
      x.lab[j] = '1/x^2'
    }
    if (x$term[i] == '-1') {
      # how to handle xi = 0?
      if (min(xi)==0)
        stop("Error: 0 is in denominator.")
      X[, j] = 1/xi
      x.lab[j] = '1/x'
    }
    if (x$term[i] == '1/sqrt') {
      # how to handle xi = 0?
      if (min(xi)==0)
        stop("Error: 0 is in denominator.")
      else if (min(xi)<0)
        stop("Error: There is negative value in sqrt term.")
      X[, j] = 1/xi^0.5
      x.lab[j] = '1/sqrt(x)'
    }
    if (x$term[i] == '-3') {
      if (min(xi)==0)
        stop("Error: 0 is in denominator.")
      X[, j] = xi^(-3)
      x.lab[j] = '1/x^3'
    }
    if (x$term[i] == 'sqrt') {
      if (min(xi)<0)
        stop("Error: There is negative value in sqrt term.")
      X[, j] = xi^0.5
      x.lab[j] = 'sqrt(x)'
    }
    if (x$term[i] == '1') {
      X[, j] = xi
      x.lab[j] = 'x'
    }
    if (x$term[i] == '2') {
      X[, j] = xi^2
      x.lab[j] = 'x^2'
    }
    if (x$term[i] == '3') {
      X[, j] = xi^3
      x.lab[j] = 'x^3'
    }
    if (x$term[i] == 'exp') {
      X[, j] = exp(xi)
      x.lab[j] = 'e^x'
    }
    if (x$term[i] == 'log') {
      if (min(xi)<0)
        stop("Error: There is negative value in sqrt term.")
      X[, j] = log(xi)
      x.lab[j] = 'log(x)'
    }
  }
  Q =t(X)%*%X/B
  dimnames(Q) = list(x.lab, x.lab)

  if(X.return) {
    Xb = cbind(xi, X)
    Xb = Xb[1:min(B, 100), ]
    tmp = sort(Xb[, 1], index.return=TRUE)
    Xb = Xb[tmp$ix, ]
    return(list(Q=Q, X=Xb))
  } else {
    return(Q)
  }
}

getQe = function(x){
  # Matrix for polynomial of x
  if (x$dist=="norm"){
    u1 = x$mu
    u2 = x$mu^2 + x$sd^2
    u3 = x$mu^3 + 3*x$mu*x$sd^2
    u4 = x$mu^4 + 6*x$mu^2*x$sd^2 + 3*x$sd^4
    u5 = x$mu^5 + 10*x$mu^3*x$sd^2 + 15*x$mu*x$sd^4
    u6 = x$mu^6 + 15*x$mu^4*x$sd^2 + 45*x$mu^2*x$sd^4+15*x$sd^6
  }
  else if (x$dist=="unif"){
    if (is.na(x$a)&&is.na(x$b)){
      a = x$mu-sqrt(3)*x$sd
      b = x$mu+sqrt(3)*x$sd
    }
    else {
      a=x$a
      b=x$b
    }
    u1 = (a+b)/2
    u2 = (a^2+a*b+b^2)/3
    u3 = (a+b)*(a^2+b^2)/4
    u4 = (a^4+a^3*b+a^2*b^2+a*b^3+b^4)/5
    u5 = (b^6-a^6)/(6*(b-a))
    u6 = (b^7-a^7)/(7*(b-a))
  }
  else {
    stop("ERROR: x shall be either normal or uniform distribution.")
  }

  Ex = c(1,u1,u2,u3,u4,u5,u6)
  d=length(x$term)+1
  x.power = matrix(c(0,as.numeric(x$term)), d, d)
  index = x.power+t(x.power)+1
  Q=matrix(Ex[index], d, d)
}

##########################################################################
# Calculate the asymptotic variance matrix of quantile regression beta hat
#
qrV=function(x, sd, tau, dist="Norm",kernel.smooth,bw,subint){
  #
  # x       Covariate
  # s       sd of error distribution
  # dist    distribution of error: Norm, Cauchy, Gamma or a vector of residual

  x.dim = length(x$term)

  if(x.dim == 1)
    Q=matrix(c(1,x$mu,x$mu,x$mu^2+x$sd^2),nrow = 2,ncol = 2)


  if(x.dim > 1) {
    # use simulation method to find Q
    if (x$method == 'sim')
      Q = getQs(x)
    else # use exact method to find Q
      Q = getQe(x)
  }

  Q=solve(Q)
  if (is.null(dist)||is.character(dist)==TRUE){
  if (is.null(dist)||dist=="Norm") {
    if (length(sd)>1)
      stop("Need a single number for standard deviation.")
      mu=uniroot(function(m) pnorm(0,m,sd)-tau,c(-1e10,1e+10))$root  
    df2 = dnorm(0,mu,sd)^2
  }
  else if (dist=="Cauchy") {
    if (length(sd)>1)
      stop("Need a single number for Cauchy scale parameter.")
    loc=uniroot(function(l) pcauchy(0,l,sd)-tau,c(-1e10,1e+10))$root
    df2 = dcauchy(0,loc,sd)^2
  }
  else if(dist=='Gamma') {
    # Error Ui has Gamma distribution with Gamma(k, 1/k), Var(Ui) = 1/k = s^2
    # log(Yi) = x%*%beta + log(Ui) - log(U_tau)
    # Vi = logUi - logU_tau, Vi has density density U_tau*exp(v)*f_Ui(U_tau*exp(v))
    # V_tau = 0
    k = 1/sd^2
    U_tau = qgamma(tau,k,scale = 1/k)
    df2 = (U_tau*dgamma(U_tau,k,scale = 1/k))^2
  }
  }
  else if(is.numeric(dist)==TRUE){
      df2=kerneldens(dist,tau,kernel.smooth,bw)$ftau
      df2=df2^2
  }
  Var = tau*(1-tau)/df2*Q
  return(Var)
}
