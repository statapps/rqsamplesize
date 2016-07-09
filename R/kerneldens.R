kerneldens<-function(dist, tau, kernel.smooth=NULL, bw=NULL){
  n=length(dist)
  m = 100
  pt = (1:m)/(m+1)
  s = sd(dist)
  ## if bw is null, then a good choice is half of Q75-Q25
  if(is.null(bw)) {
    bw = 1.5*s*n^(-0.2)
  }
  fplus=function(x) (x + abs(x))/2
  k0=function(x) dnorm(x)
  K0=function(x) pnorm(x)
  k1=function(x) dunif(x,-1/2,1/2) #also dunif(x+1/2)
  k2=function(x) (fplus(x+1)-2*fplus(x)+fplus(x-1))
  k3=function(x) ((fplus(x+3/2))^2-3*(fplus(x+1/2))^2+3*(fplus(x-1/2))^2+(fplus(x-3/2))^2)/2
  k4=function(x) ((fplus(x+2))^3-4*(fplus(x+1))^3+6*(fplus(x))^34*(fplus(x-1))^3+(fplus(x-2))^3)/6
  K1=function(x) punif(x,-1/2,1/2) #also punif(x+1/2)
  K2=function(x) ((fplus(x+1))^2-2*(fplus(x))^2+(fplus(x-1))^2)/2
  K3=function(x) ((fplus(x+3/2))^3-3*(fplus(x+1/2))^3+3*(fplus(x-1/2))^3-(fplus(x-3/2))^3)/6
  K4=function(x) ((fplus(x+2))^4-4*(fplus(x+1))^4+6*(fplus(x))^4-4*(fplus(x-1))^4+(fplus(x-2))^4)/24
  x = quantile(dist, pt)
  z=x%*%t(rep(1,n))-rep(1,m)%*%t(dist)
  if(is.null(kernel.smooth)) {
    fsmooth=k0(z/bw)%*%rep(1,n)/(bw*n)
    Fsmooth=K0(z/bw)%*%rep(1,n)/n
  }

  if (kernel.smooth==1){
    fsmooth=k1(z/bw)%*%rep(1,n)/(bw*n)
    Fsmooth=K1(z/bw)%*%rep(1,n)/n
  }
  else if (kernel.smooth==2){
    fsmooth=k2(z/bw)%*%rep(1,n)/(bw*n)
    Fsmooth=K2(z/bw)%*%rep(1,n)/n
  }
  else if (kernel.smooth==3){
    fsmooth=k3(z/bw)%*%rep(1,n)/(bw*n)
    Fsmooth=K3(z/bw)%*%rep(1,n)/n
  }
  else if (kernel.smooth==4){
    fsmooth=k4(z/bw)%*%rep(1,n)/(bw*n)
    Fsmooth=K4(z/bw)%*%rep(1,n)/n
  }
  ftau = fsmooth[sum(Fsmooth <= tau)]
  kerinfo=list(x = x, fx = fsmooth, Fx=Fsmooth, ftau = ftau)
  return(kerinfo)
}
#m = 100
#rs = rchisq(m, 16)
#ki = kerneldens(rs, 0.5, kernel.smooth=2)
#hist(rs)
#plot(ki$x, ki$fx, type = 'l')
#lines(ki$x, dchisq(ki$x, 16))
