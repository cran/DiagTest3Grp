Normal.VUS <-
function(x,y,z,p=0,q=0,alpha=0.05,subdivisions=50000,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05,...)
  {

    ######################################################################################################################################################
    #######This function calculate under normal assumtpion, the whole VUS (by setting p=q=0) over the whole domain (0,1)*(0,1)*(0,1)######################
    #######             or partial VUS (specify p and q) requriring minumum spec (p) and sens (q) by only integrating the domain of (p,1)*(0,1)*(q,1)
    ######## and their variance and CI
    ########Also calculate the sample size based on the information from the input data such that the VUS is estimated within 5%
    
    ####1. Input:
    ###(1)x,y,z:vectors are the test marker measurements for the 3 ordinal group (D-,D0,D+)
    ### NOTE: we assume mean values of the three groups increase with severity of disease from D- to D0 and D+,if reverse, simply take negated x,y,z as inputs
    ###(2)p: minimum desired specificity,0<=p<1
    ###(3)q: minimum desired sensitivity,
    ####(4) alpha: calculate (1-alpha)% Confidence interval
    ####(5) subdivisions: # of subintervals for integration using adaptive quadrature
    ####(6)lam.minus, lam0,lam.plus: for sample size calculation, the expected proportion of samples in the D-, D0 and D+ group, which can be equal or not
    ####(7)typeIerror:type I error rate for sample size calculation,default= 0.05, give 95% CI
    ####(8)margin: for sample size calculation, margin of error on the VUS estimates,,default=0.05. The normal (1-typeIerror)% CI is (VUS-Z_typeIerror*SE(VUS),VUS-Z_typeIerror*SE(VUS)), the sample size calculation will be calculated such that Z_typeIerror*SE(VUS)=margin
    
    ####(9)..., other arguments used in the R function integrate() can be passed along, such as, abs.tol,rel.tol,stop.on.error etc
    
    ###Outpout: partial VUS , a value
    ######################################################################################################################################################

    #####the functions used for integration 
    f0 <- function(s,a,b,c,d,p,q)
      {
        ###integrate for s over (-Inf, Inf) given a,b,c,d to obtain VUS
        (pnorm(a*s-b)*pnorm(-c*s+d)-p*pnorm(-c*s+d)-q*pnorm(a*s-b)+p*q)*dnorm(s)     
      }

    f1 <- function(s,a,b,c,d,q)
      {
        ##to obtain derivation of VUS w.r.t a
        (dnorm(a*s-b)*pnorm(-c*s+d)-q*dnorm(a*s-b))*dnorm(s)*s
      }


    f2 <- function(s,a,b,c,d,q)
      {
        ##to obtain derivation of VUS w.r.t b
        -(dnorm(a*s-b)*pnorm(-c*s+d)+q*dnorm(a*s-b))*dnorm(s)
      }

    f3 <- function(s,a,b,c,d,p)
      {
        ##to obtain derivation of VUS w.r.t c
        (-pnorm(a*s-b)*dnorm(-c*s+d)+p*dnorm(-c*s+d))*dnorm(s)*s
      }

    f4 <- function(s,a,b,c,d,p)
      {
        ##to obtain derivation of VUS w.r.t c
        (pnorm(a*s-b)*dnorm(-c*s+d)-p*dnorm(-c*s+d))*dnorm(s)
      }

    ###remove NAs
    x <- na.exclude(x)
    y <- na.exclude(y)
    z <- na.exclude(z)

    ####Sample size
    n.minus <- length(x)
    n0 <- length(y)
    n.plus <- length(z)
    
    ##normal sample means
    mu.minus <- mean(x,na.rm=TRUE)
    mu0 <- mean(y,na.rm=TRUE)
    mu.plus <- mean(z,na.rm=TRUE)

    if(!(mu.minus<=mu0 && mu0<=mu.plus)) stop("The means of x, y, z are not in increasing order!")

    ##normal sample SDs
    s.minus <- sd(x,na.rm=TRUE)
    s0 <- sd(y,na.rm=TRUE)
    s.plus <- sd(z,na.rm=TRUE)
    
    ###data summary
    dat.summary <- data.frame(n=c(n.minus,n0,n.plus),mu=c(mu.minus,mu0,mu.plus),sd=c(s.minus,s0,s.plus),row.names=c("D-","D0","D+"))
    
    ###reparametrize the normal sample means and SDs by a,b,c,d
    a <- s0/s.minus
    b <- (mu.minus-mu0)/s.minus
    c <- s0/s.plus
    d <- (mu.plus-mu0)/s.plus


    lower <- (qnorm(p)+b)/a#lower=-Inf when p=q=0
    upper <- (d-qnorm(q))/c
    
    VUS <- integrate(f0,a=a,b=b,c=c,d=d,p=p,q=q,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
  
   
    ##########Variance of partial VUS by equation (14) in Chengjie Xiong et al, 2006 statistics in medicine
    ####(1)variance of a,b,c,d
    a.var <- 0.5*a^2*(1/n0+1/n.minus)
    b.var <- 0.5*b^2/n.minus+a^2/n0+1/n.minus
    c.var <- 0.5*c^2*(1/n0+1/n.plus)
    d.var <- 0.5*d^2/n.plus+c^2/n0+1/n.plus

    ###(2) covariance between a,b,c,d
    
    ###The correct one
    #a.b.cov <- 0.5*a*b/(n.minus-1)
    #a.c.cov <- 0.5*a*c/(n0-1)
    #b.d.cov <- a*c/n0
    #c.d.cov <- 0.5*c*d/(n.plus-1)
    
    ###the asymptotic one
    a.b.cov <- 0.5*a*b/n.minus
    a.c.cov <- 0.5*a*c/n0
    b.d.cov <- a*c/n0
    c.d.cov <- 0.5*c*d/n.plus

    ###(3)derivations of VUS w.r.t a,b,c,d (denoted by V.a,V.b,V.c,V.d)
    V.a <- integrate(f1,a=a,b=b,c=c,d=d,q=q,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    V.b <- integrate(f2,a=a,b=b,c=c,d=d,q=q,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    V.c <- integrate(f3,a=a,b=b,c=c,d=d,p=p,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    V.d <- integrate(f4,a=a,b=b,c=c,d=d,p=p,lower=lower,upper=upper,subdivisions=subdivisions,...)$value

    ###(4)Variance    
    var0 <- V.a^2*a.var+V.b^2*b.var+V.c^2*c.var+V.d^2*d.var+2*V.a*V.b*a.b.cov+2*V.a*V.c*a.c.cov+2*V.b*V.d*b.d.cov+2*V.c*V.d*c.d.cov
    
    ###(5)95% CI
    z0 <- qnorm(alpha/2,lower=F)
    d0 <- z0*sqrt(var0)
    CI <- c(VUS-d0,VUS+d0)
    names(CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
    
    ###Sample size calculation
    
    #lam.minus <- n.minus/total.n
    #lam0 <- n0/total.n
    #lam.plus <- n.plus/total.n
        
    Mv <- 0.5*V.a^2*a^2*(1+lam0/lam.minus)+V.b^2*(a^2+0.5*b^2*lam0/lam.minus+lam0/lam.minus)+0.5*V.c^2*c^2*(1+lam0/lam.plus)+V.d^2*(c^2+0.5*d^2*lam0/lam.plus+lam0/lam.plus)+V.a*V.b*a*b*lam0/lam.minus+V.a*V.c*a*c+2*V.b*V.d*a*c+V.c*V.d*c*d*lam0/lam.plus

    z0 <- qnorm(typeIerror/2,lower=F)
    
    sampleSize <- ceiling(z0^2*Mv/margin^2)
    
    
    #return(list(dat=list(x=x,y=y,z=z),dat.summary=dat.summary, estimate=VUS,variance=var0,CI=CI,sampleSize=sampleSize,partialDeriv=c(a=a,b=b,c=c,d=d,V.a=V.a,V.b=V.b,V.c=V.c,V.d=V.d,a.var=a.var,b.var=b.var,c.var=c.var,d.var=d.var,a.b.cov=a.b.cov,a.c.cov=a.c.cov,b.d.cov=b.d.cov,c.d.cov=c.d.cov))
    return(list(dat=list(x=x,y=y,z=z),dat.summary=dat.summary, estimate=VUS,variance=var0,CI=CI,sampleSize=sampleSize,partialDeriv=data.frame(a=a,b=b,c=c,d=d,V.a=V.a,V.b=V.b,V.c=V.c,V.d=V.d)))
    
  }

