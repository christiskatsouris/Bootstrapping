# Bootstrap Applications

Learning Objectives: In this teaching page we present some key applications of the bootstrap resampling methodology for statistical inference purposes in time series regression models. 

## Bootstrap for Non-Linear Models

Simulate a data generating process that corresponds to a stationary GARCH model using the code below:

```R

# Garch parameters
alpha <- 0.07
omega <- 0.5
beta  <- 0.6

# Garch specification
fgspec <- garchSpec( model=list(alpha=alpha, omega=omega, beta=beta) )
fgsim  <- garchSim( spec=fgspec, extended=TRUE, n=100 )
series <- fgsim

garch.values<-series$garch
sigma.values<-series$sigma
eps.values<-series$eps

ts.plot(garch.values)
ts.plot(sigma.values)
ts.plot(eps.values)

est.garch <- garchFit(formula = ~ garch(1,1), data = series[1:n], include.mean = FALSE, trace = F)
est.par   <- est.garch@fit$par[1:3]
innov     <- series[1:n]/est.garch@sigma.t

ts.plot(innov)

```

### Exercise 1  

Using the bootstrap procedure below, implemented in R, consider an appropriate test statistic Tn of your choice for which you can employ as a structural-break detector when testing for the presence of parameter instability in a Garch (1,1) model.  Simulate B = 1,000 DGP based on sample sizes n = {250, 500} and obtain separately the empirical size (under the null hypothesis of no parameter instability) and the power function (under the alternative hypothesis) of your test statistic. Report the empirical size, empirical power and provide a detailed description of how you estimate the critical values of the test statistic. 

```R

method2 <- function(innov=innov,B=B,n=n, n.s=n.s, est.par=est.par)
{#begin of function
  
  t.stat <- rep(NA, times = B)
  # Store the test statistics calculated from the B time bootstrap
  
  for (i in 1:B)
  {#Begin of bootstrap step
    
    # Bootstrap Innovation
    e.tilde <- sample(innov,size=(n+n.s), replace = T)
    # Recursive Bootstrap Sample
    # Initial Values: s.boot=unconditional variance estimate
    
    s.boot <- rep(NA, times = n+n.s)             # Bootstrap sigma
    x.boot <- rep(NA, times = n+n.s)             # Bootstrap x
    
    if ( (est.par[2] + est.par[3]) < 1 )
      s.boot[1] <- est.par[1] / ( 1 - est.par[2] - est.par[3] )
    if ( (est.par[2] + est.par[3]) >= 1 )
      s.boot[1] <- var(innov)
   
   x.boot[1] <- sqrt(s.boot[1])*e.tilde[1]
   for (j in 2:(n+n.s))
    {
      s.boot[j] <- est.par[1] + est.par[2]*(x.boot[j-1])^2 + est.par[3]*s.boot[j-1]
      x.boot[j] <- sqrt(s.boot[j])*e.tilde[j]
    }
    
    x.boot <- x.boot[(n.s+1):(n.s+n)]
    t.stat[i] <- T_hat(x.boot)
    
    #T_hat is the function for the test statistic for no volatility shift
  }#end of bootstrap step
  
  return(t.stat)
  
}#End of function

```
