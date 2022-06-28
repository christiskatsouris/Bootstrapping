# Applications of Bootstrap Method

Learning Objectives: In this teaching page we present some key applications of the bootstrap resampling method for statistical inference purposes in time series regression models. In particular, the use of resampling in statistical inference goes back to Tukey (1958) and Efron (1979) while Freedman (1981) extended this idea to regression models. Therefore, it is crucial to understand the main idea of this useful tool for obtaining robust estimates of model parameters and test statistics. There are of course modelling environments especially under certain econometric conditions that a standard Bootstrap approach is not suitable, but this is beyond the scope of this teaching page.

# I. Bootstrap for Linear Time Series Regression Models

The main idea behind bootstrapping is to repeatedly draw samples with replacement from the data, compute the statistic of interest and generate the sampling distribution of the statistic. Furthermore, bootstrapping is easily programmed both in R and Matlab. In particular, the build-in function 'bootstrap' in R can be employed for constructing bootstrapped test statistics.


## Example 1





## Assignment 1  


## References

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals. Statistical science, 11(3), 189-228.

Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap. Journal of the American Statistical association, 89(428), 1303-1313.


# II. Bootstrap for Conditional Heteroscedastic Time Series Regression Models

## Example 2

Simulate a data generating process that corresponds to a stationary GARCH model using the code below:

```R

install.packages("fGarch")
library("fGarch")

# Garch parameters
alpha <- 0.07
omega <- 0.5
beta  <- 0.6

# Garch specification
fgspec <- garchSpec( model=list(alpha=alpha, omega=omega, beta=beta) )
fgsim  <- garchSim( spec=fgspec, extended=TRUE, n=100 )
series <- fgsim

# Series extraction
garch.values <- series$garch
sigma.values <- series$sigma
eps.values   <- series$eps

# Plotting sequences
ts.plot(garch.values)
ts.plot(sigma.values)
ts.plot(eps.values)

est.garch <- garchFit(formula = ~ garch(1,1), data = series[1:n], include.mean = FALSE, trace = F)
est.par   <- est.garch@fit$par[1:3]
innov     <- series[1:n]/est.garch@sigma.t

ts.plot(innov)

```

## Assignment 2  

Using the bootstrap procedure below, implemented in R, consider an appropriate test statistic Tn of your choice for which you can employ as a structural-break detector when testing for the presence of parameter instability in a Garch (1,1) model.  Simulate B = 1,000 DGP based on sample sizes n = {250, 500} and obtain separately the empirical size (under the null hypothesis of no parameter instability) and the power function (under the alternative hypothesis) of your test statistic. Report the empirical size, empirical power and provide a detailed description of how you estimate the critical values of the test statistic. 

```R

set.seed(1234)

bootstrap.step <- function( innov = innov, B = B, n = n, n.s = n.s, est.par = est.par )
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
    
    x.boot    <- x.boot[(n.s+1):(n.s+n)]
    t.stat[i] <- T_hat(x.boot)
    
    #T_hat is the function for the test statistic (under the null hypothesis there is no volatility shift)
  }#end of bootstrap step
  
  return(t.stat)
  
}#End of function

```

## References

- Berkes, István, et al. "Sequential change-point detection in GARCH (p, q) models." Econometric theory 20.6 (2004): 1140-1167.

- Horváth, Lajos, Zhenya Liu, and Shanglin Lu. "Sequential monitoring of changes in dynamic linear models, applied to the US housing market." Econometric Theory (2021): 1-64.


# Further Reading

Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. Chapman and Hall, New York.

Godfrey, L. (2009). Bootstrap tests for regression models. Springer.
