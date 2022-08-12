# Applications of Bootstrap Method

Learning Objectives: In this teaching page we present some key applications of the bootstrap resampling method for statistical inference purposes in time series regression models. In particular, the use of resampling in statistical inference goes back to Tukey (1958) and Efron (1979) while Freedman (1981) extended this idea to regression models. Therefore, it is crucial to understand the main idea of this useful tool for obtaining robust estimates of model parameters and test statistics. There are of course modelling environments especially under certain econometric conditions that a standard Bootstrap approach is not suitable, but this is beyond the scope of this teaching page.

# I. Bootstrap for Linear Time Series Regression Models

The main idea behind bootstrapping is to repeatedly draw samples with replacement from the data, compute the statistic of interest and generate the sampling distribution of the statistic. Furthermore, bootstrapping is easily programmed both in R and Matlab. In particular, the build-in function 'bootstrap' in R can be employed for constructing bootstrapped test statistics.


## Example 1

```R

# install the package Bootstrap
install.package("bootstrap")
library(bootstrap)

# to generate B=1000 bootstrap samples using the boostrap function we use

# number of boostrap
nboot = 1000

# call the bootstrap function
boot.mean = bootstrap( 1: dim(dat)[1], nboot, mean.diff, dat )

```

## Example 2

```R

# install the package coin
install.packages("coin")
library(coin)

tstat <- function(data)
{
  x    <- data[1:nEng, 2]
  y    <- data[(nEng+1):nrow(data), 2]
  tobj <- t.test(x, y)
  t    <- tobj$statistic 
  return(t)
}

bootstat <- function(data, indices) 
{
  d  <- data[indices,] # allows boot to select sample
  t  <- tstat(d)
  return(t)
}

b   <- boot(data=mice2, statistic=bootstat, R=nBoot)
p   <- length(which(b$t > b$t0)) / nBoot

```


## Assignment 1  


## References

- DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals. Statistical science, 11(3), 189-228.

- Efron, B. (1979). Computers and the theory of statistics: thinking the unthinkable. SIAM review, 21(4), 460-480.

- Freedman, D. A. (1981). Bootstrapping regression models. The Annals of Statistics, 9(6), 1218-1228.

- Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap. Journal of the American Statistical association, 89(428), 1303-1313.

- MacKinnon, J. G. (2002). Bootstrap inference in econometrics. Canadian Journal of Economics/Revue canadienne d'économique, 35(4), 615-645.

# II. Bootstrap for Conditional Heteroscedastic Time Series Regressions

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

Using the bootstrap procedure below, implemented in R, consider an appropriate test statistic Tn of your choice for which you can employ as a structural-break detector when testing for the presence of parameter instability in a Garch (1,1) model.  Simulate B = 1,000 DGP based on sample sizes $n = {250, 500}$ and obtain separately the empirical size (under the null hypothesis of no parameter instability) and the power function (under the alternative hypothesis) of your test statistic. Report the empirical size, empirical power and provide a detailed description of how you estimate the critical values of the test statistic. 

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

- Paparoditis, Efstathios, and Dimitris N. Politis. "Resampling and subsampling for financial time series." Handbook of financial time series. Springer, Berlin, Heidelberg, 2009. 983-999.
 
# III. Stochastic Processes Simulation Examples

## [A.] Classical Brownian Bridge limiting processes

We begin by considering as an example, the simulation procedure to obtain critical values from a normalized Brownian Bridge which consists of functionals of Brownian motions. 



## [B.] Nonlinear Stochastic Processes (Advanced Topics)


Consider the Cox-Ingersoll-Ross (CIR) process  which is a nonlinear stochastic process given by the following expression 



```R

# Cox-Ingersoll-Ross (CIR) process 



```


## References

- Cox, J. C., Ingersoll Jr, J. E., & Ross, S. A. (1985). An intertemporal general equilibrium model of asset prices. Econometrica: Journal of the Econometric Society, 363-384.

- Wang, J. (1993). A model of intertemporal asset prices under asymmetric information. The Review of Economic Studies, 60(2), 249-282.

# Reading List

$\textbf{[1]}$ Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. Chapman and Hall, New York.

$\textbf{[2]}$ Godfrey, L. (2009). Bootstrap tests for regression models. Springer.


# Disclaimer

The author (Christis G. Katsouris) declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved.

Any errors or omissions are the responsibility of the author.

# Acknowledgments

The author greatfully acknowledges financial support from Graduate Teaching Assistantships at the School of Economic, Social and Political Sciences of the University of Southampton as well as funding from Research Grants of various interdisciplinary Centers of research excellence based at the University of Cyprus (UCY) as well as at University College London (UCL). Furthermore, the author gratefully acknowledges financial support from the Vice-Chancellor's PhD Scholarship of the University of Southampton, for the duration of the academic years 2018 to 2021.

# Historical Background

$\textbf{David Freedman}$ (5 March 1938 – 17 October 2008) was Professor of Statistics at the University of California, Berkeley. He was a distinguished mathematical statistician whose wide-ranging research included the analysis of martingale inequalities, Markov processes, de Finetti's theorem, consistency of Bayes estimators, sampling, the bootstrap, and procedures for testing and evaluating models. He published extensively on methods for causal inference and the behavior of standard statistical models under non-standard conditions – for example, how regression models behave when fitted to data from randomized experiments. Freedman also wrote widely on the application—and misapplication—of statistics in the social sciences, including epidemiology, public policy, and law (Source: Wikipedia).
