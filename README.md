# Applications of Bootstrap Method

Applications of statistical theory and methods have positive spillover effects in the job market, especially in technology companies, financial institutions and other organizations. In particular, Applied data science has a broad scope of data science which includes researching new applications of data science and creating new methods or functions for faster retrieval of data and analysis. Applied data scientists have higher and deep technical knowledge of how data science and its methods work as compared to data scientists. Similarly, our role as econometricians and statisticians is to improve current understanding of existing methodologies and procedures.

In this teaching page we present some key applications of the bootstrap resampling method for statistical inference purposes in time series regression models. In particular, the use of resampling in statistical inference goes back to Tukey (1958) and Efron (1979) while Freedman (1981) extended this idea to regression models. Therefore, it is crucial to understand the main idea of this useful tool for obtaining robust estimates of model parameters and test statistics. There are of course modelling environments especially under certain econometric conditions that a standard Bootstrap approach is not suitable, but this is beyond the scope of this teaching page. A related introductory book chapter on Bootstrap and related methods, can be found in Chapter 21 of Davidson, R. and MacKinnon, J. G. (1993).

# [A]. Bootstrap for Linear Time Series Regression Models

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

## Example 3

```R

# Another simple bootstrapping example in R
# my.bootstrap function below implements bootstrapping with replacement

my.bootstrap <- function(x, nboot, theta, ...)
{
  data   <- matrix(sample(x, size = length(x) * nboot, replace = T), nrow = nboot)
  answer <- apply(data, 1, theta, ...)
  return( answer )
}

# loads data
data(faithful) 

# Extract the data from R
observed <- as.matrix( faithful$waiting )
hist(observed,prob=T,breaks=15) 

# Obtain 1000 bootstrap replicates of the sample mean
bootvals <- my.bootstrap( observed, 1000, mean ) 

# Histogram of the bootstrap replicates
hist(bootvals,prob=T) 

# Add the sample mean
abline(v = mean(observed),lwd=3) 

# Add estimate of bootstrap density
lines(density(bootvals,kernel="gaussian")) 

```

## Example 4 (Regression-based Bootstrap Application)

Consider the following regression-based bootstrap approach for resampling from the conditional distribution.

```R

install.packages("bootstrap")
library(bootstrap)

theta <- function(z, n1)
{
  answer <- mean( z[1:n1]) - mean( z[(n1+1):length(z)] )
  return( answer )
}

treated     <- c(94,197,16,38,99,141,23)
non.treated <- c(52,104,146,10,51,30,40,27,46)

bootvals <- bootstrap2(treated, non.treated, 1000, theta, n1 = 7)
hist(bootvals)

```

## Example 5

Consider the implementation of the following Residual and Wild Bootstrap Resampling Techniques.

```R

install.packages("evd")
install.packages("boot")
install.packages("resample")

library(evd)
library(boot)
library(resample)

Seed <- 14
set.seed(Seed)

n    <- 100
y    <- rnorm(n) #randomly generated response
x    <- rnorm(n) #randomly generated predictor

formula  <- y~x

ResidObj <- residual.boot(y~x, B=100, seed=Seed) # perform the residual bootstrap
WildObj  <- wild.boot( y~x, B=100, seed=Seed)    # perform the wild bootstrap

# residual bootstrap 95% CI for slope parameter (percentile method)
quantile(ResidObj$bootEstParam[,2], probs=c(.025, .975))

# bootstrap 95% CI for slope parameter (percentile method)
quantile(WildObj$bootEstParam[,2], probs=c(.025, .975))

```

## Remarks:

Notice that under the presence of nuisance parameters in regression models then resampling methods such as the bootstrap can be employed for estimation and inference purposes (such as in the case of threshold cointegration model e.g., either within a multivariate setting or in dynamic or non-dynamic panel data regression models). Although the mechanism for obtaining bootstrap random sequences remains the same as in the case of the standard bootstrap distribution approach, in regression models we also need to obtain the corresponding bootstrapped model estimates in order to obtain asymptotic approximations for nuisance parameters (e.g., such as the unknown threshold variable). A framework for estimation and testing in dynamic panel models with nonlinearities is proposed by Hansen, B. E. (1999).   

## Assignment 1  

Using the estimation and testing procedure proposed in the paper of Hansen, B. E. (1999) implement your own R code procedure with a suitable bootstrap reampling method for obtaining bootstrap model coefficient estimates, bootstrap estimation of the nuisance threshold variable as well as correspondoing confidence intervals for the Regime-Specific Threshold Regression model. Your estimations can be based on the 'invest' dataset (which includes the variables i,q,c,d) or a different suitable dataset of your choice which based on a related economic theory it can support the hypothesis of the presence of threshold effects.      

In particular, you can use the following part of the code for the estimation step of the model specification. Your assignment should include the functions for the bootstrap procedure, a detailed description of the econometric model and main intuition as well as indicative output from the simulation study such as graphical representation of the bootstrap distribution of the threshold variable as well as computational aspects such as executation time and the computational complexity based on the procedure proposed in the paper of Hansen, B. E. (1999).

```R

# Threshold Regression Model Estimation Step for Non-Dynamic Panel Data (see,  Hansen, B. E. (1999))

# function to obtain an estimate for the SSE measure
sse_calc <- function(y,x)
{
   e  <- y-x%*%qr.solve(x,y)
  out <- t(e)%*%e 
  return(out)
} 

thr_sse <- function(y,q,r)
{
        nq <- nrow(q)
        sse <- matrix(c(0),nq,1)
        for (qi in 1:nq)
         {
            if (r[1]==0)
            {
               rr <- q[qi] 
            }
            else 
            { rr <- rbind(r,q[qi])
            }
            rr <- as.matrix(sort(rr))
            xx <- cbind(xt,ct)
            
            for (j in 1:nrow(rr))
            {
                d <- (thresh < rr[j])
                xx <- cbind(xx,tr(cf*d))
            }
            sse[qi] <- sse_calc(y,xx)
        }
        sse
}

# Reference: Hansen, B. E. (1999). "Threshold effects in non-dynamic panels: Estimation, testing, and inference".

```

## References

- DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals. Statistical science, 11(3), 189-228.
- Efron, B. (1979). Computers and the theory of statistics: thinking the unthinkable. SIAM review, 21(4), 460-480.
- Freedman, D. A. (1981). Bootstrapping regression models. The Annals of Statistics, 9(6), 1218-1228.
- Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap. Journal of the American Statistical association, 89(428), 1303-1313.
- MacKinnon, J. G. (2002). Bootstrap inference in econometrics. Canadian Journal of Economics/Revue canadienne d'économique, 35(4), 615-645.
- Hansen, B. E. (1999). Threshold effects in non-dynamic panels: Estimation, testing, and inference. Journal of econometrics, 93(2), 345-368.

# [B]. Bootstrap for Conditional Heteroscedastic Time Series Regressions

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

Using the bootstrap procedure below, implemented in R, consider an appropriate test statistic Tn of your choice for which you can employ as a structural-break detector when testing for the presence of parameter instability in a Garch (1,1) model.  Simulate $B = 1,000$ DGPs based on sample sizes $n = {250, 500}$ and obtain separately the empirical size (under the null hypothesis of no parameter instability) and the power function (under the alternative hypothesis) of your test statistic. Report the empirical size, empirical power and provide a detailed description of how you estimate the critical values of the test statistic. 

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

- Horváth, L., Liu, Z., & Lu, S. (2022). Sequential monitoring of changes in dynamic linear models, applied to the US housing market. Econometric Theory, 38(2), 209-272. 
- Berkes, I., Gombay, E., Horváth, L., & Kokoszka, P. (2004). Sequential change-point detection in GARCH (p, q) models. Econometric theory, 20(6), 1140-1167.
- Corradi, V., & Iglesias, E. M. (2008). Bootstrap refinements for QML estimators of the GARCH (1, 1) parameters. Journal of Econometrics, 144(2), 500-510.
- Paparoditis, E., & Politis, D. N. (2009). Resampling and subsampling for financial time series. In Handbook of financial time series (pp. 983-999). Springer, Berlin, Heidelberg.

# [C]. Block Bootstrap Techniques

In this section, we demonstrate some useful examples related to Block Bootstrap techniques. 

## Example 3

Consider the following coding procedure in R which checks for the block-length for the implementation of the block bootstrap resampling approach. 

```R

set.seed(1234)

# Check block-length of subsamples
  if (is.null(sub_sample)) 
  {
    m <- round(n^(1 / 5) * n^(1 / k))
  } 
  else if (sub_sample >= n) 
  {
    stop("sub_sample must be less than series length")
  } 
  else
  {
    m <- round(sub_sample)
  }

```

Then, we can call the function 'tsbootstrap' in order to implement the block-bootstrap resampling method. 

```R
# Bootstrap variance of whole series
boot_temp <- tseries::tsbootstrap(series,
                                  statistic = stats::var,
                                  type = "block",
                                    nb = nb,
                                     b = l_star,
                                     m = bofb )
                                     
>  boot_temp

Call:
tseries::tsbootstrap(x = series, nb = nb, statistic = stats::var, 
    m = bofb, b = l_star, type = "block")

Resampled Statistic(s):
  original       bias std. error 
   1.28219   -0.01976    0.13802 

# Save updated variance of whole series
v_star <- mean(boot_temp$statistic)
        
```

# [D]. Stochastic Processes Simulation Examples

## [I]. Classical Brownian Bridge limiting processes

We begin by considering as an example, the simulation procedure to obtain critical values from a normalized Brownian Bridge which consists of functionals of Brownian motions. 

<p align="center">
  
<img src="https://github.com/christiskatsouris/Bootstrapping/blob/main/data/BM1.jpg" width="750"/>

</p>  

<p align="center">
  
<img src="https://github.com/christiskatsouris/Bootstrapping/blob/main/data/BM2.jpg" width="750"/>

</p>  

## [II]. Nonlinear Stochastic Processes (Advanced Topics)


Consider the Cox-Ingersoll-Ross (CIR) process  which is a nonlinear stochastic process given by the following expression 



```R

# Cox-Ingersoll-Ross (CIR) process 



```


## References

- Cox, J. C., Ingersoll Jr, J. E., & Ross, S. A. (1985). An intertemporal general equilibrium model of asset prices. Econometrica: Journal of the Econometric Society, 363-384.
- Cox, J. C., Ingersoll Jr, J. E., & Ross, S. A. (1985). A theory of the term structure of interest rates. Econometrica, 53(2), 385-407.
- Richard, S. F. (1978). An arbitrage model of the term structure of interest rates. Journal of Financial Economics, 6(1), 33-57.
- Wang, J. (1993). A model of intertemporal asset prices under asymmetric information. The Review of Economic Studies, 60(2), 249-282.


 
# Reading List

$\textbf{[1]}$ Davidson, R., & MacKinnon, J. G. (1993). Estimation and inference in econometrics (Vol. 63). New York: Oxford.

$\textbf{[2]}$ Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. Chapman and Hall, New York.

$\textbf{[3]}$ Godfrey, L. (2009). Bootstrap tests for regression models. Springer.


# Disclaimer

The author (Christis G. Katsouris) declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved.

Any errors or omissions are the responsibility of the author.

# Acknowledgments

The author (Christis G. Katsouris) greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

Furthermore, the author greatfully acknowledges financial support from Graduate Teaching Assistantships at the School of Economic, Social and Political Sciences of the University of Southampton as well as funding from Research Grants of various interdisciplinary Centers of research excellence based at the University of Cyprus (UCY) as well as at University College London (UCL). Furthermore, the author gratefully acknowledges financial support from the Vice-Chancellor's PhD Scholarship of the University of Southampton, for the duration of the academic years 2018 to 2021.

# Historical Background

$\textbf{David Freedman}$ (5 March 1938 – 17 October 2008) was Professor of Statistics at the University of California, Berkeley. He was a distinguished mathematical statistician whose wide-ranging research included the analysis of martingale inequalities, Markov processes, de Finetti's theorem, consistency of Bayes estimators, sampling, the bootstrap, and procedures for testing and evaluating models. He published extensively on methods for causal inference and the behavior of standard statistical models under non-standard conditions – for example, how regression models behave when fitted to data from randomized experiments. Freedman also wrote widely on the application—and misapplication—of statistics in the social sciences, including epidemiology, public policy, and law (Source: Wikipedia).
