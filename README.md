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
