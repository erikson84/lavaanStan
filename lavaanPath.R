library(rstan)
library(lavaan)
source('utils.R')

set.seed(1234)
X <- rnorm(100)
M <- 0.5*X + rnorm(100)
Y <- 0.7*M + rnorm(100)
Data <- data.frame(Y = Y, M = M, X = X)
model <- ' # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
fit <- sem(model, data = Data, meanstructure = T)
summary(fit)
stanFit <- stan('lavaanPath.stan', data=buildDataList(fit, Data),
                iter = 1000, warmup=500, chains=4, thin=2, control = list(adapt_delta=0.85),
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
