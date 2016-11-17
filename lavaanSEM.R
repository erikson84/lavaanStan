library(rstan)
library(lavaan)
source('utils.R')

data("PoliticalDemocracy")
dados <- PoliticalDemocracy[, c(9:11, 1:8)]

model <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
   # residual covariances
     y1 ~~ y5
     y2 ~~ y4 + y6
     y3 ~~ y7
     y4 ~~ y8
     y6 ~~ y8
'
fit <- sem(model,
           data=PoliticalDemocracy, meanstructure = T, fixed.x=F)

stanFit <- stan('lavaanSEM.stan', data=buildDataList(fit, dados, rep(1, dim(dados)[1])),
                iter = 2000, warmup=500, chains=4, thin=3, control = list(adapt_delta=0.85),
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
