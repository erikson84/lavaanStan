library(rstan)
library(lavaan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('utils.R')

data("PoliticalDemocracy")
dados <- PoliticalDemocracy
dados$x4 <- rbinom(75, 1, 0.5)
model <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60 + x4
   # residual covariances
     y1 ~~ y5
     y2 ~~ y4 + y6
     y3 ~~ y7
     y4 ~~ y8
     y6 ~~ y8
'
fit <- sem(model,
           data=dados, meanstructure = T, fixed.x=T)

stanFit <- stan('lavaanSEM.stan', data=buildDataList(fit, dados, rep(1, dim(dados)[1])),
                iter = 1000, warmup=500, chains=4,
                init=lapply(1:4, function(x) initf(fit)))
