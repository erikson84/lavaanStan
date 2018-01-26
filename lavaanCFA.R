library(rstan)
library(lavaan)
source('utils.R')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data("HolzingerSwineford1939")
dados <- HolzingerSwineford1939[, 7:15]

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9
'

fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure = T, group='school',
           group.equal=c('loadings', 'intercepts'))
#summary(fit)

stanFit <- stan('lavaanCFA.stan', data=buildDataList(fit, scale(dados), as.numeric(HolzingerSwineford1939$school)),
                iter = 1000, warmup=500, chains=4,
                init=lapply(1:4, function(x) initf(fit)))

