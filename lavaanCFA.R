library(rstan)
library(lavaan)
source('utils.R')
data("HolzingerSwineford1939")
dados <- HolzingerSwineford1939[, 7:15]

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9
'

fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure = T)

stanFit <- stan('lavaanCFA.stan', data=buildDataList(fit, dados, as.numeric(HolzingerSwineford1939$school)),
                iter = 1000, warmup=500, chains=4, thin=2, control = list(adapt_delta=0.85),
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
