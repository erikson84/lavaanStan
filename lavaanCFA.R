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

fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure = T, group='school')

stanFit <- stan('lavaanCFA.stan', data=buildDataList(fit, dados, as.numeric(HolzingerSwineford1939$school)),
                iter = 2000, warmup=1000, chains=4,
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
