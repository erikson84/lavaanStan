require(lavaan)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('utils.R')

lavaanStan <- function(
  ..., # Lavaan arguments
  model.type, # Type os model: SEM or CFA
  iter = 1000, # Number of MCMC iterations
  warmup=500,  # Number of warmup iterations
  chains=4     # Number of chains
) {
  fit.func <- switch(model.type,
                    cfa = lavaan::cfa,
                    sem = lavaan::sem,
                    stop('model.type must be either sem or cfa'))
  fit.stan <- switch(model.type,
                     cfa = 'lavaanCFA.stan',
                     sem = 'lavaanSEM.stan',
                     stop('model.type must be either sem or cfa'))
  lav.options <- list(...)
  fit <- fit.func(..., meanstructure=T)
  
  data <- lav.options$data
  if ('group' %in% names(lav.options)){
    group <- as.numeric(factor(data[, lav.options['group']]))
    data <- subset(data, select=-c(lav.options['group']))
  } else {
    group <- rep(1, nrow(data))
  }

  stanFit <- stan(fit.stan, data=buildDataList(fit, data, group),
                  iter = iter, warmup=warmup, chains=chains,
                  init=lapply(1:chains, function(x) initf(fit)))
  stanFit
}