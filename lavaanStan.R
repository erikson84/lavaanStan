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
    group <- as.numeric(factor(data[, lav.options$group]))
    data <- subset(data, select=-which(names(data) == lav.options$group))
  } else {
    group <- rep(1, nrow(data))
  }
  
  if ('std.ov' %in% names(lav.options)) {
    if (lav.options[['std.ov']])
      data <- scale(data)
  }

  stanFit <- stan(fit.stan, data=buildDataList(fit, data, group),
                  iter = iter, warmup=warmup, chains=chains,
                  init=lapply(1:chains, function(x) initf(fit)))
  attr(stanFit, 'data') <- data
  attr(stanFit, 'group') <- group
  stanFit
}

modelMatrix <- function(fit, what) {
  sum.func <- switch(what,
                     est=mean,
                     se=sd,
                     stop("Only 'est' (posterior means) and 'se' (posterior standard deviations) options are available."))
  pars <- extract(fit)
  out <- list()
  for (g in 1:length(unique(attr(fit, 'group')))) {
    out[[g]] <- list()
    for (par in names(pars)[grepl('^[A-Z]+[a-z]+$', names(pars))]) {
      if (length(dim(pars[[par]])) == 4) {
        out[[g]][[par]] <- apply(pars[[par]][,g,,], c(2, 3), sum.func)
      } else if (length(dim(pars[[par]])) == 3) {
        out[[g]][[par]] <- as.matrix(apply(pars[[par]][,g,], 2, sum.func))
      } else {
        out[[g]][[par]] <- sum.func(pars[[par]])
      }
      
    }
  }
  out
}

