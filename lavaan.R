library(rstan)
library(lavaan)

data("HolzingerSwineford1939")
dados <- HolzingerSwineford1939[, 7:15]

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9'

fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure = T, group='school', group.equal = 'loadings')
matrices <- lavInspect(fit, 'est')

lambda <- matrices$lambda
lambda <- ifelse(lambda != 0 & lambda != 1, 357, lambda)
lambda[2,1] <- 'est'

psi <- matrices$psi
psi <- ifelse(psi != 0 & psi != 1, 357, psi)

theta <- matrices$theta
theta <- ifelse(theta != 0 & theta != 1, 357, theta)

nu <- matrices$nu
nu <- as.vector(ifelse(abs(nu) > 1e-10 & nu != 1, 357, round(nu)))
if (all(nu == 0)){
  dados <- scale(dados)
}

alpha <- matrices$alpha
alpha <- as.vector(ifelse(alpha != 0 & alpha != 1, 357, alpha))

dataList = list(X=dados, N=nrow(dados), K=ncol(dados), F=ncol(psi), sample_cov=cov(dados),
            Lambda_const=lambda, Psi_const=psi, Theta_const=theta, Nu_const=nu, Alpha_const=alpha)

initf <- function(lavaanFit) {
  inits <- lavInspect(lavaanFit, 'est')
  
  Psi_cor <- inits$psi
  Psi_tau <- sqrt(diag(Psi_cor))
  Psi_cor <- t(chol(cov2cor(matrix(Psi_cor, nrow(Psi_cor), ncol(Psi_cor)))))
  
  Theta_cor <- inits$theta
  Theta_tau <- sqrt(diag(Theta_cor))
  Theta_cor <- t(chol(cov2cor(matrix(Theta_cor, nrow(Theta_cor), ncol(Theta_cor)))))
  
  Alpha <- as.vector(inits$alpha)
  Nu <- as.vector(inits$nu)
  Lambda <- inits$lambda
  lambda_pos <- rep(1, dim(Lambda)[2])
  eta <- predict(lavaanFit)
  
  
  list(Lambda_full=Lambda, Nu_full=Nu, Alpha_full=Alpha,
       Theta_cor=Theta_cor, Theta_tau=Theta_tau,
       Psi_cor=Psi_cor, Psi_tau=Psi_tau, 
       eta=eta, lambda_pos=lambda_pos)
}


stanFit <- stan('lavaanCFA.stan', data=dataList, iter = 800,
                warmup=300, chains=4,
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
