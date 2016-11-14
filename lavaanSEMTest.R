library(rstan)
library(lavaan)

data("PoliticalDemocracy")
dados <- PoliticalDemocracy[, c(9:11, 1:8)]

modelSim <- '
   # latent variables
     ind60 =~ 1*x1 + 0.9*x2 + 0.8*x3
     dem60 =~ 1*y1 + 0.9*y2 + 0.8*y3 + 0.7*y4
     dem65 =~ 1*y5 + 0.9*y6 + 0.8*y7 + 0.7*y8
   # regressions
     dem60 ~ 0.8*ind60
     dem65 ~ 0.8*ind60 + 0.8*dem60
   # residual covariances
    x1 ~~ 0.5*x1 
    x2 ~~ 0.3*x2
    x3 ~~ 0.8*x3
    y1 ~~ 1*y1
    y2 ~~ 0.5*y2
    y3 ~~ 1.5*y3
    y4 ~~ 0.5*y4
    y5 ~~ 1*y5
    y6 ~~ 0.8*y6
    y7 ~~ 0.4*y7
    y8 ~~ 0.2*y8
'
model <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
'
dados <- simulateData(modelSim)
fit <- sem(model,
           data=dados, meanstructure = T, fixed.x=F, conditional.x = F)

matrices <- lavInspect(fit, 'est')

lambda <- matrices$lambda
lambda <- ifelse(lambda != 0 & lambda != 1, 357, lambda)

phi <- matrices$psi
phi <- ifelse(phi != 0 & phi != 1, 357, phi)

psi <- matrices$theta
psi <- ifelse(psi != 0 & psi != 1, 357, psi)

nu <- matrices$nu
nu <- as.vector(ifelse(abs(nu) > 1e-10 & nu != 1, 357, round(nu)))
if (all(nu == 0)){
  dados <- scale(dados)
}

alpha <- matrices$alpha
alpha <- as.vector(ifelse(alpha != 0 & alpha != 1, 357, alpha))

beta <- matrices$beta
beta <- ifelse(beta != 0 & beta != 1, 357, beta)

dataList = list(X=dados, N=nrow(dados), K=ncol(dados), F=ncol(phi), sample_cov=cov(dados), Beta_const=beta,
                Lambda_const=lambda, Phi_const=phi, Psi_const=psi, Nu_const=nu, Alpha_const=alpha)

initf <- function(lavaanFit) {
  inits <- lavInspect(lavaanFit, 'est')
  
  Phi_cor <- inits$psi
  Phi_tau <- pcauchy(sqrt(diag(Phi_cor)), 0, 3)
  Phi_cor <- t(chol(cov2cor(matrix(Phi_cor, nrow(Phi_cor), ncol(Phi_cor)))))
  
  Psi_cor <- inits$theta
  Psi_tau <- pcauchy(sqrt(diag(Psi_cor)), 0, 3)
  Psi_cor <- t(chol(cov2cor(matrix(Psi_cor, nrow(Psi_cor), ncol(Psi_cor)))))
  
  Alpha <- as.vector(inits$alpha)
  Nu <- as.vector(inits$nu)
  Lambda <- inits$lambda
  lambda_pos <- rep(1, dim(Lambda)[2])
  Beta <- inits$beta
  phi <- scale(predict(lavaanFit))
  
  
  list(Lambda_full=Lambda, Nu_full=Nu, Alpha_full=Alpha,
       Psi_cor=Psi_cor, Psi_unif=Psi_tau,
       Phi_cor=Phi_cor, Phi_unif=Phi_tau, 
       phi=phi, lambda_pos=lambda_pos, 
       Beta_full=Beta)
}


stanFitMatt <- stan('lavaanSemMatt.stan', data=dataList, iter = 1000, warmup=500, chains=4,
                control = list(adapt_delta=0.8),
                #pars=c('Nu', 'Lambda', 'Beta', 'Alpha', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
