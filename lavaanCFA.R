library(rstan)
library(lavaan)

data("HolzingerSwineford1939")
dados <- HolzingerSwineford1939[, 7:15]

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9'

fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure = T, group='school',
           group.equal = c('loadings', 'intercepts'))

buildConstMatrices <- function(fit) {
  pars <- lavMatrixRepresentation(parTable(fit), add=T)
  constMat <- list()
  for (m in unique(pars$mat)){
    parTemp <- with(pars[pars$mat==m,], cbind(row, col, group, ustart))
    # Making sure we build the full matrix instead of lower triangular
    if (m %in% c('psi', 'theta')) {
      parTemp <- rbind(parTemp, parTemp[, c(2, 1, 3, 4)])
    }
    parM <- array(0, apply(parTemp[, 1:3], 2, max))
    for (idx in 1:(nrow(parTemp))) {
      sel <- parTemp[idx, 1:3]
      parM[sel[1], sel[2], sel[3]] <- parTemp[idx, 4]
    }
    outM <- reshape2::melt(parM)
    names(outM) <- c('row', 'col', 'group', 'value')
    outM <- outM[!is.na(outM$value), ]
    constMat[[m]] <- outM
  }
  constMat
}

buildEqualConst <- function(fit) {
  pars <- lavMatrixRepresentation(parTable(fit), add=T)
  equal <- pars[pars$op == '==', c('lhs', 'rhs')]
  tempEq <- data.frame(array(NA, c(dim(equal)[1], 7)))
  names(tempEq) <- c('group1', 'row1', 'col1', 'group2', 'row2', 'col2', 'mat')
  if (nrow(equal) != 0){
    for (c in 1:(dim(equal)[1])) {
      idx1 <- which(pars$plabel == equal[c, 1])
      idx2 <- which(pars$plabel == equal[c, 2])
      tempEq[c, ] <- c( pars[idx1, 'group'], pars[idx1, 'row'], pars[idx1, 'col'],
                        pars[idx2, 'group'], pars[idx2, 'row'], pars[idx2, 'col'], pars[idx1, 'mat'])
    }
  }
  
  
  matEq <- list()
  for (m in unique(tempEq[, 7])){
    matEq[[m]] <- tempEq[tempEq[, 7] == m, 1:6]
  }
  matEq
}

buildDataList <- function(fit, data, group) {
  const <- buildConstMatrices(fit)
  equal <- buildEqualConst(fit)
  pars <- lavMatrixRepresentation(parTable(fit), add = T)
  out <- list()
  # TODO Reorder data columns as needed by the model matrices
  # TODO Reorder data rows so each group is separated, scale as appropriate
  matNames <- names(const)
  matNames <- matNames[matNames != '']
  for (m in matNames) {
    out[[paste(m, 'N', sep='')]] <- nrow(const[[m]])
    out[[paste(m, 'EqN', sep='')]] <- ifelse(!is.null(equal[[m]]), nrow(equal[[m]]), 0)
    if (m %in% c('nu', 'alpha')){
      out[[paste(m, 'Par', sep='')]] <- const[[m]][, c('group', 'row')]
      if (!is.null(equal[[m]])){
        out[[paste(m, 'Equal', sep='')]] <- equal[[m]][, c(1, 2, 4, 5)]
      } else {
        out[[paste(m, 'Equal', sep='')]] <- array(0, c(0, 4))
      }
    } else {
      out[[paste(m, 'Par', sep='')]] <- const[[m]][, c('group', 'row', 'col')]
      if (!is.null(equal[[m]])){
        out[[paste(m, 'Equal', sep='')]] <- equal[[m]]
      } else {
        out[[paste(m, 'Equal', sep='')]] <- array(0, c(0, 6))
      }
    }
    out[[paste(m, 'Const', sep='')]] <- const[[m]][['value']]
  }
  out[['X']] <- data
  out[['N']] <- nrow(data)
  out[['K']] <- ncol(data)
  out[['G']] <- max(pars$group)
  out[['F']] <- max(out[['lambdaPar']][, 3])
  groupOut <- c(1)
  for (g in 2:length(group)){
    if (group[g] != group[g-1] & out[['G']] != 1){
      groupOut <- c(groupOut, g)
    }
  }
  out[['group']] <- c(groupOut, nrow(data))

  out[['sample_cov']] <- cov(data)
  out
}

initf <- function(lavaanFit) {
  inits <- lavInspect(lavaanFit, 'est')
  groups <- max(parTable(lavaanFit)$group)
  if (groups == 1){
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
  } else {
    out <- list()
    for (g in 1:groups){
      Psi_cor <- inits[[g]]$psi
      Psi_tau <- sqrt(diag(Psi_cor))
      Psi_cor <- t(chol(cov2cor(matrix(Psi_cor, nrow(Psi_cor), ncol(Psi_cor)))))
      out$Psi_cor <- abind::abind(out$Psi_cor, Psi_cor, along=3)
      if (g == groups){
        out$Psi_cor <- aperm(out$Psi_cor, c(3, 1, 2))
      }
      out$Psi_tau <- rbind(out$Psi_tau, Psi_tau+runif(1, -0.05, 0.05))
      
      
      Theta_cor <- inits[[g]]$theta
      Theta_tau <- sqrt(diag(Theta_cor))
      Theta_cor <- t(chol(cov2cor(matrix(Theta_cor, nrow(Theta_cor), ncol(Theta_cor)))))
      out$Theta_cor <- abind::abind(out$Theta_cor, Theta_cor, along=3)
      if (g == groups){
        out$Theta_cor <- aperm(out$Theta_cor, c(3, 1, 2))
      }
      out$Theta_tau <- rbind(out$Theta_tau, Theta_tau+runif(1, -0.05, 0.05))
      
      Alpha <- as.vector(inits[[g]]$alpha)
      out$Alpha_full <- rbind(out$Alpha_full, Alpha+runif(1, -0.1, 0.1))
      
      Nu <- as.vector(inits[[g]]$nu)
      out$Nu_full <- rbind(out$Nu_full, Nu+runif(1, -0.1, 0.1))
      
      Lambda <- inits[[g]]$lambda
      out$Lambda_full <- abind::abind(out$Lambda_full, Lambda+runif(1, -0.1, 0.1), along=3)
      if (g == groups){
        out$Lambda_full <- aperm(out$Lambda_full, c(3, 1, 2))
      }
      out$lambda_pos <- rbind(out$lambda_pos, rep(1, dim(Lambda)[2]))
    }
  }
  
  out$eta <- do.call(rbind, predict(fit))
  out
}

buildDataList(fit, dados, as.numeric(HolzingerSwineford1939$school))
stanFit <- stan('lavaanCFA.stan', data=buildDataList(fit, dados, as.numeric(HolzingerSwineford1939$school)),
                iter = 1000, warmup=250, chains=4, thin=3,
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                init=lapply(1:4, function(x) initf(fit)))
