buildMatrices <- function(fit) {
  pars <- lavMatrixRepresentation(parTable(fit), add.attributes = T)
  mat <- list()
  for (m in attr(pars, 'mmNames')[[1]]){
    parTemp <- as.matrix(pars[pars$mat==m, c('row', 'col', 'group', 'ustart', 'free')])
    labs <- pars[pars$mat==m, 'label']
    colnames(parTemp) <- c('row', 'col', 'group', 'value', 'free')
    parTemp[, 'free'] <- as.numeric(parTemp[, 'free'] > 0)
    
    if (sum(parTemp[, 'free']) > 0){
      if (m %in% c('theta', 'psi')){
        parTemp <- parTemp[order(parTemp[, 'group'],
                                 parTemp[, 'row'] != parTemp[, 'col']), ]
        parTemp[parTemp[, 'free'] > 0 &
                  (parTemp[, 'row'] == parTemp[, 'col']), 'free'] <- 
          1:sum(parTemp[parTemp[, 'row'] == parTemp[, 'col'], 'free'])
        parTemp[parTemp[, 'free'] > 0 &
                  (parTemp[, 'row'] != parTemp[, 'col']), 'free'] <- 
          1:sum(parTemp[parTemp[, 'row'] != parTemp[, 'col'], 'free'])
      } else {
        parTemp[parTemp[, 'free'] > 0, 'free'] <- 1:sum(parTemp[, 'free'])
      }
    }
    
    parTemp[, 'value'] <- ifelse(is.na(parTemp[, 'value']), 0, parTemp[, 'value'])
    for (l in 1:nrow(parTemp)){
      if (labs[l] != '') {
        parTemp[labs == labs[l], 'free'] <- 
          parTemp[labs == labs[l], 'free'][1]
      }
    }
    mat[[m]] <- parTemp
  }
  mat
}


buildDataList <- function(fit, data, group=rep(1, nrow(data))) {
  const <- buildMatrices(fit)
  pars <- lavMatrixRepresentation(parTable(fit), add.attributes = T)
  out <- list()
  # TODO Reorder data columns as needed by the model matrices
  # TODO Reorder data rows so each group is separated, scale as appropriate
  matNames <- names(const)
  matNames <- matNames[matNames != '']
  for (m in matNames) {
    out[[paste0(m, 'N')]] <- nrow(const[[m]])
    if (m %in% c('nu', 'alpha')){
      out[[paste0(m, 'Par')]] <- const[[m]][, c('group', 'row')]
      out[[paste0(m, 'Free')]] <- const[[m]][, c('free')]
    } else {
      out[[paste0(m, 'Par')]] <- const[[m]][, c('group', 'row', 'col')]
      out[[paste0(m, 'Free')]] <- const[[m]][, c('free')]
    }
    out[[paste0(m, 'Const')]] <- const[[m]][, 'value']
    if (m %in% c('theta', 'psi')) {
      out[[paste0(m, 'DiagN')]] <- length(which(const[[m]][, 'row'] == const[[m]][, 'col']))
      out[[paste0(m, 'OffDiagN')]] <- length(which(const[[m]][, 'row'] != const[[m]][, 'col']))
      out[[paste0(m, 'Diag')]] <- as.array(which(const[[m]][, 'row'] == const[[m]][, 'col']))
      out[[paste0(m, 'OffDiag')]] <- as.array(which(const[[m]][, 'row'] != const[[m]][, 'col']))
    }
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
  out[['group']] <- c(groupOut, nrow(data)+1)
  samp_cov <- list()
  for (g in 1:out$G){
    if (out$G == 1){
      subG <- 1:out$N
    } else { subG <- group==g }
    samp_cov[[g]] <- cov(data[subG, ])
  }
  samp_cov <- do.call(abind::abind, c(samp_cov,along=3))
  out[['sample_cov']] <- aperm(samp_cov, c(3, 1, 2))
  out
}

initf <- function(fit) {
  pars <- lavMatrixRepresentation(parTable(fit), add.attributes = T)
  parsStd <- standardizedSolution(fit)
  out <- list(
    Lambda_full = as.array(unique(round(pars$est[pars$mat == 'lambda' & pars$free > 0], digits=4))),
    Nu_full = as.array(unique(round(pars$est[pars$mat == 'nu' & pars$free > 0], digits=4))),
    Alpha_full = as.array(unique(round(pars$est[pars$mat == 'alpha' & pars$free > 0], digits=4))),
    Theta_tau = as.array(unique(round(sqrt(pars$est[pars$mat == 'theta' & pars$free > 0 & (pars$lhs == pars$rhs)]), digits=4))),
    Theta_cor = as.array(unique(round(parsStd$est.std[
      apply(parsStd[, 1:3], 1, paste, collapse='') %in%
      apply(pars[pars$mat == 'theta' & pars$free > 0 & (pars$lhs != pars$rhs), 2:4], 1, paste, collapse='')
    ]), digits=4)),
    Psi_tau = as.array(unique(round(sqrt(pars$est[pars$mat == 'psi' & pars$free > 0 & (pars$lhs == pars$rhs)]), digits=4))),
    Psi_cor = as.array(unique(round(parsStd$est.std[
      apply(parsStd[, 1:3], 1, paste, collapse='') %in%
        apply(pars[pars$mat == 'psi' & pars$free > 0 & (pars$lhs != pars$rhs), 2:4], 1, paste, collapse='')
      ], digits=4))),
    Beta_full = as.array(unique(round(pars$est[pars$mat == 'beta' & pars$free > 0], digits=4)))
  )
  out
}
