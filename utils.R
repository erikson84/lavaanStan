buildMatrices <- function(fit) {
  pars <- lavMatrixRepresentation(parTable(fit), add.attributes = T)
  mat <- list()
  for (m in attr(pars, 'mmNames')[[1]]){
    parTemp <- pars[pars$mat==m, c('row', 'col', 'group', 'ustart', 'free', 'label')]
    colnames(parTemp) <- c('row', 'col', 'group', 'value', 'free', 'label')
    parTemp[, 'free'] <- as.numeric(parTemp[, 'free'] > 0)
    if (m %in% c('theta', 'psi')){
      parTemp <- parTemp[order(parTemp[, 'group'],
                               parTemp[, 'row'] != parTemp[, 'col']), ]
      parTemp[parTemp[, 'free'] > 0 &
                parTemp[, 'row'] == parTemp[, 'col'], 'free'] <- 
        1:sum(parTemp[parTemp[, 'row'] == parTemp[, 'col'], 'free'])
      parTemp[parTemp[, 'free'] > 0 &
                parTemp[, 'row'] != parTemp[, 'col'], 'free'] <- 
        1:sum(parTemp[parTemp[, 'row'] != parTemp[, 'col'], 'free'])
    } else {
      parTemp[parTemp[, 'free'] > 0, 'free'] <- 1:sum(parTemp[, 'free'])
    }
    parTemp[, 'value'] <- ifelse(is.na(parTemp[, 'value']), 0, parTemp[, 'value'])
    for (l in 1:nrow(parTemp)){
      if (parTemp$label[l] != '') {
        parTemp$free[parTemp$label == parTemp$label[l]] <- 
          parTemp$free[parTemp$label == parTemp$label[l]][1]
      }
    }
    mat[[m]] <- parTemp[, -6]
  }
  mat
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
    # Making sure we keep cov. matrices symmetric
    if (m %in% c('psi', 'theta')) {
      mInv <- matEq[[m]][, c(1, 3, 2, 4, 6, 5)]
      names(mInv) <- names(matEq[[m]])
      matFull <- rbind(matEq[[m]], mInv)
      matEq[[m]] <- matFull[!duplicated(matFull), ]
    }
  }
  matEq
}

buildDataList <- function(fit, data, group=rep(1, nrow(data))) {
  const <- buildMatrices(fit)
  equal <- buildEqualConst(fit)
  pars <- lavMatrixRepresentation(parTable(fit), add.attributes = T)
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
      out[[paste(m, 'Free', sep='')]] <- const[[m]][, c('free')]
      if (!is.null(equal[[m]])){
        out[[paste(m, 'Equal', sep='')]] <- equal[[m]][, c(1, 2, 4, 5)]
      } else {
        out[[paste(m, 'Equal', sep='')]] <- array(0, c(0, 4))
      }
    } else {
      out[[paste(m, 'Par', sep='')]] <- const[[m]][, c('group', 'row', 'col')]
      out[[paste(m, 'Free', sep='')]] <- const[[m]][, c('free')]
      if (!is.null(equal[[m]])){
        out[[paste(m, 'Equal', sep='')]] <- equal[[m]]
      } else {
        out[[paste(m, 'Equal', sep='')]] <- array(0, c(0, 6))
      }
    }
    out[[paste(m, 'Const', sep='')]] <- const[[m]][, 'value']
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
    Lambda_full = as.array(pars$est[pars$mat == 'lambda' & pars$free > 0]),
    Nu_full = as.array(pars$est[pars$mat == 'nu' & pars$free > 0]),
    Alpha_full = as.array(pars$est[pars$mat == 'alpha' & pars$free > 0]),
    Theta_tau = as.array(sqrt(pars$est[pars$mat == 'theta' & pars$free > 0 & (pars$lhs == pars$rhs)])),
    Theta_cor = as.array(parsStd$est.std[
      apply(parsStd[, 1:3], 1, paste, collapse='') %in%
      apply(pars[pars$mat == 'theta' & pars$free > 0 & (pars$lhs != pars$rhs), 2:4], 1, paste, collapse='')
    ]),
    Psi_tau = as.array(sqrt(pars$est[pars$mat == 'psi' & pars$free > 0 & (pars$lhs == pars$rhs)])),
    Psi_cor = as.array(parsStd$est.std[
      apply(parsStd[, 1:3], 1, paste, collapse='') %in%
        apply(pars[pars$mat == 'psi' & pars$free > 0 & (pars$lhs != pars$rhs), 2:4], 1, paste, collapse='')
      ]),
    Beta_full = as.array(pars$est[pars$mat == 'beta' & pars$free > 0])
  )
  out
}
