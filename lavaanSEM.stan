functions {
  matrix cov(matrix X) {
    int N;
    int K;
    vector[rows(X)] ones;
    vector[cols(X)] means;
    
    N = rows(X);
    K = cols(X);

    for (k in 1:K){
      means[k] = mean(X[, k]);
    }
    
    for (n in 1:N){
      ones[n] = 1.0;
    }
    return ((X - ones * means')' * (X - ones * means')) / (N - 1);
  }
  
  real Fml(matrix pop_cov, matrix samp_cov){
    int K;
    K = rows(pop_cov);
    return log_determinant(pop_cov) + trace(samp_cov * inverse(pop_cov)) - log_determinant(samp_cov) - K;
  }
}

data {
  // Basic dimension for the model
  int<lower=1> G;
  int<lower=1> N;
  int<lower=3> K;
  int<lower=1> F;
  // The data matrix as vectors
  vector[K] X[N];
  // Vector indicating were each group starts
  int<lower=1> group[G+1];
  
  // Constrained parameters as sparse matrices
  // Lambda
  int lambdaN;
  int lambdaEqN;
  int lambdaPar[lambdaN, 3];
  vector[lambdaN] lambdaConst;
  // Equality constraints
  int lambdaEqual[lambdaEqN, 6];
  
  // Theta
  int thetaN;
  int thetaEqN;
  int thetaPar[thetaN, 3];
  vector[thetaN] thetaConst;
  // Equality constraints
  int thetaEqual[thetaEqN, 6];
  
  // Beta
  int betaN;
  int betaEqN;
  int betaPar[betaN, 3];
  vector[betaN] betaConst;
  // Equality constraints
  int betaEqual[betaEqN, 6];
  
  // Psi
  int psiN;
  int psiEqN;
  int psiPar[psiN, 3];
  vector[psiN] psiConst;
  // Equality constraints
  int psiEqual[psiEqN, 6];
  
  
  // Nu
  int nuN;
  int nuEqN;
  int nuPar[nuN, 2];
  vector[nuN] nuConst;
  // Equality constraints
  int nuEqual[nuEqN, 4];
  
  // Alpha
  int alphaN;
  int alphaEqN;
  int alphaPar[alphaN, 2];
  vector[alphaN] alphaConst;
  // Equality constraints
  int alphaEqual[alphaEqN, 4];
  
  // Sample covariance matrix (for PPC)
  matrix[K, K] sample_cov[G];
}

transformed data {
  vector[F] one;
  matrix[F, F] I;
  for (i in 1:F) one[i] = 1.0;
  I = diag_matrix(one);
  
}

parameters {
  matrix[K, F] Lambda_full[G];
  //vector[F] eta[N];

  matrix[F, F] Beta_full[G];

  cholesky_factor_corr[K] Theta_cor[G];
  vector<lower=0>[K] Theta_tau[G];
  
  cholesky_factor_corr[F] Psi_cor[G];
  vector<lower=0>[F] Psi_tau[G];
  
  vector[K] Nu_full[G];
  vector[F] Alpha_full[G];
}

transformed parameters {
  matrix[K, F] Lambda[G];
  matrix[F, F] Beta[G];
  cov_matrix[K] Theta[G];
  cov_matrix[F] Psi[G];
  vector[K] Nu[G];
  vector[F] Alpha[G];
  
  
  for (g in 1:G){
    Lambda[g] = Lambda_full[g];
    Beta[g] = Beta_full[g];
    Theta[g] = quad_form_diag(Theta_cor[g] * Theta_cor[g]', Theta_tau[g]);
    Psi[g] = quad_form_diag(Psi_cor[g] * Psi_cor[g]', Psi_tau[g]);
    Nu[g] = Nu_full[g];
    Alpha[g] = Alpha_full[g];
  }
  
  // Set fixed value constraints
  
  if (lambdaN > 0) for (i in 1:lambdaN){
    Lambda[lambdaPar[i, 1], lambdaPar[i, 2], lambdaPar[i, 3]] = lambdaConst[i];
  }
  
  if (betaN > 0) for (i in 1:betaN){
    Beta[betaPar[i, 1], betaPar[i, 2], betaPar[i, 3]] = betaConst[i];
  }
  
  if (thetaN > 0) for (i in 1:thetaN){
    Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 3]] = thetaConst[i];
  }
  if (psiN > 0) for (i in 1:psiN){
    Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 3]] = psiConst[i];
  }
  
  if (nuN > 0) for (i in 1:nuN){
    Nu[nuPar[i, 1], nuPar[i, 2]] = nuConst[i];
  }
  if (alphaN > 0) for (i in 1:alphaN){
    Alpha[alphaPar[i, 1], alphaPar[i, 2]] = alphaConst[i];
  }
  
  // Set equality constraints
  
  if (lambdaEqN > 0) for (i in 1:lambdaEqN){
    Lambda[lambdaEqual[i, 1], lambdaEqual[i, 2], lambdaEqual[i, 3]] = Lambda[lambdaEqual[i, 4], lambdaEqual[i, 5], lambdaEqual[i, 6]];
  }
  
  if (betaEqN > 0) for (i in 1:betaEqN){
    Beta[betaEqual[i, 1], betaEqual[i, 2], betaEqual[i, 3]] = Beta[betaEqual[i, 4], betaEqual[i, 5], betaEqual[i, 6]];
  }
  
  if (thetaEqN > 0) for (i in 1:thetaEqN){
    Theta[thetaEqual[i, 1], thetaEqual[i, 2], thetaEqual[i, 3]] = Theta[thetaEqual[i, 4], thetaEqual[i, 5], thetaEqual[i, 6]];
  }
  if (psiEqN > 0) for (i in 1:psiEqN){
    Psi[psiEqual[i, 1], psiEqual[i, 2], psiEqual[i, 3]] = Psi[psiEqual[i, 4], psiEqual[i, 5], psiEqual[i, 6]];
  }
  
  if (nuEqN > 0) for (i in 1:nuEqN){
    Nu[nuEqual[i, 1], nuEqual[i, 2]] = Nu[nuEqual[i, 3], nuEqual[i, 4]];
  }
  if (alphaEqN > 0) for (i in 1:alphaEqN){
    Alpha[alphaEqual[i, 1], alphaEqual[i, 2]] = Alpha[alphaEqual[i, 3], alphaEqual[i, 4]];
  }
}

model {
  //vector[K + F] Y[N];
  vector[K] mu[G];
  matrix[K, K] Full_matrix[G];
  matrix[F, F] Pi[G];
  matrix[F, F] covPsi[G];
  
  //for (n in 1:N){
  //  Y[n, 1:K] = X[n];
  //  Y[n, (K+1):(K+F)] = eta[n];
  //}
  
  for (g in 1:G){
    Pi[g] = inverse(I - Beta[g]);
    covPsi[g] = Pi[g] * Psi[g] * Pi[g]';
    
    Full_matrix[g, 1:K, 1:K] = Lambda[g] * covPsi[g] * Lambda[g]' + Theta[g];
    //Full_matrix[g, (K+1):(K+F), 1:K] = covPsi[g] * Lambda[g]';
    //Full_matrix[g, 1:K, (K+1):(K+F)] = Lambda[g] * covPsi[g];
    //Full_matrix[g, (K+1):(K+F), (K+1):(K+F)] = covPsi[g];
  
    mu[g, 1:K] = Nu[g] + Lambda[g] * (Pi[g] * Alpha[g]);
    //mu[g, (K+1):(K+F)] = Pi[g] * Alpha[g];
    
    X[(group[g]):(group[g+1] - 1)] ~ multi_normal(mu[g], Full_matrix[g]);  
    
    to_vector(Lambda_full[g]) ~ normal(0, 3.0);
    to_vector(Beta_full[g]) ~ normal(0, 3.0);
    Psi_cor[g] ~ lkj_corr_cholesky(2.0);
    Psi_tau[g] ~ cauchy(0, 3.0);
    Theta_cor[g] ~ lkj_corr_cholesky(2.0);
    Theta_tau[g] ~ cauchy(0, 3.0);
    Nu_full[g] ~ normal(0, 10.0);
    Alpha_full[g] ~ normal(0, 10.0);
  
  }
  
}

generated quantities {
  real Chi_sample;
  real Chi_sim;
  real PPP;
  vector[N] log_lik;
  vector[F] eta[N];
  
    for (g in 1:G){
    matrix[K, F] lowerLeft;
    matrix[F, K] upperRight;
    matrix[K, K] Full_matrix;
    matrix[F, F] Pi;
    matrix[F, F] covPsi;
    
    Pi = inverse(I - Beta[g]);
    covPsi = Pi * Psi[g] * Pi';
  //  Full_matrix[g, (K+1):(K+F), 1:K] = Psi[g] * Lambda[g]';
    upperRight = covPsi * Lambda[g]';
  //  Full_matrix[g, 1:K, (K+1):(K+F)] = Lambda[g] * Psi[g];
    lowerLeft = Lambda[g] * covPsi;
  //  Full_matrix[g, (K+1):(K+F), (K+1):(K+F)] = Psi[g];
    Full_matrix[1:K, 1:K] = inverse(Lambda[g] * covPsi * Lambda[g]' + Theta[g]);
    for (n in group[g]:(group[g+1]-1))
      eta[n] = multi_normal_rng(Alpha[g] + upperRight * (Full_matrix) * (X[n] - Nu[g]), 
      covPsi - upperRight * (Full_matrix) * lowerLeft);
  }
  
  
  {
    real Fml_group[G];
    real Fml_sim_group[G];
    
    for (g in 1:G){
    
      matrix[(group[g+1] - group[g]), K] sim_data;
      matrix[K, K] sim_cov;
      matrix[K, K] pop_cov;
      
      pop_cov = Lambda[g] * (inverse(I - Beta[g]) * Psi[g] * inverse(I - Beta[g])') * Lambda[g]' + Theta[g];
      
      for (n in 1:(group[g+1] - group[g])){
        sim_data[n, 1:K] = multi_normal_rng(Nu[g] + Lambda[g] * (inverse(I - Beta[g]) * Alpha[g]), pop_cov)';
        log_lik[n + (group[g] - 1)] = multi_normal_lpdf(X[n + (group[g] - 1)]|Nu[g] + Lambda[g] * (inverse(I - Beta[g]) * Alpha[g]), pop_cov);
      }
      sim_cov = cov(sim_data);
      
      Fml_group[g] = Fml(pop_cov, sample_cov[g]);
      Fml_sim_group[g] = Fml(pop_cov, sim_cov);
      
    }
    Chi_sample = (N - 1) * sum(Fml_group);
    Chi_sim = (N - 1) * sum(Fml_sim_group);
    PPP = Chi_sample < Chi_sim;
  }
}
