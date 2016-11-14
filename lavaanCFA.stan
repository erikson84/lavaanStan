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
  matrix[K, K] sample_cov;
}

parameters {
  matrix[K, F] Lambda_full[G];
  // Positive lambda for standardized models
  vector<lower=0>[F] lambda_pos[G];
  vector[F] eta[N];

  cholesky_factor_corr[K] Theta_cor[G];
  vector<lower=0>[K] Theta_tau[G];
  cholesky_factor_corr[F] Psi_cor[G];
  vector<lower=0>[F] Psi_tau[G];
  vector[K] Nu_full[G];
  vector[F] Alpha_full[G];
}

transformed parameters {
  matrix[K, F] Lambda[G];
  cov_matrix[K] Theta[G];
  cov_matrix[F] Psi[G];
  vector[K] Nu[G];
  vector[F] Alpha[G];
  
  
  for (g in 1:G){
    Lambda[g] = Lambda_full[g];
    Theta[g] = quad_form_diag(Theta_cor[g] * Theta_cor[g]', Theta_tau[g]);
    Psi[g] = quad_form_diag(Psi_cor[g] * Psi_cor[g]', Psi_tau[g]);
    Nu[g] = Nu_full[g];
    Alpha[g] = Alpha_full[g];
  }
  
  // Set fixed value constraints
  
  if (lambdaN > 0) for (i in 1:lambdaN){
    Lambda[lambdaPar[i, 1], lambdaPar[i, 2], lambdaPar[i, 3]] = lambdaConst[i];
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
  vector[K + F] Y[N];
  vector[K + F] mu[G];
  matrix[K + F, K + F] Full_matrix[G];
  
  for (n in 1:N){
    Y[n, 1:K] = X[n];
    Y[n, (K+1):(K+F)] = eta[n];
  }
  
  for (g in 1:G){
    Full_matrix[g, 1:K, 1:K] = Lambda[g] * Psi[g] * Lambda[g]' + Theta[g];
    Full_matrix[g, (K+1):(K+F), 1:K] = Psi[g] * Lambda[g]';
    Full_matrix[g, 1:K, (K+1):(K+F)] = Lambda[g] * Psi[g];
    Full_matrix[g, (K+1):(K+F), (K+1):(K+F)] = Psi[g];
  
    mu[g, 1:K] = Nu[g] + Lambda[g] * Alpha[g];
    mu[g, (K+1):(K+F)] = Alpha[g];
    
    Y[(group[g]):(group[g+1])] ~ multi_normal(mu[g], Full_matrix[g]);  
    
    to_vector(Lambda_full[g]) ~ normal(0, 3.0);
    lambda_pos[g] ~ normal(0, 3.0);
    Psi_cor[g] ~ lkj_corr_cholesky(2.0);
    Psi_tau[g] ~ cauchy(0, 3.0);
    Theta_cor[g] ~ lkj_corr_cholesky(2.0);
    Theta_tau[g] ~ cauchy(0, 3.0);
    Nu_full[g] ~ normal(0, 10.0);
    Alpha_full[g] ~ normal(0, 10.0);
  
  }
  
}

// generated quantities {
//   real Chi_sample;
//   real Chi_sim;
//   real PPP;
//   {
//     matrix[N, K] sim_data;
//     matrix[K, K] sim_cov;
//     matrix[K, K] cov;
//     vector[N] ones;
//     vector[K] means;
//     
//     
//     cov = Lambda * Psi * Lambda' + Theta;
//     for (n in 1:N){
//       sim_data[n, 1:K] = multi_normal_rng(Nu + Lambda * Alpha, cov)';
//       ones[n] = 1.0;
//     }
//     
//     for (k in 1:K){
//         means[k] = mean(sim_data[, k]);
//     }
//     
//     sim_cov = ((sim_data - ones * means')' * (sim_data - ones * means')) / (N - 1);
//     Chi_sample = (N-1)*(log_determinant(cov) + trace(sample_cov * inverse(cov)) - log_determinant(sample_cov) - K);
//     Chi_sim = (N-1)*(log_determinant(cov) + trace(sim_cov * inverse(cov)) - log_determinant(sim_cov) - K);
//     PPP = Chi_sample < Chi_sim;
//   }
//   
// 
//   
// }
