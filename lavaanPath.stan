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
  matrix[F, F] Beta_full[G];
  
  cholesky_factor_corr[F] Psi_cor[G];
  vector<lower=0>[F] Psi_tau[G];
  
  vector[F] Alpha_full[G];
}

transformed parameters {
  matrix[F, F] Beta[G];
  cov_matrix[F] Psi[G];
  vector[F] Alpha[G];
  
  
  for (g in 1:G){
    Beta[g] = Beta_full[g];
    Psi[g] = quad_form_diag(Psi_cor[g] * Psi_cor[g]', Psi_tau[g]);
    Alpha[g] = Alpha_full[g];
  }
  
  // Set fixed value constraints
  
  if (betaN > 0) for (i in 1:betaN){
    Beta[betaPar[i, 1], betaPar[i, 2], betaPar[i, 3]] = betaConst[i];
  }

  if (psiN > 0) for (i in 1:psiN){
    Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 3]] = psiConst[i];
  }
  
  if (alphaN > 0) for (i in 1:alphaN){
    Alpha[alphaPar[i, 1], alphaPar[i, 2]] = alphaConst[i];
  }
  
  // Set equality constraints
  
  if (betaEqN > 0) for (i in 1:betaEqN){
    Beta[betaEqual[i, 1], betaEqual[i, 2], betaEqual[i, 3]] = Beta[betaEqual[i, 4], betaEqual[i, 5], betaEqual[i, 6]];
  }
  
  if (psiEqN > 0) for (i in 1:psiEqN){
    Psi[psiEqual[i, 1], psiEqual[i, 2], psiEqual[i, 3]] = Psi[psiEqual[i, 4], psiEqual[i, 5], psiEqual[i, 6]];
  }
  
  if (alphaEqN > 0) for (i in 1:alphaEqN){
    Alpha[alphaEqual[i, 1], alphaEqual[i, 2]] = Alpha[alphaEqual[i, 3], alphaEqual[i, 4]];
  }
}

model {
  vector[K] mu[G];
  matrix[F, F] Pi[G];
  matrix[F, F] covPsi[G];
  
  for (g in 1:G){
    Pi[g] = inverse(I - Beta[g]);
    covPsi[g] = Pi[g] * Psi[g] * Pi[g]';
    

    mu[g, 1:K] = (Pi[g] * Alpha[g]);

    X[(group[g]):(group[g+1] - 1)] ~ multi_normal(mu[g], covPsi[g]);  
    
    to_vector(Beta_full[g]) ~ normal(0, 3.0);
    Psi_cor[g] ~ lkj_corr_cholesky(2.0);
    Psi_tau[g] ~ cauchy(0, 3.0);
    Alpha_full[g] ~ normal(0, 10.0);
  
  }
  
}

generated quantities {
  real Chi_sample;
  real Chi_sim;
  real PPP;
  vector[N] log_lik;
  {
    real Fml_group[G];
    real Fml_sim_group[G];
    
    for (g in 1:G){
    
      matrix[(group[g+1] - group[g]), K] sim_data;
      matrix[K, K] sim_cov;
      matrix[K, K] pop_cov;
      
      pop_cov = (inverse(I - Beta[g]) * Psi[g] * inverse(I - Beta[g])');
      
      for (n in 1:(group[g+1] - group[g])){
        sim_data[n, 1:K] = multi_normal_rng((inverse(I - Beta[g])*Alpha[g]), pop_cov)';
        log_lik[n + (group[g] - 1)] = multi_normal_lpdf(X[n + (group[g] - 1)] | (inverse(I - Beta[g])*Alpha[g]), pop_cov);
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
