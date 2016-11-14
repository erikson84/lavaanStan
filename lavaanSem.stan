data {
  // Basic dimension for the model
  int<lower=1> N;
  int<lower=3> K;
  int<lower=1> F;
  // The data matrix as vectors
  vector[K] X[N];
  
  // Constrained and free parameters matrices
  matrix[K, F] Lambda_const;
  matrix[F, F] Beta_const;
  matrix[K, K] Psi_const;
  matrix[F, F] Phi_const;
  vector[K] Nu_const;
  vector[F] Alpha_const;
  
  // Sample covariance matrix (for PPC)
  matrix[K, K] sample_cov;
}

transformed data {
  vector[F] one;
  matrix[F, F] I;
  for (i in 1:F) one[i] = 1.0;
  I = diag_matrix(one);
  
}

parameters {
  matrix[K, F] Lambda_full;
  matrix[F, F] Beta_full;
  // Positive lambda for standardized models
  vector<lower=0>[F] lambda_pos;
  vector[F] phi[N];

  cholesky_factor_corr[K] Psi_cor;
  //vector<lower=0>[K] Psi_tau;
  vector<lower=0, upper=pi()/2>[K] Psi_unif;
  cholesky_factor_corr[F] Phi_cor;
  //vector<lower=0>[F] Phi_tau;
  vector<lower=0, upper=pi()/2>[F] Phi_unif;
  vector[K] Nu_full;
  vector[F] Alpha_full;
}

transformed parameters {
  matrix[K, F] Lambda;
  matrix[F, F] Beta;
  cov_matrix[K] Psi;
  cov_matrix[F] PHI;
  vector[K] Nu;
  vector[F] Alpha;
  vector[K] Psi_tau;
  vector[F] Phi_tau;
  
  for (k in 1:K) Psi_tau[k] = 3.0 * tan(Psi_unif[k]);
  for (f in 1:F) Phi_tau[f] = 3.0 * tan(Phi_unif[f]);
  
  Psi = quad_form_diag(Psi_cor * Psi_cor', Psi_tau);
  PHI = quad_form_diag(Phi_cor * Phi_cor', Phi_tau);
  
  for (k in 1:K){
    // 357 is just an easy placeholder to indicate free parameters.
    if (Nu_const[k] != 357){
      Nu[k] = Nu_const[k];
    } else {
      Nu[k] = Nu_full[k];
    }
    
    for (f in 1:F){
      if (Lambda_const[k, f] != 357){
        Lambda[k, f] = Lambda_const[k, f];
      } else {
        if (k == 1) {
          Lambda[k, f] = lambda_pos[f];
        } else if (Lambda_const[k-1, f] == 0){
          Lambda[k, f] = lambda_pos[f];
        } else {
          Lambda[k, f] = Lambda_full[k, f];
        }
      }
      
      if (Alpha_const[f] != 357){
        Alpha[f] = Alpha_const[f];
      } else {
        Alpha[f] = Alpha_full[f];
      }
      
      for (f1 in 1:F){
        if (Phi_const[f, f1] != 357){
          PHI[f, f1] = Phi_const[f, f1];
        }
        if (Beta_const[f, f1] != 357){
          Beta[f, f1] = Beta_const[f, f1];
        } else {
          Beta[f, f1] = Beta_full[f, f1];
        }
      }
    }
    
    for (k1 in 1:K){
      if (Psi_const[k, k1] != 357){
        Psi[k, k1] = Psi_const[k, k1];
      }
    }
    
  }
}

model {
  vector[K + F] Y[N];
  vector[K + F] mu;
  matrix[K + F, K + F] Full_matrix;
  matrix[F, F] Pi;
  matrix[F, F] covPhi;
  
  for (n in 1:N){
    Y[n, 1:K] = X[n];
    Y[n, (K+1):(K+F)] = phi[n];
  }
  Pi = inverse(I - Beta);
  covPhi = (Pi * PHI * Pi');
  Full_matrix[1:K, 1:K] = Lambda * covPhi * Lambda' + Psi;
  Full_matrix[(K+1):(K+F), 1:K] = covPhi * Lambda';
  Full_matrix[1:K, (K+1):(K+F)] = Lambda * covPhi;
  Full_matrix[(K+1):(K+F), (K+1):(K+F)] = covPhi;
  
  mu[1:K] = Nu + Lambda * (Pi * Alpha);
  mu[(K+1):(K+F)] = Pi * Alpha;
  
  Y ~ multi_normal(mu, Full_matrix);
  
  to_vector(Lambda_full) ~ normal(0, 3.0);
  to_vector(Beta_full) ~ normal(0, 3.0);
  lambda_pos ~ normal(0, 3.0);
  Phi_cor ~ lkj_corr_cholesky(2.0);
  //Phi_tau ~ cauchy(0, 3.0);
  Psi_cor ~ lkj_corr_cholesky(2.0);
  //Psi_tau ~ cauchy(0, 3.0);
  Nu_full ~ normal(0, 10.0);
  Alpha_full ~ normal(0, 10.0);
  
}

generated quantities {
  real Chi_sample;
  real Chi_sim;
  real PPP;
  {
    matrix[N, K] sim_data;
    matrix[K, K] sim_cov;
    matrix[K, K] cov;
    vector[N] ones;
    vector[K] means;
    
    
    cov = Lambda * PHI * Lambda' + Psi;
    for (n in 1:N){
      sim_data[n, 1:K] = multi_normal_rng(Nu + Lambda * Alpha,
      Lambda * (inverse(I - Beta) * PHI * inverse(I - Beta)') * Lambda' + Psi)';
      ones[n] = 1.0;
    }
    
    for (k in 1:K){
        means[k] = mean(sim_data[, k]);
    }
    
    sim_cov = ((sim_data - ones * means')' * (sim_data - ones * means')) / (N - 1);
    Chi_sample = (N-1)*(log_determinant(cov) + trace(sample_cov * inverse(cov)) - log_determinant(sample_cov) - K);
    Chi_sim = (N-1)*(log_determinant(cov) + trace(sim_cov * inverse(cov)) - log_determinant(sim_cov) - K);
    PPP = Chi_sample < Chi_sim;
  }
  

  
}
