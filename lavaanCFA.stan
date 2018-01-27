functions {
  // Simple function to compute the covariance matrix, mainly for PPC
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
  // Compute log likelihood
  real Fml(matrix pop_cov, matrix samp_cov){
    int K;
    K = rows(pop_cov);
    return log_determinant(pop_cov) + trace(samp_cov * inverse(pop_cov)) - log_determinant(samp_cov) - K;
  }
  
}

data {
  // Basic dimension for the model
  int<lower=1> G; // Number of independent groups (for multiple groups analyses)
  int<lower=1> N; // Number of individuals
  int<lower=3> K; // Number of indicator variables
  int<lower=1> F; // Number of factors
  // The data matrix as vectors
  vector[K] X[N];
  // Vector indicating were each group starts
  int<lower=1> group[G+1];
  
  // Constrained parameters as sparse matrices
  // Lambda (factor loadings)
  int lambdaN;
  int lambdaFree[lambdaN];
  int lambdaPar[lambdaN, 3];
  vector[lambdaN] lambdaConst;
  
  // Theta (residual variances)
  int thetaN;
  int thetaDiagN;
  int thetaOffDiagN;
  int thetaDiag[thetaDiagN];
  int thetaOffDiag[thetaOffDiagN];
  int thetaFree[thetaN];
  int thetaPar[thetaN, 3];
  vector[thetaN] thetaConst;
  
  // Psi (factor covariance)
  int psiN;
  int psiDiagN;
  int psiOffDiagN;
  int psiDiag[psiDiagN];
  int psiOffDiag[psiOffDiagN];
  int psiFree[psiN];
  int psiPar[psiN, 3];
  vector[psiN] psiConst;
  
  // Nu (indicator variables means)
  int nuN;
  int nuFree[nuN];
  int nuPar[nuN, 2];
  vector[nuN] nuConst;
  
  // Alpha (factor means)
  int alphaN;
  int alphaFree[alphaN];
  int alphaPar[alphaN, 2];
  vector[alphaN] alphaConst;
  
  // Sample covariance matrix (for PPC)
  matrix[K, K] sample_cov[G];
}

parameters {
  vector[max(lambdaFree)] Lambda_full;

  vector<lower=0>[size(thetaDiag) == 0 ? 0 : max(thetaFree[thetaDiag])] Theta_tau;
  vector<lower=-1, upper=1>[size(thetaOffDiag) == 0 ? 0 : max(thetaFree[thetaOffDiag])] Theta_cor;

  vector<lower=0>[size(psiDiag) == 0 ? 0 : max(psiFree[psiDiag])] Psi_tau;
  vector<lower=-1, upper=1>[size(psiOffDiag) == 0 ? 0 : max(psiFree[psiOffDiag])] Psi_cor;
  
  vector[max(nuFree)] Nu_full;
  vector[max(alphaFree)] Alpha_full;
}

transformed parameters {
  matrix[max(lambdaPar[, 2]), max(lambdaPar[, 3])] Lambda[G];
  matrix[max(thetaPar[, 2]), max(thetaPar[, 3])] Theta[G];
  matrix[max(psiPar[, 2]), max(psiPar[, 3])] Psi[G];
  vector[max(nuPar[, 2])] Nu[G];
  vector[max(alphaPar[, 2])] Alpha[G];
  
  // Set matrix values
  
  for (g in 1:G){
    for (v1 in 1:max(lambdaPar[, 2])){
     for (v2 in 1:max(lambdaPar[, 3])){
       Lambda[g, v1, v2] = 0.0;
     }
    }
    for (v1 in 1:max(thetaPar[, 2])){
     for (v2 in 1:max(thetaPar[, 3])){
       Theta[g, v1, v2] = 0.0;
      }
    }
    for (v1 in 1:max(psiPar[, 2])){
     for (v2 in 1:max(psiPar[, 3])){
       Psi[g, v1, v2] = 0.0;
      }
    }
    for (v1 in 1:max(nuPar[, 2])){
       Nu[g, v1] = 0.0;
    }
    for (v1 in 1:max(alphaPar[, 2])){
       Alpha[g, v1] = 0.0;
    }
  }  
  for (i in 1:lambdaN){
    if (lambdaFree[i] != 0){
      Lambda[lambdaPar[i, 1], lambdaPar[i, 2], lambdaPar[i, 3]] = Lambda_full[lambdaFree[i]];
    } else {
      Lambda[lambdaPar[i, 1], lambdaPar[i, 2], lambdaPar[i, 3]] = lambdaConst[i];
    }
  }

  for (i in 1:thetaN){
    if (thetaFree[i] != 0){
      if (thetaPar[i, 2] == thetaPar[i, 3]){
        Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 3]] = 
        Theta_tau[thetaFree[i]]*Theta_tau[thetaFree[i]];
      } else {
        
        Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 3]] = 
        Theta_cor[thetaFree[i]]*
        sqrt(Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 2]])*
        sqrt(Theta[thetaPar[i, 1], thetaPar[i, 3], thetaPar[i, 3]]);
        
        Theta[thetaPar[i, 1], thetaPar[i, 3], thetaPar[i, 2]] = 
        Theta_cor[thetaFree[i]]*
        sqrt(Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 2]])*
        sqrt(Theta[thetaPar[i, 1], thetaPar[i, 3], thetaPar[i, 3]]);
      }
    } else {
      if (thetaPar[i, 2] == thetaPar[i, 3]){
        Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 3]] = thetaConst[i];
      } else {
        Theta[thetaPar[i, 1], thetaPar[i, 2], thetaPar[i, 3]] = thetaConst[i];
        Theta[thetaPar[i, 1], thetaPar[i, 3], thetaPar[i, 2]] = thetaConst[i];
      }
    }
  }
  
  for (i in 1:psiN){
    if (psiFree[i] != 0){
      if (psiPar[i, 2] == psiPar[i, 3]){
        Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 3]] = Psi_tau[psiFree[i]]*Psi_tau[psiFree[i]];
      } else {
        Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 3]] = 
        Psi_cor[psiFree[i]]*
        sqrt(Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 2]])*
        sqrt(Psi[psiPar[i, 1], psiPar[i, 3], psiPar[i, 3]]);
        
        Psi[psiPar[i, 1], psiPar[i, 3], psiPar[i, 2]] = 
        Psi_cor[psiFree[i]]*
        sqrt(Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 2]])*
        sqrt(Psi[psiPar[i, 1], psiPar[i, 3], psiPar[i, 3]]);
      }
    } else {
      if (psiPar[i, 2] == psiPar[i, 3]){
        Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 3]] = psiConst[i];
      } else {
        Psi[psiPar[i, 1], psiPar[i, 2], psiPar[i, 3]] = psiConst[i];
        Psi[psiPar[i, 1], psiPar[i, 3], psiPar[i, 2]] = psiConst[i];
      }
    }
  }
  
  for (i in 1:nuN){
    if (nuFree[i]){
      Nu[nuPar[i, 1], nuPar[i, 2]] = Nu_full[nuFree[i]];
    } else {
      Nu[nuPar[i, 1], nuPar[i, 2]] = nuConst[i];
    }
  }
  for (i in 1:alphaN){
    if (alphaFree[i]){
      Alpha[alphaPar[i, 1], alphaPar[i, 2]] = Alpha_full[alphaFree[i]];
    } else {
      Alpha[alphaPar[i, 1], alphaPar[i, 2]] = alphaConst[i];
    }
  }

}

model {
  vector[K] mu[G];
  matrix[K, K] Full_matrix[G];
  for (g in 1:G){
    Full_matrix[g, 1:K, 1:K] = Lambda[g] * Psi[g] * Lambda[g]' + Theta[g];

    mu[g, 1:K] = Nu[g] + Lambda[g] * Alpha[g];
    
    X[(group[g]):(group[g+1] - 1)] ~ multi_normal(mu[g], Full_matrix[g]);  
  
  }
  
    Lambda_full ~ normal(0, 3.0);
    //Psi_cor[g] ~ normal(0, 0.5);
    Psi_tau ~ cauchy(0, 3.0);
    //Theta_cor[g] ~ normal(0, 0.5);
    Theta_tau ~ cauchy(0, 3.0);
    Nu_full ~ normal(0, 10.0);
    Alpha_full ~ normal(0, 10.0);
  
}

generated quantities {
  real Chi_sample;
  real Chi_sim;
  real PPP;
  vector[N] log_lik;
  vector[max(psiPar[, 2])] eta[N];
  
  for (g in 1:G){
    matrix[K, max(psiPar[, 2])] lowerLeft;
    matrix[max(psiPar[, 2]), K] upperRight;
    matrix[K, K] Full_matrix;
    upperRight = Psi[g] * Lambda[g]';
    lowerLeft = Lambda[g] * Psi[g];
    Full_matrix[1:K, 1:K] = inverse(Lambda[g] * Psi[g] * Lambda[g]' + Theta[g]);
    for (n in group[g]:(group[g+1]-1))
      eta[n] = multi_normal_rng(Alpha[g] + upperRight * (Full_matrix) * (X[n] - Nu[g]), 
      Psi[g] - upperRight * (Full_matrix) * lowerLeft);
  }
  
 
  {
    real Fml_group[G];
    real Fml_sim_group[G];
    
    for (g in 1:G){
    
      matrix[(group[g+1] - group[g]), K] sim_data;
      matrix[K, K] sim_cov;
      matrix[K, K] pop_cov;
      
      pop_cov = Lambda[g] * Psi[g] * Lambda[g]' + Theta[g];
      
      for (n in 1:(group[g+1] - group[g])){
        sim_data[n, 1:K] = multi_normal_rng(Nu[g] + Lambda[g] * Alpha[g], pop_cov)';
        log_lik[n + (group[g] - 1)] = multi_normal_lpdf(X[n + (group[g] - 1)]|Nu[g] + Lambda[g] * Alpha[g], pop_cov);
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
