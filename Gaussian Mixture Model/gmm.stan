data {
  int<lower=1> N_data;
  array[N_data] real y;
}

parameters {
  simplex[3] psi; // Mixture weights sum to 1
  real<lower=0,upper=0.5> X; // left peak position
  real<lower=1> N; // effective number of alleles
}

transformed parameters {
  // Means of the three Gaussian peaks
  real mu_left = X;
  real mu_right = 1 - X;
  real mu_middle = 0.5;

  // Variances (p*(1-p)/N) for each peak
  real var_left = X*(1-X)/N;
  real var_right = X*(1-X)/N;
  real var_middle = 0.25/N; // 0.5*(1-0.5)=0.25

  // Standard deviations
  real sigma_left = sqrt(var_left);
  real sigma_right = sqrt(var_right);
  real sigma_middle = sqrt(var_middle);
}

model {
  psi ~ dirichlet([5.0, 10.0, 5.0]); // quarter half quarter weighting with tuning for variance
  X ~ beta(1,4);
  N ~ gamma(2, 0.1);

  for (n in 1:N_data) {
    array[3] real lps;
    lps[1] = log(psi[1]) + normal_lpdf(y[n] | mu_left, sigma_left);
    lps[2] = log(psi[2]) + normal_lpdf(y[n] | mu_middle, sigma_middle);
    lps[3] = log(psi[3]) + normal_lpdf(y[n] | mu_right, sigma_right);

    target += log_sum_exp(lps);
  }
}

generated quantities {
  real phi = -log1m(2 * X) / 2;
}