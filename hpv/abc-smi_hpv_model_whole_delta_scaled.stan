
// saved as plummer_hpv.stan

data {
  int<lower=0> n_obs; // number of observations

  int<lower=0> nhpv[n_obs];
  int<lower=0> Npart[n_obs];
  int<lower=0> ncases[n_obs];
  real<lower=0> Npop[n_obs];
  int<lower=0> delta[n_obs];
}

parameters {
  real theta1;
  real<lower=0> theta2;
  real<lower=0,upper=1> phi[n_obs];
}

transformed parameters {
  real<lower=0> mu[n_obs];
  for (i in 1:n_obs) {
    mu[i] = (Npop[i]/1000)* exp( theta1 + phi[i] * theta2);
  }
}

model {
  // The likelihood
  
  
  for (i in 1:n_obs) {
    int nbrLwr = ncases[i] - delta[i] - 1;
    int nbrUpr = ncases[i] + delta[i];
    real logUprCDF;
    
    
    logUprCDF = poisson_lcdf( nbrUpr | mu[i] );
    
    if(nbrLwr < 1){
      target += (logUprCDF - log( delta[i]*2 + 1 ));
    }else{
      real logLwrCDF = poisson_lcdf( nbrLwr | mu[i] );
      real logCdfDiff = log(exp(logUprCDF - logLwrCDF) - 1) + logLwrCDF;
      target += (logCdfDiff - log( delta[i]*2 + 1 )); // Corrected 21 Feb 2022
    }
    
    target += binomial_lpmf( nhpv[i] | Npart[i], phi[i] );
    
  }
}
