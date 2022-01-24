# This function is from the Carnoma and Nicholls (2020)
mcmc_hpv <- function( HPV, # HPV data
                      
                      hpv_model_whole_stan,
                      hpv_model_poisson_stan,
                      
                      # SMI degree of influence for each module
                      eta_pois = 1,
                      eta_binom = 1,
                      
                      # Number of iterations
                      n_iter = 10000, # main chain
                      n_warmup = 1000,
                      n_chains_mcmc = 4, # Number of chains
                      n_iter_sub = 200, # Subchain
                      
                      out_file_rda=NULL,
                      n_cores=1 ) {
  
  # Data in stan format #
  HPV_data_stan <- list( n_obs=nrow(HPV),
                         nhpv = HPV$nhpv,
                         Npart = HPV$Npart,
                         ncases = HPV$ncases,
                         Npop = HPV$Npop,
                         eta_pois=eta_pois,
                         eta_binom=eta_binom,
                         phi=NULL )
  
  ### Power likelihood ###
  cat(paste("Stage 1 starts.", "\n", sep=""), append = T, 
      file = paste(out_file_rda, "progress.txt", sep = "_"))
  stan_fit <- rstan::sampling( hpv_model_whole_stan,
                               data = HPV_data_stan,
                               iter = n_iter,
                               warmup = n_warmup,
                               chains = n_chains_mcmc,
                               cores = n_cores,
                               # pars=c("phi","theta1","theta2"),
                               show_messages=FALSE )
  hpv_mcmc_pow <- rstan::extract(stan_fit)
  
  ### Multiple imputation ###
  # Conventional Bayesian inference on poisson model, conditional on imputed phi #
  
  cat(paste("Stage 2 starts.", "\n", sep=""), append = T, 
      file = paste(out_file_rda, "progress.txt", sep = "_"))
  
  HPV_data_stan$eta_pois = 1
  HPV_data_stan$eta_binom = 0
  
  imp_i=1
  hpv_mcmc_stage2 <- foreach::foreach( imp_i = 1:nrow(hpv_mcmc_pow$phi),
                                       .combine = rbind ) %dorng% {
                                         # imp_i=1
                                         
                                         # Select one imputed value of phi
                                         HPV_data_stan$phi = hpv_mcmc_pow$phi[imp_i, ]
                                         
                                         # Sample the poisson module conditional on such imputed phi
                                         stan_fit <- rstan::sampling( hpv_model_poisson_stan,
                                                                      data = HPV_data_stan,
                                                                      iter = n_iter_sub,
                                                                      chains = 1,
                                                                      pars=c("theta1","theta2"),
                                                                      show_messages=FALSE )
                                         stan_fit_mcmc <- rstan::extract(stan_fit)
                                         
                                         # Return only the last sample
                                         n_aux <- length(stan_fit_mcmc$theta1)
                                         c(stan_fit_mcmc$theta1[n_aux], stan_fit_mcmc$theta2[n_aux])
                                       }
  
  colnames(hpv_mcmc_pow$phi) <- paste("phi_", 1:nrow(HPV), sep = "")
  colnames(hpv_mcmc_stage2) <- paste("theta", 1:2, sep = "")
  
  hpv_mcmc_smi <- data.frame( hpv_mcmc_pow$phi,
                              hpv_mcmc_stage2,
                              theta_tilde_1 = hpv_mcmc_pow$theta1,
                              theta_tilde_2 = hpv_mcmc_pow$theta2 )
  rownames(hpv_mcmc_smi) = NULL
  
  # Save results #
  if( !is.null(out_file_rda) ){
    save( hpv_mcmc_smi, file=out_file_rda )
  }
  
  return( hpv_mcmc_smi )
  
}

library(rstan)
library(doRNG)

## Load data
inputPath = args = commandArgs(trailingOnly=TRUE)
input.df = read.table(file = inputPath, header = T, as.is = T, sep = "\t")


## Set up Monte Carlo parameters
n_iter = as.numeric(input.df$value[which(input.df$param == "n_iter")])# main chain
n_warmup = as.numeric(input.df$value[which(input.df$param == "n_warmup")])
n_chains_mcmc = as.numeric(input.df$value[which(input.df$param == "n_chains_mcmc")]) # Number of chains
n_iter_sub = as.numeric(input.df$value[which(input.df$param == "n_iter_sub")]) # Subchain
n_cores = as.numeric(input.df$value[which(input.df$param == "n_cores")])
eta_pois = as.numeric(input.df$value[which(input.df$param == "eta_pois")])
out_file_rda = input.df$value[which(input.df$param == "out_file_rda")]
hpv_path = input.df$value[which(input.df$param == "hpv_data_path")]
hpv_model_poisson_stan_path = input.df$value[which(input.df$param == "hpv_model_poisson_stan_path")]
hpv_model_whole_stan_path = input.df$value[which(input.df$param == "hpv_model_whole_stan_path")]

data = read.csv(file = hpv_path, header = T)

hpv_model_whole_stan = stan_model(file = hpv_model_whole_stan_path )
hpv_model_poisson_stan = rstan::stan_model(file = hpv_model_poisson_stan_path )

hpv_mcmc_smi =  mcmc_hpv( HPV = data,
                          hpv_model_whole_stan = hpv_model_whole_stan,
                          hpv_model_poisson_stan = hpv_model_poisson_stan,
                          
                          # Number of iterations
                          n_iter = n_iter, # main chain
                          n_warmup = n_warmup,
                          n_chains_mcmc = n_chains_mcmc, # Number of chains
                          n_iter_sub = n_iter_sub, # Subchain
                          
                          # Cut rate
                          eta_pois = eta_pois,
                          eta_binom = 1,
                          out_file_rda = out_file_rda,
                          n_cores = n_cores )