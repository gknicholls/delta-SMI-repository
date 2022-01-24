## Load data from the HPV.csv file
HPV = read.csv(file=file.choose(), header = T)

library(rstan)

## Set up Monte Carlo parameters
n_iter = 10000# main chain
n_warmup = 1000
n_chains_mcmc = 4 # Number of chains
n_iter_sub = 200 # Subchain

out_file_rda=NULL
n_cores=1 

## Create models in stan
hpv_model_whole_stan_path = "hpv_model_whole.stan" # Assuming hpv_model_whole.stan is in the working directory.
hpv_model_whole_stan = stan_model(file = hpv_model_whole_stan_path)
hpv_model_poisson_stan_path = "hpv_model_poisson.stan" # Assuming hpv_model_poisson.stan is in the working directory.
hpv_model_poisson_stan = stan_model(file = hpv_model_poisson_stan_path)

###############################################
################ Full Bayes  ##################
###############################################


## Set up data list for stan 
## with eta = 1 for both Poisson and binomial modules.
HPV_bayes_data_stan <- list( n_obs=nrow(HPV),
                       nhpv = HPV$nhpv,
                       Npart = HPV$Npart,
                       ncases = HPV$ncases,
                       Npop = HPV$Npop,
                       eta_pois=1,
                       eta_binom=1,
                       phi=NULL)

## Obtain estimates of full Bayes model via stan
hpv_bayes_fit = sampling( hpv_model_whole_stan,
                          data = HPV_bayes_data_stan,
                          iter = n_iter,
                          warmup = n_warmup,
                          chains = n_chains_mcmc,
                          cores = n_cores,
                          show_messages = FALSE )

## Extract the posteriors of full Bayes estimates
hpv_bayes_param = rstan::extract(hpv_bayes_fit)

## Save posteriors of full Bayes estimates to the working directory.
save(hpv_bayes_param, file = "hpv_bayes.RData") 

###############################################
################ Cut Model  ###################
###############################################

## Set up data list for stan 
## with eta = 0 for the Poisson module
## and eta = 1 for binomial module.
HPV_cut_data_stan <- list( n_obs=nrow(HPV),
                           nhpv = HPV$nhpv,
                           Npart = HPV$Npart,
                           ncases = HPV$ncases,
                           Npop = HPV$Npop,
                           eta_pois = 0,
                           eta_binom = 1,
                           phi = NULL)

# Stage 1: Obtain phi estimates under the binomial model only
phiCutMat = matrix(nrow = 36000, ncol = 13)
for(paramIndex in 1:13){
  phiCutMat[, paramIndex] = rbeta(nrow(phiCutMat),
                                  shape1 = HPV$nhpv[paramIndex] + 1, 
                                  shape2 = HPV$Npart[paramIndex] - HPV$nhpv[paramIndex] + 1)  
}

# Stage 2: Obtain theta estimates conditioned 
# on the phi estimates in stage 1.
hpv_theta_cut = matrix(nrow = nrow(phiCutMat), ncol = 2)
for(index in 1:nrow(phiCutMat)){
  
  ## Obtain estimates of cut model via stan
  HPV_cut_data_stan$phi = phiCutMat[index,]
  hpv_cut_stage2_stan_fit = rstan::sampling( hpv_model_poisson_stan,
                                             data = HPV_cut_data_stan,
                                             iter = n_iter_sub,
                                             chains = 1,
                                             pars=c("theta1","theta2"),
                                             show_messages=FALSE )
  hpv_stage2_cut = rstan::extract(hpv_cut_stage2_stan_fit)
  s2N = length(hpv_stage2_cut$theta1)
  hpv_theta_cut[index,] = c(hpv_stage2_cut$theta1[s2N],
                            hpv_stage2_cut$theta2[s2N])
}

## Extract the posteriors of cut estimates
colnames(phiCutMat) = paste("phi", c(1:ncol(phiCutMat)), sep="_")
colnames(hpv_theta_cut) = c("theta_1", "theta_2")
hpv_cut_param =  cbind(phiCutMat, hpv_theta_cut)

## Save posteriors of cut estimates to the working directory.
save(hpv_cut_param, file = "hpv_cut.RData")