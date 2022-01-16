library(rstan)

## Load data
inputPath = args = commandArgs(trailingOnly=TRUE)
input.df = read.table(file = inputPath, header = T, as.is = T, sep = "\t")


## Set up Monte Carlo parameters
### Number of iteration in the main chain
n_iter = as.numeric(input.df$value[which(input.df$param == "n_iter")])
### Number of steps as burn-in for sampling() in stan
n_warmup = as.numeric(input.df$value[which(input.df$param == "n_warmup")])
### Length of main MCMC chain in stage 1 of delta-SMI
n_chains_mcmc = as.numeric(input.df$value[which(input.df$param == "n_chains_mcmc")])
### Length of main MCMC chain in stage 1 of delta-SMI
n_iter_sub = as.numeric(input.df$value[which(input.df$param == "n_iter_sub")])
### Path of HCV data
hpv_path = input.df$value[which(input.df$param == "hpv_data_path")]
### Path of output file
out_file_rda = input.df$value[which(input.df$param == "out_file_rda")]
### Number of cores
n_cores = as.numeric(input.df$value[which(input.df$param == "n_cores")])
### delta value used for delta-SMI
delta = as.numeric(input.df$value[which(input.df$param == "delta")])
### Path to the stan script for the HPV model used for the subchains in stage 2 of delta-SMI
hpv_model_poisson_stan_path = input.df$value[which(input.df$param == "hpv_model_poisson_stan_path")]
### Path to the stan script for the HPV model used for the main chain 
### in stage 1 of delta-SMI allowing delta to by scaled by sqrt(y_i)
hpv_model_whole_abc_smi_delta_scaled_stan_path = 
  input.df$value[which(input.df$param == "hpv_model_whole_abc_smi_delta_scaled_stan_path")]
### Initial theta1 value (Bayesian posterior mean is recommended)
initTheta1 = as.numeric(input.df$value[which(input.df$param == "initTheta1")])
### Initial theta2 value (Bayesian posterior mean is recommended)
initTheta2 = as.numeric(input.df$value[which(input.df$param == "initTheta2")])
### Initial phi1--13 values (Bayesian posterior mean is recommended)
initPhi = as.numeric(trimws(unlist(strsplit(input.df$value[which(input.df$param == "initPhi")], split = ","))))

## Load HPV data
HPV = read.csv(file = hpv_path, header = T)


## Create initial value lists
init_param_list_temp = list(theta1 = initTheta1, theta2 = initTheta2, phi = initPhi)
init_param = list()
for(index in 1:n_chains_mcmc){
  init_param[[index]] = init_param_list_temp
}


## Create models in stan
hpv_model_poisson_stan = stan_model(file = hpv_model_poisson_stan_path)
hpv_model_whole_abc_smi_delta_scaled_stan = stan_model(file = hpv_model_whole_abc_smi_delta_scaled_stan_path)

## Scaled delta
delta_scaled = delta*round(sqrt(HPV$ncases))

## Create data list for stan
hpv_abc_smi_data_stan = list( n_obs=nrow(HPV),
                              nhpv = HPV$nhpv,
                              Npart = HPV$Npart,
                              ncases = HPV$ncases,
                              Npop = HPV$Npop,
                              delta = delta_scaled,
                              phi = NULL)

## Stage 1
cat("ABC-SMI stage 1 begins.\n", file = paste(out_file_rda,"_progress.txt", sep=""))

hpv_abc_smi_stan_stage1 = sampling( hpv_model_whole_abc_smi_delta_scaled_stan,
                                    data = hpv_abc_smi_data_stan,
                                    iter = n_iter,
                                    warmup = n_warmup,
                                    chains = n_chains_mcmc,
                                    cores = n_cores,
                                    init = init_param,
                                    show_messages = FALSE)

hpv_abc_smi_stage1_param = rstan::extract(hpv_abc_smi_stan_stage1)

cat("ABC-SMI stage 1 completed.\n", append = T, file = paste(out_file_rda,"_progress.txt", sep=""))


## Stage 2
cat("ABC-SMI stage 2 begins.\n", append = T, file = paste(out_file_rda,"_progress.txt", sep=""))

hpv_abc_smi_data_stan_temp = hpv_abc_smi_data_stan
hpv_theta_abc_smi = matrix(nrow = nrow(hpv_abc_smi_stage1_param$phi), ncol = 2)
colnames(hpv_theta_abc_smi) = c("theta1", "theta2")
for(index in 1:nrow(hpv_theta_abc_smi)){
#for(index in 1:1000){
  
  # Conditioning on phi estimates obtained in stage 1
  hpv_abc_smi_data_stan_temp$phi = hpv_abc_smi_stage1_param $phi[index,]
  
  # Obtain theta estimates from stage 2 of ABC-SMI 
  hpv_abc_smi_stage2_stan = rstan::sampling( hpv_model_poisson_stan,
                                             data = hpv_abc_smi_data_stan_temp,
                                             iter = n_iter_sub,
                                             chains = 1,
                                             pars = c("theta1", "theta2"),
                                             show_messages = FALSE )
  
  # Extract last theta estimates from stage 2 of ABC-SMI 
  hpv_abc_smi_stage2_param = rstan::extract(hpv_abc_smi_stage2_stan)
  s2N = length(hpv_abc_smi_stage2_param$theta1)
  hpv_theta_abc_smi[index,] = c(hpv_abc_smi_stage2_param$theta1[s2N],
                                hpv_abc_smi_stage2_param$theta2[s2N])
  if(index%%1000 == 0){
    cat(paste(index,"\n", sep=""), append = T, file = paste(out_file_rda,"_progress.txt", sep=""))
  }
}


## Create output file
hpv_abc_smi_output.df = data.frame(hpv_abc_smi_stage1_param$phi, 
                                   theta1_tilde = hpv_abc_smi_stage1_param$theta1,
                                   theta2_tilde = hpv_abc_smi_stage1_param$theta2,
                                   hpv_theta_abc_smi)
names(hpv_abc_smi_output.df)[1:ncol(hpv_abc_smi_stage1_param$phi)] = 
  paste("phi", 1:ncol(hpv_abc_smi_stage1_param$phi), sep = "")


## Write to output file
cat(paste("Output file to", out_file_rda), append = T, file = paste(out_file_rda,"_progress.txt", sep=""))
save(hpv_abc_smi_output.df, 
     hpv_abc_smi_stan_stage1, 
     file = out_file_rda)