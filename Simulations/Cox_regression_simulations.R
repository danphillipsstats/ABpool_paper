### Pilot simulation study
# Load packages
t_init <- Sys.time()

project_directory <- getwd()
simulation_directory <- file.path(project_directory, "Simulations")
program_directory <- simulation_directory
output_directory <- file.path(simulation_directory, "Output")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
name <- "v20_cox_regression_nobs_propmis_20_scen_Bayesimp_cens20_B5000"

require(mice)
require(parallel)
require(survival)
require(abind)

exponential_model = rstan::stan_model(file=paste0(program_directory,'/Normal_covariate_exponential_model.stan'))

# Function for skewness
skewness <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  sum(((x - m) / s)^3) * n / ((n - 1) * (n - 2))
}

# Function for kurtosis
kurtosis <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  num <- n * (n + 1) * sum(((x - m) / s)^4)
  denom <- (n - 1) * (n - 2) * (n - 3)
  adj <- 3 * (n - 1)^2 / ((n - 2) * (n - 3))
  (num / denom) - adj
}
#######
# Simulating normal data
alpha <- 0.05
power <- 0.8

# Calculates the power in a simple linear regression
# for sample size n, effect, effect beta
# using Monte-Carlo integration with B_small samples
# effect_size is then sqrt(beta^2/(1+beta^2))
# Y = X*X_eff + Z*Z_eff + epsilon, where X and epsilon N(0,1) independent, Z binary 1, 0.
sim_cox_pow <- function(n, X_eff, Z_eff, basehaz = 1, B_small = 10000){
  pow.small <- numeric(B_small)
  for (i in 1:B_small){
    X <- rnorm(n) # continuous covariate
    Z <- c(rep(1,n/2),rep(0,n/2)) # binary covariate
    eta <- X*X_eff+Z*Z_eff
    # Simulate event times
    exponential_times <- (-log(runif(n=n))/(basehaz*exp(eta)))
    # Censor at cutoff
    cutoff <- exponential_times[order(exponential_times)[ceiling(n*0.8)]] # Censor all individuals at the 80th percentile of event times
    eventtime <- pmin(exponential_times,cutoff)
    status <- exponential_times<=cutoff
    dat <- data.frame("X" = X,
                      "Z" = Z,
                      "eventtime" = eventtime,
                      "status" = status)
    cox.out.full <- survival::coxph(Surv(eventtime, status) ~ X + Z, data = dat)
    pow.small[i] <- confint(cox.out.full)["X",1]>0 | confint(cox.out.full)["X",2]<0
  }
  return(mean(pow.small))
}
# candidate_params <- list(list("X_eff" = seq(from=0.85,0.87,by=0.005), "n" = 20),
#                          list("X_eff" = seq(from=0.471,0.477,by=0.001), "n" = 50),
#                          list("X_eff" = seq(from=0.322,0.328,by=0.001), "n" = 100),
#                          list("X_eff" = seq(from=0.222,0.227,by=0.001), "n" = 200))
# power_vec <- function(params, B){
#   print(params)
#   pow <- numeric(length(params$X_eff))
#   n <- params$n
#   if(!is.null(params$B)){B <- params$B} else{B <- 10000}
#   Z_eff <- 0.2
#   for (i in 1:length(params$X_eff)){
#     X_eff <- params$X_eff[i]
#     pow[i] <- sim_cox_pow(n,Z_eff = Z_eff,X_eff = X_eff,B=B)
#   }
#   return(list("power" =pow,
#               "param" = params))
# }
# power_list2 <- lapply(candidate_params,power_vec)
# plot_power_list <- function(power_list){
#   for (i in 1:length(power_list)){
#     power_i <- power_list[[i]]$power
#     X_eff_i <- power_list[[i]]$param$X_eff
#     plot(power_i~X_eff_i)
#     lm_i <- lm(power_i~X_eff_i)
#     abline(lm_i)
#     abline(h=0.8,lty=2)
#     X_eff_pred <- (0.8-lm_i$coefficients[1])/lm_i$coefficients[2]
#     abline(v=X_eff_pred)
#     print(paste0(i,": n = ",power_list[[i]]$param$n,", predicted X_eff = ", round(X_eff_pred,5)))
#   }
# }
# plot_power_list(power_list2)

n_obss <- c(20,50,100,200) # 20,50,100,200,500
print(n_obss)
betas <- c(0.865,0.475,0.325,0.225) # 0.872,0.475,0.324,0.225,0.141
length(betas)==length(n_obss)
prop_misses <- c(1-10/11,1-10/14,1-10/20,1-10/35,1-10/100)
parameters_list <- list()
{count <- 0
  for (i in 1:length(betas)){
    for (j in 1:length(prop_misses)){
      count <- count + 1
      parameters_list[[count]] <- list("n_obs" = n_obss[i], "X_eff" = betas[i], "prop_mis" = prop_misses[j], "Z_eff" = 0.2,
                                       "basehaz" = 1)
    }
  }
}
paste0("Number of scenarios: ",length(parameters_list))
paste0("Largest scenario: ",parameters_list[length(parameters_list)])
paste0("Smallest scenario: ",parameters_list[1])

# Check which scenarios in overall study are included in this one
n_obss_simstudy <- c(10,20,50,100,200) # ,500
prop_misses <- c(1-10/11,1-10/14,1-10/20,1-10/35,1-10/100)
nprop_misses <- length(prop_misses)
nobs_groups <- match(n_obss,n_obss_simstudy)
scen <- rep(1:nprop_misses,length(nobs_groups)) + rep(5*(nobs_groups-1),each=nprop_misses)

# Reverse order so simulation which requires largest memory runs first, so if it crashes it does so early.
parameters_list <- rev(parameters_list)

# Random seeds for snowfall
set.seed(45678)                         # master seed for reproducibility
n_scenarios <- length(parameters_list)
scenario_seeds <- sample.int(1e9, n_scenarios)  # reproducible list of scenario seeds
saveRDS(scenario_seeds,paste0(output_directory,"/","scenario_seeds_",name,".RDS"))

sample_ABpool_norm <- function(out,coeff_name){
  B <- length(out)
  MLE <- numeric(B)
  var_vec <- numeric(B)
  samples <- numeric(B)
  MLE <- sapply(out,function(x) x$coefficients[coeff_name])
  var_vec <- sapply(out,function(x) vcov(x)[coeff_name,coeff_name])
  samples <- rnorm(n=B,mean = MLE,sd=sqrt(var_vec))
  return(samples)
}

sample_ABpool_t_Cox <- function(out,coeff_name, df=NULL){
  B <- length(out)
  MLE <- numeric(B)
  var_vec <- numeric(B)
  samples <- numeric(B)
  if (is.null(df)) {df <- out[[1]]$nevent-length(out[[1]]$coefficients)} 
  MLE <- sapply(out,function(x) x$coefficients[coeff_name])
  var_vec <- sapply(out,function(x) vcov(x)[coeff_name,coeff_name])
  samples <- rt(n=B,df=df)*sqrt(var_vec)+MLE
  return(samples)
}

# Simulation study
m <- 10000
methnames <- c("full","Bayes_100","Bayes_1000","Bayes_10000","CC_ML",
               "Rubin_5","Rubin_20","Rubin_100","Rubin_1000","Rubin_norm_1000",
               "ABpool_norm_100","ABpool_norm_1000","ABpool_norm_10000",
               "ABpool_t_100","ABpool_t_1000","ABpool_t_10000")
results <- data.frame("2.5 %" = rep(0,length(methnames)),
                      "97.5 %"= rep(0,length(methnames)),
                      "CI_width"= rep(0,length(methnames)),
                      "power"= rep(0,length(methnames)),
                      "coverage"= rep(0,length(methnames)),
                      "miss_low"= rep(0,length(methnames)),
                      row.names = methnames)
colnames(results) <- c("2.5 %","97.5 %","CI_width","power","coverage","miss_low")
results_i <- results
nworkers <- length(parallelly::availableWorkers())
# Run in parallel using mclapply or something.

cox_regression_sims <- function(B, parameters_in){
  n_obs <- parameters_in$n_obs
  X_eff <- parameters_in$X_eff
  Z_eff <- parameters_in$Z_eff
  shape <- parameters_in$shape
  scale <- parameters_in$scale
  cutoff <- parameters_in$cutoff
  prop_mis <- parameters_in$prop_mis
  basehaz <- parameters_in$basehaz
  n <- round(n_obs/(1-prop_mis))
  parameters_in$n <- n
  
  runtimes <- matrix(nrow=B,ncol=8)
  lb_95_CI <- matrix(nrow=B,ncol=nrow(results_i), dimnames = list(1:B,rownames(results_i))) # Note results_i is slightly naughtily defined outside the function
  ub_95_CI <- matrix(nrow=B,ncol=nrow(results_i), dimnames = list(1:B,rownames(results_i)))
  CI_width <- matrix(nrow=B,ncol=nrow(results_i), dimnames = list(1:B,rownames(results_i)))
  
  other_info_types <- c("ess_bulk","ess_tail","fmi","df.com","df.barnard","df_rubin_norm",
                        "diffmedian_Bayes","bias_Bayes","bias_MI",
                        "meanbias_X_mis","bias_beta_Z",
                        "skew_MLE","kurt_MLE","Cov_MLE_U","Var_U","Cov_MLEsq_U",
                        "var_par","var_Rubin","skew_term_MLE","skew_term_Cov","skew_calc",
                        "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc",
                        "skew_par","kurt_par","kurt_rubin",
                        "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                        "CF_approx_t","CI_diff_t","divergences")
  other_info <- matrix(nrow=B,ncol=length(other_info_types))
  colnames(other_info) <- other_info_types
  
  key_methods <- c("Bayes_10000","Rubin_1000","Rubin_norm_1000","ABpool_t_10000","ABpool_norm_10000")
  key_results <- c("2.5 %","97.5 %","power","coverage")
  sim_results_array <- array(dim=c(B,length(key_methods),length(key_results)), dimnames=list(NULL,key_methods,key_results))
  
  count <- 1
  error_count <- 0
  CF_NA_count <- 0
  while (count <= B){
    t1 <- Sys.time()
    
    print(count)
    X_true <- rnorm(n) # continuous covariate
    if (n%%2==0){
      Z <- c(rep(1,n/2),rep(0,n/2)) # binary covariate (even n)
    } else if (n%%2==1){
      Z <- c(rep(1,(n+1)/2),rep(0,(n-1)/2)) # binary covariate (odd n)
    }
    eta <- X_true*X_eff+Z*Z_eff
    # Simulate event times
    exponential_times <- (-log(runif(n=n))/(basehaz*exp(eta)))
    # Censor at cutoff
    cutoff <- exponential_times[order(exponential_times)[ceiling(n*0.8)]] # Censor all individuals at the 80th percentile of event times
    eventtime <- pmin(exponential_times,cutoff)
    status <- exponential_times<=cutoff
    data_full <- data.frame("X" = X_true,
                            "Z" = Z,
                            "eventtime" = eventtime,
                            "status" = status)
    cox.out.full <- tryCatch(survival::coxph(Surv(eventtime, status) ~ X + Z, data = data_full),
                             error = function(e) e)
    if(!inherits(cox.out.full, "error")){
      X_obs <- X_true
      mis <- sample(n,n*prop_mis)
      X_obs[mis] <- NA
      obs <- which(!is.na(X_obs))
      # Fully observed data analysis
      data_i <- data.frame("X" = X_obs,
                           "Z" = Z,
                           "eventtime" = eventtime,
                           "status" = status)
      t2 <- Sys.time()
      runtimes[count,1] <- difftime(t2,t1, units = "secs")
      
      data_stan <- list("t"=eventtime,
                        "event"=status,
                        "X_obs"=X_true[obs],
                        "Z"=Z,
                        "obs"=array(obs),
                        "mis"=array(mis),
                        "n"=length(eventtime),
                        "n_obs"=length(obs),
                        "n_mis"=length(mis))
      warmup_iterations = 1e3
      sampling_iterations = m/4 + warmup_iterations #best to use 1e3 or higher
      
      # Stan sampling step - run in parallel with 4 cores (can be done on a local computer but will take a very long time with this many iterations)
      t_stan_1 <- Sys.time()
      stan_obj = rstan::sampling(
        object = exponential_model,
        data = data_stan,
        chains = 4,
        cores = 1,
        iter = sampling_iterations,
        warmup = warmup_iterations,
        control = list(adapt_delta = 0.97, max_treedepth = 15),
        refresh = sampling_iterations/10)
      t_stan_2 <- Sys.time()
      difftime(t_stan_2,t_stan_1)
      
      # Count number of divergences
      sp <- rstan::get_sampler_params(stan_obj, inc_warmup = FALSE)
      divergences_per_chain <- sapply(sp, function(chain) {
        sum(chain[, "divergent__"])
      })
      total_divergences <- sum(divergences_per_chain)
      rm(sp)
      
      stan_out <- rstan::extract(stan_obj,permuted=T)
      meanbias_X_mis <- mean(sweep(stan_out$X_mis,2,X_true[mis],FUN="-"))
      bias_betaZ <- mean(stan_out$beta_Z) - Z_eff
      # bias_beta0 <- mean(stan_out$beta_0) - bias_betaZ/2 - 0 # Correct for centering in stan
      # bias_sigma <- mean(stan_out$sigma) - 1
      stan_out <- stan_out[c("beta_X","X_mis")] # Keep only beta_X and X_mis
      rm(stan_obj)
      t2 <- Sys.time()
      runtimes[count,1] <- difftime(t2,t1, units = "secs")
      
      ess_b <- rstan::ess_bulk(stan_out$X_mis[,1])
      ess_t <- rstan::ess_tail(stan_out$X_mis[,1])
      
      # imp_data_i <- tryCatch(mice(data_i, predictorMatrix = pred_mat,
      #                             m = m, method = "norm", # Note m is slightly naughtily defined outside the function
      #                             printFlag = F, maxit=1),
      #                        error = function(e) e) #maxit=1 because there is only one variable with missing data and the posterior can be sampled from exactly due to conjugate priors.
      
      {
        t3 <- Sys.time()
        runtimes[count,2] <- difftime(t3,t2, units = "secs")
        
        cox.imp <- as.list(numeric(m))
        for (k in 1:m){
          X <- X_true
          X[mis] <- stan_out$X_mis[k,]
          cox.imp[[k]] <- survival::coxph(Surv(eventtime, status) ~ X + Z)
        }
        Bayes_beta <- stan_out$beta_X
        rm(stan_out)
        gc()
        
        t4 <- Sys.time()
        runtimes[count,3] <- difftime(t4,t3, units = "secs")
        
        cox.imp5 <- cox.imp[round(seq(1,m,length.out=5))]
        cox.imp20 <- cox.imp[round(seq(1,m,length.out=20))]
        cox.imp100 <- cox.imp[round(seq(1,m,length.out=100))]
        cox.imp1000 <- cox.imp[round(seq(1,m,length.out=1000))]
        
        t5 <- Sys.time()
        runtimes[count,4] <- difftime(t5,t4, units = "secs")
        
        # Rubin's rules
        pool.imp5 <- pool(cox.imp5)
        pool.imp20 <- pool(cox.imp20)
        pool.imp100 <- pool(cox.imp100)
        
        t6 <- Sys.time()
        runtimes[count,5] <- difftime(t6,t5, units = "secs")
        
        pool.imp1000 <- pool(cox.imp1000)
        rm(cox.imp1000)
        
        t7 <- Sys.time()
        runtimes[count,6] <- difftime(t7,t6, units = "secs")
        # Sampling method
        ABpool_samples_norm <- sample_ABpool_norm(out=cox.imp,coeff_name = "X")
        ABpool_samples_t <- sample_ABpool_t_Cox(out=cox.imp,coeff_name = "X")
        
        # Skewness and kurtosis
        MLEs <- sapply(cox.imp,function(x) x$coefficients["X"])
        Us <- sapply(cox.imp,function(x) vcov(x)["X","X"])
        rm(cox.imp)
        
        t8 <- Sys.time()
        runtimes[count,7] <- difftime(t8,t7, units = "secs")
        
        # Complete case analysis
        data_obs <- data.frame("eventtime" = eventtime[obs],
                               "status" = status[obs],
                               "X" = X_true[obs],
                               "Z" = Z[obs])
        cox.out <- survival::coxph(Surv(eventtime, status) ~ X + Z, data = data_obs)
        
        df_rubin_norm <- (pool.imp1000$m-1)/pool.imp1000$pooled$fmi[1]^2
        df_rubin_t <- pool.imp1000$pooled$df[1]
        sd_rubin <- sqrt(pool.imp1000$pooled$t[1])
        mean_rubin <- pool.imp1000$pooled$estimate[1]
        
        results_i[,1:2] <- rbind(confint(cox.out.full)[1,],
                                 quantile(Bayes_beta[1:100],c(0.025,0.975)),
                                 quantile(Bayes_beta[1:1000],c(0.025,0.975)),
                                 quantile(Bayes_beta[1:10000],c(0.025,0.975)),
                                 confint(cox.out)[1,],
                                 summary(pool.imp5, conf.int=T)[1,7:8],
                                 summary(pool.imp20, conf.int=T)[1,7:8],
                                 summary(pool.imp100, conf.int=T)[1,7:8],
                                 summary(pool.imp1000, conf.int=T)[1,7:8],
                                 mean_rubin + c(-1,1)*qt(0.975,df=df_rubin_norm)*sd_rubin,
                                 quantile(ABpool_samples_norm[1:100],c(0.025,0.975)),
                                 quantile(ABpool_samples_norm[1:1000],c(0.025,0.975)),
                                 quantile(ABpool_samples_norm[1:10000],c(0.025,0.975)),
                                 quantile(ABpool_samples_t[1:100],c(0.025,0.975)),
                                 quantile(ABpool_samples_t[1:1000],c(0.025,0.975)),
                                 quantile(ABpool_samples_t[1:10000],c(0.025,0.975)))
        results_i[,"CI_width"] <- results_i[,"97.5 %"] - results_i[,"2.5 %"]
        results_i[,"power"] <- results_i[,"97.5 %"] < 0 | results_i[,"2.5 %"] > 0
        results_i[,"coverage"] <- results_i[,"97.5 %"] > X_eff & results_i[,"2.5 %"] < X_eff
        results_i[,"miss_low"] <- results_i[,"97.5 %"] < X_eff
        
        if(any(is.na(c(as.matrix(results_i))))) {error_count <- error_count + 1; next} # Skip the iterations if any results are NA
        
        results <- results + results_i
        lb_95_CI[count,] <- results_i[,"2.5 %"]
        ub_95_CI[count,] <- results_i[,"97.5 %"]
        CI_width[count,] <- results_i[,"CI_width"]
        
        sim_results_array[count,,] <- as.matrix(results_i[key_methods,key_results])
        
        # Calculate quantities for other info
        b_m <- var(MLEs)
        u_m <- mean(Us)
        nu <- pool.imp1000$pooled$dfcom[1]
        df_rubin <- pool.imp1000$pooled$df[1]
        nu_nu_2 <- nu/(nu-2)
        total_var <- b_m + nu_nu_2*u_m
        lambda_nu <- exp(log(b_m)-log(total_var))
        z_a <- qnorm(0.025)
        Cov_MLEsq_U <- cov((MLEs-mean(MLEs))^2,Us)
        sd_rubin <- sqrt((pool.imp1000$pooled$ubar[1]+pool.imp1000$pooled$b[1])*df_rubin/(df_rubin-2))
        sd_par <- sd(ABpool_samples_t)
        skew_par <- skewness(ABpool_samples_t)
        kurt_par <- kurtosis(ABpool_samples_t)
        CF_t_terms <- c(z_a*(sd_par-sd_rubin),
                        skew_par*sd_par*(z_a^2-1)/6,
                        (kurt_par*sd_par - sd_rubin*6/(df_rubin-4))*(z_a^3-3*z_a)/24,
                        -skew_par^2*sd_par*(2*z_a^3-5*z_a)/36)
        skew_calc_terms <- c(skewness(MLEs)*lambda_nu^(3/2),
                             3*nu_nu_2*cov(MLEs,Us)/total_var^(3/2))
        kurt_calc_terms <- c(kurtosis(MLEs)*lambda_nu^2,
                             6*(1-lambda_nu)^2/(nu-4),
                             3*nu^2/((nu-2)*(nu-4))*var(Us)/total_var^2,
                             6*nu_nu_2*Cov_MLEsq_U/total_var^2)
        
        CI_diff_t <- quantile(ABpool_samples_t,c(0.025))-summary(pool.imp1000, conf.int=T)[1,7]
        if (df_rubin<=4 | nu <= 4){CF_t_terms <- rep(NaN,4); kurt_calc_terms <- rep(NaN,4); skew_calc_terms <- rep(NaN,3); CI_diff_t <- NA; CF_NA_count <- CF_NA_count + 1} # If df < 4 set these terms to NA and count
        
        other_info[count,] <- c(ess_b,
                                ess_t,
                                pool.imp1000$pooled$fmi[1],
                                nu,
                                df_rubin,
                                df_rubin_norm,
                                median(Bayes_beta)-X_eff,
                                mean(Bayes_beta)-X_eff,
                                mean(MLEs)-X_eff,
                                meanbias_X_mis,
                                bias_betaZ,
                                skewness(MLEs),
                                kurtosis(MLEs), # excess kurtosis
                                cov(MLEs,Us),
                                var(Us),
                                Cov_MLEsq_U,
                                sd_par^2,
                                sd_rubin^2,
                                skew_calc_terms[1],
                                skew_calc_terms[2],
                                sum(skew_calc_terms),
                                kurt_calc_terms[1],
                                kurt_calc_terms[2],
                                kurt_calc_terms[3],
                                kurt_calc_terms[4],
                                sum(kurt_calc_terms),
                                skew_par,
                                kurt_par,
                                6/(df_rubin-4),
                                CF_t_terms[1],
                                CF_t_terms[2],
                                CF_t_terms[3],
                                CF_t_terms[4],
                                sum(CF_t_terms),
                                CI_diff_t,
                                total_divergences)
        
        t9 <- Sys.time()
        runtimes[count,8] <- difftime(t9,t8, units = "secs")
        count <- count + 1
      } 
    gc()
  } else {error_count <- error_count + 1}} 
  results <- results/B
  lb_95_CI_diff_Bayes <- abs(sweep(lb_95_CI[,key_methods],1,lb_95_CI[,"Bayes_10000"],FUN="-"))
  ub_95_CI_diff_Bayes <- abs(sweep(ub_95_CI[,key_methods],1,ub_95_CI[,"Bayes_10000"],FUN="-"))
  CI_total_diff_Bayes <- lb_95_CI_diff_Bayes + ub_95_CI_diff_Bayes
  return(list("results" = results, "CI_total_diff_Bayes" = CI_total_diff_Bayes, "lb_95_CI" = lb_95_CI, "ub_95_CI" = ub_95_CI, "CI_width" = CI_width, "parameters" = parameters_in, "other_info" = other_info,  "lb_95_CI_diff_Bayes" = lb_95_CI_diff_Bayes, "ub_95_CI_diff_Bayes" = ub_95_CI_diff_Bayes, "runtimes" = runtimes, "B" = B, "error_count" = error_count, "CF_NA_count" = CF_NA_count, "sim_results_array" = sim_results_array))
}
# # Testing the time it will take
# out_test <- cox_regression_sims(2,parameters_list[[length(parameters_list)]])
# parameters_in <- parameters_list[[rev(1:20)[5]]]

B <- 5000
nwork <- nworkers
out_list <- list()
for (i in 1:length(parameters_list)){
  parameters_in <- parameters_list[[i]]
  B_small <- floor(B/nwork)
  B_seq <- rep(list(B_small),nwork)
  B_diff <- Reduce("+",B_seq)
  if (B-B_diff!=0){B_seq[1:(B-B_diff)] <- B_small+1} # Splits B into a vector of length nwork which sums to B
  
  # Begin a cluster in snowfall package
  cluster <- snowfall::sfInit(nwork,type="SOCK", parallel=T)
  
  # Load the relevant objects in the parallel workspaces - including 
  snowfall::sfExport("B_seq","m","results","results_i",
                     "sample_ABpool_norm","sample_ABpool_t_Cox","exponential_model",
                     "skewness","kurtosis")
  
  # Load libraries
  snowfall::sfLibrary(mice)
  snowfall::sfLibrary(parallel)
  snowfall::sfLibrary(survival)
  snowfall::sfLibrary(abind)
  
  # Set up the random number generation
  snowfall::sfClusterSetupRNG(seed=scenario_seeds[i])
  print(paste("Number of workers:",nwork))
  
  # Runs the function on each element of the list in parallel
  out_sublist <- snowfall::sfLapply(x = B_seq, fun = cox_regression_sims, parameters_in = parameters_in)
  snowfall::sfStop()
  
  my_quantile <- function(x, na.rm = FALSE ){ c("mean"=mean(x),quantile(x,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975), na.rm = na.rm)) } # Get quantiles and mean
  # Bind the list
  out <- list()
  out$results <- Reduce("+",lapply(out_sublist,FUN=function(x){x$results*x$B}))/Reduce("+",lapply(out_sublist,FUN=function(x){x$B}))
  out$CI_total_diff_Bayes <-  t(apply( do.call(rbind,lapply(out_sublist,FUN=function(x){x$CI_total_diff_Bayes})) ,2,my_quantile)) # The do.call and things inside create an rbind matrix of the CI_total_diff_Bayes matrices. The apply creates a smaller matrix of the quantiles of the columns.
  out$lb_95_CI_diff_Bayes <-  t(apply( do.call(rbind,lapply(out_sublist,FUN=function(x){x$lb_95_CI_diff_Bayes})) ,2,my_quantile)) # The do.call and things inside create an rbind matrix of the lb_95_CI_diff_Bayes matrices. The apply creates a smaller matrix of the quantiles of the columns.
  out$ub_95_CI_diff_Bayes <-  t(apply( do.call(rbind,lapply(out_sublist,FUN=function(x){x$ub_95_CI_diff_Bayes})) ,2,my_quantile)) # The do.call and things inside create an rbind matrix of the ub_95_CI_diff_Bayes matrices. The apply creates a smaller matrix of the quantiles of the columns.
  out$lb_95_CI <-  t(apply( do.call(rbind,lapply(out_sublist,FUN=function(x){x$lb_95_CI})) ,2,my_quantile)) # The do.call and things inside create an rbind matrix of the lb_95_CI matrices. The apply creates a smaller matrix of the quantiles of the columns.
  out$ub_95_CI <-  t(apply( do.call(rbind,lapply(out_sublist,FUN=function(x){x$ub_95_CI})) ,2,my_quantile))
  out$CI_width <-  t(apply( do.call(rbind,lapply(out_sublist,FUN=function(x){x$CI_width})) ,2,my_quantile))
  out$parameters <- out_sublist[[1]]$parameters
  out$sim_results_array <- do.call(abind,c(lapply(out_sublist,FUN=function(x){x$sim_results_array}),along=1))
  out$other_info <- do.call(rbind,lapply(out_sublist,FUN=function(x){x$other_info}))
  out$other_info_summary <- rbind(colMeans(out$other_info,na.rm=TRUE),apply(out$other_info,2,my_quantile,na.rm=TRUE))
  runtimes <-  do.call(rbind,lapply(out_sublist,FUN=function(x){x$runtimes}))
  out$runtimes <- colMeans(runtimes)
  out$runtimes_lb_95 <- apply(runtimes,2,quantile,0.025)
  out$runtimes_ub_95 <- apply(runtimes,2,quantile,0.975)
  out$error_count <- Reduce("+",lapply(out_sublist,FUN=function(x){x$error_count}))
  out$CF_NA_count <- Reduce("+",lapply(out_sublist,FUN=function(x){x$CF_NA_count}))
  out$B <- Reduce("+",lapply(out_sublist,FUN=function(x){x$B}))
  out$stan_model <- exponential_model
  rm(runtimes,out_sublist)
  
  out_list[[i]] <- out
}

out_list <- rev(out_list) # un-reverse the simulation scenario order.
names(out_list) <- scen # Name the scenarios

saveRDS(out_list,paste0(output_directory,"/",name,".RDS"))
t_end <- Sys.time()
print(difftime(t_end,t_init))