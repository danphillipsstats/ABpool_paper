# Data

# The data is available in the supplementary material of the paper:
# Plummer, M. Cuts in Bayesian graphical models. Stat Comput 25, 37–43 (2015). 
# https://doi.org/10.1007/s11222-014-9503-z
# The data can be extracted from the file using the following links: 
# https://static-content.springer.com/esm/art%3A10.1007%2Fs11222-014-9503-z/MediaObjects/11222_2014_9503_MOESM2_ESM.rnw
# The data is defined in Section 1, at the bottom of page 1 of the supplementary material.
# The vectors defining the variables nhpv, Npart, ncases, and Npop can be copied directly
# from the code given in the supplementary material above.

#########################
# Load packages
require(mice)
require(mvtnorm)
require(moments)
require(rstan)
require(bayesplot)
require(devtools)
devtools::install_github("danphillipsstats/abpool")
require(abpool)

HPV_directory <- paste0(getwd(),"/HPV")
output_directory <- paste0(getwd(),"/HPV/Output")
plot_directory <- paste0(output_directory,"/Plots")

outcome_cols <- palette("Okabe-Ito")[c(7,6)]

# Load functions
source(paste0(getwd(),"/functions.R"))

###################################################
# Latent variable imputation
# For the model from Plummer (2015)

# The data is available in the supplementary material of the paper:
# Plummer, M. Cuts in Bayesian graphical models. Stat Comput 25, 37–43 (2015). 
# https://doi.org/10.1007/s11222-014-9503-z
# The data can be extracted from the file using the following links: 
# https://static-content.springer.com/esm/art%3A10.1007%2Fs11222-014-9503-z/MediaObjects/11222_2014_9503_MOESM2_ESM.rnw
# The data is defined in Section 1, at the bottom of page 1 of the supplementary material.
# The vectors defining the variables nhpv, Npart, ncases, and Npop can be copied directly
# from the code given in the supplementary material above.

# The data must first be extracted from the above link, 
# defining the variables nhpv, Npart, ncases, and Npop.
HPV_cancer_data <- data.frame("Z" = nhpv,
                          "N" = Npart,
                          "Y" = ncases,
                          "Time" = Npop)
HPV_cancer_data$Time <- HPV_cancer_data$Time*1e-3 # Convert time to smaller value so theta_0 is a good scale

###
# Sample from Stan model.
data_hpv <- list("n" = nrow(HPV_cancer_data),
                 "Z" = HPV_cancer_data[,"Z"],
                 "N" = HPV_cancer_data[,"N"],
                 "Y" = HPV_cancer_data[,"Y"],
                 "t" = HPV_cancer_data[,"Time"])

t1 <- Sys.time()
hpv_impute_model = rstan::stan_model(file=paste0(HPV_directory,'/HPV_cancer_cut.stan'))
t2 <- Sys.time()
time_Plummer_stan_compile_impute <- difftime(t2,t1)
t1 <- Sys.time()
stage2_model = rstan::stan_model(file=paste0(HPV_directory,'/HPV_cancer_stage_two.stan'))
t2 <- Sys.time()
time_Plummer_stan_compile_BayesMI <- difftime(t2,t1)

set.seed(1)

m <- 10000
warmup_iterations = 1e3
nchains <- 4
sampling_iterations = m/nchains + warmup_iterations #best to use 1e3 or higher
t1 <- Sys.time()
impute_out = rstan::sampling(
  object = hpv_impute_model,
  data = data_hpv,
  chains = nchains,
  cores = nchains,
  iter = sampling_iterations,
  warmup = warmup_iterations,
  refresh = sampling_iterations/10, #show an update @ each %10
  seed = 1000)
t2 <- Sys.time()
time_Plummer_impute <- difftime(t2,t1)

impute_summary <- rstan::summary(
  impute_out,
  pars = "phi"
)$summary

impute_summary[, c("n_eff", "Rhat")] # n_eff > m so draws are independent

# pairs(impute_out, pars = c("phi[1]","phi[13]","lp__"))
par(mfrow=c(1,1))

rstan::traceplot(impute_out, pars=c("phi[1]","phi[13]"))

extract_impute_out <- rstan::extract(impute_out, pars = "phi", permuted = FALSE)
extract_impute_out <- apply(extract_impute_out,3,c)

range(apply(extract_impute_out,2,skewness))

## Stage 2
# MLE based methods
data_stage2 <- data.frame("Y" = data_hpv$Y,
                          "t" = data_hpv$t,
                          "phi" = numeric(data_hpv$n))
times <- data_hpv$t
Y <- data_hpv$Y
glm_theta <- matrix(nrow=m,ncol=2)
mean_mu <- numeric(m)
n <- data_hpv$n

set.seed(1)

t1 <- Sys.time()
glms <- list()
for (i in 1:m){
  phi_i <- extract_impute_out[i,]
  X_mat <- model.matrix(Y~phi_i)
  glm_i <- glm(Y~phi_i,offset=log(times),family=poisson())
  mean_mu[i] <- mean(predict.glm(glm_i,type="response"))
  glm_theta[i,] <- glm_i$coeff
  glms[[i]] <- glm_i
}
t2 <- Sys.time()
(time_Plummer_MLE <- difftime(t2,t1))
qqnorm(glm_theta[,1], main = "MLEs for theta_1 after imputation")
qqline(glm_theta[,1])
qqnorm(glm_theta[,2], main = "MLEs for theta_2 after imputation")
qqline(glm_theta[,2])
# MLEs are not normally distributed

set.seed(1)

# Diagnostic plot
HPV_abpool_200 <- abpool::abpool(glms[1:200], parameters = 2, dfcom = Inf)

png(paste0(plot_directory,"/Plummer_QQ_diagnostic.png"),width = 1400,height=1000, pointsize=23)
abpool::qqplot_rubin_abpool(HPV_abpool_200)
dev.off()

# Approximate Bayesian multiple imputation - dfcom = inf
require(mvtnorm)
sample_ABpool_norm_Plummer <- function(index,model.output.list=glms){
  model.output <- model.output.list[[index]]
  MLE <- model.output$coeff
  var_mat <- vcov(model.output)
  return(mvtnorm::rmvnorm(n=1,mean=MLE,sigma = var_mat))
}

set.seed(1)

t1 <- Sys.time()
ABpool_samples_norm <- t(vapply(1:m,FUN=sample_ABpool_norm_Plummer,FUN.VALUE=numeric(2)))
t2 <- Sys.time()
(ABpool_thetas_time <- difftime(t2,t1))
qqnorm(ABpool_samples_norm[,1])
qqline(ABpool_samples_norm[,1])
qqnorm(ABpool_samples_norm[,2])
qqline(ABpool_samples_norm[,2])

# Approximate Bayesian multiple imputation - dfcom = n - p
sample_ABpool_t_Plummer <- function(index,model.output.list=glms, n = glms[[1]]$data$n){
  model.output <- model.output.list[[index]]
  MLE <- model.output$coeff
  var_mat <- vcov(model.output)
  dfcom <- nobs(model.output) - length(MLE)
  return(mvtnorm::rmvt(n=1, sigma = var_mat, df = dfcom)+MLE) 
}
t1 <- Sys.time()
ABpool_samples_t <- t(vapply(1:m,FUN=sample_ABpool_t_Plummer,FUN.VALUE=numeric(2)))
t2 <- Sys.time()
(time_Plummer_ABpool <- difftime(t2,t1))
qqnorm(ABpool_samples_t[,1])
qqline(ABpool_samples_t[,1])
qqnorm(ABpool_samples_t[,2])
qqline(ABpool_samples_t[,2])

# Rubin's rules
t1 <- Sys.time()
glm_var <- t(vapply(glms,vcov, FUN.VALUE=numeric(4)))
rubin_norm <- mice::pool(glms, dfcom=Inf)
rubin_t <- mice::pool(glms, dfcom=glms[[1]]$df.residual) # Note mice::pool(glms) gives the same result.
t2 <- Sys.time()
(time_Plummer_Rubin <- difftime(t2,t1))
rubin_t$pooled$estimate
rubin_t$pooled$t
# df
rubin_t$pooled$df # 0.42, 0.18
# fmi
rubin_t$pooled$b/(rubin_t$pooled$ubar*rubin_t$pooled$dfcom/(rubin_t$pooled$dfcom-2)+rubin_t$pooled$b) # 0.95, 0.98
# U_bar
rubin_t$pooled$ubar # 0.000090, 0.13
# b
rubin_t$pooled$b # 0.020, 6.3

## Bayesian multiple imputation
# Check how many iterations are needed

set.seed(1)

data_stage2 <- list("n" = nrow(HPV_cancer_data),
                    "Y" = HPV_cancer_data[,"Y"],
                    "t" = HPV_cancer_data[,"Time"],
                    "phi" = numeric(n))
data_stage2$phi <- extract_impute_out[1,]
stage_2_out = rstan::sampling(
  object = stage2_model,
  data = data_stage2,
  chains = nchains,
  cores = nchains,
  iter = sampling_iterations,
  warmup = warmup_iterations,
  refresh = sampling_iterations/10, #show an update @ each %10
  seed = 1000 + i)

rstan::traceplot(stage_2_out,pars=c("theta_1_std","theta_2_std"), inc_warmup=TRUE,window=c(1,500))
rstan::traceplot(stage_2_out,pars=c("theta_1","theta_2"), inc_warmup=TRUE,window=c(1,500))
# Looks like 100 iterations is enough, say 200 to be safe.
BayesMI_thetas <- matrix(nrow=m,ncol=2)
warmup_iterations <- 200
sampling_iterations <- 1 + warmup_iterations
nchains <- 1
t1 <- Sys.time()
for (i in 1:m){
  data_stage2$phi <- extract_impute_out[i,]
  stage_2_out = rstan::sampling(
    object = stage2_model,
    data = data_stage2,
    chains = nchains,
    cores = nchains,
    iter = sampling_iterations,
    warmup = warmup_iterations,
    refresh = 0, #show an update @ each %10
    seed = 1000)
  stage_2_samples <- rstan::extract(stage_2_out, permuted = FALSE)
  BayesMI_thetas[i,] <- stage_2_samples[1,1,c("theta_1","theta_2")]
}
t2 <- Sys.time()
saveRDS(BayesMI_thetas,paste0(output_directory,"/Plummer_BayesMI_thetas.RDS"))
BayesMI_thetas <- readRDS(paste0(output_directory,"/Plummer_BayesMI_thetas.RDS"))
(time_Plummer_BayesMI <- difftime(t2,t1))

CIs_theta_1 <- rbind(quantile(BayesMI_thetas[,1],c(0.025,0.975)),
                     quantile(ABpool_samples_norm[,1],c(0.025,0.975)),
                     quantile(ABpool_samples_t[,1],c(0.025,0.975)),
                     rubin_norm$pooled$estimate[1] + c(-1,1)*qt(0.975,rubin_norm$pooled$df[1])*sqrt(rubin_norm$pooled$t[1]),
                     rubin_t$pooled$estimate[1] + c(-1,1)*qt(p=0.975, df=rubin_t$pooled$df[1])*sqrt(rubin_t$pooled$t[1]))
results_theta_1 <- cbind("Posterior mean"=c(mean(BayesMI_thetas[,1]),
                                            mean(ABpool_samples_norm[,1]),
                                            mean(ABpool_samples_t[,1]),
                                            rubin_norm$pooled$estimate[1],
                                            rubin_t$pooled$estimate[1]),
                         "Posterior median"=c(median(BayesMI_thetas[,1]),
                                              median(ABpool_samples_norm[,1]),
                                              median(ABpool_samples_t[,1]),
                                              rubin_norm$pooled$estimate[1],
                                              rubin_t$pooled$estimate[1]),
                         CIs_theta_1)
rownames(results_theta_1) <- c("Bayes_MI","ABpool_norm","ABpool_t","Rubin_norm","Rubin_t")
write.csv(round(results_theta_1,2),paste0(output_directory,"/Plummer_results_theta_1.csv"),row.names=T)

CIs_theta_2 <- rbind(quantile(BayesMI_thetas[,2],c(0.025,0.975)),
                     quantile(ABpool_samples_norm[,2],c(0.025,0.975)),
                     quantile(ABpool_samples_t[,2],c(0.025,0.975)),
                     rubin_norm$pooled$estimate[2] + c(-1,1)*qt(0.975,rubin_norm$pooled$df[2])*sqrt(rubin_norm$pooled$t[2]),
                     rubin_t$pooled$estimate[2] + c(-1,1)*qt(p=0.975, df=rubin_t$pooled$df[2])*sqrt(rubin_t$pooled$t[2]))
results_theta_2 <- cbind("Posterior mean"=c(mean(BayesMI_thetas[,2]),
                                            mean(ABpool_samples_norm[,2]),
                                            mean(ABpool_samples_t[,2]),
                                            rubin_norm$pooled$estimate[2],
                                            rubin_t$pooled$estimate[2]),
                         "Posterior median"=c(median(BayesMI_thetas[,2]),
                                              median(ABpool_samples_norm[,2]),
                                              median(ABpool_samples_t[,2]),
                                              rubin_norm$pooled$estimate[2],
                                              rubin_t$pooled$estimate[2]),
                         CIs_theta_2)
rownames(results_theta_2) <- c("Bayes_MI","ABpool_norm","ABpool_t","Rubin_norm","Rubin_t")
write.csv(round(results_theta_2,1),paste0(output_directory,"/Plummer_results_theta_2.csv"),row.names=T)

other_info <- round(rbind("Bayespool estimate" = colMeans(BayesMI_thetas),
                          "ABpool & Rubin estimate" = rubin_norm$pooled$estimate, # Equal to the ABpool estimates
                          "dfcom" = rubin_norm$pooled$dfcom,
                          "fmi" = rubin_norm$pooled$fmi,
                          "df" = rubin_norm$pooled$df,
                          "MLEs skewness" = apply(glm_theta, 2, skewness),
                          "MLEs kurtosis" = apply(glm_theta, 2, kurtosis),
                          "Posterior skewness" = apply(BayesMI_thetas, 2, skewness),
                          "Posterior kurtosis" = apply(BayesMI_thetas, 2, kurtosis)),2)
colnames(other_info) <- c("theta_1","theta_2")
write.csv(other_info,paste0(output_directory,"/Plummer_other_info.csv"),row.names=T)

impute_time <- time_Plummer_impute + time_Plummer_stan_compile_impute 
ABpool_time <- time_Plummer_ABpool + time_Plummer_MLE 
Rubin_time <- time_Plummer_Rubin + time_Plummer_MLE 
BayesMI_time <- time_Plummer_BayesMI + time_Plummer_stan_compile_BayesMI 
as.numeric(BayesMI_time,units="secs")/as.numeric(ABpool_time,units="secs")
as.numeric(BayesMI_time+impute_time,units="secs")/as.numeric(ABpool_time+impute_time,units="secs")

times_Plummer <- mget(c("time_Plummer_stan_compile_impute",
                        "time_Plummer_impute",
                        "time_Plummer_MLE",
                        "time_Plummer_ABpool",
                        "time_Plummer_Rubin",
                        "time_Plummer_stan_compile_BayesMI",
                        "time_Plummer_BayesMI"))
saveRDS(times_Plummer,paste0(output_directory,"/Plummer_runtimes.RDS"))
times_Plummer <- readRDS(paste0(output_directory,"/Plummer_runtimes.RDS"))


# Plot of distribution for theta_1
png(paste0(plot_directory,"/Plummer_theta_1.png"),width = 1400,height=1000, pointsize=23)
plot(density(BayesMI_thetas[,1]), main = "Distribution of theta_1",xlab="theta_1",lwd=3,ylim=c(0,3.2))
lines(density(ABpool_samples_norm[,1]), col= outcome_cols[1],lwd=2)
lines(density(ABpool_samples_t[,1]), col= outcome_cols[1], lty=2, lwd=2)
xx <- seq(min(ABpool_samples_norm[,1])-1,max(ABpool_samples_norm[,1])+1,0.0001)
xn <- (xx-rubin_norm$pooled$estimate[1])/sqrt(rubin_norm$pooled$t[1])
xt <- (xx-rubin_t$pooled$estimate[1])/sqrt(rubin_t$pooled$t[1])
tt <- dt(xt,df=rubin_t$pooled$df[1])/sqrt(rubin_t$pooled$t[1]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
nn <- dt(xn,df=rubin_norm$pooled$df[1])/sqrt(rubin_norm$pooled$t[1])
lines(nn~xx, col = outcome_cols[2],lwd=2)
lines(tt~xx, col = outcome_cols[2],lty=2 ,lwd=2)
legend(x=-2.4,y=2.0,
       legend = c("MI_Bayes","ABpool_norm","ABpool_t","Rubin_norm","Rubin_t"),
       col = c("black",rep(outcome_cols,each=2)),
       lty = c(1,1,2,1,2),lwd=c(3,2,2,2,2))
dev.off()

# Similar to  above for theta_2
png(paste0(plot_directory,"/Plummer_theta_2.png"),width = 1400,height=1000, pointsize=23)
plot(density(BayesMI_thetas[,2]), main = "Distribution of theta_2", ylim = c(0,0.18),xlab="theta_2",lwd=3)
lines(density(ABpool_samples_norm[,2]), col= outcome_cols[1],lwd=2)
lines(density(ABpool_samples_t[,2]), col= outcome_cols[1], lty=2,lwd=2)
xx <- seq(min(ABpool_samples_norm[,2])-1,max(ABpool_samples_norm[,2])+1,0.0001)
xn <- (xx-rubin_norm$pooled$estimate[2])/sqrt(rubin_norm$pooled$t[2])
xt <- (xx-rubin_t$pooled$estimate[2])/sqrt(rubin_t$pooled$t[2])
tt <- dt(xt,df=rubin_t$pooled$df[2])/sqrt(rubin_t$pooled$t[2]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
nn <- dt(xn,df=rubin_norm$pooled$df[2])/sqrt(rubin_norm$pooled$t[2])
lines(nn~xx, col = outcome_cols[2],lwd=2)
lines(tt~xx, col = outcome_cols[2],lty=2,lwd=2)
legend(x=22,y=0.11,
       legend = c("MI_Bayes","ABpool_norm","ABpool_t","Rubin_norm","Rubin_t"),
       col = c("black",rep(outcome_cols,each=2)),
       lty = c(1,1,2,1,2),lwd=c(3,2,2,2,2))
dev.off()

# Cornish-Fisher approximation - normal version

MLEs <- glm_theta
Us <- lapply(glms,function(x) vcov(x))

skew_kurt_CF_calc_norm <- skew_kurt_CF_calc(MLEs, Us, rubin_norm, ABpool_samples_norm, term = 2)
write.csv(skew_kurt_CF_calc_norm,paste0(output_directory,"/Plummer_calculations_norm_info_theta_2.csv"),row.names=T)

# Large diagnostic plot
set.seed(1)
HPV_abpool <- abpool::abpool(glms, parameters = 2, dfcom = Inf)

png(paste0(plot_directory,"/Plummer_QQ_diagnostic_large.png"),width = 1400,height=1000, pointsize=23)
abpool::qqplot_rubin_abpool(HPV_abpool)
dev.off()
