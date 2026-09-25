##############################################################################
# Real data examples
##############################################################################
require(mice)
require(mvtnorm)
require(moments)
require(mixpoissonreg)
require(rstan)
require(bayesplot)
require(glmmTMB)
require(rjags)

output_directory <- paste0(getwd(),"/HPV/Output")
plot_directory <- paste0(output_directory,"/Plots")

outcome_cols <- palette("Okabe-Ito")[c(7,6)]

# Functions
bootstrap.cox <- function(data, formula, B=1000){
  cox.main <- coxph(formula, data)
  effects <- matrix(NA,B,length(coef(cox.main)))
  for (i in 1:B){
    data_boot <- data[sample(1:nrow(data),nrow(data),TRUE),]
    cox.boot <- coxph(formula, data=data_boot)
    effects[i,] <- coef(cox.boot)
  }
  return(effects)
}
get.ci.mice.pred.i <- function(i, pool.out){
  pool.out$pooled$estimate[i] + c(-1,1)*qt(p=0.975, df=pool.out$pooled$df[i])*sqrt(pool.out$pooled$t[i])
}
get.ci.mice <- function(pool.out){
  out <- vapply(1:nrow(pool.out$pooled),get.ci.mice.pred.i,numeric(2),pool.out=pool.out)
  colnames(out) <- pool.out$pooled[,1]
  return(out)
}
sample_ABpool_norm <- function(index,mice.with.list){
  out <- mice.with.list$analyses[[index]]
  MLE <- out$coefficients
  var_mat <- vcov(out)
  return(mvtnorm::rmvnorm(n=1,mean=MLE,sigma = var_mat))
}
sample_ABpool_t_cox <- function(index,mice.with.list){
  out <- mice.with.list$analyses[[index]]
  MLE <- out$coefficients
  df <- out$nevent - length(MLE)
  var_mat <- vcov(out)
  return(mvtnorm::rmvt(n=1, sigma = var_mat, df = df)+MLE)
}
nmis <- function(x)sum(is.na(x))
propmis <- function(x)sum(is.na(x))/length(x)

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
skew_kurt_CF_calc <- function(MLEs,Us,rubin_pool,ABpool_samples,distr=c("t","norm"),alpha=0.025,term=1){
  if (is.null(dim(MLEs))){
    MLEs <- matrix(MLEs,ncol=1);
    Us <- matrix(Us,ncol=1);
    ABpool_samples <- matrix(ABpool_samples,ncol=1);
  }
  
  distr <- match.arg(distr)
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
  
  if (distr == "t"){
    b_m <- var(MLEs[,term])
    u_m <- mean(Us[,term])
    nu <- rubin_pool$pooled$dfcom[term]
    df_rubin <- rubin_pool$pooled$df[term]
    nu_nu_2 <- nu/(nu-2)
    total_var <- b_m + nu_nu_2*u_m
    lambda_nu <- exp(log(b_m)-log(total_var))
    z_a <- qnorm(alpha)
    Cov_MLEsq_U <- cov((MLEs[,term]-mean(MLEs[,term]))^2,Us[,term])
    sd_rubin <- sqrt((rubin_pool$pooled$ubar[term]+rubin_pool$pooled$b[term])*df_rubin/(df_rubin-2))
    sd_par <- sd(ABpool_samples[,term])
    skew_par <- skewness(ABpool_samples[,term])
    kurt_par <- kurtosis(ABpool_samples[,term])
    CF_t_terms <- c(z_a*(sd_par-sd_rubin),
                    skew_par*sd_par*(z_a^2-1)/6,
                    (kurt_par*sd_par - sd_rubin*6/(df_rubin-4))*(z_a^3-3*z_a)/24,
                    -skew_par^2*sd_par*(2*z_a^3-5*z_a)/36)
    skew_calc_terms <- c(skewness(MLEs[,term])*lambda_nu^(3/2),
                         3*nu_nu_2*cov(MLEs[,term],Us[,term])/total_var^(3/2))
    kurt_calc_terms <- c(kurtosis(MLEs[,term])*lambda_nu^2,
                         6*(1-lambda_nu)^2/(nu-4),
                         3*nu^2/((nu-2)*(nu-4))*var(Us[,term])/total_var^2,
                         6*nu_nu_2*Cov_MLEsq_U/total_var^2)
    calculations_info <- c(skew_calc_terms,
                           sum(skew_calc_terms),
                           skew_par,
                           kurt_calc_terms,
                           sum(kurt_calc_terms),
                           kurt_par,
                           CF_t_terms,
                           sum(CF_t_terms),
                           quantile(ABpool_samples[,term],c(alpha)) - (rubin_pool$pooled$estimate[term] + qt(alpha,(rubin_pool$m-1)/rubin_pool$pooled$fmi[term]^2)*sqrt(rubin_pool$pooled$t[term])))
    names(calculations_info) <- c("skew_term_MLE","skew_term_Cov","skew_calc",
                                  "skew_par",
                                  "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc","kurt_par",
                                  "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                                  "CF_approx_t","CI_diff_t")
  }
  if (distr == "norm"){
    b_m <- var(MLEs[,term])
    u_m <- mean(Us[,term])
    nu <- rubin_pool$pooled$dfcom[term]
    df_rubin <- rubin_pool$pooled$df[term]
    total_var <- b_m + u_m
    lambda <- exp(log(b_m)-log(total_var))
    z_a <- qnorm(alpha)
    Cov_MLEsq_U <- cov((MLEs[,term]-mean(MLEs[,term]))^2,Us[,term])
    sd_rubin <- sqrt((rubin_pool$pooled$ubar[term]+rubin_pool$pooled$b[term])*df_rubin/(df_rubin-2))
    sd_par <- sd(ABpool_samples[,term])
    skew_par <- skewness(ABpool_samples[,term])
    kurt_par <- kurtosis(ABpool_samples[,term])
    CF_norm_terms <- c(0,
                       skew_par*sd_par*(z_a^2-1)/6,
                       (kurt_par*sd_par)*(z_a^3-3*z_a)/24,
                       -skew_par^2*sd_par*(2*z_a^3-5*z_a)/36)
    skew_calc_terms <- c(skewness(MLEs[,term])*lambda^(3/2),
                         3*cov(MLEs[,term],Us[,term])/total_var^(3/2))
    kurt_calc_terms <- c(kurtosis(MLEs[,term])*lambda^2,
                         0,
                         3*var(Us[,term])/total_var^2,
                         6*Cov_MLEsq_U/total_var^2)
    calculations_info <- c(skew_calc_terms,
                           sum(skew_calc_terms),
                           skew_par,
                           kurt_calc_terms,
                           sum(kurt_calc_terms),
                           kurt_par,
                           CF_norm_terms,
                           sum(CF_norm_terms),
                           quantile(ABpool_samples,c(alpha)) - (rubin_pool$pooled$estimate[term] + qt(alpha,(rubin_pool$m-1)/rubin_pool$pooled$fmi[term]^2)*sqrt(rubin_pool$pooled$t[term])))
    names(calculations_info) <- c("skew_term_MLE","skew_term_Cov","skew_calc",
                                  "skew_par",
                                  "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc","kurt_par",
                                  "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                                  "CF_approx_norm","CI_diff_norm")
  }
  return(calculations_info)
}

###################################################
# Latent variable imputation
# For the model from Plummer (2015)

# The data can be extracted from the following file: 
# https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs11222-014-9503-z/MediaObjects/11222_2014_9503_MOESM1_ESM.pdf
# This is available as the supplementary material of the paper:
# Plummer, M. Cuts in Bayesian graphical models. Stat Comput 25, 37–43 (2015). 
# https://doi.org/10.1007/s11222-014-9503-z

# First download and save the data, then
data_directory <- paste0(dirname(dirname(getwd())),"/2_Data_for_Analysis")
HPV_cancer_data <- read.csv(paste0(data_directory,"/HPV_cancer.csv"))
HPV_cancer_data$Time <- HPV_cancer_data$Time*1e-3 # Convert time to smaller value so theta_0 is a good scale

###
# Sample from Stan model.
data_hpv <- list("n" = nrow(HPV_cancer_data),
                 "Z" = HPV_cancer_data[,"Z"],
                 "N" = HPV_cancer_data[,"N"],
                 "Y" = HPV_cancer_data[,"Y"],
                 "t" = HPV_cancer_data[,"Time"])
wd <- getwd()
# hpv_full_model = rstan::stan_model(file=paste0(wd,'/HPV_cancer_full.stan'))
t1 <- Sys.time()
hpv_impute_model = rstan::stan_model(file=paste0(wd,'/HPV_cancer_cut.stan'))
t2 <- Sys.time()
time_Plummer_stan_compile_impute <- difftime(t2,t1)
t1 <- Sys.time()
stage2_model = rstan::stan_model(file=paste0(wd,'/HPV_cancer_stage_two.stan'))
t2 <- Sys.time()
time_Plummer_stan_compile_BayesMI <- difftime(t2,t1)

m <- 10000
warmup_iterations = 1e3
nchains <- 4
sampling_iterations = m/nchains + warmup_iterations #best to use 1e3 or higher

# full_out = rstan::sampling(
#   object = hpv_full_model,
#   data = data_hpv,
#   chains = nchains,
#   cores = nchains,
#   iter = sampling_iterations,
#   warmup = warmup_iterations,
#   refresh = sampling_iterations/10, #show an update @ each %10
#   seed = 1000)

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

s_impute_out <- rstan::summary(impute_out)$summary
s_impute_phi_means <- s_impute_out[substring(rownames(s_impute_out),1,3)=="phi","mean"]

# pairs(full_out, pars = c("phi[1]","phi[13]","theta_1","theta_2","lp__"))

pairs(impute_out, pars = c("phi[1]","phi[13]","lp__"))
par(mfrow=c(1,1))

rstan::traceplot(impute_out, pars=c("phi[1]","phi[13]"))

# s_full_out <- rstan::summary(full_out)$summary
# s_full_eta_means <- s_full_out[substring(rownames(s_full_out),1,3)=="eta","mean"]
# # cbind(s_impute_eta_means,s_full_eta_means)
# extract_full_out <- rstan::extract(full_out)

round(s_impute_out[,c("mean","sd","2.5%","25%","75%","97.5%")],2)
# round(s_full_out[,c("mean","sd","2.5%","25%","75%","97.5%")],2)

extract_impute_out <- rstan::extract(impute_out)

range(apply(extract_impute_out$phi,2,skewness))

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
  phi_i <- extract_impute_out$phi[i,]
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

# Approximate Bayesian multiple imputation
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

# Approximate Bayesian multiple imputation - t distributed
sample_ABpool_t_Plummer <- function(index,model.output.list=glms, n = glms[[1]]$data$n){
  model.output <- model.output.list[[index]]
  MLE <- model.output$coeff
  var_mat <- vcov(model.output)
  return(mvtnorm::rmvt(n=1, sigma = var_mat, df = n-2)+MLE) 
}
set.seed(1)
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
rubin <- mice::pool(glms, dfcom=glms[[1]]$df.residual) # Note mice::pool(glms) gives the same result.
rubin_small <- mice::pool(glms, dfcom=glms[[1]]$df.residual)
t2 <- Sys.time()
(time_Plummer_Rubin <- difftime(t2,t1))
rubin$pooled$estimate
rubin$pooled$t
# df
rubin$pooled$df # 0.43, 0.18 (no wonder the CIs are stupidly wide!)
# fmi
rubin$pooled$fmi # 0.98, 0.99
# U_bar
rubin$pooled$ubar # 0.000090, 0.12
# b
rubin$pooled$b # 0.019, 6.3

## Bayesian multiple imputation
# Check how many iterations are needed

data_stage2 <- list("n" = nrow(HPV_cancer_data),
                    "Y" = HPV_cancer_data[,"Y"],
                    "t" = HPV_cancer_data[,"Time"],
                    "phi" = numeric(n))
data_stage2$phi <- extract_impute_out$phi[1,]
stage_2_out = rstan::sampling(
  object = stage2_model,
  data = data_stage2,
  chains = nchains,
  cores = nchains,
  iter = sampling_iterations,
  warmup = warmup_iterations,
  refresh = sampling_iterations/10, #show an update @ each %10
  seed = 1000)
stage_2_samples <- rstan::extract(stage_2_out)
# stage_2_samples_mat <- rstan::extract(stage_2_out,permuted=F)
rstan::traceplot(stage_2_out,pars=c("theta_1_std","theta_2_std"), inc_warmup=TRUE,window=c(1,500))
rstan::traceplot(stage_2_out,pars=c("theta_1","theta_2"), inc_warmup=TRUE,window=c(1,500))
# Looks like 100 iterations is enough, say 200 to be safe.
BayesMI_thetas <- matrix(nrow=m,ncol=2)
warmup_iterations <- 200
sampling_iterations <- 1 + warmup_iterations
nchains <- 1
t1 <- Sys.time()
for (i in 1:m){
  print(i)
  data_stage2$phi <- extract_impute_out$phi[i,]
  inits_list <- list(list(
    theta_1_std = rev(stage_2_samples$theta_1_std)[1],
    theta_2_std = rev(stage_2_samples$theta_2_std)[1]
  ))
  stage_2_out = rstan::sampling(
    object = stage2_model,
    data = data_stage2,
    chains = nchains,
    cores = nchains,
    iter = sampling_iterations,
    init = inits_list,
    warmup = warmup_iterations,
    refresh = 0, #show an update @ each %10
    seed = 1000)
  stage_2_samples <- rstan::extract(stage_2_out)
  BayesMI_thetas[i,1] <- stage_2_samples$theta_1[1]
  BayesMI_thetas[i,2] <- stage_2_samples$theta_2[1]
}
t2 <- Sys.time()
saveRDS(BayesMI_thetas,paste0(output_directory,"/Plummer_BayesMI_thetas.RDS"))
BayesMI_thetas <- readRDS(paste0(output_directory,"/Plummer_BayesMI_thetas.RDS"))
(time_Plummer_BayesMI <- difftime(t2,t1))

png(paste0(plot_directory,"/Plummer_QQ.png"),width = 1400,height=1000, pointsize=23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(4.1,4.1,0.1,1))
qqnorm(glm_theta[,1], main = "MLEs for theta_1 after imputation")
qqline(glm_theta[,1])
qqnorm(glm_theta[,2], main = "MLEs for theta_2 after imputation")
qqline(glm_theta[,2])
qqnorm(BayesMI_thetas[,1], main = "Posterior for theta_1 after imputation")
qqline(BayesMI_thetas[,1])
qqnorm(BayesMI_thetas[,2], main = "Posterior for theta_2 after imputation")
qqline(BayesMI_thetas[,2])
dev.off()

subsample <- ABpool_samples_norm[1:200,2]
png(paste0(plot_directory,"/Plummer_QQ_2.png"),width = 1400,height=1000, pointsize=23)
p <- ppoints(length(subsample))
q_Rubin <- rubin$pooled$estimate[2] + sqrt(rubin$pooled$t[2])*qt(p,df=(rubin$m-1)/rubin$pooled$fmi[1]^2)
qqplot(q_Rubin,subsample, xlab = "Rubin_norm theoretical quantiles", ylab = "ABpool_norm samples")
abline(0,1)
dev.off()


CIs_theta_1 <- rbind(quantile(BayesMI_thetas[,1],c(0.025,0.975)),
                     quantile(ABpool_samples_t[,1],c(0.025,0.975)),
                     quantile(ABpool_samples_norm[,1],c(0.025,0.975)),
                     rubin$pooled$estimate[1] + c(-1,1)*qt(p=0.975, df=rubin$pooled$df[1])*sqrt(rubin$pooled$t[1]),
                     rubin$pooled$estimate[1] + c(-1,1)*qt(0.975,(rubin$m-1)/rubin$pooled$fmi[1]^2)*sqrt(rubin$pooled$t[1]))
results_theta_1 <- cbind("Posterior mean"=c(mean(BayesMI_thetas[,1]),
                                            mean(ABpool_samples_t[,1]),
                                            mean(ABpool_samples_norm[,1]),
                                            rubin$pooled$estimate[1],
                                            rubin$pooled$estimate[1]),
                         "Posterior median"=c(median(BayesMI_thetas[,1]),
                                              median(ABpool_samples_t[,1]),
                                              median(ABpool_samples_norm[,1]),
                                              rubin$pooled$estimate[1],
                                              rubin$pooled$estimate[1]),
                         CIs_theta_1)
rownames(results_theta_1) <- c("Bayespool","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm")
write.csv(round(results_theta_1,2),paste0(output_directory,"/Plummer_results_theta_1.csv"),row.names=T)

CIs_theta_2 <- rbind(quantile(BayesMI_thetas[,2],c(0.025,0.975)),
                     quantile(ABpool_samples_t[,2],c(0.025,0.975)),
                     quantile(ABpool_samples_norm[,2],c(0.025,0.975)),
                     rubin$pooled$estimate[2] + c(-1,1)*qt(p=0.975, df=rubin$pooled$df[2])*sqrt(rubin$pooled$t[2]),
                     rubin$pooled$estimate[2] + c(-1,1)*qt(0.975,(rubin$m-1)/rubin$pooled$fmi[2]^2)*sqrt(rubin$pooled$t[2]))
results_theta_2 <- cbind("Posterior mean"=c(mean(BayesMI_thetas[,2]),
                                            mean(ABpool_samples_t[,2]),
                                            mean(ABpool_samples_norm[,2]),
                                            rubin$pooled$estimate[2],
                                            rubin$pooled$estimate[2]),
                         "Posterior median"=c(median(BayesMI_thetas[,2]),
                                              median(ABpool_samples_t[,2]),
                                              median(ABpool_samples_norm[,2]),
                                              rubin$pooled$estimate[2],
                                              rubin$pooled$estimate[2]),
                         CIs_theta_2)
rownames(results_theta_2) <- c("Bayespool","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm")
write.csv(round(results_theta_2,2),paste0(output_directory,"/Plummer_results_theta_2.csv"),row.names=T)

other_info <- round(rbind("Bayespool estimate" = colMeans(BayesMI_thetas),
                          "ABpool & Rubin estimate" = rubin$pooled$estimate, # Equal to the ABpool estimates
                          "dfcom" = rubin$pooled$dfcom,
                          "fmi" = rubin$pooled$fmi,
                          "df" = rubin$pooled$df,
                          "MLEs skewness" = skewness(glm_theta),
                          "MLEs kurtosis" = kurtosis(glm_theta),
                          "Posterior skewness" = skewness(BayesMI_thetas),
                          "Posterior kurtosis" = kurtosis(BayesMI_thetas)),2)
colnames(other_info) <- c("theta_1","theta_2")
write.csv(other_info,paste0(output_directory,"/Plummer_other_info.csv"),row.names=T)

time_Plummer_impute + time_Plummer_stan_compile_impute # 1 min 27s (1min 16s to compile, 11s to run)
time_Plummer_ABpool + time_Plummer_MLE # 28 seconds
time_Plummer_Rubin + time_Plummer_MLE # 55 seconds (note using m =10000 unnecessary)
time_Plummer_BayesMI + time_Plummer_stan_compile_BayesMI #  7 mins 50 s

times_Plummer <- mget(c("time_Plummer_stan_compile_impute",
                        "time_Plummer_impute",
                        "time_Plummer_MLE",
                        "time_Plummer_ABpool",
                        "time_Plummer_Rubin",
                        "time_Plummer_stan_compile_BayesMI",
                        "time_Plummer_BayesMI"))
saveRDS(times_Plummer,paste0(output_directory,"/Plummer_runtimes.RDS"))

# Plot of distribution for theta_1
png(paste0(plot_directory,"/Plummer_theta_1.png"),width = 1400,height=1000, pointsize=23)
plot(density(BayesMI_thetas[,1]), main = "Distribution of theta_1",xlab="theta_1",lwd=3,ylim=c(0,3.2))
lines(density(ABpool_samples_t[,1]), col= outcome_cols[1],lwd=2)
lines(density(ABpool_samples_norm[,1]), col= outcome_cols[1], lty=2,lwd=2)
xx <- seq(min(ABpool_samples_norm[,1])-1,max(ABpool_samples_norm[,1])+1,0.0001)
xt <- (xx-rubin$pooled$estimate[1])/sqrt(rubin$pooled$t[1])
tt <- dt(xt,df=rubin$pooled$df[1])/sqrt(rubin$pooled$t[1]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(rubin$m-1)/rubin$pooled$fmi[1]^2)/sqrt(rubin$pooled$t[1])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
legend(x=-2.4,y=2.0,
       legend = c("MI_Bayes","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm"),
       col = c("black",rep(outcome_cols,each=2)),
       lty = c(1,1,2,1,2),lwd=c(3,2,2,2,2))
dev.off()

# Similar to  above for theta_2
png(paste0(plot_directory,"/Plummer_theta_2.png"),width = 1400,height=1000, pointsize=23)
plot(density(BayesMI_thetas[,2]), main = "Distribution of theta_2", ylim = c(0,0.18),xlab="theta_2",lwd=3)
lines(density(ABpool_samples_t[,2]), col= outcome_cols[1],lwd=2)
lines(density(ABpool_samples_norm[,2]), col= outcome_cols[1], lty=2,lwd=2)
xx <- seq(min(ABpool_samples_norm[,2])-1,max(ABpool_samples_norm[,2])+1,0.0001)
xt <- (xx-rubin$pooled$estimate[2])/sqrt(rubin$pooled$t[2])
tt <- dt(xt,df=rubin$pooled$df[2])/sqrt(rubin$pooled$t[2]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(rubin$m-1)/rubin$pooled$fmi[2]^2)/sqrt(rubin$pooled$t[2])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
legend(x=22,y=0.11,
       legend = c("MI_Bayes","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm"),
       col = c("black",rep(outcome_cols,each=2)),
       lty = c(1,1,2,1,2),lwd=c(3,2,2,2,2))
dev.off()

# Cornish-Fisher approximation - normal version

# skew_kurt_CF_calc(glm_theta[,2],sapply(glms,function(x) vcov(x)[2,2]),rubin,ABpool_samples_norm,alpha=0.025,term=1,distr="norm")
MLEs <- glm_theta[,2]
Us <- sapply(glms,function(x) vcov(x)[2,2])

b_m <- var(MLEs)
u_m <- mean(Us)
nu <- rubin$pooled$dfcom[2]
df_rubin <- rubin$pooled$df[2]
total_var <- b_m + u_m
lambda <- exp(log(b_m)-log(total_var))
z_a <- qnorm(0.025)
Cov_MLEsq_U <- cov((MLEs-mean(MLEs))^2,Us)
sd_rubin <- sqrt((rubin$pooled$ubar[2]+rubin$pooled$b[2])*df_rubin/(df_rubin-2))
sd_par <- sd(ABpool_samples_norm[,2])
skew_par <- skewness(ABpool_samples_norm[,2])
kurt_par <- kurtosis(ABpool_samples_norm[,2])
CF_norm_terms <- c(0,
                   skew_par*sd_par*(z_a^2-1)/6,
                   (kurt_par*sd_par)*(z_a^3-3*z_a)/24,
                   -skew_par^2*sd_par*(2*z_a^3-5*z_a)/36)
skew_calc_terms <- c(skewness(MLEs)*lambda^(3/2),
                     3*cov(MLEs,Us)/total_var^(3/2))
kurt_calc_terms <- c(kurtosis(MLEs)*lambda^2,
                     0,
                     3*var(Us)/total_var^2,
                     6*Cov_MLEsq_U/total_var^2)
calculations_info <- c(skew_calc_terms,
                       sum(skew_calc_terms),
                       skew_par,
                       kurt_calc_terms,
                       sum(kurt_calc_terms),
                       kurt_par,
                       CF_norm_terms,
                       sum(CF_norm_terms),
                       quantile(ABpool_samples_norm[,2],c(0.025)) - (rubin$pooled$estimate[2] + c(-1)*qt(0.975,(rubin$m-1)/rubin$pooled$fmi[2]^2)*sqrt(rubin$pooled$t[2])))
names(calculations_info) <- c("skew_term_MLE","skew_term_Cov","skew_calc",
                              "skew_par",
                              "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc","kurt_par",
                              "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                              "CF_approx_norm","CI_diff_norm")
write.csv(calculations_info,paste0(output_directory,"/Plummer_calculations_norm_info_theta_2.csv"),row.names=T)
calculations_info <- read.csv(paste0(output_directory,"/Plummer_calculations_norm_info_theta_2.csv"))
calc_names <- calculations_info[,1]
calculations_info <- calculations_info[,2]
names(calculations_info) <- calc_names

