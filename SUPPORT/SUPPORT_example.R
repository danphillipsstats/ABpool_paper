##############################################################################
# SUPPORT example
##############################################################################
require(mice)
require(mvtnorm)
require(moments)
require(abpool)
require(survival)
require(Hmisc)
require(survminer)

output_directory <- paste0(getwd(),"/SUPPORT/Output")
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
get.ci.mice.pred.i <- function(i, pooled.out){
  pooled.out$estimate[i] + c(-1,1)*qt(p=0.975, df=pooled.out$df[i])*sqrt(pooled.out$t[i])
}
get.ci.mice <- function(pooled.out){
  out <- vapply(1:nrow(pooled.out),get.ci.mice.pred.i,numeric(2),pooled.out=pooled.out)
  colnames(out) <- pooled.out[,1]
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

# Function for excess kurtosis
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
# Cornish-Fisher approximation to the differences between ABpool and Rubin's
# lower alpha quantiles. Assumes the posterior means agree, since the ABpool sample mean
# can be replaced by the mean of the estimates, matching Rubin's mean and with a lower variance
# The function assumes the ABpool samples and rubin_pooled were generated using the same
# complete-data degrees of freedom dfcom. rubin_pooled <- mice::pool(...)$pooled
skew_kurt_CF_calc <- function(estimates,variances,rubin_pooled,ABpool_samples,alpha=0.025,term=NULL){
  if (!is.numeric(estimates)){stop("Estimates must be a numeric vector or matrix, with each row corresponding to each imputation")}
  if (!is.null(dim(variances))){stop("Variances must be a numeric vector or list of matrices, with entry in the vector/list corresponding to the variance(-covariance matrix) for each imputation")}
  if (is.numeric(variances)){
    # Single parameter
    if(length(variances)!=length(estimates)){stop("Estimates and variances must be vectors of the same length, with each entry corresponding to one imputation. Or a matrix and a list of matrices respectively.")}
    if (is.null(term)) term <- 1
  } else if (is.list(variances)){
    # Multi parameter
    # Checks
    if (length(dim(variances[[1]]))!=2 || dim(variances[[1]])[1]!=dim(variances[[1]])[2]){stop("Where variances is supplied as a list, each entry must be a matrix giving the variance-covariance matrix for that imputation.")}
    if (length(dim(estimates))!=2){stop("Where variances is supplied as a list, estimates must be a matrix")}
    if (dim(estimates)[1]!=length(variances)){stop("The number of rows in estimates must match the length of the variances list")}
    if (dim(estimates)[2]!=dim(variances[[1]])[1]){stop("The number of columns in estimates must match the dimension of each entry of variances")}
    if (is.null(term)){stop("term must be supplied where multiple parameters are present")}
    # Select the correct term
    estimates <- estimates[,term]
    variances <- vapply(variances,function(x) x[term,term], numeric(1))
    ABpool_samples <- ABpool_samples[,term]
  }
  
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
  
  rownames(rubin_pooled) <- rubin_pooled$term
  b_m <- var(estimates)
  u_m <- mean(variances)
  nu <- rubin_pooled[term,"dfcom"]
  df_rubin <- rubin_pooled[term,"df"]
  if (nu <= 4 || df_rubin <= 4) {
    stop("Degrees of freedom must exceed 4 for kurtosis calculations.")
  }
  if (!is.infinite(nu)) nu_nu_2 <- nu/(nu-2) else nu_nu_2 <- 1
  total_var <- b_m + nu_nu_2*u_m
  lambda_nu <- b_m/total_var
  z_a <- qnorm(alpha)
  Cov_MLEsq_U <- cov((estimates-mean(estimates))^2,variances)
  sd_rubin <- sqrt(rubin_pooled[term,"t"]*df_rubin/(df_rubin-2))
  sd_par <- sd(ABpool_samples)
  skew_par <- skewness(ABpool_samples)
  kurt_par <- kurtosis(ABpool_samples)
  CF_terms <- c(z_a*(sd_par-sd_rubin),
                skew_par*sd_par*(z_a^2-1)/6,
                (kurt_par*sd_par - sd_rubin*6/(df_rubin-4))*(z_a^3-3*z_a)/24,
                -skew_par^2*sd_par*(2*z_a^3-5*z_a)/36)
  skew_calc_terms <- c(skewness(estimates)*lambda_nu^(3/2),
                       3*nu_nu_2*cov(estimates,variances)/total_var^(3/2))
  if (!is.infinite(nu)) kurt_calc_2 <- 6*(1-lambda_nu)^2/(nu-4) else kurt_calc_2 <- 0
  if (!is.infinite(nu)) kurt_nu_frac <- nu^2/((nu-2)*(nu-4)) else kurt_nu_frac <- 1
  kurt_calc_terms <- c(kurtosis(estimates)*lambda_nu^2,
                       kurt_calc_2,
                       3*kurt_nu_frac*var(variances)/total_var^2,
                       6*nu_nu_2*Cov_MLEsq_U/total_var^2)
  calculations_info <- c(skew_calc_terms,
                         sum(skew_calc_terms),
                         skew_par,
                         kurt_calc_terms,
                         sum(kurt_calc_terms),
                         kurt_par,
                         CF_terms,
                         sum(CF_terms),
                         quantile(ABpool_samples,c(alpha)) - (rubin_pooled[term,"estimate"] + qt(alpha,df_rubin)*sqrt(rubin_pooled[term,"t"])))
  names(calculations_info) <- c("skew_term_MLE","skew_term_Cov","skew_calc",
                                "skew_par",
                                "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc","kurt_par",
                                "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                                "CF_approx","CI_diff")
  return(calculations_info)
}
###################################################
# SUPPORT

getHdata(support) # getHdata(support, "all") to view data dictionary
summary(support)
apply(support,2,function(vec) sum(is.na(vec)))/nrow(support)*100
table(support$death,support$dzgroup,is.na(support$bili))
support_filled_in <- support
support_filled_in$alb[which(is.na(support_filled_in$alb))] <- 3.5
support_filled_in$pafi[which(is.na(support_filled_in$pafi))] <- 333.3
support_filled_in$bili[which(is.na(support_filled_in$bili))] <- 1.01
support_filled_in$crea[which(is.na(support_filled_in$crea))] <- 1.01
support_filled_in$bun[which(is.na(support_filled_in$bun))] <- 6.51
support_filled_in$wblc[which(is.na(support_filled_in$wblc))] <- 9
support_filled_in$urine[which(is.na(support_filled_in$urine))] <- 2502
support.overall.formula <- Surv(d.time, death) ~ age + num.co + scoma + alb + bili + pafi + crea + meanbp + hrt + resp + temp + adlsc + dzgroup
support.overall.cox <- coxph(support.overall.formula,data=support_filled_in)
summary(support.overall.cox)
confint(support.overall.cox)
# Chosen by backward deletion from the full model above, in a complete-case analysis (any individuals with a missing value for any variable included in the current model is excluded.)
support.formula <- Surv(d.time, death) ~ age + num.co + scoma + alb + bili + pafi + crea + meanbp + hrt + resp + temp + adlsc + dzgroup
support.vs.cox <- coxph(support.formula,data=support_filled_in)
summary(support.vs.cox)
#####
skewness(support$bili,na.rm=T); skewness(support$crea,na.rm=T); kurtosis(support$bili,na.rm=T); kurtosis(support$crea,na.rm=T)
cox.support.form.char <- "Surv(d.time, death) ~ age + num.co + scoma + alb + bili + pafi + crea + meanbp + hrt + resp + temp + adlsc + dzgroup"
(support.cox <- coxph(as.formula(cox.support.form.char),data=support))
summary(support.cox)
key_parameters <- c("alb", "bili", "pafi", "crea")

npred.support <- length(coef(support.cox))

support.data <- support[,all.vars(as.formula(cox.support.form.char))]
apply(support.data,2,function(vec) sum(is.na(vec)))/nrow(support.data)*100 # 37.8% for alb, 29.7% for bili, 25.3% for pafi, 0.3% for crea, all others fully observed
support.data.obs <- support.data[complete.cases(support.data),]

#########################################################################
m <- 100000
# The following took hours to run on my laptop. The output can be loaded below to avoid the need to rerun.
{
  set.seed(1)
  t1 <- Sys.time()
  support.data$H0_NA <- mice::nelsonaalen(support.data, timevar = "d.time", statusvar = "death")
  pred <- mice::make.predictorMatrix(support.data)
  # choose methods: do not impute death or H_NA; impute others with pmm
  method <- make.method(support.data)
  for (v in names(method)) {
    method[v] <- "pmm"
  }
  pred[, "d.time"] <- 0 # do not predict using d.time
  # use H0_NA and death to predict other vars
  pred[, "H0_NA"] <- 1
  pred[, "death"] <- 1
  pred["H0_NA", ] <- 0
  pred["death", ] <- 0
  # Note since we're using predictive mean matching, there is only 80 odd values they might realistically choose.
  # May be a computationally better alternative to impute say m = 1000 and then use J = 100 Cox models for each one. But we don't do this.
  # The imputation model takes hours to run
  imp_data_support <- mice(support.data, m = m, method = method, predictorMatrix = pred)
  t2 <- Sys.time()
  (pmm100000 <- difftime(t2,t1))
  saveRDS(imp_data_support,paste0(output_directory,"/imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
  imp_data_support <- readRDS(paste0(output_directory,"/imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
  # The Cox model took about 30 mins to run on my laptop.
  # I haven't uploaded the output onto GitHub as it is too large (about 3GB).
  set.seed(1)
  t1 <- Sys.time()
  support.cox.imp <- with(imp_data_support,coxph(as.formula(cox.support.form.char))) # Took about 30 mins to run
  t2 <- Sys.time()
  (cox100000 <- difftime(t2,t1))
  saveRDS(support.cox.imp,paste0(output_directory,"/cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
  support.cox.imp <- readRDS(paste0(output_directory,"/cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
  
  set.seed(1)
  support.cox.abpool.norm100000 <- abpool::abpool(support.cox.imp, parameters = key_parameters, dfcom = Inf)
  # QQ diagnostic plot comparing ABpool and Rubin's rules from all 100,000 imputations
  png(paste0(plot_directory,"/SUPPORT_QQplot_large.png"),width = 1400,height=1000, pointsize=23)
  par(mfrow = c(2,2))
  abpool::qqplot_rubin_abpool(support.cox.abpool.norm100000)
  dev.off()
  support.cox.abpool.norm <- support.cox.abpool.norm100000$samples
  rm(support.cox.abpool.norm100000)
  
  # Save 200 imputations for later QQ diagnostic
  support.cox.imp200 <- support.cox.imp
  support.cox.imp200$analyses <- support.cox.imp$analyses[1:200]
  
  support.estimates <- t(vapply(1:m, function(i) coef(support.cox.imp$analyses[[i]])[key_parameters],numeric(length(key_parameters))))
  support.variances <- lapply(1:m, function(i) vcov(support.cox.imp$analyses[[i]])[key_parameters,key_parameters])
  
  set.seed(1)
  t1 <- Sys.time()
  support.cox.pool.t <- mice::pool(support.cox.imp)
  support.cox.pooled.t <- support.cox.pool.t$pooled # Since this is too large and we don't use it
  rownames(support.cox.pooled.t) <- support.cox.pooled.t$term
  t2 <- Sys.time()
  (rubint100000 <- difftime(t2,t1))
  rm(support.cox.pool.t)
  
  t1 <- Sys.time()
  support.cox.pool.norm <- mice::pool(support.cox.imp, dfcom = Inf)
  support.cox.pooled.norm <- support.cox.pool.norm$pooled # Since this is too large and we don't use it
  rownames(support.cox.pooled.norm) <- support.cox.pooled.norm$term
  t2 <- Sys.time()
  (rubin100000 <- difftime(t2,t1))
  rm(support.cox.pool.norm)
  
  set.seed(1000)
  dfcom_finite <- support.cox.imp$analyses[[1]]$nevent - length(support.cox.imp$analyses[[1]]$coefficients)
  t1 <- Sys.time()
  support.cox.abpool.t <- abpool::abpool(support.cox.imp, parameters = key_parameters, dfcom = dfcom_finite)$samples
  t2 <- Sys.time()
  (abpoolt100000 <- difftime(t2,t1))
  
  saveRDS(mget(c("support.cox.pooled.norm","support.cox.pooled.t","support.cox.abpool.norm","support.cox.abpool.t",
                 "support.cox.imp200", "support.estimates", "support.variances")),
          paste0(output_directory,"/SUPPORT_output_objects.RDS"))
}
########################################################################
# Load output
output_objects <- readRDS(paste0(output_directory,"/SUPPORT_output_objects.RDS"))
mgcv::ls.size(output_objects)
list2env(output_objects, envir = .GlobalEnv)
rm(output_objects)

set.seed(1)
for (par_name in key_parameters){
  qqnorm(support.estimates[,par_name], main = par_name)
  qqline(support.estimates[,par_name])
}

round(get.ci.mice(support.cox.pooled.norm),4)
round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4)
round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4)

# Get skewness and fmi
fmis <- support.cox.pooled.norm$fmi
names(fmis) <- names(coef(support.cox))
fmis[key_parameters] # alb and bili have non-negligible fmi
support.Us <- t(vapply(1:m, function(i) diag(support.variances[[i]]),numeric(ncol(support.estimates))))
(skew_MLE <- apply(support.estimates,2,skewness)) # crea and bili have considerable skewness (-1.1,-0.5)
(skew_par <- apply(support.cox.abpool.norm,2,skewness)) # but not for the posterior(-0.1,-0.1)
(kurt_MLE <- apply(support.estimates,2,kurtosis))
(kurt_par <- apply(support.cox.abpool.norm,2,kurtosis)) # All parameters have small excess kurtosis
(skew_U <- apply(support.Us,2,skewness)) # crea and bili have considerable skewness (-1.1,-0.5)
(kurt_U <- apply(support.Us,2,kurtosis))
rm(support.Us)

set.seed(1)
support.cox.abpool.norm200 <- abpool::abpool(support.cox.imp200, parameters = key_parameters, dfcom = Inf)
# QQ diagnostic plot comparing ABpool and Rubin's rules from 200 imputations
png(paste0(plot_directory,"/SUPPORT_QQplot_200.png"),width = 1400,height=1000, pointsize=23)
par(mfrow = c(2,2))
abpool::qqplot_rubin_abpool(support.cox.abpool.norm200)
dev.off()

# Some MLEs noticeably non-normal
png(paste0(plot_directory,"/SUPPORT_QQplot_MLEs.png"),width = 1400,height=800, pointsize=23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(0.1,0.1,0.0,1))
for (i in 1:length(key_parameters)){
  vars <- key_parameters[i]
  qqnorm(support.estimates[,vars], main = key_parameters[i],cex=1)
  qqline(support.estimates[,vars])
}
dev.off()

support.cox.pooled.norm[key_parameters,"fmi"]
full_names <- c("albumin","bilirubin","PaO2/FiO2","creatinine")

png(paste0(plot_directory,"/SUPPORT_compare_post_missing.png"),width = 1400,height=1000, pointsize=23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(4.1,4.1,0.1,1))
for (i in 1:4){
  plot(density(support.cox.abpool.t[,key_parameters[i]]),
       ylim = c(0,max(density(support.cox.abpool.t[,key_parameters[i]])$y)*1.2),
       main = paste0("Distribution of ",full_names[i]," effect"),
       xlab= paste0("Effect due to ",full_names[i]), col=NA)
  xx <- seq(min(support.cox.abpool.norm[,key_parameters[i]])-1,max(support.cox.abpool.norm[,key_parameters[i]])+1,0.0001)
  ss <- dnorm(xx,coef(support.overall.cox)[key_parameters[i]],sqrt(vcov(support.overall.cox)[key_parameters[i],key_parameters[i]])) # Fill-in single imputation
  lines(ss~xx, col = "grey",lty=2,lwd=2)
  lines(density(support.cox.abpool.t[,key_parameters[i]]), col= outcome_cols[1],lwd=2)
  lines(density(support.cox.abpool.norm[,key_parameters[i]]), col= outcome_cols[1], lty=2,lwd=2)
  xt <- (xx-support.cox.pooled.t[key_parameters[i],"estimate"])/sqrt(support.cox.pooled.t[key_parameters[i],"t"])
  tt <- dt(xt,df=support.cox.pooled.t[key_parameters[i],"df"])/sqrt(support.cox.pooled.t[key_parameters[i],"t"]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
  lines(tt~xx, col = outcome_cols[2],lwd=2)
  xn <- (xx-support.cox.pooled.norm[key_parameters[i],"estimate"])/sqrt(support.cox.pooled.norm[key_parameters[i],"t"])
  nn <- dt(xn,df=(support.cox.pooled.norm[key_parameters[i],"m"]-1)/support.cox.pooled.norm[key_parameters[i],"fmi"]^2)/sqrt(support.cox.pooled.norm[key_parameters[i],"t"])
  lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
}
dev.off()

i <- 2
png(paste0(plot_directory,"/SUPPORT_compare_post_bili.png"),width = 1400,height=1000, pointsize=23)
plot(density(support.cox.abpool.t[,key_parameters[i]]), 
     main = paste0("Distribution of ",full_names[i]," effect"),
     xlab= paste0("Effect due to ",full_names[i]), col=NA,ylim=c(0,46))
xx <- seq(min(support.cox.abpool.norm[,key_parameters[i]])-1,max(support.cox.abpool.norm[,key_parameters[i]])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[key_parameters[i]],sqrt(vcov(support.overall.cox)[key_parameters[i],key_parameters[i]])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abpool.t[,key_parameters[i]]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abpool.norm[,key_parameters[i]]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pooled.t[key_parameters[i],"estimate"])/sqrt(support.cox.pooled.t[key_parameters[i],"t"])
tt <- dt(xt,df=support.cox.pooled.t[key_parameters[i],"df"])/sqrt(support.cox.pooled.t[key_parameters[i],"t"]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
xn <- (xx-support.cox.pooled.norm[key_parameters[i],"estimate"])/sqrt(support.cox.pooled.norm[key_parameters[i],"t"])
nn <- dt(xn,df=(support.cox.pooled.norm[key_parameters[i],"m"]-1)/support.cox.pooled.norm[key_parameters[i],"fmi"]^2)/sqrt(support.cox.pooled.norm[key_parameters[i],"t"])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
dev.off()

posterior_summaries <- cbind("term" = rep(c("alb","bili","pafi","crea"),5),
                             "Method" = rep(c("Single_imp","Rubin_t","Rubin_norm","ABpool_t","ABpool_norm"),each=4),
                             "Mean" = c(coef(support.overall.cox)[key_parameters],
                                        support.cox.pooled.t[key_parameters,"estimate"],
                                        support.cox.pooled.norm[key_parameters,"estimate"],
                                        apply(support.cox.abpool.t,2,mean)[key_parameters],
                                        apply(support.cox.abpool.norm,2,mean)[key_parameters]),
                             "Median" = c(coef(support.overall.cox)[key_parameters],
                                          support.cox.pooled.t[key_parameters,"estimate"],
                                          support.cox.pooled.norm[key_parameters,"estimate"],
                                          apply(support.cox.abpool.t,2,median)[key_parameters],
                                          apply(support.cox.abpool.norm,2,median)[key_parameters]),
                             "2.5%" = c(t(confint(support.overall.cox))[1,key_parameters],
                                        round(get.ci.mice(support.cox.pooled.t),4)[1,key_parameters],
                                        round(get.ci.mice(support.cox.pooled.norm),4)[1,key_parameters],
                                        round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4)[1,key_parameters],
                                        round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4)[1,key_parameters]),
                             "97.5%" = c(t(confint(support.overall.cox))[2,key_parameters],
                                         round(get.ci.mice(support.cox.pooled.t),4)[2,key_parameters],
                                         round(get.ci.mice(support.cox.pooled.norm),4)[2,key_parameters],
                                         round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4)[2,key_parameters],
                                         round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4)[2,key_parameters])
)
posterior_summaries[,c("Method","term")] <- as.character(posterior_summaries[,c("Method","term")])
posterior_summaries[,c("Mean","Median","2.5%","97.5%")] <- format(round(as.numeric(as.matrix(posterior_summaries[,c("Mean","Median","2.5%","97.5%")])),4),4)
posterior_summaries <- posterior_summaries[c(which(posterior_summaries[,"term"]=="alb"),
                                             which(posterior_summaries[,"term"]=="bili"),
                                             which(posterior_summaries[,"term"]=="pafi"),
                                             which(posterior_summaries[,"term"]=="crea")),]

missing_var_outputs <- data.frame("dfcom" = support.cox.pooled.norm[key_parameters,"dfcom"],
                                  "df_Rubin_norm" = support.cox.pooled.norm[key_parameters,"df"],
                                  "df_Rubin_t" = support.cox.pooled.t[key_parameters,"df"],
                                  "fmi" = fmis[key_parameters],
                                  "skew_MLE" = skew_MLE[key_parameters],
                                  "kurt_MLE" = kurt_MLE[key_parameters],
                                  "skew_post" = skew_par[key_parameters],
                                  "kurt_post" = kurt_par[key_parameters],
                                  "skew_U"   = skew_U[key_parameters],
                                  "kurt_U"   = kurt_U[key_parameters])

outputs <- list("Posterior_summaries" = posterior_summaries,
                "missing_var_outputs" = missing_var_outputs,
                "CI_Rubin_t" = round(get.ci.mice(support.cox.pooled.t),4),
                "CI_Rubin_norm" = round(get.ci.mice(support.cox.pooled.norm),4),
                "CI_ABpool_t" = round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4),
                "CI_ABpool_norm" = round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4),
                "estimates_Rubin" = support.cox.pooled.norm$estimate,
                "estimates_ABpool_t" = apply(support.cox.abpool.t,2,mean),
                "estimates_ABpool_norm" = apply(support.cox.abpool.norm,2,mean))

saveRDS(outputs,paste0(output_directory,"/SUPPORT_outputs_results.RDS"))
outputs <- readRDS(paste0(output_directory,"/SUPPORT_outputs_results.RDS"))

write.csv(outputs$Posterior_summaries,
          paste0(output_directory,"/SUPPORT_results.csv"),
          row.names = FALSE)

# Create LaTeX table
df <- read.csv(paste0(output_directory,"/SUPPORT_results.csv"),
               check.names = FALSE,
               colClasses = "character")

# Escape underscores for LaTeX
df$Method <- gsub("_", "\\_", df$Method, fixed = TRUE)

# Only display each parameter name on its first row
df$term[duplicated(df$term)] <- ""

# Build LaTeX table rows
rows <- apply(df, 1, function(x) {
  paste0("\t\t", paste(x, collapse = " & "), " \\\\")
})

# Add horizontal lines after each parameter group
group_ends <- which(c(df$term[-1] != "", TRUE))

rows[group_ends] <- paste0(rows[group_ends], "\n\t\t\\hline")

# Build complete LaTeX table
latex <- paste0(
  "\\begin{supptable}\n",
  "\t\\captionof{supptable}{\\textbf{Posterior summaries for SUPPORT}}\n",
  "\t\\centering\n",
  "\t\\begin{tabular}{|cc|cccc|}\n",
  "\t\t\\hline\n",
  "\t\t& Method & Mean & Median & 2.5\\% & 97.5\\% \\\\\n",
  "\t\t\\hline\n",
  paste(rows, collapse = "\n"),
  "\n\t\\end{tabular}\\flushleft\n",
  "\tPosterior mean, median and 95\\% uncertainty interval for four of the parameters in the Cox model on the SUPPORT dataset. We compare the single imputation model recommended by the authors (Single\\_imp) with multiple imputation pooled by Rubin's rules and Approximate Bayesian pooling (ABpool), each with a Gaussian (norm) and $t$ assumption.\n",
  "\t\\label{tab:SUPPORT}\n",
  "\\end{supptable}"
)

# Print LaTeX code to console
cat(latex)

write.csv(format(round(outputs$missing_var_outputs,2),2),
          paste0(output_directory,"/SUPPORT_supplementary_info.csv"),
          row.names = TRUE)

# times100000 <- list(pmm100000,cox100000,rubin100000,abpoolt100000,abpoolnorm100000)
# saveRDS(times100000,paste0(output_directory,"/SUPPORT_times100000.RDS"))
# times100000 <- readRDS(paste0(output_directory,"/SUPPORT_times10000.RDS"))

# Cornish-Fisher
# Test Cornish-Fisher where average U is used to sample ABpool
CF_calcs <- rbind(skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[1],alpha=0.025),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[2],alpha=0.025),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[3],alpha=0.025),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[4],alpha=0.025),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[1],alpha=0.975),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[2],alpha=0.975),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[3],alpha=0.975),
                  skew_kurt_CF_calc(support.estimates,support.variances,support.cox.pooled.norm,support.cox.abpool.norm,term=key_parameters[4],alpha=0.975))
rownames(CF_calcs) <- paste0(rep(key_parameters,2),rep(c("_2.5%","_97.5%"),each=4))

write.csv(round(CF_calcs,6),
          paste0(output_directory,"/SUPPORT_skew_kurt_CornishFisher.csv"),
          row.names = TRUE)
CF_calcs <- read.csv(paste0(output_directory,"/SUPPORT_skew_kurt_CornishFisher.csv"),row.names = 1)
