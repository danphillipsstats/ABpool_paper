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
# Cornish-Fisher approximation to the differences between ABpool and Rubin's
# lower alpha quantiles. Assumes the posterior means agree, since the ABpool sample mean
# can be replaced by the mean of the estimates, matching Rubin's mean and with a lower variance
# The function assumes the ABpool samples and rubin_pool were generated using the same
# complete-data degrees of freedom dfcom.
skew_kurt_CF_calc <- function(estimates,variances,rubin_pool,ABpool_samples,alpha=0.025,term=NULL){
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
  
  b_m <- var(estimates)
  u_m <- mean(variances)
  nu <- rubin_pool$pooled$dfcom[term]
  df_rubin <- rubin_pool$pooled$df[term]
  if (nu <= 4 || df_rubin <= 4) {
    stop("Degrees of freedom must exceed 4 for kurtosis calculations.")
  }
  if (!is.infinite(nu)) nu_nu_2 <- nu/(nu-2) else nu_nu_2 <- 1
  total_var <- b_m + nu_nu_2*u_m
  lambda_nu <- b_m/total_var
  z_a <- qnorm(alpha)
  Cov_MLEsq_U <- cov((estimates-mean(estimates))^2,variances)
  sd_rubin <- sqrt(rubin_pool$pooled$t[term]*df_rubin/(df_rubin-2))
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
                         quantile(ABpool_samples,c(alpha)) - (rubin_pool$pooled$estimate[term] + qt(alpha,df_rubin)*sqrt(rubin_pool$pooled$t[term])))
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

set.seed(1)
# Note since we're using predictive mean matching, there is only 80 odd values they might realistically choose.
# May be computationally better to impute say m = 1000 and then use J = 100 Cox models for each one
m <- 100000
# The following takes a long time - the output can instead be loaded below.
{
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
  imp_data_support <- mice(support.data, m = m, method = method, predictorMatrix = pred)
  t2 <- Sys.time()
  (pmm100000 <- difftime(t2,t1))
  saveRDS(imp_data_support,paste0(output_directory,"/imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
}
imp_data_support <- readRDS(paste0(output_directory,"/imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
t1 <- Sys.time()
support.cox.imp <- with(imp_data_support,coxph(as.formula(cox.support.form.char)))
t2 <- Sys.time()
(cox100000 <- difftime(t2,t1))
saveRDS(support.cox.imp,paste0(output_directory,"/cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
support.cox.imp <- readRDS(paste0(output_directory,"/cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))

coef_list <- list()
for (j in 1:npred.support){
  coef_list[[names(coef(support.cox))[j]]] <- sapply(1:m, function(i) coef(support.cox.imp$analyses[[i]])[j])
  qqnorm(coef_list[[j]], main = names(coef(support.cox))[j])
  qqline(coef_list[[j]])
}
t1 <- Sys.time()
support.cox.pool.t <- mice::pool(support.cox.imp)
t2 <- Sys.time()
(rubint100000 <- difftime(t2,t1))

t1 <- Sys.time()
support.cox.pool.norm <- mice::pool(support.cox.imp, dfcom = Inf)
t2 <- Sys.time()
(rubin100000 <- difftime(t2,t1))

set.seed(1)
t1 <- Sys.time()
support.cox.abpool.norm <- t(vapply(1:m,FUN=sample_ABpool_norm,FUN.VALUE=numeric(length(coef(support.cox))), "mice.with.list"=support.cox.imp))
t2 <- Sys.time()
(abminorm100000 <- difftime(t2,t1))
t1 <- Sys.time()
support.cox.abpool.t <- t(vapply(1:m,FUN=sample_ABpool_t_cox,FUN.VALUE=numeric(length(coef(support.cox))), "mice.with.list"=support.cox.imp))
t2 <- Sys.time()
(abpoolt100000 <- difftime(t2,t1))

saveRDS(mget(c("support.cox.pool","support.cox.abpool.norm","support.cox.abpool.t")),
        paste0(output_directory,"/SUPPORT_output_objects.RDS"))
output_objects <- readRDS(paste0(output_directory,"/SUPPORT_output_objects.RDS"))
support.cox.pool <- output_objects$support.cox.pool
support.cox.abpool.norm <- output_objects$support.cox.abpool.norm
support.cox.abpool.t <- output_objects$support.cox.abpool.t
rm(output_objects)

round(get.ci.mice(support.cox.pool),4)
round(apply(support.cox.abmi.norm,2,quantile,c(0.025,0.975)),4)
round(apply(support.cox.abmi.t,2,quantile,c(0.025,0.975)),4)

# Get skewness and fmi
fmis <- support.cox.pool.norm$pooled$fmi
names(fmis) <- names(coef(support.cox))
fmis # alb and bili have non-negligible fmi
support.MLEs <- t(vapply(1:m, function(i) coef(support.cox.imp$analyses[[i]]),numeric(length(coef(support.cox.imp$analyses[[1]])))))
support.Us <- t(vapply(1:m, function(i) diag(vcov(support.cox.imp$analyses[[i]])),numeric(length(coef(support.cox.imp$analyses[[1]])))))
(skew_MLE <- apply(support.MLEs,2,skewness)) # crea and bili have considerable skewness (-1.1,-0.5)
(skew_par <- apply(support.cox.abpool.norm,2,skewness)) # but not for the posterior(-0.1,-0.1)
(kurt_MLE <- apply(support.MLEs,2,kurtosis))
(kurt_par <- apply(support.cox.abpool.norm,2,kurtosis)) # All parameters have small excess kurtosis
(skew_U <- apply(support.Us,2,skewness)) # crea and bili have considerable skewness (-1.1,-0.5)
(kurt_U <- apply(support.Us,2,kurtosis))

support.cox.imp200 <- support.cox.imp
support.cox.imp200$analyses <- support.cox.imp$analyses[1:200]
support.cox.abpool.norm200 <- abpool::abpool(support.cox.imp200, parameters = key_parameters, dfcom = Inf)
# QQ diagnostic plot comparing ABpool and Rubin's rules from 200 imputations
png(paste0(plot_directory,"/SUPPORT_QQplot_200.png"),width = 1400,height=1000, pointsize=23)
par(mfrow = c(2,2))
abpool::qqplot_rubin_abpool(support.cox.abpool.norm200)
dev.off()

support.cox.abpool.norm100000 <- abpool::abpool(support.cox.imp, parameters = key_parameters, dfcom = Inf)
# QQ diagnostic plot comparing ABpool and Rubin's rules from all 100,000 imputations
png(paste0(plot_directory,"/SUPPORT_QQplot_large.png"),width = 1400,height=1000, pointsize=23)
par(mfrow = c(2,2))
abpool::qqplot_rubin_abpool(support.cox.abpool.norm100000)
dev.off()

# Some MLEs noticeably non-normal
png(paste0(plot_directory,"/SUPPORT_QQplot_MLEs.png"),width = 1400,height=800, pointsize=23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(0.1,0.1,0.0,1))
for (i in 1:length(key_parameters)){
  vars <- key_parameters[i]
  qqnorm(support.MLEs[,vars], main = key_parameters[i],cex=1)
  qqline(support.MLEs[,vars])
}
dev.off()

alb_term <- which(support.cox.pool.norm$pooled$term=="alb")
bili_term <- which(support.cox.pool.norm$pooled$term=="bili")
pafi_term <- which(support.cox.pool.norm$pooled$term=="pafi")
crea_term <- which(support.cox.pool.norm$pooled$term=="crea")

support.cox.pool.norm$pooled$fmi[c(alb_term,bili_term,pafi_term,crea_term)]

png(paste0(plot_directory,"/SUPPORT_compare_post_missing.png"),width = 1400,height=1000, pointsize=23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(4.1,4.1,0.1,1))
plot(density(support.cox.abpool.t[,"alb"]), main = "Distribution of albumin effect",xlab="Effect due to albumin", col=NA,ylim=c(0,6))
xx <- seq(min(support.cox.abpool.norm[,"alb"])-1,max(support.cox.abpool.norm[,"alb"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[alb_term],sqrt(vcov(support.overall.cox)[alb_term,alb_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abpool.t[,"alb"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abpool.norm[,"alb"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool.t$pooled$estimate[alb_term])/sqrt(support.cox.pool.t$pooled$t[alb_term])
tt <- dt(xt,df=support.cox.pool.t$pooled$df[alb_term])/sqrt(support.cox.pool.t$pooled$t[alb_term]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
xn <- (xx-support.cox.pool.norm$pooled$estimate[alb_term])/sqrt(support.cox.pool.norm$pooled$t[alb_term])
nn <- dt(xn,df=(support.cox.pool.norm$m-1)/support.cox.pool.norm$pooled$fmi[alb_term]^2)/sqrt(support.cox.pool.norm$pooled$t[alb_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)

plot(density(support.cox.abpool.t[,"bili"]), main = "Distribution of bilirubin effect",xlab="Effect due to bilirubin", col=NA,ylim=c(0,46))
xx <- seq(min(support.cox.abpool.norm[,"bili"])-1,max(support.cox.abpool.norm[,"bili"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[bili_term],sqrt(vcov(support.overall.cox)[bili_term,bili_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abpool.t[,"bili"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abpool.norm[,"bili"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool.t$pooled$estimate[bili_term])/sqrt(support.cox.pool.t$pooled$t[bili_term])
tt <- dt(xt,df=support.cox.pool.t$pooled$df[bili_term])/sqrt(support.cox.pool.t$pooled$t[bili_term]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
xn <- (xx-support.cox.pool.norm$pooled$estimate[bili_term])/sqrt(support.cox.pool.norm$pooled$t[bili_term])
nn <- dt(xn,df=(support.cox.pool.norm$m-1)/support.cox.pool.norm$pooled$fmi[bili_term]^2)/sqrt(support.cox.pool.norm$pooled$t[bili_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)

plot(density(support.cox.abpool.t[,"pafi"]), main = "Distribution of PaO2/FiO2 ratio effect",xlab="Effect due to PaO2/FiO2 ratio", col=NA)
xx <- seq(min(support.cox.abpool.norm[,"pafi"])-1,max(support.cox.abpool.norm[,"pafi"])+1,0.000001)
ss <- dnorm(xx,coef(support.overall.cox)[pafi_term],sqrt(vcov(support.overall.cox)[pafi_term,pafi_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abpool.t[,"pafi"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abpool.norm[,"pafi"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool.t$pooled$estimate[pafi_term])/sqrt(support.cox.pool.t$pooled$t[pafi_term])
tt <- dt(xt,df=support.cox.pool.t$pooled$df[pafi_term])/sqrt(support.cox.pool.t$pooled$t[pafi_term]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
xn <- (xx-support.cox.pool.norm$pooled$estimate[pafi_term])/sqrt(support.cox.pool.norm$pooled$t[pafi_term])
nn <- dt(xn,df=(support.cox.pool.norm$m-1)/support.cox.pool.norm$pooled$fmi[pafi_term]^2)/sqrt(support.cox.pool.norm$pooled$t[pafi_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)

plot(density(support.cox.abpool.norm[,"crea"]), main = "Distribution of creatinine effect",xlab="Effect due to creatinine", col=NA,ylim=c(0,19))
xx <- seq(min(support.cox.abpool.norm[,"crea"])-1,max(support.cox.abpool.norm[,"crea"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[crea_term],sqrt(vcov(support.overall.cox)[crea_term,crea_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abpool.t[,"crea"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abpool.norm[,"crea"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool.t$pooled$estimate[crea_term])/sqrt(support.cox.pool.t$pooled$t[crea_term])
tt <- dt(xt,df=support.cox.pool.t$pooled$df[crea_term])/sqrt(support.cox.pool.t$pooled$t[crea_term]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
xn <- (xx-support.cox.pool.norm$pooled$estimate[crea_term])/sqrt(support.cox.pool.norm$pooled$t[crea_term])
nn <- dt(xn,df=(support.cox.pool.norm$m-1)/support.cox.pool.norm$pooled$fmi[crea_term]^2)/sqrt(support.cox.pool.norm$pooled$t[crea_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
dev.off()

png(paste0(plot_directory,"/SUPPORT_compare_post_bili.png"),width = 1400,height=1000, pointsize=23)
plot(density(support.cox.abpool.t[,"bili"]), main = "Distribution of bilirubin effect",xlab="Effect due to bilirubin", col=NA,ylim=c(0,46))
xx <- seq(min(support.cox.abpool.norm[,"bili"])-1,max(support.cox.abpool.norm[,"bili"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[bili_term],sqrt(vcov(support.overall.cox)[bili_term,bili_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abpool.t[,"bili"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abpool.norm[,"bili"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool.t$pooled$estimate[bili_term])/sqrt(support.cox.pool.t$pooled$t[bili_term])
tt <- dt(xt,df=support.cox.pool.t$pooled$df[bili_term])/sqrt(support.cox.pool.t$pooled$t[bili_term]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
xn <- (xx-support.cox.pool.norm$pooled$estimate[bili_term])/sqrt(support.cox.pool.norm$pooled$t[bili_term])
nn <- dt(xn,df=(support.cox.pool.norm$m-1)/support.cox.pool.norm$pooled$fmi[bili_term]^2)/sqrt(support.cox.pool.norm$pooled$t[bili_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
legend(x=0,y=20,
       legend = c("Single_imp","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm"),
       col = c("grey",rep(outcome_cols,each=2)),
       lty = c(2,1,2,1,2),lwd=c(2,2,2,2,2))
dev.off()
missing_terms <- c(alb_term,bili_term,pafi_term,crea_term)
names(missing_terms) <- c("alb","bili","pafi","crea")
posterior_summaries <- cbind("term" = rep(c("alb","bili","pafi","crea"),5),
                             "Method" = rep(c("Single_imp","Rubin_t","Rubin_norm","ABpool_t","ABpool_norm"),each=4),
                             "Mean" = c(coef(support.overall.cox)[missing_terms],
                                        support.cox.pool.t$pooled$estimate[missing_terms],
                                        support.cox.pool.norm$pooled$estimate[missing_terms],
                                        apply(support.cox.abpool.t,2,mean)[missing_terms],
                                        apply(support.cox.abpool.norm,2,mean)[missing_terms]),
                             "Median" = c(coef(support.overall.cox)[missing_terms],
                                          support.cox.pool.t$pooled$estimate[missing_terms],
                                          support.cox.pool$pooled$estimate[missing_terms],
                                          apply(support.cox.abpool.t,2,median)[missing_terms],
                                          apply(support.cox.abpool.norm,2,median)[missing_terms]),
                             "2.5%" = c(t(confint(support.overall.cox))[1,missing_terms],
                                        round(get.ci.mice(support.cox.pool.t),4)[1,missing_terms],
                                        round(get.ci.mice(support.cox.pool.norm),4)[1,missing_terms],
                                        round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4)[1,missing_terms],
                                        round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4)[1,missing_terms]),
                             "97.5%" = c(t(confint(support.overall.cox))[2,missing_terms],
                                         round(get.ci.mice(support.cox.pool.norm),4)[2,missing_terms],
                                         CIs_Rubin_norm[2,missing_terms],
                                         round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4)[2,missing_terms],
                                         round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4)[2,missing_terms])
)
posterior_summaries[,c("Method","term")] <- as.character(posterior_summaries[,c("Method","term")])
posterior_summaries[,c("Mean","Median","2.5%","97.5%")] <- format(round(as.numeric(as.matrix(posterior_summaries[,c("Mean","Median","2.5%","97.5%")])),4),4)
posterior_summaries <- posterior_summaries[c(which(posterior_summaries[,"term"]=="alb"),
                                             which(posterior_summaries[,"term"]=="bili"),
                                             which(posterior_summaries[,"term"]=="pafi"),
                                             which(posterior_summaries[,"term"]=="crea")),]

missing_var_outputs <- data.frame("dfcom" = support.cox.pool.norm$pooled[missing_terms,"dfcom"],
                                  "df_Rubin_norm" = support.cox.pool.norm$pooled[missing_terms,"df"],
                                  "df_Rubin_t" = support.cox.pool.t$pooled[missing_terms,"df"],
                                  "fmi" = fmis[missing_terms],
                                  "skew_MLE" = skew_MLE[missing_terms],
                                  "kurt_MLE" = kurt_MLE[missing_terms],
                                  "skew_post" = skew_par[missing_terms],
                                  "kurt_post" = kurt_par[missing_terms],
                                  "skew_U"   = skew_U[missing_terms],
                                  "kurt_U"   = kurt_U[missing_terms])

outputs <- list("Posterior_summaries" = posterior_summaries,
                "missing_var_outputs" = missing_var_outputs,
                "CI_Rubin_t" = round(get.ci.mice(support.cox.pool.t),4),
                "CI_Rubin_norm" = round(get.ci.mice(support.cox.pool.norm),4),
                "CI_ABpool_t" = round(apply(support.cox.abpool.t,2,quantile,c(0.025,0.975)),4),
                "CI_ABpool_norm" = round(apply(support.cox.abpool.norm,2,quantile,c(0.025,0.975)),4),
                "estimates_Rubin" = support.cox.pool.norm$pooled$estimate,
                "estimates_ABpool_t" = apply(support.cox.abpool.t,2,mean),
                "estimates_ABpool_norm" = apply(support.cox.abpool.norm,2,mean))

saveRDS(outputs,paste0(output_directory,"/SUPPORT_outputs_results.RDS"))
outputs <- readRDS(paste0(output_directory,"/SUPPORT_outputs_results.RDS"))

write.csv(outputs$Posterior_summaries,
          paste0(output_directory,"/SUPPORT_results.csv"),
          row.names = FALSE)

write.csv(format(round(outputs$missing_var_outputs,2),2),
          paste0(output_directory,"/SUPPORT_supplementary_info.csv"),
          row.names = TRUE)

# times100000 <- list(pmm100000,cox100000,rubin100000,abmit100000,abminorm100000)
# saveRDS(times100000,paste0(output_directory,"/SUPPORT_times100000.RDS"))
# times100000 <- readRDS(paste0(output_directory,"/SUPPORT_times10000.RDS"))

# Cornish-Fisher
# Test Cornish-Fisher where average U is used to sample ABMI

CF_calcs <- rbind(skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[1],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[2],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[3],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[4],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[1],alpha=0.975),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[2],alpha=0.975),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[3],alpha=0.975),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool.norm,support.cox.abpool.norm,term=missing_terms[4],alpha=0.975))
rownames(CF_calcs) <- paste0(rep(names(missing_terms),2),rep(c("_2.5%","_97.5%"),each=4))

write.csv(round(CF_calcs,6),
          paste0(output_directory,"/SUPPORT_skew_kurt_CornishFisher.csv"),
          row.names = TRUE)
CF_calcs <- read.csv(paste0(output_directory,"/SUPPORT_skew_kurt_CornishFisher.csv"),row.names = 1)
