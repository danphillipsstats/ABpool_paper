##############################################################################
# SUPPORT example
##############################################################################
require(mice)
require(mvtnorm)
require(moments)
require(mixpoissonreg)
require(rstan)
require(bayesplot)
require(glmmTMB)
require(rjags)

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
# Cox model
###################################################
require(survival)
require(Hmisc)
require(survminer)
###################################################
# SUPPORT
# Could be promising

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
# Mention in supps - Chosen by backward deletion from the full model above, in a complete-case analysis (any individuals with a missing value for any variable included in the current model is excluded.)
support.formula <- Surv(d.time, death) ~ age + num.co + scoma + alb + bili + pafi + crea + meanbp + hrt + resp + temp + adlsc + dzgroup
support.vs.cox <- coxph(support.formula,data=support_filled_in)
summary(support.vs.cox)
#####
# support - bili: fmi=0.4, skew(thetahat) = -0.5, skew(theta) = -0.1, kurtosis(theta) = 0, theta not obviously non-normal
#           bili: no difference in CIs that I'm convinced is more than randomness (not noticeable anyway)
skewness(support$bili,na.rm=T); skewness(support$crea,na.rm=T); kurtosis(support$bili,na.rm=T); kurtosis(support$crea,na.rm=T)
cox.support.form.char <- "Surv(d.time, death) ~ age + num.co + scoma + alb + bili + pafi + crea + meanbp + hrt + resp + temp + adlsc + dzgroup"
(support.cox <- coxph(as.formula(cox.support.form.char),data=support))
summary(support.cox)

npred.support <- length(coef(support.cox))

support.data <- support[,all.vars(as.formula(cox.support.form.char))]
apply(support.data,2,function(vec) sum(is.na(vec)))/nrow(support.data)*100 # 37.8% for alb, 29.7% for bili, 25.3% for pafi, 0.3% for crea, all others fully observed
support.data.obs <- support.data[complete.cases(support.data),]


# Create dummy variables for each dzgroup to prepate for schoenfeld residuals
for(group in levels(support$dzgroup)[-1]){
  support[,group] <- as.numeric(support$dzgroup==group)
}
support <- data.frame(support)
# Create formula using the dummy variables
cox.support.form.dummy <- paste0(strsplit(cox.support.form.char,"dzgroup")[[1]][1],paste(make.names(levels(support$dzgroup))[-1],collapse=" + "))
(support.cox.dummy <- coxph(as.formula(cox.support.form.dummy),data=support))
# Schoenfeld residuals - testing time-constant effect
czph <- cox.zph(support.cox.dummy, transform="identity")
png(paste0(plot_directory,"/SUPPORT_schoenfeld_res_cc.png"),width = 1400,height=800, pointsize=23)
par(mfrow=c(4,5), mar = c(2.5,2.5,2.5,0.1),oma=c(1.1,1.1,0.1,1))
for (i in 1:length(support.cox.dummy$coefficients)){
  plot(czph[i],resid=F, main = names(support.cox.dummy$coefficients)[i], ylim=range(czph$y[,i]),col=NA)
  points(czph$x,czph$y[,i], col="grey",cex=0.5)
  abline(h=0,col="grey",lty=3)
  par(new=TRUE)
  plot(czph[i],resid=F, ylim=range(czph$y[,i]))
  
}
dev.off()
rownames(czph$table)

# Martingale residuals - testing functional form
cont_preds <- setdiff(colnames(support.data),c("d.time","death","dzgroup"))
n_cont_pred <- length(cont_preds) # = ncol(support.data)-3
png(paste0(plot_directory,"/SUPPORT_martingale_res_si.png"),width = 1400,height=800, pointsize=23)
par(mfrow=c(3,4), mar = c(4.1,2.5,1.5,0.1),oma=c(1.1,1.1,0.1,1))
for (i in 1:n_cont_pred){
  print(i)
  support.data_i <- support_filled_in[,c("d.time","death",cont_preds[-i], "dzgroup")] # Data excluding ith covariate
  cox_i <- coxph(Surv(d.time, death)~.,data=support.data_i)
  mart_resid <- residuals(cox_i, type="martingale")
  data_i <- data.frame("covariate" = support_filled_in[,cont_preds[i]],
                       "martingale_residual" = mart_resid)
  data_i <- data_i[order(data_i$covariate),]
  plot(martingale_residual~covariate,data=data_i,
       xlab = cont_preds[i],
       ylab = "Martingale residuals", cex=0.5,col="grey")
  loess_fit <- loess(martingale_residual~covariate,data=data_i)
  pred <- predict(loess_fit,se=TRUE)
  lines(data_i$covariate, pred$fit)
  lines(data_i$covariate, pred$fit + 2 * pred$se.fit, lty = 2)
  lines(data_i$covariate, pred$fit - 2 * pred$se.fit, lty = 2)
  # lines(lowess(mart_resid~support_filled_in[,cont_preds[i]]))
}
dev.off()
# Martingale residuals when log(bilirubin) is used instead. Either look reasonable.
{i <- 3
  support.data_i <- support_filled_in[,c("d.time","death",cont_preds[-which(cont_preds=="bili")], "dzgroup")] # Data excluding ith covariate
  cox_i <- coxph(Surv(d.time, death)~.,data=support.data_i)
  mart_resid <- residuals(cox_i, type="martingale")
  data_i <- data.frame("covariate" = log(support_filled_in[,cont_preds[which(cont_preds=="bili")]]),
                       "martingale_residual" = mart_resid)
  data_i <- data_i[order(data_i$covariate),]
  plot(martingale_residual~covariate,data=data_i,
       xlab = "log(bili)",
       ylab = "Martingale residuals", cex=0.5,col="grey")
  loess_fit <- loess(martingale_residual~covariate,data=data_i)
  pred <- predict(loess_fit,se=TRUE)
  lines(data_i$covariate, pred$fit)
  lines(data_i$covariate, pred$fit + 2 * pred$se.fit, lty = 2)
  lines(data_i$covariate, pred$fit - 2 * pred$se.fit, lty = 2)}

# Note since we're using predictive mean matching, there is only 80 odd values they might realistically choose.
# May be computationally better to impute say m = 1000 and then use J = 100 Cox models for each one
m <- 100000
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
imp_data_support <- readRDS(paste0(output_directory,"/imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
t1 <- Sys.time()
support.cox.imp <- with(imp_data_support,coxph(as.formula(cox.support.form.char)))
t2 <- Sys.time()
(cox100000 <- difftime(t2,t1))
saveRDS(support.cox.imp,paste0(output_directory,"/cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))
support.cox.imp <- readRDS(paste0(output_directory,"/cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS"))

head(support.cox.imp$analyses)
# Plot observed and imputed data density
thin <- 10
imp_1 <- c(as.matrix(round(imp_data_support$imp$alb[,seq(from=thin,to=dim(imp_data_support$imp$alb)[2],by=thin)],1)))
png(paste0(plot_directory,"/SUPPORT_imputed_observed_data_density_alb.png"),width = 1400,height=800, pointsize=23)
plot(density(support.data$alb,bw=0.5,na.rm=T),col=outcome_cols[2], xlab = "alb",main="",lwd=3) # blue observed data
lines(density(imp_1,bw=0.5),col=outcome_cols[1],lwd=3) # red imputed data
dev.off()

thin <- 10
imp_1 <- c(as.matrix(round(imp_data_support$imp$bili[,seq(from=thin,to=dim(imp_data_support$imp$bili)[2],by=thin)],1)))
png(paste0(plot_directory,"/SUPPORT_imputed_observed_data_density_bili.png"),width = 1400,height=800, pointsize=23)
plot(density(imp_1,bw=0.5),col=outcome_cols[1], xlab = "bili",main="",lwd=3) # red imputed data
lines(density(support.data$bili,bw=0.5,na.rm=T),col=outcome_cols[2],lwd=3) # blue observed data
dev.off()

thin <- 10
imp_1 <- c(as.matrix(round(imp_data_support$imp$pafi[,seq(from=thin,to=dim(imp_data_support$imp$pafi)[2],by=thin)],1)))
png(paste0(plot_directory,"/SUPPORT_imputed_observed_data_density_pafi.png"),width = 1400,height=800, pointsize=23)
plot(density(imp_1,bw=0.5),col=outcome_cols[1], xlab = "pafi",main="",lwd=3) # red imputed data
lines(density(support.data$pafi,bw=0.5,na.rm=T),col=outcome_cols[2],lwd=3) # blue observed data
dev.off()
rm(imp_data_support)

coef_list <- list()
for (j in 1:npred.support){
  coef_list[[names(coef(support.cox))[j]]] <- sapply(1:m, function(i) coef(support.cox.imp$analyses[[i]])[j])
  qqnorm(coef_list[[j]], main = names(coef(support.cox))[j])
  qqline(coef_list[[j]])
}
t1 <- Sys.time()
support.cox.pool <- mice::pool(support.cox.imp)
t2 <- Sys.time()
(rubin100000 <- difftime(t2,t1))

t1 <- Sys.time()
support.cox.abmi.norm <- t(vapply(1:m,FUN=sample_ABpool_norm,FUN.VALUE=numeric(length(coef(support.cox))), "mice.with.list"=support.cox.imp))
t2 <- Sys.time()
(abminorm100000 <- difftime(t2,t1))
t1 <- Sys.time()
support.cox.abmi.t <- t(vapply(1:m,FUN=sample_ABpool_t_cox,FUN.VALUE=numeric(length(coef(support.cox))), "mice.with.list"=support.cox.imp))
t2 <- Sys.time()
(abmit100000 <- difftime(t2,t1))
colnames(support.cox.abmi.norm) <- names(coef_list); colnames(support.cox.abmi.t) <- names(coef_list)

saveRDS(mget(c("support.cox.pool","support.cox.abmi.norm","support.cox.abmi.t")),
        paste0(output_directory,"/SUPPORT_output_objects.RDS"))
output_objects <- readRDS(paste0(output_directory,"/SUPPORT_output_objects.RDS"))
support.cox.pool <- output_objects$support.cox.pool
support.cox.abmi.norm <- output_objects$support.cox.abmi.norm
support.cox.abmi.t <- output_objects$support.cox.abmi.t
rm(output_objects)

round(get.ci.mice(support.cox.pool),4)
CIs_Rubin_norm <- round(rbind(support.cox.pool$pooled$estimate + c(-1,1)[1]*qt(0.975,(support.cox.pool$m-1)/support.cox.pool$pooled$fmi^2)*sqrt(support.cox.pool$pooled$t),
                              support.cox.pool$pooled$estimate + c(-1,1)[2]*qt(0.975,(support.cox.pool$m-1)/support.cox.pool$pooled$fmi^2)*sqrt(support.cox.pool$pooled$t)),4) # Normal estimate
colnames(CIs_Rubin_norm) <- support.cox.pool$pooled$term
round(apply(support.cox.abmi.norm,2,quantile,c(0.025,0.975)),4)
round(apply(support.cox.abmi.t,2,quantile,c(0.025,0.975)),4)

# Get skewness and fmi
fmis <- support.cox.pool$pooled$fmi
names(fmis) <- names(coef(support.cox))
fmis # alb and bili have decent sized fmi
support.MLEs <- t(vapply(1:m, function(i) coef(support.cox.imp$analyses[[i]]),numeric(length(coef(support.cox.imp$analyses[[1]])))))
support.Us <- t(vapply(1:m, function(i) diag(vcov(support.cox.imp$analyses[[i]])),numeric(length(coef(support.cox.imp$analyses[[1]])))))
(skew_MLE <- apply(support.MLEs,2,skewness)) # crea and bili have considerable skewness (-1.1,-0.5)
(skew_par <- apply(support.cox.abmi.t,2,skewness)) # but not for the posterior(-0.1,-0.1)
(kurt_MLE <- apply(support.MLEs,2,kurtosis))
(kurt_par <- apply(support.cox.abmi.t,2,kurtosis)) # All parameters have small excess kurtosis
(skew_U <- apply(support.Us,2,skewness)) # crea and bili have considerable skewness (-1.1,-0.5)
(kurt_U <- apply(support.Us,2,kurtosis))

# QQ plot of MLEs and posterior samples from just 100 imputations
thin <- 200
missing_vars <- c("alb","bili","pafi","crea")
png(paste0(plot_directory,"/SUPPORT_QQplot_missing200.png"),width = 1400,height=1000, pointsize=23)
# define a 4-row layout: title, 4 plots, title, 4 plots
lay <- matrix(c(
  1, 1, 1, 1,   # row title for top row (spans 4 columns)
  2, 3, 4, 5,   # 4 plots (top row)
  6, 6, 6, 6,   # row title for bottom row (spans 4 columns)
  7, 8, 9, 10   # 4 plots (bottom row)
), nrow = 4, byrow = TRUE)

layout(lay, heights = c(1, 6, 1, 6))

# top row title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.4, "Q-Q plot of MLEs", cex = 1.4, font = 2)

par(mar = c(4, 4, 2, 1))
for (i in 1:length(missing_vars)){
  vars <- missing_vars[i]
  qqnorm(support.MLEs[1:thin,vars], main = missing_vars[i],cex=0.4)
  qqline(support.MLEs[1:thin,vars])
}

# bottom row title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.4, "Q-Q plot of posterior samples", cex = 1.4, font = 2)

par(mar = c(4, 4, 2, 1))
for (i in 1:length(missing_vars)){
  vars <- missing_vars[i]
  qqnorm(support.cox.abmi.t[1:thin,vars], main = missing_vars[i],cex=0.4)
  qqline(support.cox.abmi.t[1:thin,vars])
}
dev.off()

varnames <- names(coef_list)
varnames <- sub("^dzgroup","",varnames) # removes the substring dzgroup from the start of the name
# Some MLEs noticeably non-normal
png(paste0(plot_directory,"/SUPPORT_QQplot_MLEs.png"),width = 1400,height=800, pointsize=23)
par(mfrow=c(4,5), mar = c(2.5,2.5,2.5,0.1),oma=c(0.1,0.1,0.0,1))
for (i in 1:length(names(coef_list))){
  vars <- names(coef_list)[i]
  qqnorm(support.MLEs[,vars], main = varnames[i],cex=0.4)
  qqline(support.MLEs[,vars])
}
dev.off()
# but not really the posterior samples
png(paste0(plot_directory,"/SUPPORT_QQplot_posterior_abmi_t.png"),width = 1400,height=800, pointsize=23)
par(mfrow=c(4,5), mar = c(2.5,2.5,2.5,0.1),oma=c(0.1,0.1,0.0,1))
for (i in 1:length(names(coef_list))){
  vars <- names(coef_list)[i]
  qqnorm(support.cox.abmi.t[,vars], main = varnames[i],cex=0.4)
  qqline(support.cox.abmi.t[,vars])
}
dev.off()

alb_term <- which(support.cox.pool$pooled$term=="alb")
bili_term <- which(support.cox.pool$pooled$term=="bili")
pafi_term <- which(support.cox.pool$pooled$term=="pafi")
crea_term <- which(support.cox.pool$pooled$term=="crea")

png(paste0(plot_directory,"/SUPPORT_compare_post_missing.png"),width = 1400,height=1000, pointsize=23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(4.1,4.1,0.1,1))
plot(density(support.cox.abmi.t[,"alb"]), main = "Distribution of albumin effect",xlab="Effect due to albumin", col=NA,ylim=c(0,6))
xx <- seq(min(support.cox.abmi.norm[,"alb"])-1,max(support.cox.abmi.norm[,"alb"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[alb_term],sqrt(vcov(support.overall.cox)[alb_term,alb_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abmi.t[,"alb"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abmi.norm[,"alb"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool$pooled$estimate[alb_term])/sqrt(support.cox.pool$pooled$t[alb_term])
tt <- dt(xt,df=support.cox.pool$pooled$df[alb_term])/sqrt(support.cox.pool$pooled$t[alb_term]) # Rubin's t - Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(support.cox.pool$m-1)/support.cox.pool$pooled$fmi[alb_term]^2)/sqrt(support.cox.pool$pooled$t[alb_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)

plot(density(support.cox.abmi.t[,"bili"]), main = "Distribution of bilirubin effect",xlab="Effect due to bilirubin", col=NA,ylim=c(0,46))
xx <- seq(min(support.cox.abmi.norm[,"bili"])-1,max(support.cox.abmi.norm[,"bili"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[bili_term],sqrt(vcov(support.overall.cox)[bili_term,bili_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abmi.t[,"bili"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abmi.norm[,"bili"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool$pooled$estimate[bili_term])/sqrt(support.cox.pool$pooled$t[bili_term])
tt <- dt(xt,df=support.cox.pool$pooled$df[bili_term])/sqrt(support.cox.pool$pooled$t[bili_term]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(support.cox.pool$m-1)/support.cox.pool$pooled$fmi[bili_term]^2)/sqrt(support.cox.pool$pooled$t[bili_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)

plot(density(support.cox.abmi.t[,"pafi"]), main = "Distribution of PaO2/FiO2 ratio effect",xlab="Effect due to PaO2/FiO2 ratio", col=NA)
xx <- seq(min(support.cox.abmi.norm[,"pafi"])-1,max(support.cox.abmi.norm[,"pafi"])+1,0.000001)
ss <- dnorm(xx,coef(support.overall.cox)[pafi_term],sqrt(vcov(support.overall.cox)[pafi_term,pafi_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abmi.t[,"pafi"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abmi.norm[,"pafi"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool$pooled$estimate[pafi_term])/sqrt(support.cox.pool$pooled$t[pafi_term])
tt <- dt(xt,df=support.cox.pool$pooled$df[pafi_term])/sqrt(support.cox.pool$pooled$t[pafi_term]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(support.cox.pool$m-1)/support.cox.pool$pooled$fmi[pafi_term]^2)/sqrt(support.cox.pool$pooled$t[pafi_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)

plot(density(support.cox.abmi.norm[,"crea"]), main = "Distribution of creatinine effect",xlab="Effect due to creatinine", col=NA,ylim=c(0,19))
xx <- seq(min(support.cox.abmi.norm[,"crea"])-1,max(support.cox.abmi.norm[,"crea"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[crea_term],sqrt(vcov(support.overall.cox)[crea_term,crea_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abmi.t[,"crea"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abmi.norm[,"crea"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool$pooled$estimate[crea_term])/sqrt(support.cox.pool$pooled$t[crea_term])
tt <- dt(xt,df=support.cox.pool$pooled$df[crea_term])/sqrt(support.cox.pool$pooled$t[crea_term]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(support.cox.pool$m-1)/support.cox.pool$pooled$fmi[crea_term]^2)/sqrt(support.cox.pool$pooled$t[crea_term])
lines(nn~xx, col = outcome_cols[2],lty=2,lwd=2)
dev.off()

png(paste0(plot_directory,"/SUPPORT_compare_post_bili.png"),width = 1400,height=1000, pointsize=23)
plot(density(support.cox.abmi.t[,"bili"]), main = "Distribution of bilirubin effect",xlab="Effect due to bilirubin", col=NA,ylim=c(0,46))
xx <- seq(min(support.cox.abmi.norm[,"bili"])-1,max(support.cox.abmi.norm[,"bili"])+1,0.0001)
ss <- dnorm(xx,coef(support.overall.cox)[bili_term],sqrt(vcov(support.overall.cox)[bili_term,bili_term])) # Fill-in single imputation
lines(ss~xx, col = "grey",lty=2,lwd=2)
lines(density(support.cox.abmi.t[,"bili"]), col= outcome_cols[1],lwd=2)
lines(density(support.cox.abmi.norm[,"bili"]), col= outcome_cols[1], lty=2,lwd=2)
xt <- (xx-support.cox.pool$pooled$estimate[bili_term])/sqrt(support.cox.pool$pooled$t[bili_term])
tt <- dt(xt,df=support.cox.pool$pooled$df[bili_term])/sqrt(support.cox.pool$pooled$t[bili_term]) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
lines(tt~xx, col = outcome_cols[2],lwd=2)
nn <- dt(xt,df=(support.cox.pool$m-1)/support.cox.pool$pooled$fmi[bili_term]^2)/sqrt(support.cox.pool$pooled$t[bili_term])
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
                                        support.cox.pool$pooled$estimate[missing_terms],
                                        support.cox.pool$pooled$estimate[missing_terms],
                                        apply(support.cox.abmi.t,2,mean)[missing_terms],
                                        apply(support.cox.abmi.norm,2,mean)[missing_terms]),
                             "Median" = c(coef(support.overall.cox)[missing_terms],
                                          support.cox.pool$pooled$estimate[missing_terms],
                                          support.cox.pool$pooled$estimate[missing_terms],
                                          apply(support.cox.abmi.t,2,median)[missing_terms],
                                          apply(support.cox.abmi.norm,2,median)[missing_terms]),
                             "2.5%" = c(t(confint(support.overall.cox))[1,missing_terms],
                                        round(get.ci.mice(support.cox.pool),4)[1,missing_terms],
                                        CIs_Rubin_norm[1,missing_terms],
                                        round(apply(support.cox.abmi.t,2,quantile,c(0.025,0.975)),4)[1,missing_terms],
                                        round(apply(support.cox.abmi.norm,2,quantile,c(0.025,0.975)),4)[1,missing_terms]),
                             "97.5%" = c(t(confint(support.overall.cox))[2,missing_terms],
                                         round(get.ci.mice(support.cox.pool),4)[2,missing_terms],
                                         CIs_Rubin_norm[2,missing_terms],
                                         round(apply(support.cox.abmi.t,2,quantile,c(0.025,0.975)),4)[2,missing_terms],
                                         round(apply(support.cox.abmi.norm,2,quantile,c(0.025,0.975)),4)[2,missing_terms])
)
posterior_summaries[,c("Method","term")] <- as.character(posterior_summaries[,c("Method","term")])
posterior_summaries[,c("Mean","Median","2.5%","97.5%")] <- format(round(as.numeric(as.matrix(posterior_summaries[,c("Mean","Median","2.5%","97.5%")])),4),4)
posterior_summaries <- posterior_summaries[c(which(posterior_summaries[,"term"]=="alb"),
                                             which(posterior_summaries[,"term"]=="bili"),
                                             which(posterior_summaries[,"term"]=="pafi"),
                                             which(posterior_summaries[,"term"]=="crea")),]

missing_var_outputs <- data.frame("dfcom" = support.cox.pool$pooled[missing_terms,"dfcom"],
                                  "df" = support.cox.pool$pooled[missing_terms,"df"],
                                  "fmi" = fmis[missing_terms],
                                  "skew_MLE" = skew_MLE[missing_terms],
                                  "kurt_MLE" = kurt_MLE[missing_terms],
                                  "skew_post" = skew_par[missing_terms],
                                  "kurt_post" = kurt_par[missing_terms],
                                  "skew_U"   = skew_U[missing_terms],
                                  "kurt_U"   = kurt_U[missing_terms])

outputs <- list("Posterior_summaries" = posterior_summaries,
                "missing_var_outputs" = missing_var_outputs,
                "CI_Rubin_t" = round(get.ci.mice(support.cox.pool),4),
                "CI_Rubin_norm" = CIs_Rubin_norm,
                "CI_ABpool_t" = round(apply(support.cox.abmi.t,2,quantile,c(0.025,0.975)),4),
                "CI_ABpool_norm" = round(apply(support.cox.abmi.norm,2,quantile,c(0.025,0.975)),4),
                "estimates_Rubin" = support.cox.pool$pooled$estimate,
                "estimates_ABpool_t" = apply(support.cox.abmi.t,2,mean),
                "estimates_ABpool_norm" = apply(support.cox.abmi.norm,2,mean))

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

CF_calcs <- rbind(skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[1],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[2],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[3],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[4],alpha=0.025),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[1],alpha=0.975),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[2],alpha=0.975),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[3],alpha=0.975),
                  skew_kurt_CF_calc(support.MLEs,support.Us,support.cox.pool,support.cox.abmi.t,term=missing_terms[4],alpha=0.975))
rownames(CF_calcs) <- paste0(rep(names(missing_terms),2),rep(c("_2.5%","_97.5%"),each=4))

write.csv(round(CF_calcs,6),
          paste0(output_directory,"/SUPPORT_skew_kurt_CornishFisher.csv"),
          row.names = TRUE)
CF_calcs <- read.csv(paste0(output_directory,"/SUPPORT_skew_kurt_CornishFisher.csv"),row.names = 1)
