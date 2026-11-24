# Functions
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

