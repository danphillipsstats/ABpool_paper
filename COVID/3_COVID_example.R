##############
# This file uses the output from the Cox model in Cox_infection_model.R
# to calculate the plots, figures and results needed for the paper.

require(mitools)
require(RColorBrewer)
require(scales)
require(moments)
require(mvtnorm)

outcome_cols <- palette("Okabe-Ito")[c(2,7,3,6)]
names(outcome_cols) <- c("ABpool_t","ABpool_norm","Rubin_t","Rubin_norm")

# Create a list to store the time taken to run different parts of the file
runtimes <- list()
t_init <- Sys.time()
print(t_init)
runtimes$start <- t_init

# Set what the outcome for this analysis will be 
# (pos = all COVID-19 infections, prim = primary symptomatic COVID-19 infections)
# we use prim

# Set which type of antibodies we are using
antibody_type <- c("S","neuts")[1]
long_type <- c("7_day_incubator","0_day_incubator")[1]
# The name of the longitudinal model
long_file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_t0PB28"
cox_model_name <- "prim_site_parallel_simple_7inc"
outcome_type <- "symptomatic"
outcome_title <- "Primary symptomatic cases"

model_name <- "submission"
if(!(substring(long_file_name,1,5)=="neuts") == (antibody_type=="neuts")){warning("Your long_file_name and antibody_type do not match. You may be running a neuts analysis on S data or vice versa. Please double check.")}

# Directories
# Define the directories in which files are saved - these may need adjusting for other users.
directory <- getwd() 
directory_project <- paste0(directory,"/COVID")
output_directory <- paste0(directory_project,"/4_Output")
plot_directory <- paste0(directory_project,"/4_Output/Plots")

#load packages
library(rstan)
library(survival)
library(parallel)

#######
# Load user-written functions
source(paste0(directory_project,"/functions.R"))
################################################################################
# Two-stage model, second stage
n.outcomes <- 2
out <- list()
cox_model_name_symptomatic <- paste0("_",cox_model_name,"_",model_name)
cox_stage_two_symptomatic <- readRDS(paste0(output_directory,"/Cox_infection_model_long_",long_file_name,"_",cox_model_name,".RDS"))
# Read the output from Posterior_estimates.R
out <- readRDS(paste0(output_directory,"/Posterior_estimates_cox_infection_model_long_",long_file_name,"_",cox_model_name,"_",model_name,".RDS"))
# Save the entries which don't depend on outcome separately
data_long <- out$data_long
pred_antibody_values_28_hist <- out$pred_antibody_values_28_hist
pred_antibody_values_182_hist <- out$pred_antibody_values_182_hist
tens_BAU <- out$tens_BAU
tens_BAU2 <- tens_BAU[which(tens_BAU>1)]
antibody_quantiles <- out$antibody_quantiles
VE_vs_time_omicron <- out$VE_vs_time_omicron
n <- nrow(cox_stage_two_symptomatic$cox_model_mat)
n_pred <- ncol(cox_stage_two_symptomatic$cox_model_mat)

# Longitudinal model
# For an output which is a list, first element being (possibly thinned) extracted samples
# the second element being the data
# Take a subsample of size
nsamples <- dim(cox_stage_two_symptomatic$cox_model_pred)[1]
nchains <- dim(cox_stage_two_symptomatic$cox_model_pred)[2]
long_out <- readRDS(paste0(output_directory,"/Correlates_long_",long_file_name,"_noreffects.RDS"))
# File is a list, first entry is the Stan output, second entry the inputted data.
data_long <- long_out[[2]]
rm(long_out)
runtimes$a_array <- Sys.time();

###############
# Read joint_correlates and long_correlates
# Set their factors to be correct
source(paste0(directory_project,"/read_data_set_factors.R"))
# Create long_correlates to be the long_correlates dataset with the appropriate antibodies and log antibodies
long_correlates <- long_correlates_S
long_correlates$antibody <- long_correlates$S.PPD
long_correlates$log_antibody <- log(long_correlates$antibody)
# Limit of detection
LOD <- 33
# Conversion factor from arbitrary units AU/mL in our data to international unit BAU/mL
conv_factor <- 0.00645
print(LOD*conv_factor)
antibody_name <- "anti-spike IgG"
Antibody_name <- "Anti-spike IgG"
antibody_units <- "BAU/mL"
# Axis limits for plots of log antibodies
# Hard-coding this in so can more easily compare models.
log_antibody_lims <- c(-0.092,14.8)
# log_antibody_lims <- c(min(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]),max(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]))
# Marking outliers
long_correlates$outlier <- F
# Exclude the exceptionally small results as they appear to be outliers and this appears to be preventing good mixing in the model (rho = -1)
long_correlates$outlier[which(long_correlates$log_antibody<7 & long_correlates$As_vaccinated_arm_2=="ChAdOx1")] <- T
# Exclude the exceptionally large result at PB90
long_correlates$outlier[which((long_correlates$antibody_time>50 &
                                 long_correlates$antibody_time<130 &
                                 long_correlates$log_antibody>13) &
                                long_correlates$As_vaccinated_arm_2=="ChAdOx1")] <- T

## Left-censoring below the limit of detection
# Check how many antibody observations are less than, or equal to, the limit of detection, and treat them as left-censored
nantibody_LOD <- sum(long_correlates$log_antibody==log(LOD))
nantibody_leq_LOD <- sum(long_correlates$log_antibody<=log(LOD))
if(nantibody_LOD<nantibody_leq_LOD & nantibody_LOD>0){suppressWarnings(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
if(nantibody_LOD==nantibody_leq_LOD & nantibody_LOD>0){warning(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
long_correlates$left_censored <- FALSE
long_correlates$left_censored[which(long_correlates$antibody<LOD)] <- TRUE
long_correlates$antibody[which(long_correlates$left_censored)] <- LOD
long_correlates$log_antibody[which(long_correlates$left_censored)] <- log(LOD)
long_correlates <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),] 

# Defining variables
long_correlates <- long_correlates[,c("redcap_data_access_group","As_vaccinated_arm_2","sc_repeat_pid","record_id","sc_age","age_group","cor2dose_positive","cor2dose_primary","cor2dose_outcome_pos_prim","antibody","log_antibody","antibody_time","ll","timepoint","visit2","left_censored","outlier")]

# times
# Final antibody time we will report
amax_time <- 209-7*(long_type=="7_day_incubator") # Minus 7 days if we're using the 7_day_incubator method
# Earliest antibody time we will report
amin_time <- 28
adsb <- amin_time:amax_time # Days since boost (ds_t0+21)
if(long_type=="7_day_incubator") {
  t0 <- 21 + data_long$t_center
  # VE time is 7 days after antibody time due to assuming a 7 day gap between exposure and reported infection
  VEdsb <- adsb + 7 # Days since boost (ds_t0+21)
} else if(long_type=="0_day_incubator") {t0 <- 14 + data_long$t_center; VEdsb <- adsb} # Days since boost (ds_t0+21)
# Number of days from boost to day 0 in the model (i.e. 28)
ds_t0 <- adsb-t0 # Days since day 0 in the model for the antibodies
nts <- length(ds_t0)
time_units <- c("Days","Weeks","Months")[1]
# div is the Amount to divide days by to get the time units we're using (1 if days, 7 if weeks, 365/12 if months)
if (time_units=="Days"){div <- 1}
if (time_units=="Weeks"){div <- 7}
if (time_units=="Months"){div <- 365/12}
atsb <- adsb/div # time since boost in the appropriate units
VEtsb <- VEdsb/div # time since boost in the appropriate units

# The data for vaccinated individuals only
joint_correlates_ChAd <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),]

n_event <- sum(joint_correlates$cor2dose_primary_ind==1)

# Generate the antibody values needed to measure VE vs antibody
quants <- c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
nquants <- length(quants)
# Convert BAU/ml antibody levels of 1,10,...,10^5 to log AU/ml units
tens_BAU <- 10^seq(-1,7,1)
tens <- log(10^seq(-1,7,1)/conv_factor)
# Exclude those less than the smallest observed value or greater than the largest value
which_tens <- which(tens>min(long_correlates$log_antibody) & tens<max(long_correlates$log_antibody))
tens <- tens[which_tens]
tens_BAU <- tens_BAU[which_tens]

# Set the antibody values at which to test VE vs antibody
# Observe log antibody levels at an even sequence from min antibodies to max antibodies,
# and at the quantiles of the observed antibody levels, and at the powers of 10
log_antibody_values <- c(seq(from=min(long_correlates$log_antibody), to=max(long_correlates$log_antibody), by = 0.1),quantile(long_correlates$log_antibody,c(0,quants,1)))
log_antibody_values <- log_antibody_values[!duplicated(log_antibody_values)]
antibody_values <- exp(log_antibody_values)
BAU_antibodies <- antibody_values*conv_factor
BAU_antibodies <- c(BAU_antibodies,tens_BAU)
log_antibody_values <- c(log_antibody_values,tens)
log_antibody_values <- log_antibody_values[!duplicated(log_antibody_values)] # Remove duplicated values
log_antibody_values <- log_antibody_values[order(log_antibody_values)]
antibody_values <- exp(log_antibody_values)
BAU_antibodies <- BAU_antibodies[!duplicated(BAU_antibodies)]
BAU_antibodies <- BAU_antibodies[order(BAU_antibodies)]
nant_val <- length(log_antibody_values)

time_axis_lab <- paste0(time_units," since second dose vaccination")
antibody_lab <- paste0(Antibody_name," levels (",antibody_units,")")
time_axis_marks <- c(0,30,60,90,120,150,180,210) # Choose the marks for the x axis in weeks
long_time_axis_marks <- c(time_axis_marks,240,270,300,330,360)
# Place the marks for the y-axis in powers of 10 BAU/mL converted to the log AU/mL antibody scale
antibody_ax_mark_loc <- log(10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=1)/conv_factor)
# Write marks for the y-axis in powers of 10 BAU/mL
antibody_ax_mark <- 10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=1)
VE_lab <- "Vaccine efficacy (%)"
VE_lims <- c(-50,100)
VE_axis_marks <- c(seq(-60,100,20))
VE_axis_ticks <- c(seq(-60,100,10))
# Define colours
blue_for_histogram <- rgb(0.583,0.830,0.956,alpha=0.6)
green_for_histogram <- rgb(0.875,0.933,0.600,alpha=0.6)

method_colours <- palette("Okabe-Ito")[c(2,7)]
method_colours <- colorRampPalette(brewer.pal(5,"Spectral"))(5)
method_colours <- method_colours[c(1,length(method_colours))]

#########################################
# Calculate using symptomatic cases

# Complete data degrees of freedom
df.com <- n_event-n_pred
# Generate correlated gamma and zeta from ABpool
set.seed(1234)
Z0_norm <- matrix(rnorm(n=nsamples*nchains),nrow=nsamples,ncol=nchains) # Standard normal random variables
Z1_norm <- matrix(rnorm(n=nsamples*nchains),nrow=nsamples,ncol=nchains)
Chisq_sample <- matrix(rchisq(n=nsamples*nchains,df = df.com),nrow=nsamples,ncol=nchains)
Z0_t <- Z0_norm*sqrt(df.com/Chisq_sample)
Z1_t <- Z1_norm*sqrt(df.com/Chisq_sample)
# Calculate the correlation
cor_gamma_zeta <- cox_stage_two_symptomatic$cox_model_var[,,"antibody.As_vaccinated_arm_2ChAdOx1"]/sqrt(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"]*cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])
# Transform the standard normal random variables to normal random variables with correct mean, and covariances
gammas_ABpool_norm <- cox_stage_two_symptomatic$cox_model_pred[,,"antibody"] + Z0_norm*sqrt(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"])
zetas_ABpool_norm <- cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"] + (Z0_norm*cor_gamma_zeta + Z1_norm*sqrt(1-cor_gamma_zeta^2))*sqrt(cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])
# Transform the t random variables to t variables with correct mean, and covariances
gammas_ABpool_t <- cox_stage_two_symptomatic$cox_model_pred[,,"antibody"] + Z0_t*sqrt(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"])
zetas_ABpool_t <- cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"] + (Z0_t*cor_gamma_zeta + Z1_t*sqrt(1-cor_gamma_zeta^2))*sqrt(cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])

# Overall summary - norm
VE_samples_norm <- 100*(1 - exp(outer(c(gammas_ABpool_norm),antibody_values) + c(zetas_ABpool_norm))) # Given as a percentage
VE_summary_norm <- cbind(colMeans(VE_samples_norm),t(apply(VE_samples_norm,2,quantile,c(0.5,0.025,0.975))))
rm(VE_samples_norm)
rownames(VE_summary_norm) <- BAU_antibodies
colnames(VE_summary_norm) <- c("mean","median","lower95CI","upper95CI")

VE_samples_norm <- 100*(1 - exp(gammas_ABpool_norm*antibody_values[length(antibody_values)] + zetas_ABpool_norm)) # Given as a percentage

# Overall summary - t
VE_samples_t <- 100*(1 - exp(outer(c(gammas_ABpool_t),antibody_values) + c(zetas_ABpool_t))) # Given as a percentage
VE_summary_t <- cbind(colMeans(VE_samples_t),t(apply(VE_samples_t,2,quantile,c(0.5,0.025,0.975))))
rm(VE_samples_t)
rownames(VE_summary_t) <- BAU_antibodies
colnames(VE_summary_t) <- c("mean","median","lower95CI","upper95CI")

VE_samples_t <- 100*(1 - exp(gammas_ABpool_t*antibody_values[length(antibody_values)] + zetas_ABpool_t)) # Given as a percentage

#######
# Compare for antibody coefficient
skew_MLE <- skewness(c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]))
skew_par <- skewness(c(gammas_ABpool_norm))
kurt_MLE <- kurtosis(c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]))
kurt_par <- kurtosis(c(gammas_ABpool_norm))

var_between <- var(c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]))
var_within <- mean(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"]) 
skewness(c(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"])) # these are even more skewed than the estimates
lambda <- var_between/(var_within+var_between)
c(skew_par,skew_MLE*lambda^(3/2))
c(kurt_par,kurt_MLE*lambda^2)
var_U <- skewness(c(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"])) 
# The approximations of skew and kurtosis from Cornish-Fisher with fixed variance are bad, perhaps due to the variances U_i varying

# Calculate the CIs from Rubin's rules and ABpool
antibody_effect_Rubin_t <- rubin(c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]),
                                 c(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"]),
                                 df.com)
antibody_effect_Rubin_norm <- rubin(c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]),
                                    c(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"]),
                                    Inf)
c(antibody_effect_Rubin_t$estimate,antibody_effect_Rubin_t$ci)
c(antibody_effect_Rubin_norm$estimate,antibody_effect_Rubin_norm$ci)
quantile(c(gammas_ABpool_norm),c(0.5,0.025,0.975))
# Quite different CIs
k <- length(c(gammas_ABpool_norm))
p <- ppoints(k)
rubin_q <- function(p,rubin_out){
  rubin_out$estimate + qt(p, rubin_out$df) * sqrt(rubin_out$variance)
}
q_rubin_gamma <- rubin_q(ppoints(k),antibody_effect_Rubin_norm)
# Rubin's rules for zeta
vaccine_effect_Rubin_norm <- rubin(c(cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"]),
                                   c(cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"]),
                                   Inf)
q_rubin_zeta <- rubin_q(ppoints(k),vaccine_effect_Rubin_norm)


png(paste0(plot_directory,"/COVID_QQ_plot.png"),width = 1400,height=1000, pointsize = 23)
par(mfrow=c(2,2), mar = c(2.5,2.5,2.5,0.1),oma=c(4.1,4.1,0.1,1))
qqnorm(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"], main = "MLEs for gamma after imputation")
qqline(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"])
qqnorm(cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"], main = "MLEs for zeta after imputation")
qqline(cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"])
qqplot(x = q_rubin_gamma, y = gammas_ABpool_norm[order(gammas_ABpool_norm)],
       xlab = "Quantiles from Rubin's t-approximation",
       ylab = "Posterior samples from ABpool",
       main = "Diagnostic - gamma")
qqline(y = gammas_ABpool_norm[order(gammas_ABpool_norm)], 
       distribution = function(p) rubin_q(p,antibody_effect_Rubin_norm))
qqplot(x = q_rubin_zeta, y = zetas_ABpool_norm[order(zetas_ABpool_norm)],
       xlab = "Quantiles from Rubin's t-approximation",
       ylab = "Posterior samples from ABpool",
       main = "Diagnostic - zeta")
qqline(y = zetas_ABpool_norm[order(zetas_ABpool_norm)], 
       distribution = function(p) rubin_q(p,vaccine_effect_Rubin_norm))
dev.off()

methnames <- c("ABpool_norm","ABpool_t","Rubin_norm","Rubin_t")
png(paste0(plot_directory,"/COVID_gamma_distribution_compare_Rubin.png"),width = 1400,height=1000, pointsize = 23)
plot(density(gammas_ABpool_norm))
plot(density(gammas_ABpool_t), main = "Distribution of gamma (antibody effect)",xlab="theta_1",lwd=3, col= NA)
lines(density(gammas_ABpool_norm), col= outcome_cols["ABpool_norm"],lwd=2)
lines(density(gammas_ABpool_t), col= outcome_cols["ABpool_t"], lty=2,lwd=2)
xx <- seq(min(gammas_ABpool_norm)-1,max(gammas_ABpool_norm)+1,length.out=500000)
xt <- (xx-antibody_effect_Rubin_t$estimate)/sqrt(antibody_effect_Rubin_t$variance)
tt <- dt(xt,df=antibody_effect_Rubin_t$df)/sqrt(antibody_effect_Rubin_t$variance) # Have to divide by the sd to correct the density as we are using the location-scale t distribution
nn <- dt(xt,df=antibody_effect_Rubin_norm$df)/sqrt(antibody_effect_Rubin_norm$variance)
lines(nn~xx, col = outcome_cols["Rubin_norm"],lwd=2)
lines(tt~xx, col = outcome_cols["Rubin_t"],lty=2,lwd=2)
xn <- (xx-antibody_effect_Rubin_norm$estimate)/sqrt(antibody_effect_Rubin_norm$variance)
legend(x=-3.5e-04,y=4000,
       legend = methnames,
       col = outcome_cols[methnames],
       lty = c(1,2,1,2),lwd=c(2,2,2,2))
dev.off()

#########
# Rubin's rules
ci_antibody_values_Rubin_norm <- matrix(nrow=length(antibody_values),ncol=3)
ci_antibody_values_Rubin_t <- matrix(nrow=length(antibody_values),ncol=3)
for (i in 1:length(antibody_values)){
  means_vec <- c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]*antibody_values[i]+cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"])
  ci_antibody_values_Rubin_norm[i,1] <- mean(means_vec);ci_antibody_values_Rubin_t[i,1] <- mean(means_vec);
  vcov_est <- c(
    cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"]*(antibody_values[i])^2+
      2*cox_stage_two_symptomatic$cox_model_var[,,"antibody.As_vaccinated_arm_2ChAdOx1"]*antibody_values[i]+
      cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])
  # Calculate confidence interval using Barnard and Rubin's approach
  df.com <- n_event-n_pred
  rubin_pool_VE_i_norm <- rubin(means_vec,vcov_est,dfcom = Inf)
  rubin_pool_VE_i_t <- rubin(means_vec,vcov_est,dfcom = df.com)
  ci_antibody_values_Rubin_norm[i,c(2,3)] <- rubin_pool_VE_i_norm$ci
  ci_antibody_values_Rubin_t[i,c(2,3)] <- rubin_pool_VE_i_t$ci
}
VE_summary_rubin_t <- (1-exp(ci_antibody_values_Rubin_t[,c(1,3,2)]))*100
colnames(VE_summary_rubin_t) <- c("mean","lower95CI","upper95CI")
rownames(VE_summary_rubin_t) <- BAU_antibodies
VE_summary_rubin_norm <- (1-exp(ci_antibody_values_Rubin_norm[,c(1,3,2)]))*100
colnames(VE_summary_rubin_norm) <- c("mean","lower95CI","upper95CI")
rownames(VE_summary_rubin_norm) <- BAU_antibodies

png(paste0(plot_directory,"/COVID_VE_vs_antibody_Rubins_rules_symptomatic_compare.png"),width = 1400,height=1000, pointsize = 23)
par(mar=c(5.1,4.1,1.1,1.1))
plot(VE_summary_t[,"lower95CI"]~BAU_antibodies, lty=2, log="x", col="NA",
     type="l", xlab = antibody_lab, yaxt="n", ylim=VE_lims, xaxt="n",
     ylab = VE_lab, main = "Vaccine efficacy against COVID-19")
abline(v = 10^seq(-1,4,0.5),col="grey95")
abline(h = seq(-60,100,10),col="grey95")
abline(v = 10^seq(-1,4,1),col="grey90")
abline(h = seq(-60,100,20),col="grey90")
axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
log10_antibodies <- log10(BAU_antibodies)
min_anti <- floor(min(log10_antibodies))
max_anti <- ceiling(max(log10_antibodies))
axis(1,at=10^seq(0,max_anti,1),lab=gsub(" ","",format(10^seq(0,max_anti,1),scientific=F),fixed=T))

# Add a histogram of predicted antibody values at end of study day 209 and of observed antibody values
plot(pred_antibody_values_182_hist,add=TRUE,border=green_for_histogram,col=green_for_histogram)
plot(pred_antibody_values_28_hist,add=TRUE,border=green_for_histogram,col=blue_for_histogram)
# Plot rubin's rules
polygon(c(BAU_antibodies,rev(BAU_antibodies)), c(VE_summary_rubin_t[,"upper95CI"], rev(VE_summary_rubin_t[,"lower95CI"])),
        col = alpha(outcome_cols["Rubin_t"],0.3),lty=0)
polygon(c(BAU_antibodies,rev(BAU_antibodies)), c(VE_summary_rubin_norm[,"upper95CI"], rev(VE_summary_rubin_norm[,"lower95CI"])),
        col = alpha(outcome_cols["Rubin_norm"],0.3),lty=0)

polygon(c(BAU_antibodies,rev(BAU_antibodies)), c(VE_summary_t[,"upper95CI"], rev(VE_summary_t[,"lower95CI"])),
        col = alpha(outcome_cols["ABpool_t"],0.3),lty=0)

lines(VE_summary_rubin_t[,"mean"]~BAU_antibodies, col = outcome_cols["Rubin_t"], lwd=2)
lines(VE_summary_rubin_t[,"lower95CI"]~BAU_antibodies,lty=2, col = outcome_cols["Rubin_t"], lwd=1.5)
lines(VE_summary_rubin_t[,"upper95CI"]~BAU_antibodies,lty=2, col = outcome_cols["Rubin_t"], lwd=1.5)

lines(VE_summary_rubin_norm[,"mean"]~BAU_antibodies, col = outcome_cols["Rubin_norm"], lwd=2)
lines(VE_summary_rubin_norm[,"lower95CI"]~BAU_antibodies,lty=2, col = outcome_cols["Rubin_norm"], lwd=1.5)
lines(VE_summary_rubin_norm[,"upper95CI"]~BAU_antibodies,lty=2, col = outcome_cols["Rubin_norm"], lwd=1.5)

lines(VE_summary_t[,"median"]~BAU_antibodies, col = outcome_cols["ABpool_t"], lwd=2)
lines(VE_summary_t[,"lower95CI"]~BAU_antibodies,lty=2, col = outcome_cols["ABpool_t"], lwd=1.5)
lines(VE_summary_t[,"upper95CI"]~BAU_antibodies,lty=2, col = outcome_cols["ABpool_t"], lwd=1.5)

lines(VE_summary_norm[,"median"]~BAU_antibodies, col = outcome_cols["ABpool_norm"], lwd=2)
lines(VE_summary_norm[,"lower95CI"]~BAU_antibodies,lty=2, col = outcome_cols["ABpool_norm"], lwd=1.5)
lines(VE_summary_norm[,"upper95CI"]~BAU_antibodies,lty=2, col = outcome_cols["ABpool_norm"], lwd=1.5)

abline(v=LOD*conv_factor,lty=2,lwd=2)
abline(h=0,lty=2)
legend(x=10^2,y=-20,
       legend = methnames,
       col = outcome_cols[methnames],
       lty = 2,lwd=c(2,2,2,2))
dev.off()

# skew_kurt_CF_calc(MLEs = c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"]), 
#                   Us = c(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"]),
#                   rubin_pool = )

# Cornish-Fisher

MLEs <- c(cox_stage_two_symptomatic$cox_model_pred[,,"antibody"])
Us <- c(cox_stage_two_symptomatic$cox_model_var[,,"antibody.antibody"])

b_m <- var(MLEs)
u_m <- mean(Us)
nu <- df.com
nu_nu_2 <- nu/(nu-2)
df_rubin <- antibody_effect_Rubin_t$df
total_var <- b_m + u_m*nu/(nu-2)
lambda_nu <- exp(log(b_m)-log(total_var))
lambda <- exp(log(b_m)-log(total_var))
z_a <- qnorm(0.025)
Cov_MLEsq_U <- cov((MLEs-mean(MLEs))^2,Us)
sd_rubin <- sqrt(antibody_effect_Rubin_t$variance*df_rubin/(df_rubin-2))
sd_par <- sd(gammas_ABpool_t)
skew_par <- skewness(gammas_ABpool_t)
kurt_par <- kurtosis(gammas_ABpool_t)
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
calculations_info_gamma <- c(skew_calc_terms,
                             sum(skew_calc_terms),
                             skew_par,
                             kurt_calc_terms,
                             sum(kurt_calc_terms),
                             kurt_par,
                             CF_t_terms/conv_factor,
                             sum(CF_t_terms)/conv_factor,
                             (quantile(gammas_ABpool_t,c(0.025)) - antibody_effect_Rubin_t$ci[1])/conv_factor,
                             sum(CF_t_terms*c(-1,1,-1,-1))/conv_factor,
                             (quantile(gammas_ABpool_t,c(0.975)) - antibody_effect_Rubin_t$ci[2])/conv_factor)

names(calculations_info_gamma) <- c("skew_term_MLE","skew_term_Cov","skew_calc",
                                    "skew_par",
                                    "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc","kurt_par",
                                    "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                                    "CF_approx_t_lb","CI_diff_t_lb","CF_approx_t_ub","CI_diff_t_ub")

# Intercept
vaccine_effect_Rubin_norm <- rubin(c(cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"]),
                                   c(cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"]),
                                   Inf)
vaccine_effect_Rubin_t <- rubin(c(cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"]),
                                   c(cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"]),
                                   df.com)

MLEs <- c(cox_stage_two_symptomatic$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"])
Us <- c(cox_stage_two_symptomatic$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])

b_m <- var(MLEs)
u_m <- mean(Us)
nu <- df.com
nu_nu_2 <- nu/(nu-2)
df_rubin <- vaccine_effect_Rubin_t$df
total_var <- b_m + u_m*nu/(nu-2)
lambda_nu <- exp(log(b_m)-log(total_var))
z_a <- qnorm(0.025)
Cov_MLEsq_U <- cov((MLEs-mean(MLEs))^2,Us)
sd_rubin <- sqrt(vaccine_effect_Rubin_t$variance*df_rubin/(df_rubin-2))
sd_par <- sd(zetas_ABpool_t)
skew_par <- skewness(zetas_ABpool_t)
kurt_par <- kurtosis(zetas_ABpool_t)
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
calculations_info_zeta <- c(skew_calc_terms,
                            sum(skew_calc_terms),
                            skew_par,
                            kurt_calc_terms,
                            sum(kurt_calc_terms),
                            kurt_par,
                            CF_t_terms,
                            sum(CF_t_terms),
                            (quantile(zetas_ABpool_t,c(0.025)) - vaccine_effect_Rubin_t$ci[1]),
                            sum(CF_t_terms*c(-1,1,-1,-1)),
                            (quantile(zetas_ABpool_t,c(0.975)) - vaccine_effect_Rubin_t$ci[2]))

names(calculations_info_zeta) <- c("skew_term_MLE","skew_term_Cov","skew_calc",
                                   "skew_par",
                                   "kurt_term_MLE","kurt_term_df","kurt_term_varU","kurt_term_Cov","kurt_calc","kurt_par",
                                   "CF_term_var","CF_term_skew","CF_term_kurt","CF_term_skewsq",
                                   "CF_approx_t_lb","CI_diff_t_lb","CF_approx_t_ub","CI_diff_t_ub")
CF_skew_kurt <- round(rbind(calculations_info_gamma,calculations_info_zeta),5)
rownames(CF_skew_kurt) <- c("gamma","zeta")

results_summary <- matrix(nrow=4*4,ncol=6)
colnames(results_summary) <- c("Outcome","Method","Posterior mean","Posterior median","2.5%","97.5%")
results_summary[,"Outcome"] <- rep(c("Antibody effect","VE(10)","VE(100)","VE(1000)"),each=4)
results_summary[,"Method"] <- rep(c("ABpool_norm","ABpool_t","Rubin_norm","Rubin_t"),4)
results_summary[1:4,3:6] <- round(rbind(c(mean(gammas_ABpool_norm),quantile(c(gammas_ABpool_norm),c(0.5,0.025,0.975))),
                                        c(mean(gammas_ABpool_t),quantile(c(gammas_ABpool_t),c(0.5,0.025,0.975))),
                                        c(antibody_effect_Rubin_norm$estimate[c(1,1)],antibody_effect_Rubin_norm$ci),
                                        c(antibody_effect_Rubin_t$estimate[c(1,1)],antibody_effect_Rubin_t$ci))/conv_factor,5)
results_summary[c(5,9,13),3:6] <- VE_summary_norm[c("10","100","1000"),]
results_summary[1+c(5,9,13),3:6] <- VE_summary_t[c("10","100","1000"),]
results_summary[2+c(5,9,13),3:6] <- VE_summary_rubin_norm[c("10","100","1000"),c(1,1,2,3)]
results_summary[3+c(5,9,13),3:6] <- VE_summary_rubin_t[c("10","100","1000"),c(1,1,2,3)]
results_summary[1:4,3:6] <- format(as.numeric(results_summary[1:4,3:6]),nsmall=5)
results_summary[5:16,3:6] <- format(round(as.numeric(results_summary[5:16,3:6]),1),nsmall=1)
write.csv(results_summary,paste0(output_directory,"/COVID_results_summary.csv"),row.names = FALSE)

make_latex_table <- function(results_summary) {
  
  x <- as.matrix(results_summary)
  x <- trimws(x)
  
  # Escape underscores in method names
  x[, "Method"] <- gsub("_", "\\_", x[, "Method"], fixed = TRUE)
  
  # Show outcome only on the first row of each group
  group_start <- c(TRUE, x[-1, "Outcome"] != x[-nrow(x), "Outcome"])
  x[!group_start, "Outcome"] <- ""
  
  # Format table rows
  rows <- apply(x, 1, function(z) {
    paste0(paste(z, collapse = " & "), " \\\\")
  })
  
  # Add horizontal lines between outcome groups
  rows <- ifelse(group_start, paste0("\\hline\n", rows), rows)
  
  # Construct complete LaTeX table
  latex <- paste0(
    "\\captionof{supptable}{\\textbf{Uncertainty intervals for antibody effect and vaccine efficacy at 10, 100, 1000 BAU/mL}}\n",
    "\t\\centering\n",
    "\t\\begin{tabular}{|c|c|cccc|}\n",
    "\t\t\\hline\n",
    "\t\tOutcome & Method & Mean & Median & 2.5\\% & 97.5\\% \\\\\n",
    "\t\t", paste(rows, collapse = "\n\t\t"), "\n",
    "\t\t\\hline\n",
    "\t\\end{tabular}\n",
    "\t\\\\\\vspace{0.2cm}\\flushleft\n",
    "\tMean and median are the posterior mean and median respectively, as estimated by each method. ",
    "Antibody effect is the effect on the log-hazard per increase of 1 BAU/mL. ",
    "Vaccine efficacy (VE) is given in percentages. ",
    "\\label{tab:COVID_results}"
  )
  
  return(latex)
}
latex_table <- make_latex_table(results_summary)
cat(latex_table)
