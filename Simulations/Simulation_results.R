### Pilot simulation study results
# Load packages
t_init <- Sys.time()

# scales::show_col(palette("Okabe-Ito"))
scales::show_col(palette("Okabe-Ito")[c(1,4,2,7,3,6)])
outcome_cols <- palette("Okabe-Ito")[c(1,4,2,7,3,6)]
names(outcome_cols) <- c("CC","Bayes","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm")
pch_vals <- c(4,21, 23, 22, 24, 25)
names(pch_vals) <- c("CC","Bayes","ABpool_t","ABpool_norm","Rubin_t","Rubin_norm")
#######
project_directory <- getwd()
simulation_directory <- file.path(project_directory, "Simulations")
program_directory <- simulation_directory
output_directory <- file.path(simulation_directory, "Output")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
results_output_directory <- output_directory
sims_output_directory <- output_directory

results_lin <- readRDS(paste0(results_output_directory,"/v20_linear_regression_nobs_propmis_25_scen_Bayesimp_B5000.RDS"))
results_lin_lnorm <- readRDS(paste0(results_output_directory,"/v21_linear_regression_nobs_propmis_25_scen_Bayesimp_lognormalXskew2_B5000.RDS"))
results_log <- readRDS(paste0(results_output_directory,"/v20_logistic_regression_nobs_propmis_15_scen_Bayesimp_B5000.RDS"))
results_log_lnorm <- readRDS(paste0(results_output_directory,"/v21_logistic_regression_nobs_propmis_15_scen_Bayesimp_lognormalXskew2_B5000.RDS"))
results_cox <- readRDS(paste0(results_output_directory,"/v20_cox_regression_nobs_propmis_20_scen_Bayesimp_cens20_B5000.RDS"))
results_cox_lnorm <- readRDS(paste0(results_output_directory,"/v20_cox_regression_nobs_propmis_20_scen_Bayesimp_cens20_lognormalXskew2_B5000.RDS"))

other_info_types <- c("ess_bulk","ess_tail","fmi","df.com","df.barnard","df.norm","skew_MLE","kurt_MLE","skew_U","kurt_U","skew_par","kurt_par","divergences")

# lapply(list(results_lin,results_lin_lnorm,results_log,results_log_lnorm,results_cox,results_cox_lnorm),
#        function(x) lapply(x, function(y) summary(y$other_info[,"divergences"])))
# lapply(results_cox_lnorm, function(x) mean(x$other_info[,"divergences"]>5))

# summary(results_lin_lnorm[[25]]$other_info[,"fmi"])
# results_lin_lnorm[[25]]$ub_95_CI
# summary(results_lin_lnorm[[25]]$other_info[,"kurt_par"])
# vapply(results_lin_lnorm,function(x){x$other_info_summary[,c("skew_par")]},FUN.VALUE = numeric(length(results_lin_lnorm[[1]]$other_info_summary[,c("kurt_par")])))
# results_lin[[5]]$results

# Get vector of which scenarios present for each model from names of scenarios
scen_lin <- as.numeric(names(results_lin))
scen_lin_lnorm <- as.numeric(names(results_lin_lnorm))
scen_log <- as.numeric(names(results_log))
scen_log_lnorm <- as.numeric(names(results_log_lnorm))
scen_cox <- as.numeric(names(results_cox))
scen_cox_lnorm <- as.numeric(names(results_cox_lnorm))

add_group_labels <- function(group_centers = c(3, 8, 13, 18, 23),
                             group_labels = c(10, 20, 50, 100, 200),
                             vertical_offset = 0.06) {
  # Remove default axis and ticks completely
  axis(side = 1, at = NULL, labels = FALSE, tck = 0)
  
  # Compute position just below x-axis
  usr <- par("usr")            # user coordinate limits: (x1, x2, y1, y2)
  y_pos <- usr[3] - vertical_offset * (usr[4] - usr[3])
  
  # Add centered group labels
  text(x = group_centers,
       y = y_pos,
       labels = paste0("n_obs = ", group_labels),
       xpd = NA, font = 1, cex = 1.5)
}

# define group centers and labels
group_centers <- c(3, 8, 13, 18, 23)
group_labels  <- c(10, 20, 50, 100, 200)
scenarios_all <- 1:25

#####
# Coverage (all models)
coverages_lin <- sapply(results_lin,function(x){x$results[,"coverage"]})
rownames(coverages_lin) <- rownames(results_lin[[1]]$results)
colnames(coverages_lin) <- 1:ncol(coverages_lin)
coverages_lin_lnorm <- sapply(results_lin_lnorm,function(x){x$results[,"coverage"]})
rownames(coverages_lin_lnorm) <- rownames(results_lin_lnorm[[1]]$results)
colnames(coverages_lin_lnorm) <- 1:ncol(coverages_lin_lnorm)
coverages_log <- sapply(results_log,function(x){x$results[,"coverage"]})
rownames(coverages_log) <- rownames(results_log[[1]]$results)
colnames(coverages_log) <- 1:ncol(coverages_log)
coverages_log_lnorm <- sapply(results_log_lnorm,function(x){x$results[,"coverage"]})
rownames(coverages_log_lnorm) <- rownames(results_log_lnorm[[1]]$results)
colnames(coverages_log_lnorm) <- 1:ncol(coverages_log_lnorm)
coverages_cox <- sapply(results_cox,function(x){x$results[,"coverage"]})
rownames(coverages_cox) <- rownames(results_cox[[1]]$results)
colnames(coverages_cox) <- 1:ncol(coverages_cox)
coverages_cox_lnorm <- sapply(results_cox_lnorm,function(x){x$results[,"coverage"]})
rownames(coverages_cox_lnorm) <- rownames(results_cox_lnorm[[1]]$results)
colnames(coverages_cox_lnorm) <- 1:ncol(coverages_cox_lnorm)

# Get observed coverage
coverages_list <- list(coverages_lin,coverages_lin_lnorm,coverages_cox,coverages_log)
get_observed_coverages <- function(method, results_list = coverages_list){
  sapply(results_list, function(x) {sum(as.numeric(x[method,]>qbinom(p=0.025,size=5000,prob=0.95)/5000))})
} 
get_observed_coverages("Bayes_10000")
get_observed_coverages("Rubin_1000")
get_observed_coverages("ABpool_t_10000")
get_observed_coverages("Rubin_norm_1000")
get_observed_coverages("ABpool_norm_10000")
require(dplyr)
get_observed_coverages("Bayes_10000") %>% sum()
get_observed_coverages("Rubin_1000") %>% sum()
get_observed_coverages("ABpool_t_10000") %>% sum()
get_observed_coverages("Rubin_norm_1000") %>% sum()
get_observed_coverages("ABpool_norm_10000") %>% sum()

png(paste0(sims_output_directory,"/Plots/","v21_stan_combined_coverage.png"),width = 2000,height=1000, pointsize=16)
layout_mat <- matrix(c(1,2,3,
                       7,7,7,
                       4,5,6), nrow = 3, byrow = TRUE)
layout(layout_mat, heights = c(0.47, 0.06, 0.47))
par(mar=c(2.5,3.5,2.5,0.1), oma=c(0.5,0.1,0.1,0.5))
plot(scen_lin,coverages_lin["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.95,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.95)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.95)/5000,lty=2,lwd=2)
points(scen_lin, coverages_lin["CC_ML",], pch = pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin,coverages_lin["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin, coverages_lin["Rubin_norm_1000",], pch = pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin, coverages_lin["Rubin_1000",], pch = pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin, coverages_lin["ABpool_norm_10000",], pch = pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin, coverages_lin["ABpool_t_10000",], pch = pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression
plot(scen_log,coverages_log["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.95,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.95)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.95)/5000,lty=2,lwd=2)
points(scen_log,coverages_log["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log,coverages_log["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log,coverages_log["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log,coverages_log["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log,coverages_log["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log,coverages_log["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[3:5], group_labels = group_labels[3:5])
# cox regression
plot(scen_cox,coverages_cox["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.95,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.95)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.95)/5000,lty=2,lwd=2)
points(scen_cox,coverages_cox["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox,coverages_cox["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox,coverages_cox["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox,coverages_cox["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox,coverages_cox["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox,coverages_cox["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[2:5], group_labels = group_labels[2:5])
# linear regression - lognormal covariate
plot(scen_lin_lnorm,coverages_lin_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.95,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.95)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.95)/5000,lty=2,lwd=2)
points(scen_lin_lnorm,coverages_lin_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin_lnorm,coverages_lin_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin_lnorm,coverages_lin_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,coverages_lin_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin_lnorm,coverages_lin_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,coverages_lin_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression - lognormal covariate
plot(scen_log_lnorm,coverages_log_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.95,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.95)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.95)/5000,lty=2,lwd=2)
points(scen_log_lnorm,coverages_log_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log_lnorm,coverages_log_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log_lnorm,coverages_log_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log_lnorm,coverages_log_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log_lnorm,coverages_log_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log_lnorm,coverages_log_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Cox regression - lognormal covariate
plot(scen_cox_lnorm,coverages_cox_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.95,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.95)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.95)/5000,lty=2,lwd=2)
points(scen_cox_lnorm,coverages_cox_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox_lnorm,coverages_cox_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox_lnorm,coverages_cox_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,coverages_cox_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox_lnorm,coverages_cox_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,coverages_cox_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# legend
par(mar = c(0, 0, 0, 0)); plot.new()
legend(x = "center",
       legend = names(outcome_cols),
       col = outcome_cols,
       pch = pch_vals,
       pt.cex = 3,
       cex = 2,
       lty = NA,
       lwd = 3,
       ncol = length(outcome_cols),
       bty = "n",
       xpd = TRUE)
dev.off()

#####
# Power
powers_lin <- sapply(results_lin,function(x){x$results[,"power"]})
rownames(powers_lin) <- rownames(results_lin[[1]]$results)
colnames(powers_lin) <- 1:ncol(powers_lin)
powers_lin_lnorm <- sapply(results_lin_lnorm,function(x){x$results[,"power"]})
rownames(powers_lin_lnorm) <- rownames(results_lin_lnorm[[1]]$results)
colnames(powers_lin_lnorm) <- 1:ncol(powers_lin_lnorm)
powers_log <- sapply(results_log,function(x){x$results[,"power"]})
rownames(powers_log) <- rownames(results_log[[1]]$results)
colnames(powers_log) <- 1:ncol(powers_log)
powers_log_lnorm <- sapply(results_log_lnorm,function(x){x$results[,"power"]})
rownames(powers_log_lnorm) <- rownames(results_log_lnorm[[1]]$results)
colnames(powers_log_lnorm) <- 1:ncol(powers_log_lnorm)
powers_cox <- sapply(results_cox,function(x){x$results[,"power"]})
rownames(powers_cox) <- rownames(results_cox[[1]]$results)
colnames(powers_cox) <- 1:ncol(powers_cox)
powers_cox_lnorm <- sapply(results_cox_lnorm,function(x){x$results[,"power"]})
rownames(powers_cox_lnorm) <- rownames(results_cox_lnorm[[1]]$results)
colnames(powers_cox_lnorm) <- 1:ncol(powers_cox_lnorm)

png(paste0(sims_output_directory,"/Plots/","v21_stan_combined_power.png"),width = 2000,height=1000, pointsize=16)
layout_mat <- matrix(c(1,2,3,
                       7,7,7,
                       4,5,6), nrow = 3, byrow = TRUE)
layout(layout_mat, heights = c(0.47, 0.06, 0.47))
par(mar=c(2.5,3.5,2.5,0.1), oma=c(0.5,0.1,0.1,0.5))
# linear regression
plot(scen_lin,powers_lin["Bayes_10000",], ylim=c(0.6,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0.8,lwd=1.5)
points(scen_lin,powers_lin["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin,powers_lin["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin,powers_lin["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin,powers_lin["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin,powers_lin["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin,powers_lin["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression
plot(scen_log,powers_log["Bayes_10000",], ylim=c(0.6,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0.8,lwd=1.5)
points(scen_log,powers_log["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log,powers_log["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log,powers_log["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log,powers_log["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log,powers_log["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log,powers_log["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[3:5], group_labels = group_labels[3:5])
# cox regression
plot(scen_cox,powers_cox["Bayes_10000",], ylim=c(0.6,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0.8,lwd=1.5)
points(scen_cox,powers_cox["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox,powers_cox["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox,powers_cox["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox,powers_cox["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox,powers_cox["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox,powers_cox["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[2:5], group_labels = group_labels[2:5])
# linear regression - lognormal covariate
plot(scen_lin_lnorm,powers_lin_lnorm["Bayes_10000",], ylim=c(0.6,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression - lognormal covariate")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0.8,lwd=1.5)
points(scen_lin_lnorm,powers_lin_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin_lnorm,powers_lin_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin_lnorm,powers_lin_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,powers_lin_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin_lnorm,powers_lin_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,powers_lin_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Logistic regression - lognormal covariate
plot(scen_log_lnorm,powers_log_lnorm["Bayes_10000",], ylim=c(0.6,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression - lognormal covariate")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0.8,lwd=1.5)
points(scen_log_lnorm,powers_log_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log_lnorm,powers_log_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log_lnorm,powers_log_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log_lnorm,powers_log_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log_lnorm,powers_log_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log_lnorm,powers_log_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Cox regression - lognormal covariate
plot(scen_cox_lnorm,powers_cox_lnorm["Bayes_10000",], ylim=c(0.6,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression - lognormal covariate")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0.8,lwd=1.5)
points(scen_cox_lnorm,powers_cox_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox_lnorm,powers_cox_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox_lnorm,powers_cox_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,powers_cox_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox_lnorm,powers_cox_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,powers_cox_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# legend
par(mar = c(0, 0, 0, 0)); plot.new()
legend(x = "center",
       legend = names(outcome_cols),
       col = outcome_cols,
       pch = pch_vals,
       pt.cex = 3,
       cex = 2,
       lty = NA,
       lwd = 3,
       ncol = length(outcome_cols),
       bty = "n",
       xpd = TRUE)
dev.off()

#####
# Upper coverage
upper_coverages_lin <- 1 - sapply(results_lin,function(x){x$results[,"miss_low"]})
rownames(upper_coverages_lin) <- rownames(results_lin[[1]]$results)
colnames(upper_coverages_lin) <- 1:ncol(upper_coverages_lin)
upper_coverages_lin_lnorm <- 1 - sapply(results_lin_lnorm,function(x){x$results[,"miss_low"]})
rownames(upper_coverages_lin_lnorm) <- rownames(results_lin_lnorm[[1]]$results)
colnames(upper_coverages_lin_lnorm) <- 1:ncol(upper_coverages_lin_lnorm)
upper_coverages_log <- 1 - sapply(results_log,function(x){x$results[,"miss_low"]})
rownames(upper_coverages_log) <- rownames(results_log[[1]]$results)
colnames(upper_coverages_log) <- 1:ncol(upper_coverages_log)
upper_coverages_log_lnorm <- 1 - sapply(results_log_lnorm,function(x){x$results[,"miss_low"]})
rownames(upper_coverages_log_lnorm) <- rownames(results_log_lnorm[[1]]$results)
colnames(upper_coverages_log_lnorm) <- 1:ncol(upper_coverages_log_lnorm)
upper_coverages_cox <- 1 - sapply(results_cox,function(x){x$results[,"miss_low"]})
rownames(upper_coverages_cox) <- rownames(results_cox[[1]]$results)
colnames(upper_coverages_cox) <- 1:ncol(upper_coverages_cox)
upper_coverages_cox_lnorm <- 1 - sapply(results_cox_lnorm,function(x){x$results[,"miss_low"]})
rownames(upper_coverages_cox_lnorm) <- rownames(results_cox_lnorm[[1]]$results)
colnames(upper_coverages_cox_lnorm) <- 1:ncol(upper_coverages_cox_lnorm)

png(paste0(sims_output_directory,"/Plots/","v21_stan_combined_coverage_upper.png"),width = 2000,height=1000, pointsize=16)
layout_mat <- matrix(c(1,2,3,
                       7,7,7,
                       4,5,6), nrow = 3, byrow = TRUE)
layout(layout_mat, heights = c(0.47, 0.06, 0.47))
par(mar=c(2.5,3.5,2.5,0.1), oma=c(0.5,0.1,0.1,0.5))
plot(scen_lin,upper_coverages_lin["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_lin, upper_coverages_lin["CC_ML",], pch = pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin,upper_coverages_lin["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin, upper_coverages_lin["Rubin_norm_1000",], pch = pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin, upper_coverages_lin["Rubin_1000",], pch = pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin, upper_coverages_lin["ABpool_norm_10000",], pch = pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin, upper_coverages_lin["ABpool_t_10000",], pch = pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression
plot(scen_log,upper_coverages_log["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_log,upper_coverages_log["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log,upper_coverages_log["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log,upper_coverages_log["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log,upper_coverages_log["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log,upper_coverages_log["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log,upper_coverages_log["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[3:5], group_labels = group_labels[3:5])
# cox regression
plot(scen_cox,upper_coverages_cox["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_cox,upper_coverages_cox["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox,upper_coverages_cox["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox,upper_coverages_cox["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox,upper_coverages_cox["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox,upper_coverages_cox["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox,upper_coverages_cox["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[2:5], group_labels = group_labels[2:5])
# linear regression - lognormal covariate
plot(scen_lin_lnorm,upper_coverages_lin_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_lin_lnorm,upper_coverages_lin_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin_lnorm,upper_coverages_lin_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin_lnorm,upper_coverages_lin_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,upper_coverages_lin_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin_lnorm,upper_coverages_lin_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,upper_coverages_lin_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression - lognormal covariate
plot(scen_log_lnorm,upper_coverages_log_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_log_lnorm,upper_coverages_log_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log_lnorm,upper_coverages_log_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log_lnorm,upper_coverages_log_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log_lnorm,upper_coverages_log_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log_lnorm,upper_coverages_log_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log_lnorm,upper_coverages_log_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Cox regression - lognormal covariate
plot(scen_cox_lnorm,upper_coverages_cox_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_cox_lnorm,upper_coverages_cox_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox_lnorm,upper_coverages_cox_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox_lnorm,upper_coverages_cox_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,upper_coverages_cox_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox_lnorm,upper_coverages_cox_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,upper_coverages_cox_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# legend
par(mar = c(0, 0, 0, 0)); plot.new()
legend(x = "center",
       legend = names(outcome_cols),
       col = outcome_cols,
       pch = pch_vals,
       pt.cex = 3,
       cex = 2,
       lty = NA,
       lwd = 3,
       ncol = length(outcome_cols),
       bty = "n",
       xpd = TRUE)
dev.off()

#####
# Lower coverage
lower_coverages_lin <- 1 - lower_coverages_lin + coverages_lin
lower_coverages_lin_lnorm <- 1 - lower_coverages_lin_lnorm + coverages_lin_lnorm
lower_coverages_log <- 1 - lower_coverages_log + coverages_log
lower_coverages_log_lnorm <- 1 - lower_coverages_log_lnorm + coverages_log_lnorm
lower_coverages_cox <- 1 - lower_coverages_cox + coverages_cox
lower_coverages_cox_lnorm <- 1 - lower_coverages_cox_lnorm + coverages_cox_lnorm

png(paste0(sims_output_directory,"/Plots/","v21_stan_combined_coverage_lower.png"),width = 2000,height=1000, pointsize=16)
layout_mat <- matrix(c(1,2,3,
                       7,7,7,
                       4,5,6), nrow = 3, byrow = TRUE)
layout(layout_mat, heights = c(0.47, 0.06, 0.47))
par(mar=c(2.5,3.5,2.5,0.1), oma=c(0.5,0.1,0.1,0.5))
plot(scen_lin,lower_coverages_lin["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_lin, lower_coverages_lin["CC_ML",], pch = pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin,lower_coverages_lin["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin, lower_coverages_lin["Rubin_norm_1000",], pch = pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin, lower_coverages_lin["Rubin_1000",], pch = pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin, lower_coverages_lin["ABpool_norm_10000",], pch = pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin, lower_coverages_lin["ABpool_t_10000",], pch = pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression
plot(scen_log,lower_coverages_log["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_log,lower_coverages_log["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log,lower_coverages_log["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log,lower_coverages_log["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log,lower_coverages_log["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log,lower_coverages_log["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log,lower_coverages_log["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[3:5], group_labels = group_labels[3:5])
# cox regression
plot(scen_cox,lower_coverages_cox["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_cox,lower_coverages_cox["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox,lower_coverages_cox["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox,lower_coverages_cox["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox,lower_coverages_cox["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox,lower_coverages_cox["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox,lower_coverages_cox["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[2:5], group_labels = group_labels[2:5])
# linear regression - lognormal covariate
plot(scen_lin_lnorm,lower_coverages_lin_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_lin_lnorm,lower_coverages_lin_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_lin_lnorm,lower_coverages_lin_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_lin_lnorm,lower_coverages_lin_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,lower_coverages_lin_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin_lnorm,lower_coverages_lin_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,lower_coverages_lin_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression - lognormal covariate
plot(scen_log_lnorm,lower_coverages_log_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_log_lnorm,lower_coverages_log_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_log_lnorm,lower_coverages_log_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_log_lnorm,lower_coverages_log_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log_lnorm,lower_coverages_log_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log_lnorm,lower_coverages_log_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log_lnorm,lower_coverages_log_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Cox regression - lognormal covariate
plot(scen_cox_lnorm,lower_coverages_cox_lnorm["Bayes_10000",], ylim=c(0.9,1), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression - lognormal covariate")
abline(h = seq(0.9,1,0.01),col="grey95")
abline(h = seq(0.9,1,0.02),col="grey90")
abline(h=0.975,lwd=2)
abline(h=qbinom(p=0.025,size=5000,prob=0.975)/5000,lty=2,lwd=2)
abline(h=qbinom(p=0.975,size=5000,prob=0.975)/5000,lty=2,lwd=2)
points(scen_cox_lnorm,lower_coverages_cox_lnorm["CC_ML",], ylim=c(0,1), pch=pch_vals["CC"], col = outcome_cols["CC"],cex=3, lwd=3)
points(scen_cox_lnorm,lower_coverages_cox_lnorm["Bayes_10000",], ylim=c(0,1), pch=pch_vals["Bayes"], col = outcome_cols["Bayes"],cex=3, lwd=3)
points(scen_cox_lnorm,lower_coverages_cox_lnorm["Rubin_norm_1000",], ylim=c(0,1), pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,lower_coverages_cox_lnorm["Rubin_1000",], ylim=c(0,1), pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox_lnorm,lower_coverages_cox_lnorm["ABpool_norm_10000",], ylim=c(0,1), pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,lower_coverages_cox_lnorm["ABpool_t_10000",], ylim=c(0,1), pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# legend
par(mar = c(0, 0, 0, 0)); plot.new()
legend(x = "center",
       legend = names(outcome_cols),
       col = outcome_cols,
       pch = pch_vals,
       pt.cex = 3,
       cex = 2,
       lty = NA,
       lwd = 3,
       ncol = length(outcome_cols),
       bty = "n",
       xpd = TRUE)
dev.off()

# Midpoint discrepancy from Bayes
# CI midpoint
CI_midpoint_diff_Bayes_fun <- function(x){
  CI_midpoint <- ((x$ub_95_CI+x$lb_95_CI)/2)[,"mean"]
  CI_width_Bayes <- (x$ub_95_CI-x$lb_95_CI)["Bayes_10000","mean"]
  CI_midpoint_Bayes <- CI_midpoint["Bayes_10000"]
  # mean CI midpoint discrepancy from Bayes
  CI_midpoint_diff_Bayes <- (CI_midpoint-CI_midpoint_Bayes)/CI_width_Bayes
  CI_midpoint_diff_Bayes
}
CI_midpoint_diff_Bayes_fun(results_lin_lnorm[[1]])
midpoint_diff_lin <- vapply(results_lin,CI_midpoint_diff_Bayes_fun,numeric(nrow(results_lin[[1]]$results)))
midpoint_diff_lin_lnorm <- vapply(results_lin_lnorm,CI_midpoint_diff_Bayes_fun,numeric(nrow(results_lin_lnorm[[1]]$results)))
midpoint_diff_log <- vapply(results_log,CI_midpoint_diff_Bayes_fun,numeric(nrow(results_log[[1]]$results)))
midpoint_diff_log_lnorm <- vapply(results_log_lnorm,CI_midpoint_diff_Bayes_fun,numeric(nrow(results_log_lnorm[[1]]$results)))
midpoint_diff_cox <- vapply(results_cox,CI_midpoint_diff_Bayes_fun,numeric(nrow(results_cox[[1]]$results)))
midpoint_diff_cox_lnorm <- vapply(results_cox_lnorm,CI_midpoint_diff_Bayes_fun,numeric(nrow(results_cox_lnorm[[1]]$results)))
# Midpoint outcomes
mp_outcomes <- c("ABpool_t", "ABpool_norm", "Rubin_t", "Rubin_norm")
mp_outcome_cols <- outcome_cols[mp_outcomes]
mp_pch_vals <- pch_vals[mp_outcomes]

png(paste0(sims_output_directory,"/Plots/","v21_stan_combined_midpoint_diff.png"),width = 2000,height=1000, pointsize=16)
layout_mat <- matrix(c(1,2,3,
                       7,7,7,
                       4,5,6), nrow = 3, byrow = TRUE)
layout(layout_mat, heights = c(0.47, 0.06, 0.47))
par(mar=c(2.5,3.5,2.5,0.1), oma=c(0.5,0.1,0.1,0.5))
# linear regression
plot(scen_lin,midpoint_diff_lin["Bayes_10000",], ylim=c(-0.2,0.2), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0,lwd=2)
points(scen_lin,midpoint_diff_lin["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin,midpoint_diff_lin["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin,midpoint_diff_lin["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin,midpoint_diff_lin["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression
plot(scen_log,midpoint_diff_log["Bayes_10000",], ylim=c(-0.2,0.2), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0,lwd=2)
points(scen_log,midpoint_diff_log["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log,midpoint_diff_log["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log,midpoint_diff_log["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log,midpoint_diff_log["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[3:5], group_labels = group_labels[3:5])
# cox regression
plot(scen_cox,midpoint_diff_cox["Bayes_10000",], ylim=c(-0.2,0.2), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0,lwd=2)
points(scen_cox,midpoint_diff_cox["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox,midpoint_diff_cox["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox,midpoint_diff_cox["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox,midpoint_diff_cox["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[2:5], group_labels = group_labels[2:5])
# linear regression - lognormal covariate
plot(scen_lin_lnorm,midpoint_diff_lin_lnorm["Bayes_10000",], ylim=c(-0.2,0.2), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression - lognormal covariate")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0,lwd=2)
points(scen_lin_lnorm,midpoint_diff_lin_lnorm["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,midpoint_diff_lin_lnorm["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin_lnorm,midpoint_diff_lin_lnorm["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,midpoint_diff_lin_lnorm["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Logistic regression - lognormal covariate
plot(scen_log_lnorm,midpoint_diff_log_lnorm["Bayes_10000",], ylim=c(-0.2,0.2), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression - lognormal covariate")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0,lwd=2)
points(scen_log_lnorm,midpoint_diff_log_lnorm["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log_lnorm,midpoint_diff_log_lnorm["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log_lnorm,midpoint_diff_log_lnorm["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log_lnorm,midpoint_diff_log_lnorm["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Cox regression - lognormal covariate
plot(scen_cox_lnorm,midpoint_diff_cox_lnorm["Bayes_10000",], ylim=c(-0.2,0.2), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression - lognormal covariate")
abline(h = seq(0.6,1,0.01),col="grey97")
abline(h = seq(0.6,1,0.05),col="grey90")
abline(h = seq(0.6,1,0.1),col="grey70")
abline(h=0,lwd=2)
points(scen_cox_lnorm,midpoint_diff_cox_lnorm["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,midpoint_diff_cox_lnorm["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox_lnorm,midpoint_diff_cox_lnorm["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,midpoint_diff_cox_lnorm["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# legend
par(mar = c(0, 0, 0, 0)); plot.new()
legend(x = "center",
       legend = names(mp_outcome_cols),
       col = mp_outcome_cols,
       pch = mp_pch_vals,
       pt.cex = 3,
       cex = 2,
       lty = NA,
       lwd = 3,
       ncol = length(mp_outcome_cols),
       bty = "n",
       xpd = TRUE)
dev.off()

# Get observed midpoints
midpoint_diff_list <- cbind(midpoint_diff_lin,midpoint_diff_lin_lnorm,midpoint_diff_cox,midpoint_diff_log)
rowMeans(abs(midpoint_diff_list))[c("Rubin_1000","ABpool_t_10000","Rubin_norm_1000", "ABpool_norm_10000")]

# Relative CI width to Bayes
# CI width
CI_relative_width_Bayes_fun <- function(x){
  CI_width <- (x$ub_95_CI-x$lb_95_CI)[,"mean"]
  CI_width_Bayes <- CI_width["Bayes_10000"]
  CI_relative_width_Bayes <- CI_width/CI_width_Bayes
  CI_relative_width_Bayes
}
relative_width_lin <- vapply(results_lin,CI_relative_width_Bayes_fun,numeric(nrow(results_lin[[1]]$results)))
relative_width_lin_lnorm <- vapply(results_lin_lnorm,CI_relative_width_Bayes_fun,numeric(nrow(results_lin_lnorm[[1]]$results)))
relative_width_log <- vapply(results_log,CI_relative_width_Bayes_fun,numeric(nrow(results_log[[1]]$results)))
relative_width_log_lnorm <- vapply(results_log_lnorm,CI_relative_width_Bayes_fun,numeric(nrow(results_log_lnorm[[1]]$results)))
relative_width_cox <- vapply(results_cox,CI_relative_width_Bayes_fun,numeric(nrow(results_cox[[1]]$results)))
relative_width_cox_lnorm <- vapply(results_cox_lnorm,CI_relative_width_Bayes_fun,numeric(nrow(results_cox_lnorm[[1]]$results)))

png(paste0(sims_output_directory,"/Plots/","v21_stan_combined_relative_width.png"),width = 2000,height=1000, pointsize=16)
layout_mat <- matrix(c(1,2,3,
                       7,7,7,
                       4,5,6), nrow = 3, byrow = TRUE)
layout(layout_mat, heights = c(0.47, 0.06, 0.47))
par(mar=c(2.5,3.5,2.5,0.1), oma=c(0.5,0.1,0.1,0.5))
# linear regression
plot(scen_lin,relative_width_lin["Bayes_10000",], ylim=c(0.5,1.5), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression")
abline(h = seq(0.5,1.5,0.05),col="grey90")
abline(h = seq(0.5,1.5,0.1),col="grey70")
abline(h=1,lwd=2)
points(scen_lin,relative_width_lin["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin,relative_width_lin["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin,relative_width_lin["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin,relative_width_lin["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# logistic regression
plot(scen_log,relative_width_log["Bayes_10000",], ylim=c(0.5,1.5), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression")
abline(h = seq(0.5,1.5,0.05),col="grey90")
abline(h = seq(0.5,1.5,0.1),col="grey70")
abline(h=1,lwd=2)
points(scen_log,relative_width_log["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log,relative_width_log["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log,relative_width_log["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log,relative_width_log["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[3:5], group_labels = group_labels[3:5])
# cox regression
plot(scen_cox,relative_width_cox["Bayes_10000",], ylim=c(0.5,1.5), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression")
abline(h = seq(0.5,1.5,0.05),col="grey90")
abline(h = seq(0.5,1.5,0.1),col="grey70")
abline(h=1,lwd=2)
points(scen_cox,relative_width_cox["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox,relative_width_cox["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox,relative_width_cox["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox,relative_width_cox["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers[2:5], group_labels = group_labels[2:5])
# linear regression - lognormal covariate
plot(scen_lin_lnorm,relative_width_lin_lnorm["Bayes_10000",], ylim=c(0.5,1.5), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Linear regression - lognormal covariate")
abline(h = seq(0.5,1.5,0.05),col="grey90")
abline(h = seq(0.5,1.5,0.1),col="grey70")
abline(h=1,lwd=2)
points(scen_lin_lnorm,relative_width_lin_lnorm["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,relative_width_lin_lnorm["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_lin_lnorm,relative_width_lin_lnorm["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_lin_lnorm,relative_width_lin_lnorm["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Logistic regression - lognormal covariate
plot(scen_log_lnorm,relative_width_log_lnorm["Bayes_10000",], ylim=c(0.5,1.5), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Logistic regression - lognormal covariate")
abline(h = seq(0.5,1.5,0.05),col="grey90")
abline(h = seq(0.5,1.5,0.1),col="grey70")
abline(h=1,lwd=2)
points(scen_log_lnorm,relative_width_log_lnorm["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_log_lnorm,relative_width_log_lnorm["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_log_lnorm,relative_width_log_lnorm["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_log_lnorm,relative_width_log_lnorm["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# Cox regression - lognormal covariate
plot(scen_cox_lnorm,relative_width_cox_lnorm["Bayes_10000",], ylim=c(0.5,1.5), ylab = NA, xaxt="n", xlim=c(1,25), col = NA, cex.main = 1.5, cex.axis = 1.5, main = "Cox regression - lognormal covariate")
abline(h = seq(0.5,1.5,0.05),col="grey90")
abline(h = seq(0.5,1.5,0.1),col="grey70")
abline(h=1,lwd=2)
points(scen_cox_lnorm,relative_width_cox_lnorm["Rubin_norm_1000",], pch=pch_vals["Rubin_norm"], col = outcome_cols["Rubin_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,relative_width_cox_lnorm["Rubin_1000",], pch=pch_vals["Rubin_t"], col = outcome_cols["Rubin_t"],cex=3, lwd=3)
points(scen_cox_lnorm,relative_width_cox_lnorm["ABpool_norm_10000",], pch=pch_vals["ABpool_norm"], col = outcome_cols["ABpool_norm"],cex=3, lwd=3)
points(scen_cox_lnorm,relative_width_cox_lnorm["ABpool_t_10000",], pch=pch_vals["ABpool_t"], col = outcome_cols["ABpool_t"],cex=3, lwd=3)
abline(v=c(5,10,15,20)+0.5,lty=1)
add_group_labels(group_centers = group_centers, group_labels = group_labels)
# legend
par(mar = c(0, 0, 0, 0)); plot.new()
legend(x = "center",
       legend = names(mp_outcome_cols),
       col = mp_outcome_cols,
       pch = mp_pch_vals,
       pt.cex = 3,
       cex = 2,
       lty = NA,
       lwd = 3,
       ncol = length(mp_outcome_cols),
       bty = "n",
       xpd = TRUE)
dev.off()


# Results for different m
results_all <- c(results_lin,results_log,results_cox,results_lin_lnorm,results_log_lnorm,results_cox_lnorm)
average_across_all_sims <- Reduce("+",lapply(results_all,FUN=function(x){x$results}))/length(results_all)
average_across_all_sims <- average_across_all_sims[c(5,2,3,4,6:nrow(average_across_all_sims)),c("power","coverage")]
average_across_all_sims <- round(average_across_all_sims,3)
average_across_all_sims$method <- c("CC",rep("Bayes",3),rep("Rubin_t",4),"Rubin_norm",rep("ABpool_norm",3),rep("ABpool_t",3))
sampling_m <- c(100,1000,10000)
average_across_all_sims$m <- c(NA,sampling_m,5,20,100,1000,1000,sampling_m,sampling_m)
average_across_all_sims <- average_across_all_sims[,c("method","m","power","coverage")]
write.csv(average_across_all_sims,paste0(sims_output_directory,"/Simulations_average_power_coverage_v21.csv"),row.names=F)

# Main text table
library(dplyr)
library(knitr)
library(kableExtra)

powers_and_coverages <- round(rbind(cbind(t(coverages_lin[c("Rubin_1000","ABpool_t_10000"),c(1,3,5)]),
                                          t(powers_lin[c("Rubin_1000","ABpool_t_10000"),c(1,3,5)])),
                                    cbind(t(coverages_log[c("Rubin_norm_1000","ABpool_norm_10000"),c(1,3,5)]),
                                          t(powers_log[c("Rubin_norm_1000","ABpool_norm_10000"),c(1,3,5)]))),2)
table.main <- cbind("model"=c(rep("Linear regression (n_obs = 10)",3),rep("Logistic regression (n_obs = 50)",3)),"n"=c(11,20,100)*c(rep(1,3),rep(5,3)),"Missing X"=rep(c("9%","50%","90%"),2),powers_and_coverages)
colnames(table.main)[4:7] <- c("Rubin_c","ABpool_c","Rubin_p","ABpool_p")
tab <- as.data.frame(table.main, stringsAsFactors = FALSE)
tab %>%
  select(-model) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c(
      "$n$",
      "Missing $X$ (\\%)",
      "Rubin",
      "ABpool",
      "Rubin",
      "ABpool"
    ),
    row.names = FALSE,
    align = c("r", "r", "r", "r", "r", "r")
  ) %>%
  add_header_above(
    c(" " = 2, "Coverage" = 2, "Power" = 2)
  ) %>%
  pack_rows(
    "Linear regression ($n_{\\mathrm{obs}}=10$)",
    1, 3,
    escape = FALSE
  ) %>%
  pack_rows(
    "Logistic regression ($n_{\\mathrm{obs}}=50$)",
    4, 6,
    escape = FALSE
  )

# Important rows
imp_rows <- c("CC_ML","Bayes_10000","Rubin_1000","Rubin_norm_1000","ABpool_norm_10000","ABpool_t_10000")
names(results_lin[[1]])
# Explaining linear regression
results_lin[[1]]$results[imp_rows,]
results_lin[[1]]$other_info_summary
results_lin[[5]]$results[imp_rows,]
results_lin[[5]]$other_info_summary
X_eff <- results_lin[[5]]$parameters$X_eff
sim_results_scen <- results_lin[[5]]$sim_results_array
dimnames(sim_results_scen)
rubin_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==1 &
                            sim_results_scen[,"ABpool_t_10000","power"]==0)
rubin_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==1 &
                          sim_results_scen[,"ABpool_t_10000","coverage"]==0)
ABpool_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==0 &
                           sim_results_scen[,"ABpool_t_10000","coverage"]==1)
ABpool_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==0 &
                           sim_results_scen[,"ABpool_t_10000","power"]==1)
apply(sim_results_scen[results_lin[[5]]$other_info[,"bias_Bayes"]>0,,"2.5 %"],2,quantile,c(0.1,0.5,0.9))
apply(results_lin[[5]]$other_info[rubin_wins_coverage,],2,quantile,c(0.1,0.5,0.9),na.rm=T)
apply(results_lin[[5]]$other_info[rubin_wins_power,],2,quantile,c(0.1,0.5,0.9),na.rm=T)

results_lin[[10]]$results[imp_rows,]
results_lin[[10]]$other_info_summary
X_eff <- results_lin[[10]]$parameters$X_eff
sim_results_scen <- results_lin[[10]]$sim_results_array
dimnames(sim_results_scen)
rubin_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==1 &
                       sim_results_scen[,"ABpool_t_10000","power"]==0)
rubin_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==1 &
                          sim_results_scen[,"ABpool_t_10000","coverage"]==0)
ABpool_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==0 &
                           sim_results_scen[,"ABpool_t_10000","coverage"]==1)
ABpool_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==0 &
                        sim_results_scen[,"ABpool_t_10000","power"]==1)
apply(results_lin[[10]]$other_info[rubin_wins_coverage,],2,quantile,c(0.1,0.5,0.9))
apply(results_lin[[10]]$other_info[rubin_wins_power,],2,quantile,c(0.1,0.5,0.9))
sim_results_scen[rubin_wins_power,,"2.5 %"]
sim_results_scen[rubin_wins_power,,"97.5 %"]

# Explaining linear regression (lognormal covariate)
results_lin_lnorm[[1]]$results[imp_rows,]
results_lin_lnorm[[1]]$other_info_summary
results_lin_lnorm[[5]]$results[imp_rows,]
results_lin_lnorm[[5]]$other_info_summary
results_lin_lnorm[[10]]$results[imp_rows,]
results_lin_lnorm[[20]]$other_info_summary
results_lin_lnorm[[20]]$results[imp_rows,]
X_eff <- results_lin_lnorm[[20]]$parameters$X_eff
sim_results_scen <- results_lin_lnorm[[20]]$sim_results_array
dimnames(sim_results_scen)
rubin_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==1 &
                       sim_results_scen[,"ABpool_t_10000","power"]==0)
rubin_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==1 &
                          sim_results_scen[,"ABpool_t_10000","coverage"]==0)
ABpool_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==0 &
                           sim_results_scen[,"ABpool_t_10000","coverage"]==1)
ABpool_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==0 &
                        sim_results_scen[,"ABpool_t_10000","power"]==1)
sum(rubin_wins_coverage&rubin_wins_power)
apply(results_lin_lnorm[[20]]$other_info[rubin_wins_coverage,],2,quantile,c(0.1,0.5,0.9))
apply(results_lin_lnorm[[20]]$other_info[rubin_wins_power,],2,quantile,c(0.1,0.5,0.9))

results_lin_lnorm[[25]]$other_info_summary
results_lin_lnorm[[25]]$results[imp_rows,]

# Explaining logistic regression
results_log[[1]]$results[imp_rows,]
results_log[[1]]$other_info_summary
results_log[[5]]$results[imp_rows,]
results_log[[5]]$other_info_summary
X_eff <- results_log[[1]]$parameters$X_eff
sim_results_scen <- results_log[[1]]$sim_results_array
dimnames(sim_results_scen)
rubin_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==1 &
                       sim_results_scen[,"ABpool_t_10000","power"]==0)
rubin_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==1 &
                          sim_results_scen[,"ABpool_t_10000","coverage"]==0)
ABpool_wins_coverage <- (sim_results_scen[,"Rubin_1000","coverage"]==0 &
                           sim_results_scen[,"ABpool_t_10000","coverage"]==1)
ABpool_wins_power <- (sim_results_scen[,"Rubin_1000","power"]==0 &
                        sim_results_scen[,"ABpool_t_10000","power"]==1)
sum(ABpool_wins_power&ABpool_wins_coverage)
apply(results_log[[1]]$other_info[ABpool_wins_power,],2,quantile,c(0.1,0.5,0.9))
results_log[[6]]$results[imp_rows,]
results_log[[6]]$other_info_summary

# Explaining cox regression
results_cox[[1]]$other_info_summary
results_cox[[3]]$other_info_summary

results_cox[[5]]$other_info
results_cox[[5]]$other_info_summary
results_cox[[10]]$other_info_summary

X_eff <- results_cox[[5]]$parameters$X_eff
sim_results_scen <- results_cox[[5]]$sim_results_array
dimnames(sim_results_scen)
rubin_wins_power <- (sim_results_scen[,"Rubin_norm_1000","power"]==1 &
                       sim_results_scen[,"ABpool_norm_10000","power"]==0)
rubin_wins_coverage <- (sim_results_scen[,"Rubin_norm_1000","coverage"]==1 &
                          sim_results_scen[,"ABpool_norm_10000","coverage"]==0)
ABpool_wins_coverage <- (sim_results_scen[,"Rubin_norm_1000","coverage"]==0 &
                           sim_results_scen[,"ABpool_norm_10000","coverage"]==1)
ABpool_wins_power <- (sim_results_scen[,"Rubin_norm_1000","power"]==0 &
                        sim_results_scen[,"ABpool_norm_10000","power"]==1)
sum(rubin_wins_coverage&rubin_wins_power)
apply(results_cox[[5]]$other_info[rubin_wins_power,],2,quantile,c(0.1,0.5,0.9))