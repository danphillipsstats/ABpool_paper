output_directory <- paste0(getwd(),"/Output")
plot_directory <- paste0(output_directory,"/Plots")
##############
# Cornish-Fisher
# Gaussian case
# T = 1

scales::show_col(palette("Okabe-Ito"))
outcome_cols <- palette("Okabe-Ito")[c(6,7)]

Cornish_Fisher_norm_approx <- function(skew,kurt,alpha=0.025){
  z_a <- qnorm(alpha)
  CI_diff <- skew*(z_a^2-1)/6 + kurt*(z_a^3-3*z_a)/24 - skew^2*(2*z_a^3-5*z_a)/36
  return(CI_diff/abs(qnorm(alpha)))
}
sskew <- seq(-2,2,0.001)
Cornish_Fisher_norm_approx(2,0,0.025)
Cornish_Fisher_norm_approx(2,0,0.975)
CF_sskew_lb <- Cornish_Fisher_norm_approx(sskew,0,0.025)
CF_sskew_ub <- Cornish_Fisher_norm_approx(sskew,0,0.975)
kkurt <- seq(-2,2,0.001)
CF_kkurt_lb <- Cornish_Fisher_norm_approx(0,kkurt,0.025)
CF_kkurt_ub <- Cornish_Fisher_norm_approx(0,kkurt,0.975)
range(c(CF_sskew_lb,CF_sskew_ub,CF_kkurt_lb,CF_kkurt_ub))

png(paste0(plot_directory,"/Cornish_Fisher_error.png"),width = 1400,height=1000, pointsize=23)
plot(CF_sskew_ub~sskew,ylim=c(-1,1),col="NA",xlab="Skewness or kurtosis",ylab="Error")
lines(CF_sskew_ub~sskew,col=outcome_cols[1],lwd=2)
lines(CF_sskew_lb~sskew,col=outcome_cols[1],lty=2,lwd=2)
lines(CF_kkurt_ub~sskew,col=outcome_cols[2],lwd=2)
lines(CF_kkurt_lb~sskew,col=outcome_cols[2],lty=2,lwd=2)
abline(h=0,lty=1,lwd=2)
dev.off()