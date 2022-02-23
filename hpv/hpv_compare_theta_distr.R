# Load input files
## hpv_eta_smi_0.050.rda
load(file = file.choose())
hpv_eta_smi_0.05.df  = hpv_smi_mcmc_all[hpv_smi_mcmc_all$eta == 0.05,]
## load hpv_abc-smi_delta_128.RData
load(file = file.choose())
hpv_abc_smi_128.df = hpv_abc_smi_output.df
## load hpv_abc-smi_delta_scaled_1.RData
load(file = file.choose())
hpv_abc_smi_ds_1.df = hpv_abc_smi_output.df
## load hpv_bayes_param object from hpv_bayes.RData
load(file = file.choose())
## load hpv_cut_param object from hpv_cut.RData
load(file = file.choose())

# Load MASS library
library(MASS)

# Create contours for the bivariate theta posterior distribution from delta-SMI.
hpv_abc_smi_delta_128_contour_z = 
  kde2d(x = hpv_abc_smi_128.df$theta1, 
        y = hpv_abc_smi_128.df$theta2, n = 50)

# Create contours for the bivariate theta posterior distribution 
# from delta-SMI with delta scaled by sqrt(y).
hpv_abc_smi_ds_1_contour_z = 
  kde2d(x = hpv_abc_smi_ds_1.df$theta1, 
        y = hpv_abc_smi_ds_1.df$theta2, n = 50)

# Create contours for the bivariate theta posterior distribution from eta-SMI 
hpv_eta_smi_0.05_contour_z = 
  kde2d(x = hpv_eta_smi_0.05.df$theta1, 
        y = hpv_eta_smi_0.05.df$theta2, n = 50)

# Save the plot to working directory
cairo_pdf(file = "hpv_compare_theta_distr_unscaled.pdf", height = 4.5, width = 5)
windows(family='serif')
par(mar = c(5, 4, 1, 1) + 0.2, lend = 1)

# Plot the bivariate posterior theta distribution of full Bayes analysis.
plot(hpv_bayes_param$theta1[c(1:7200)*5], 
     hpv_bayes_param$theta2[c(1:7200)*5], 
     cex = 0.1, las = 1,
     xlim = c(-3, -1.2), ylim = c(7,37),
     col="#243757", xlab = expression(italic("\u03B8")[1]), 
     ylab = expression(italic("\u03B8")[2]))

# Plot the bivariate posterior theta distribution of cut analysis.
points(hpv_cut_param[c(1:7200)*5, 'theta_1'], 
       hpv_cut_param[c(1:7200)*5, 'theta_2'], 
       cex = 0.1, col = "#FAD77D")

# Contour of the bivariate posterior theta distribution of 
# delta-SMI analysis with delta (= 128) on the Poisson module.
contour(hpv_abc_smi_delta_128_contour_z, 
        drawlabels=FALSE, nlevels=7, col= "#985775", add=TRUE)

# Contour of the bivariate posterior theta distribution of 
# eta-SMI analysis with eta (= 0.05) on the Poisson module.
contour(hpv_eta_smi_0.05_contour_z, lty = "dashed",
        drawlabels=FALSE, nlevels=7, col= "#985775", add=TRUE)

# Plot legend
legend("topright", 
       c("Bayes", 
         expression(italic("\u03B4")*"-SMI"), 
         expression("("*italic("\u03B4")~"= 128)"),  
         expression(italic("\u03B7")*"-SMI"), 
         expression("("*italic("\u03B7")~"= 0.05)"), 
         "Cut"), 
       bty="n", 
       col=c("#243757", "#985775", NA, "#985775", NA,  "#FAD77D"),
       lty = c(NA, "solid", NA, "dashed", NA,  NA),
       pch = c(16, NA, NA, NA, NA, 16))

dev.off()


# Save the plot to working directory
cairo_pdf(file = "hpv_compare_theta_distr_scaled.pdf", height = 4.5, width = 5)
windows(family='serif')
par(mar = c(5, 4, 1, 1) + 0.2, lend = 1)

# Plot the bivariate posterior theta distribution of full Bayes analysis.
plot(hpv_bayes_param$theta1[c(1:7200)*5], 
     hpv_bayes_param$theta2[c(1:7200)*5], 
     cex = 0.1, las = 1,
     xlim = c(-3, -1.2), ylim = c(7,37),
     col="#243757", xlab = expression(italic("\u03B8")[1]), 
     ylab = expression(italic("\u03B8")[2]))

# Plot the bivariate posterior theta distribution of cut analysis.
points(hpv_cut_param[c(1:7200)*5, 'theta_1'], 
       hpv_cut_param[c(1:7200)*5, 'theta_2'], 
       cex = 0.1, col = "#FAD77D")

# Contour of the bivariate posterior theta distribution of 
# delta-SMI analysis with delta (= 1)*sqrt(y) on the Poisson module.
contour(hpv_abc_smi_ds_1_contour_z, 
        drawlabels=FALSE, nlevels = 7, col= "#985775", add=TRUE)

# Contour of the bivariate posterior theta distribution of 
# eta-SMI analysis with eta (= 0.05) on the Poisson module.
contour(hpv_eta_smi_0.05_contour_z, lty = "dashed",
        drawlabels=FALSE, nlevels = 7, col= "#985775", add=TRUE)

# Plot legend
legend("topright", c("Bayes", expression(italic("\u03B4")*"-SMI"), expression("("*italic("\u03B4")~"= 1)"),  
                     expression(italic("\u03B7")*"-SMI"), expression("("*italic("\u03B7")~"= 0.05)"), "Cut"), 
       bty="n", 
       col=c("#243757", "#985775", NA, "#985775",  NA, "#FAD77D"),
       lty = c(NA, "solid", NA, "dashed",  NA, NA),
       pch = c(16, NA, NA, NA, NA, 16))

dev.off()

