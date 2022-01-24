etaCurve = function(par = NULL, 
                    etaVals = NULL, 
                    etaWAIC = NULL, 
                    delta_fit = NULL, 
                    targetVals = NULL){
  deltaVals = sqrt((par[1]/etaVals - par[1]))
  eta_fit = sapply(targetVals, 
                   mapDeltaToEta,
                   spectrum = deltaVals, 
                   mappedVal = etaWAIC)
  res = sum(abs(eta_fit - delta_fit))
  return(res)
}


# Load hpv_bayes_waic_pois and hpv_bayes_waic_binom.
# Assuming hpv_bayes_elpd.RData is in the working directory.
load("hpv_bayes_elpd.RData")
# Load hpv_cut_waic_pois and hpv_cut_waic_binom.
# Assuming hpv_cut_elpd.RData is in the working directory.
load("hpv_cut_elpd.RData")

# Load hpv_abc_smi_waic_pois and hpv_abc_smi_waic_binom.
# Assuming hpv_abc_smi_elpd.RData is in the working directory.
load("hpv_abc_smi_elpd.RData")
# Load hpv_eta_smi_pois_waic_all.df and hpv_eta_smi_binom_waic_all.df
# Assuming hpv_abc_smi_elpd.RData is in the working directory.
load("hpv_eta_smi_elpd.RData")

etaPoisCurveFit = optim(par = c(0.8), etaCurve, 
                        targetVals = c(0:8),
                        etaVals = hpv_eta_smi_pois_waic_all.df$eta_pois[-1],
                        etaWAIC = hpv_eta_smi_pois_waic_all.df$elpd_waic[-1],
                        delta_fit = hpv_abc_smi_waic_pois[1:9,'elpd_waic'],
                        method="Brent", lower=c(0), upper=c(1000))
# Re-scaling eta for the Poisson module
etaToDeltaPois = sqrt((etaPoisCurveFit$par / hpv_eta_smi_pois_waic_all.df$eta_pois[-1] - 
                         etaPoisCurveFit$par))

etaBinomCurveFit = optim(par = c(0.8), etaCurve, 
                         targetVals = c(0:10),
                         etaVals = hpv_eta_smi_binom_waic_all.df$eta_pois[-1],
                         etaWAIC = hpv_eta_smi_binom_waic_all.df$elpd_waic[-1],
                         delta_fit = hpv_abc_smi_waic_binom[,'elpd_waic'],
                         method="Brent", lower=c(0), upper=c(1000))
# Re-scaling eta for the binomial module
etaToDeltaBinom = sqrt((etaBinomCurveFit$par /hpv_eta_smi_binom_waic_all.df$eta_pois[-1]  
                        - etaBinomCurveFit$par))

cairo_pdf(file = "hpv_waic_eta_vs_delta.pdf", height = 5.25, width = 10)
par(mfrow = c(1,2), mar = c(5, 4, 5, 1)+0.2, cex.axis = 0.75, lend = 2)

# Plot the ELPD (based on WAIC) for the Poisson module
## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the delta-SMI analysis
plot(c(0:10), hpv_abc_smi_waic_pois[,'elpd_waic'], las = 1,
     xlab = expression(italic("\u03B4")~"over Poisson module"), 
     ylab = "ELPD (Poisson)", xaxt="n", pch = 4)

## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the eta-SMI analysis
lines(etaToDeltaPois, hpv_eta_smi_pois_waic_all.df$elpd_waic[-1], col="red")

## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the full Bayes analysis
abline(h = hpv_bayes_waic_pois['elpd_waic'], lty = "dashed")

## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the cut analysis
abline(h = hpv_cut_waic_pois['elpd_waic'], lty = "dotdash")

## x-axis labels  on the top side.
mtext(text = expression(italic("\u03B7")~"over Poisson module"), side = 3, line = 3.5)

## x-axis tick labels on the top side (i.e., the eta values).
axis(side = 3, at = etaToDeltaPois[c(2:5, 7, 9,  12, 21)], las = 2,
     labels = formatC(hpv_eta_smi_pois_waic_all.df$eta_pois[-1][c(2:5, 7, 9,  12, 21)], 
                      digits=4, format="f", flag="0"))

## x-axis lab on the bottom side (u.e., the delta values).
axis(side = 1, at = c(0:5)*2, 
     label = c(expression(1, 2^2, 2^4, 2^6, 2^8, 2^10)))

## Plot the legend.
legend("left", 
       c("Bayes", 
         expression(italic("\u03B4")*"-SMI"),
         expression(italic("\u03B7")*"-SMI"), 
         "Cut"),
       pch = c(NA, 4, NA,  NA), bty = "n",
       lty = c("dashed", NA, "solid", "dotdash"),
       col = c("black", "black", "red", "black"))

# Plot the ELPD (based on WAIC) for the binomial module
## Plot the the ELPD (based on WAIC) 
## for the binomial module from the delta-SMI analysis
plot(c(0:10), hpv_abc_smi_waic_binom[,'elpd_waic'], ylim  = c(-48, -33.5),
     xlab = expression(italic("\u03B4")~"over Poisson module"), 
     ylab = "ELPD (Binomial)", xaxt="n", pch = 4, las = 1)

## Plot the the ELPD (based on WAIC) 
## for the binomial module from the eta-SMI analysis
lines(etaToDeltaBinom, hpv_eta_smi_binom_waic_all.df$elpd_waic[-1], col="red")

## Plot the the ELPD (based on WAIC) 
## for the binomial module from the eta-SMI analysis
abline(h = hpv_bayes_waic_binom['elpd_waic'], lty = "dashed")

## Plot the the ELPD (based on WAIC) 
## for the binomial module from the cut analysis
abline(h = hpv_cut_waic_binom['elpd_waic'], lty = "dotdash")

## x-axis labels  on the top side.
mtext(text = expression(italic("\u03B7")~"over Poisson module"), side = 3, line = 3.5)

## x-axis tick labels on the top side (i.e., the eta values).
axis(side = 3, at = etaToDeltaBinom[c(2:5, 7, 9, 12, 21)], las = 2,
     labels = formatC(hpv_eta_smi_pois_waic_all.df$eta_pois[-1][c(2:5, 7, 9,  12, 21)], 
                      digits=4, format="f", flag="0"))
## x-axis lab on the bottom side (u.e., the delta values).
axis(side = 1, at = c(0:5)*2, label = c(expression(1, 2^2, 2^4, 2^6, 2^8, 2^10)))

## Plot the legend.
legend("right", 
       c("Bayes", 
         expression(italic("\u03B4")*"-SMI"),
         expression(italic("\u03B7")*"-SMI"), 
         "Cut"),
       pch = c(NA, 4, NA,  NA), bty = "n",
       lty = c("dashed", NA, "solid", "dotdash"),
       col = c("black", "black", "red", "black"))

dev.off()

