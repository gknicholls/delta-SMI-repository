etaCurve = function(par, etaVals = NULL, etaWAIC = NULL, delta_fit = NULL, targetVals = NULL){
  deltaVals = sqrt((par[1]/etaVals - par[1])) 
  eta_fit = sapply(targetVals, 
                   mapDeltaToEta,
                   spectrum = deltaVals, 
                   mappedVal = etaWAIC)
  res = sum(abs(eta_fit - delta_fit))
  return(res)
}

etaToDelta = function(par, etaVals = NULL, etaWAIC = NULL, delta_fit = NULL, targetVals = NULL){
  deltaVals =  par[1] - (log(etaVals/(1 - etaVals)))/par[2]
  
  
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

etaPoisCurveFit = optim(par = c(1, 1), etaToDelta, 
                        targetVals = c(0:8, 8.2, 8.4, 8.6, 8.8, 9, 10),
                        etaVals = hpv_eta_smi_pois_waic_all.df$eta_pois,
                        etaWAIC = hpv_eta_smi_pois_waic_all.df$elpd_waic,
                        delta_fit = hpv_abc_smi_waic_pois[,'elpd_waic'])
# Re-scaling eta for the Poisson module

etaToDeltaPois = sqrt(etaPoisCurveFit$par[1]/hpv_eta_smi_pois_waic_all.df$eta_pois[1:18] - etaPoisCurveFit$par[1]) 
etaToDeltaPois = etaPoisCurveFit$par[1]/hpv_eta_smi_pois_waic_all.df$eta_pois[1:18]^etaPoisCurveFit$par[2]

etaToDeltaPois = etaPoisCurveFit$par[1] - 
  log(hpv_eta_smi_pois_waic_all.df$eta_pois/(1 - hpv_eta_smi_pois_waic_all.df$eta_pois))/etaPoisCurveFit$par[2]

etaBinomCurveFit = optim(par = c(3, 0.8), etaToDeltaScaled, 
                         targetVals = c(0:8, 8.2, 8.4, 8.6, 8.8, 9, 10),
                         etaVals = hpv_eta_smi_binom_waic_all.df$eta_pois,
                         etaWAIC = hpv_eta_smi_binom_waic_all.df$elpd_waic,
                         delta_fit = hpv_abc_smi_waic_binom[,'elpd_waic'])
# Re-scaling eta for the binomial module
etaToDeltaBinom = etaBinomCurveFit$par[1] - 
  log(hpv_eta_smi_binom_waic_all.df$eta_pois/(1 - hpv_eta_smi_binom_waic_all.df$eta_pois))/etaBinomCurveFit$par[2]

library(plotrix)
cairo_pdf(file = "hpv_waic_eta_vs_delta.pdf", height = 5.25, width = 10)
par(mfrow = c(1,2), mar = c(5, 4, 5, 1)+0.2, cex.axis = 0.75, lend = 2)

# Plot the ELPD (based on WAIC) for the Poisson module
## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the delta-SMI analysis
plot(c(4, 5:8, 8.2, 8.4, 8.6, 8.8, 9:10)-4, hpv_abc_smi_waic_pois[c(1,6:15),'elpd_waic'], las = 1,
     xlab = expression(italic("\u03B4")~"over Poisson module"), 
     ylab = "ELPD (Poisson)", xaxt="n", pch = 4)


## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the eta-SMI analysis
lines(etaToDeltaPois[1:16]-4, hpv_eta_smi_pois_waic_all.df$elpd_waic[1:16], col="red")
lines(x = c(0, 1), y = hpv_eta_smi_pois_waic_all.df$elpd_waic[c(19,23)], lty = "dashed", col = "red") 

## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the full Bayes analysis
abline(h = hpv_bayes_waic_pois['elpd_waic'], lty = "dashed")

## Plot the the ELPD (based on WAIC) 
## for the Poisson module from the cut analysis
abline(h = hpv_cut_waic_pois['elpd_waic'], lty = "dotdash")

## x-axis labels  on the top side.
mtext(text = expression(italic("\u03B7")~"over Poisson module"), side = 3, line = 3.5)

## x-axis tick labels on the top side (i.e., the eta values).
axis(side = 3, at = c(etaToDeltaPois[c(2,4,9,12,13, 14, 15)] - 4, 0), las = 2,
     labels = formatC(hpv_eta_smi_pois_waic_all.df$eta_pois[c(2,4,9,12,13, 14, 15,22)], 
                      digits=4, format="f", flag="0"))

## x-axis lab on the bottom side (u.e., the delta values).
axis(side = 1, at = c(0:6), 
     label = c(expression(1,  2^5, 2^6, 2^7, 2^8, 2^9, 2^10)))


axis.break(1, 0.5, style="slash") 
axis.break(3, 0.5, style="slash") 

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
plot(x = c(0:10), 
     y = hpv_abc_smi_waic_binom[-c(10:13),'elpd_waic'], ylim  = c(-48, -33.5),
     xlab = expression(italic("\u03B4")~"over Poisson module"), 
     ylab = "ELPD (Binomial)", xaxt="n", pch = 4, las = 1)

## Plot the the ELPD (based on WAIC) 
## for the binomial module from the eta-SMI analysis
lines(etaToDeltaBinom, hpv_eta_smi_binom_waic_all.df$elpd_waic, col="red")

## Plot the the ELPD (based on WAIC) 
## for the binomial module from the eta-SMI analysis
abline(h = hpv_bayes_waic_binom['elpd_waic'], lty = "dashed")

## Plot the the ELPD (based on WAIC) 
## for the binomial module from the cut analysis
abline(h = hpv_cut_waic_binom['elpd_waic'], lty = "dotdash")

## x-axis labels  on the top side.
mtext(text = expression(italic("\u03B7")~"over Poisson module"), side = 3, line = 3.5)

## x-axis tick labels on the top side (i.e., the eta values).

axis(side = 3, at =  etaToDeltaBinom[c(2,4,9,12,14, 18,21, 22)], las = 2,
     labels = formatC(hpv_eta_smi_pois_waic_all.df$eta_pois[c(2,4,9,12, 14, 18,21,22)], 
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

