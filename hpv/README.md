**delta-SMI**

To run a delta-SMI analysis on the HPV dataset, the command at the command line is

Rscript hpv_abc_smi_delta_gen.R [path to inputFile]

The template for the input file of hpv_abc_smi_delta_gen.R is abc-smi_hpv_delta_gen_input_template.txt.

---

**delta-SMI scaled by square root of y**

To run a delta-SMI analysis with delta scaled by the square root of y, the command at the command line is

Rscript hpv_abc_smi_delta_scaled_gen.R [path to inputFile]

The template for the input file of hpv_abc_smi_delta_scaled_gen.R is abc-smi_hpv_model_whole_delta_scaled_input_template.txt.

---

**eta-SMI**

To run a eta-SMI analysis on the HPV dataset, the command is

Rscript hpv_eta_smi.R [path to inputFile]

The template for the input file of hpv_eta_smi.R is hpv_eta_smi_input_template.txt.

---

**Full Bayes and cut analyses**

R commands for the full Bayes and cut analyses can be found in hpv_bayes_and_cut.R.

---

**ELPD calculation**

R commands for calculating the ELPD of the delta-SMI, eta-SMI, full Bayes and cut model analyses on the HPV dataset can be found hpv_elpd_all.Rmd.

This file also contain R commands for generating the files required for making the figures.

---

**Figures**

The R commands for figure 5 in the manuscript can be found in hpv_compare_theta_distr.R.

The R commands for figure 6 in the manuscript can be found in hpv_waic_eta_vs_delta_plot.R

Note: input files for generating the figures are already available in the hpv Github folder. See comments in hpv_compare_theta_distr.R and hpv_waic_eta_vs_delta_plot.R for details.

---

**Corrections**

The likelihood calculation has been corrected in abc-smi_hpv_model_whole.stan and abc-smi_hpv_model_whole_delta_scaled.stan files.
