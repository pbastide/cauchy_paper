# Data and scripts to reproduce Pybus et al 2012.

## `xml` file analyses

* `WNV_GTR_gamma.xml`: Reproduces the analysis from Pybus et al 2012.
   Set up using BEAST tutorial: https://beast.community/workshop_continuous_diffusion_wnv.
   Produces:
   * `WNV_GTR_gamma_tips.log`: tip traits log file from the MCMC run.
   * `WNV_GTR_gamma_MCC.tree`: MCC tree using treeannotator (cf `R_script/wnv.R`)
   * `WNV_GTR_gamma_data.txt`: jittered tip data extracted from the log (cf `R_script/wnv.R`)

* `WNV_HMC_Cauchy.xml`: Fits a "Cauchy RRW" on the jittered data `WNV_GTR_gamma_data.txt`,
  using the HMC technique from Fisher et al 2021 (www.doi.org/10.1093/sysbio/syaa056).
  Produces:
   * `WNV_HMC_Cauchy.log`: parameters log file from the MCMC run.
   * `WNV_HMC_Cauchy.trees`: tree log files from the MCMC run.
   
* `WNV_HMC_Cauchy_lat.xml`: Fits an independent "Cauchy RRW" on the lattitude trait from the jittered data `WNV_GTR_gamma_data.txt`,
  using the HMC technique from Fisher et al 2021 (www.doi.org/10.1093/sysbio/syaa056).
  Produces:
   * `WNV_HMC_Cauchy_lat.log`: parameters log file from the MCMC run.

* `WNV_HMC_Cauchy_long.xml`: Fits an independent "Cauchy RRW" on the longitude trait from the jittered data `WNV_GTR_gamma_data.txt`,
  using the HMC technique from Fisher et al 2021 (www.doi.org/10.1093/sysbio/syaa056).
  Produces:
   * `WNV_HMC_Cauchy_long.log`: parameters log file from the MCMC run.

## Data

Data downloaded from: https://beast.community/workshop_continuous_diffusion_wnv:

* `WNV.fas`: sequence data file.
* `WNV_lat_long.txt`: location data file.
* `ordered_accession_numbers.txt`: accession number for all sequences.
* `ordered_states.txt`: state of each sequence, ordered as `ordered_accession_numbers.txt`.
