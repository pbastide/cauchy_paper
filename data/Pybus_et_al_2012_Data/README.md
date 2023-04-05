# Data and scripts to reproduce Pybus et al 2012.

## `xml` file analyses

* `WNV_GTR_gamma.xml`: Reproduces the analysis from Pybus et al 2012.
   Set up using BEAST tutorial: https://beast.community/workshop_continuous_diffusion_wnv.
   Produces:
   * `WNV_GTR_gamma.log`, `WNV_GTR_gamma.trees`: log files, omitted in the git repo because of their large size
   * `WNV_GTR_gamma_tips.log`, `WNV_GTR_gamma.ops`: log files
   * `WNV_GTR_gamma_MCC.tree`: MCC tree using treeannotator (cf `R_script/wnv.R`)
   * `WNV_GTR_gamma_data.txt`: jittered tip data extracted from the log (cf `R_script/wnv.R`)

* `WNV_HMC_Cauchy.xml`: Fits a "Cauchy RRW" on the jittered data `WNV_GTR_gamma_data.txt`,
  using the HMC technique from Fisher et al 2021 (www.doi.org/10.1093/sysbio/syaa056).
  Produces:
   * `WNV_HMC_Cauchy.log`, `WNV_HMC_Cauchy.trees`: log files
   
* `WNV_HMC_Cauchy_lat.xml` and `WNV_HMC_Cauchy_long.xml`: Fits independent "Cauchy RRW" on the lattitude and longitude traits from the jittered data `WNV_GTR_gamma_data.txt`,
  using the HMC technique from Fisher et al 2021 (www.doi.org/10.1093/sysbio/syaa056).
  Produces:
   * `WNV_HMC_Cauchy_lat.log`, `WNV_HMC_Cauchy_lat.trees`, `WNV_HMC_Cauchy_long.log`, `WNV_HMC_Cauchy_long.trees`: log files

## Data

* `WNV.fas`, `WNV_lat_long.txt`: data downloaded from: https://beast.community/workshop_continuous_diffusion_wnv.
