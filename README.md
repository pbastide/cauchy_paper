# Repository to reproduce the analyses in Bastide & Didier 2023.

## Content

* `cauchy_paper.Rproj`: `R` project that can be opened with `R Studio` for convenience. 
All the `R` scripts assume that the working directory is the root directory of the project.

* `R_scripts`: `R` scripts to reproduce the simulations and plots.
Please see the header of each file for a more detailed description of each script.
  * `lizard.R`: Analysis of the lizard dataset of Mahler et al. 2013 (www.doi.org/10.1126/science.1232392).
  * `wnv.R`: Analysis of the WNV dataset of Pybus et al. 2012 (www.doi.org/10.1073/pnas.1206598109).
  * `simulations_models.R`: simulation study for model selection.
  * `simulations_parameters.R`: simulation study for parameter inference.

* `data`: data used in the empirical analyses.
  * `Mahler_et_al_2013_Data`: data from the Dryad repository associated with Mahler et al. 2013 (www.doi.org/10.5061/dryad.9g182).
  * `Pybus_et_al_2012_Data`: data, `xml` analysis scripts and results to reproduce the study of Pybus et al. 2012, using beast tutorial (https://beast.community/workshop_continuous_diffusion_wnv)

* `results`: saved results from the `R_scripts`.
  * `2023-07-25_lizards_model_selection.rda`: results from script `R_scripts/simulations_models.R`. Contains:
    * `alldat`: a list with all simulated trait values.
    * `phy`: a `phylo` object with the normalized lizard tree used in simulations.
    * `resAll`: a data frame with the results of all analyses on all datasets.
    * `params_BM`, `params_EB`, `params_OU`, `a`, `all_alpha_NIG`, `all_delta_NIG`, `all_sigma_noise`, `alpha_1_all`, `alpha_OU`, `beta_NIG`, `datestamp`, `delta_1`, `disp_cau`, `ex_kurt`, `mad_EB`, `mado`, `models`, `mu`, `Ncores`, `Nrep`, `r_EB`, `reqpckg`, `results_directory`, `sigma_BM`, `sigma_EB`, `sigma_OU`, `sim_models`, `t_12_OU`, `var_EB`, `vari`: parameter vectors used in the simulations.
    * `git_all_models`, `fit_phylolm`, `fit_phylolm_all`, `model_nig`, `sim_model_fun`: instrumental simulation and fitting functions.
  * `2023-03-23_simulations_paramters.rda`: results from script `R_scripts/simulations_parameters.R`. Contains:
    * `all_sims`: a list with all simulated trees and trait values.
    * `res`: a data frame with the results of all analyses on all datasets.
    * `all_disp_cau`, `all_init`, `all_lambda`, `all_methods`, `all_models`, `all_n_tips`, `datestamp`, `Ncores`, `Nrep`, `reqpckg`, `result_directory`, `values_root`: parameter vectors used in the simulations.
    * `fit_params`, `get_name`, `get_params`, `sim`: instrumental simulation and fitting functions.

## Requirements

* [`R`](https://cran.r-project.org/index.html). Version 4.0 or later.

* `R` package [`cauphy`](https://github.com/gilles-didier/cauphy/).
The development version of the package can be installed as:
```R
devtools::install_github("gilles-didier/cauphy")
```

* `BEAST` v1.10.5 pre-release `#dc4f12c`, installed on a mac with:
```
brew install beast --HEAD
```
