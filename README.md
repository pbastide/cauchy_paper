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
  * `2023-03-14_lizards_model_selection.rda`: results from script `R_scripts/simulations_models.R`.
  * `2023-03-23_simulations_paramters.rda`: results from script `R_scripts/simulations_parameters.R`.

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
