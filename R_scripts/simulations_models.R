# Copyright (C) {2023} {PB, GD}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Program to conduct the simulations for the assessment of the model selection.

library(here)
library(ape)
library(cauphy)
library(phylolm)
library(fBasics) # For NIG simulations

################################################################################
## File management
################################################################################

datestamp <- format(Sys.time(), "%Y-%m-%d")
results_directory <- here("results", paste0(datestamp, "_"))

################################################################################
## Lizard Tree
################################################################################

# Empirical tree
phy <- read.tree(file = here("data", "Mahler_et_al_2013_Data", "GA_Anolis_MCC.tre"))
# normalize to unit height
phy$edge.length <- phy$edge.length / max(vcv(phy))

################################################################################
## Parameters for Data Simulation
################################################################################
mu <- 0
vari <- 1.0
mado <- qnorm(3/4, sd = sqrt(vari))

# BM
sigma_BM <- sqrt(vari)
params_BM <- list(ancestral.state = mu, sigma2 = sigma_BM^2)
all.equal(qnorm(3/4, sd = sigma_BM), mado)

# Cauchy
disp_cau <- mado # qnorm(3/4) to ensure mad(x) = 1 same as Gaussian
all.equal(qcauchy(3/4, scale = disp_cau), mado)

# OU
t_12_OU <- 0.2
alpha_OU <- log(2) / t_12_OU
sigma_OU <- sqrt(2 * alpha_OU * vari)
params_OU <- list(ancestral.state = mu, sigma2 = sigma_OU^2, optimal.value = mu, alpha = alpha_OU)
all.equal(sigma_OU^2 / (2 * alpha_OU), vari) # stationary variance of 1
all.equal(qnorm(3/4, sd = sqrt(sigma_OU^2 / (2 * alpha_OU))), mado)

# EB
sigma_EB <- sigma_OU
r_EB <- -alpha_OU
params_EB <- list(ancestral.state = mu, sigma2 = sigma_EB^2, rate = r_EB)
var_EB <- sigma_EB^2 * (exp(r_EB) - 1) / (r_EB)
mad_EB <- qnorm(3/4, sd = sqrt(var_EB))

# NIG
beta_NIG <- 0
ex_kurt <- c(0.1, 0.5, 1, 2, 3, 4, 5, 10)
alpha_1_all <- 1 / ex_kurt * 3
delta_1 <- 1
a <- sapply(alpha_1_all, function(alpha_1) qnorm(3/4) / qnig(3/4, alpha = alpha_1, delta = delta_1, beta = 0, mu = 0))
all_alpha_NIG <- alpha_1_all / a
all_delta_NIG <- delta_1 * a
for (i in seq_along(alpha_1_all)) {
  stopifnot(all.equal(qnig(3/4, alpha = all_alpha_NIG[i], delta = all_delta_NIG[i], beta = 0, mu = 0), mado, tol = 1e-2))
}
all_alpha_NIG
all_delta_NIG
all_delta_NIG / all_alpha_NIG # variances
3 / (all_delta_NIG * all_alpha_NIG) # excess kurtosis
# Noise
all_sigma_noise <- c(0)

################################################################################
## Simulation Functions
################################################################################

# Models
sim_models <- c("BM", "OU", "EB", paste0("NIG", seq_along(alpha_1_all)), "Cauchy")

# Simulations NIG
model_nig <- function(x, l, alpha_NIG, beta_NIG, delta_NIG) {
  y <- rnig(1, alpha = alpha_NIG, beta = beta_NIG, delta = delta_NIG * l, mu = x)
}

sim_model_fun <- function(model, phy,
                          sigma_error) {
  
  if (model == "Cauchy") {
    dat <- rTraitCauchy(1, phy, model = "cauchy", parameters = list(disp = disp_cau, root.value = mu))
  } else if (model == "BM") {
    dat <- rTrait(1, phy, model = model, parameters = params_BM)
  } else if (model == "OU") {
    dat <- rTrait(1, phy, model = model, parameters = params_OU)
  } else if (model == "EB"){
    dat <- rTrait(1, phy, model = model, parameters = params_EB)
  } else if (grepl("NIG", model)){
    imod <- as.numeric(substr(model, nchar(model), nchar(model)))
    dat <- rTraitCont(phy, model = model_nig, root.value = mu,
                      alpha_NIG = all_alpha_NIG[imod], beta_NIG = beta_NIG, delta_NIG = all_delta_NIG[imod])
  }
  
  err <- rnorm(length(phy$tip.label), mean = 0, sd = sigma_error)
  dat <- dat + err
  
  return(dat)
}

set.seed(1289)
alldatOne <- sapply(sim_models, sim_model_fun, phy, sigma_error = 0)

library(PhylogeneticEM)
plot(params_BM(p = ncol(alldatOne)), data = t(alldatOne), phylo = phy)

################################################################################
## Fit Functions
################################################################################

models <- c("BM", "OUfixedRoot", "EB", "BMerr", "OUfixedRooterr", "EBerr", "cauchy", "cauchyLambda")

fit_phylolm <- function(dat, model) {
  if (model == "cauchy") {
    res <- fitCauchy(phy, dat, method = "fixed.root", model = model)
  } else if (model == "cauchyLambda") {
    res <- fitCauchy(phy, dat, method = "fixed.root", model = "lambda")
  } else if (model == "OUfixedRoot") {
    res <- phylolm(dat ~ 1, phy = phy, model = model, measurement_error = FALSE, upper.bound = 10)
  } else if (model == "EB") {
    res <- phylolm(dat ~ 1, phy = phy, model = model, measurement_error = FALSE, lower.bound = -10)
  } else if (model == "BM") {
    res <- phylolm(dat ~ 1, phy = phy, model = model, measurement_error = FALSE)
  } else if (model == "OUfixedRooterr") {
    res <- phylolm(dat ~ 1, phy = phy, model = "OUfixedRoot", measurement_error = TRUE, upper.bound = 10)
  } else if (model == "EBerr") {
    res <- phylolm(dat ~ 1, phy = phy, model = "EB", measurement_error = TRUE, lower.bound = -10)
  } else if (model == "BMerr") {
    res <- phylolm(dat ~ 1, phy = phy, model = "BM", measurement_error = TRUE)
  }
  if (!inherits(try, "try-error")) {
    return(list(logLik = res$logLik, AIC = res$aic))
  } else {
    return(list(logLik = -Inf, AIC = +Inf))
  }
}

fit_all_models <- function(dat, error, i) {
  res <- sapply(models, fit_phylolm, dat = dat)
  res <- as.data.frame(res)
  res$score <- rownames(res)
  res <- as.data.frame(tidyr::pivot_longer(res, cols = all_of(models)))
  res$rep <- i
  res$error <- error
  return(res)
}

fit_phylolm_all <- function(alldat, error, i) {
  res <- apply(alldat, 2, fit_all_models, error = error, i = i)
  res <- do.call(rbind, res)
  res$sim <- sub("\\.[0-9]*", "", rownames(res))
  return(res)
}

################################################################################
## Simulations
################################################################################

Nrep <- 100

set.seed(1289)
alldat <- NULL
for (sigma_error in c(0)) {
  alldat[[paste(sigma_error)]] <- NULL
  for (i in 1:Nrep) {
    ## all data
    alldat[[paste(sigma_error)]][[i]] <- sapply(sim_models, sim_model_fun, phy, sigma_error = sigma_error)
  }
}

# Check vars
mean(sapply(1:Ntip(phy), function(i) var(sapply(alldat[["0"]], function(x) x[i, 1]))))
mean(sapply(1:Ntip(phy), function(i) var(sapply(alldat[["0"]], function(x) x[i, 2]))))
mean(sapply(1:Ntip(phy), function(i) var(sapply(alldat[["0"]], function(x) x[i, 3]))))
# Check MAD
for (k in 1:ncol(alldat[["0"]][[1]])) {
  if (sim_models[k] == "EB") {
    stopifnot(all.equal(mean(sapply(1:100, function(i) mad(sapply(alldat[["0"]], function(x) x[i, k]), constant = 1))),
                        mad_EB, tol = 1e-2))
  } else {
    stopifnot(all.equal(mean(sapply(1:100, function(i) mad(sapply(alldat[["0"]], function(x) x[i, k]), constant = 1))),
                        mado, tol = 1e-1))
  }
}

################################################################################
## Fits
################################################################################

library(doParallel)
## Required packages to pass to all nodes
reqpckg <- c("cauphy", "phylolm")
## Number of cores
Ncores <- 6
## Register nodes
cl <- makeCluster(Ncores, outfile = "")
registerDoParallel(cl)

resAll <- foreach (sigma_error = all_sigma_noise) %:%
  foreach(i = 1:Nrep, .packages = reqpckg, .verbose = FALSE) %dopar% {
    res <- fit_phylolm_all(alldat[[paste0(sigma_error)]][[i]], error = sigma_error, i)
    return(res)
  }

## Stop the parallel cluster
stopCluster(cl)

resAll <- lapply(resAll, function(x) do.call(rbind, x))
resAll <- do.call(rbind, resAll)

save.image(file = paste0(results_directory, "lizards_model_selection.rda"))


################################################################################
## Plots
################################################################################

load(file = paste0(results_directory, "lizards_model_selection.rda"))

library(dplyr)
library(tidyr)

resAll$name <- sub("err", "", resAll$name)

sel <- subset(resAll, score == "AIC") %>% group_by(rep, error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
sel$sim <- factor(sel$sim, levels = sim_models)


library(colorspace)
colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), 0:2/10),
             lighten(rep(colorblindr:::palette_OkabeIto[6], 2), 0:1/10))
# colPlot[3] <- darken(colPlot[3], 0.3)

library(ggplot2)
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  theme_minimal() + 
  facet_wrap(vars(error), nrow = 2) + 
  scale_fill_manual(values = colPlot)
