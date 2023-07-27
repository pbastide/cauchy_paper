# cauchy_paper by PB, GD
#   
# To the extent possible under law, the person who associated CC0 with
# cauchy_paper has waived all copyright and related or neighboring rights
# to cauchy_paper.
# 
# You should have received a copy of the CC0 legal code along with this
# work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

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
# Root ancestral value (fixed for all simulations)
mu <- 0
# Reference variance of the BM process, and of the Gaussian tip additional error.
all_vars <- data.frame(var_BM = c(1.0, 5.0, 10.0, 50.0, 100.0, rep(1.0, 4)),
                       var_error = c(rep(0.0, 5), c(0.01, 0.1, 1.0, 10.0)))
# Excess kurtosis for the NIG
ex_kurt <- c(1, 3, 5, 10)

#' @title Get simulation parameters
#'
#' @description
#' Compute the parameters of the various processes so that all have the same MAD
#' (except for the EB, see Bastide & Didier, 2023).
#' 
#' @param mu the root ancestral value.
#' @param vari the reference variance of the BM process.
#' 
#' @return list of all simulation parameters
#' 
get_simulation_parameters <- function(mu, vari) {
  ## MAD value corresponding to a BM with variance vari
  mado <- qnorm(3/4, sd = sqrt(vari))
  
  # BM
  sigma_BM <- sqrt(vari)
  params_BM <- list(ancestral.state = mu, sigma2 = sigma_BM^2)
  stopifnot(isTRUE(all.equal(qnorm(3/4, sd = sigma_BM), mado)))
  
  # Cauchy
  disp_cau <- mado # qnorm(3/4) to ensure mad(x) = 1 same as Gaussian
  params_CP <- list(disp = disp_cau, root.value = mu)
  stopifnot(isTRUE(all.equal(qcauchy(3/4, scale = disp_cau), mado)))
  
  # OU
  t_12_OU <- 0.2
  alpha_OU <- log(2) / t_12_OU
  sigma_OU <- sqrt(2 * alpha_OU * vari)
  params_OU <- list(ancestral.state = mu, sigma2 = sigma_OU^2, optimal.value = mu, alpha = alpha_OU)
  stopifnot(isTRUE(all.equal(sigma_OU^2 / (2 * alpha_OU), vari))) # stationary variance of vari
  stopifnot(isTRUE(all.equal(qnorm(3/4, sd = sqrt(sigma_OU^2 / (2 * alpha_OU))), mado)))
  
  # EB
  sigma_EB <- sigma_OU
  r_EB <- -alpha_OU
  params_EB <- list(ancestral.state = mu, sigma2 = sigma_EB^2, rate = r_EB)
  var_EB <- sigma_EB^2 * (exp(r_EB) - 1) / (r_EB)
  mad_EB <- qnorm(3/4, sd = sqrt(var_EB))
  
  # NIG
  beta_NIG <- 0
  alpha_1_all <- 1 / ex_kurt * 3
  delta_1 <- 1
  a <- sapply(alpha_1_all, function(alpha_1) qnorm(3/4) / qnig(3/4, alpha = alpha_1, delta = delta_1, beta = 0, mu = 0))
  all_alpha_NIG <-  rep(0,length(ex_kurt))
  all_delta_NIG <- rep(0,length(ex_kurt))
  mad_NIG = rep(0,length(ex_kurt))
  for (i in seq_along(ex_kurt)) {
  	sca = optimize(function(x) abs(qnorm(3/4, sd = sigma_BM) - qnig(3/4, alpha = 3/(ex_kurt[i]*a[i]*x), delta = a[i]*x, beta = 0, mu = 0)), c(0,1000), tol = 0.00001)
  	all_alpha_NIG[i] = 3/(ex_kurt[i]*a[i]*sca$minimum)
  	all_delta_NIG[i] = a[i]*sca$minimum
  	mad_NIG[i] = qnig(3/4, alpha = all_alpha_NIG[i], delta = all_delta_NIG[i], beta = 0, mu = 0)
  	stopifnot(isTRUE(all.equal(mad_NIG[i], mado, tol = 1e-5)))
  }
  stopifnot(isTRUE(all.equal(3 / (all_delta_NIG * all_alpha_NIG), ex_kurt)))
  
  return(list(params_BM = params_BM,
              params_OU = params_OU,
              params_EB = params_EB,
              params_CP = params_CP,
              all_alpha_NIG = all_alpha_NIG,
              beta_NIG = beta_NIG,
              all_delta_NIG = all_delta_NIG))
}

################################################################################
## Simulation Functions
################################################################################

# Models
sim_models <- c("BM", "OU", "EB", paste0("NIG", seq_along(ex_kurt)), "Cauchy")

#' @title Simulate a NIG increment
#' 
#' @param x the value at the beginning of the branch
#' @param l the length of the branch
#' @param alpha_NIG alpha value of the NIG
#' @param beta_NIG beta value of the NIG
#' @param delta_NIG delta value of the NIG
#' 
#' @return simulated value at the end of the branch
#' 
model_nig <- function(x, l, alpha_NIG, beta_NIG, delta_NIG) {
  y <- rnig(1, alpha = alpha_NIG, beta = beta_NIG, delta = delta_NIG * l, mu = x)
}

#' @title Simulate process of a given model
#' 
#' @param model the model to simulate
#' @param phy the fixed phylogeny
#' @param var_BM the variance of the reference BM process
#' @param sigma2_error the variance of the extra Gaussian noise
#' 
#' @return simulated traits at the tips of the tree
#' 
sim_model_fun <- function(model, phy, var_BM, var_error) {
  
  all_params_process <-  get_simulation_parameters(mu = mu, vari = var_BM)
  
  if (model == "Cauchy") {
    dat <- rTraitCauchy(1, phy, model = "cauchy", parameters = all_params_process$params_CP)
  } else if (model == "BM") {
    dat <- rTrait(1, phy, model = model, parameters = all_params_process$params_BM)
  } else if (model == "OU") {
    dat <- rTrait(1, phy, model = model, parameters = all_params_process$params_OU)
  } else if (model == "EB"){
    dat <- rTrait(1, phy, model = model, parameters = all_params_process$params_EB)
  } else if (grepl("NIG", model)){
    imod <- as.numeric(substr(model, nchar(model), nchar(model)))
    dat <- rTraitCont(phy, model = model_nig, root.value = mu,
                      alpha_NIG = all_params_process$all_alpha_NIG[imod],
                      beta_NIG = all_params_process$beta_NIG,
                      delta_NIG = all_params_process$all_delta_NIG[imod])
  }
  
  err <- rnorm(length(phy$tip.label), mean = 0, sd = sqrt(var_error))
  dat <- dat + err
  
  return(dat)
}

## Test one simulation
set.seed(1289)
alldatOne <- sapply(sim_models, sim_model_fun, phy, var_BM = 1, var_error = 0)

library(PhylogeneticEM)
plot(params_BM(p = ncol(alldatOne)), data = t(alldatOne), phylo = phy)

################################################################################
## Fit Functions
################################################################################

models <- c("BM", "OUfixedRoot", "EB", "BMerr", "OUfixedRooterr", "EBerr", "cauchy", "cauchyLambda")
curBound <- 10

#' @title Fit data with a model
#' 
#' @param dat tip trait data
#' @param model the model to fit
#' 
#' @return a list, with the likelihood and the AIC values of the fitted model.
#' 
fit_phylolm <- function(dat, model) {
  if (model == "cauchy") {
    res <- fitCauchy(phy, dat, method = "fixed.root", model = model)
  } else if (model == "cauchyLambda") {
    res <- fitCauchy(phy, dat, method = "fixed.root", model = "lambda")
  } else if (model == "OUfixedRoot") {
    res <- phylolm(dat ~ 1, phy = phy, model = model, measurement_error = FALSE, upper.bound = curBound)
  } else if (model == "EB") {
    res <- phylolm(dat ~ 1, phy = phy, model = model, measurement_error = FALSE, lower.bound = -curBound)
  } else if (model == "BM") {
    res <- phylolm(dat ~ 1, phy = phy, model = model, measurement_error = FALSE)
  } else if (model == "OUfixedRooterr") {
    res <- phylolm(dat ~ 1, phy = phy, model = "OUfixedRoot", measurement_error = TRUE, upper.bound = curBound)
  } else if (model == "EBerr") {
    res <- phylolm(dat ~ 1, phy = phy, model = "EB", measurement_error = TRUE, lower.bound = -curBound)
  } else if (model == "BMerr") {
    res <- phylolm(dat ~ 1, phy = phy, model = "BM", measurement_error = TRUE)
  }
  if (!inherits(try, "try-error")) {
    return(list(logLik = res$logLik, AIC = res$aic))
  } else {
    return(list(logLik = -Inf, AIC = +Inf))
  }
}

#' @title Fit data with all models
#' 
#' @param dat tip trait data
#' @param var_BM the BM variance used for the simulation (for book keeping)
#' @param var_error the error variance used for the simulation (for book keeping)
#' @param i the replication number
#' 
#' @return a data frame, with the likelihoods and AIC values of the fit of all
#' models on the data.
#' 
fit_all_models <- function(dat, var_BM, var_error, i) {
  res <- sapply(models, fit_phylolm, dat = dat)
  res <- as.data.frame(res)
  res$score <- rownames(res)
  res <- as.data.frame(tidyr::pivot_longer(res, cols = all_of(models)))
  res$rep <- i
  res$var_BM <- var_BM
  res$var_error <- var_error
  return(res)
}

#' @title Fit all data with all models
#' 
#' @param alldat tip trait data generated by all models for a given replicate
#' @param var_BM the BM variance used for the simulation (for book keeping)
#' @param var_error the error variance used for the simulation (for book keeping)
#' @param i the replication number
#' 
#' @return a data frame, with the likelihoods and AIC values of the fit of all
#' models on the data.
#' 
fit_phylolm_all <- function(alldat, var_BM, var_error, i) {
  res <- apply(alldat, 2, fit_all_models, var_BM = var_BM, var_error = var_error, i = i)
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
for (pp in seq_len(nrow(all_vars))) {
  var_BM <- all_vars[pp, 1]
  var_error <- all_vars[pp, 2]
  for (i in 1:Nrep) {
    ## all data
    alldat <- rbind(alldat, data.frame(dat = sapply(sim_models, sim_model_fun, phy, var_BM = var_BM, var_error = var_error),
                                       var_BM = var_BM,
                                       var_error = var_error,
                                       rep = i))
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

resAll <- foreach (pp = seq_len(nrow(all_vars))) %:%
  foreach(i = 1:Nrep, .packages = reqpckg, .verbose = FALSE) %dopar% {
    var_BM <- all_vars[pp, 1]
    var_error <- all_vars[pp, 2]
    dd <- alldat[alldat$var_BM == var_BM & alldat$var_error == var_error & alldat$rep == i, seq_along(sim_models)]
    rownames(dd) <- sub("[0-9]*$", "", rownames(dd))
    res <- fit_phylolm_all(dd, var_BM, var_error, i)
    return(res)
  }

## Stop the parallel cluster
stopCluster(cl)

resAll <- lapply(resAll, function(x) do.call(rbind, x))
resAll <- do.call(rbind, resAll)
resAll$sim <- sub("\\.[0-9]*$", "", resAll$sim)
resAll$sim <- sub("^dat", "", resAll$sim)

save.image(file = paste0(results_directory, "lizards_model_selection.rda")) 

################################################################################
## Plot from main text : var_BM = 1, var_error = 0
################################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

load(paste0(results_directory, "lizards_model_selection.rda")) 

# Merge error and no error for gaussian models
resAll$name <- sub("err", "", resAll$name)
# Names
for (i in seq_along(ex_kurt)) {
  resAll$sim <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), resAll$sim)
  sim_models <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), sim_models)
}
resAll$sim <- sub("Cauchy", "CP", resAll$sim)
sim_models <- sub("Cauchy", "CP", sim_models)
models <- c("BM", "OUfixedRoot", "EB", "cauchy", "cauchyLambda")
models_pretty <- c("BM", "OU", "EB", "CP", "$\\lambda$ CP")

# Subsets Format
resAll <- subset(resAll, var_error == 0 & var_BM == 1)
resAll <- subset(resAll, name != "cauchyLambda")
sel <- subset(resAll, score == "AIC") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)

# colors
# colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), c(-0.8, -0.3, 0.2)),
#              lighten(rep(colorblindr:::palette_OkabeIto[1], 2), c(0.6, 0.8)))
colPlot <- c("#002418", "#026D4E", "#3AB589", "#FFD6AD", "#FFEBD9")

# Plot
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  theme_minimal() + 
  xlab("Simulation Model") + ylab("Count") + 
  scale_fill_manual(values = colPlot, name = "Selected Model:") + 
  theme(#legend.direction="vertical",
    # panel.grid = element_blank(),
    legend.position="top",
    legend.justification = c("left", "top"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.key.size = unit(12, 'pt'),
    text = element_text(size = 9))
# guides(fill = guide_legend(ncol = 4, byrow = TRUE))
colorblindr::cvd_grid()

################################################################################
## Plot : var_BM changes, var_error = 0
################################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

load(paste0(results_directory, "lizards_model_selection.rda")) 

# Merge error and no error for gaussian models
resAll$name <- sub("err", "", resAll$name)
# Names
for (i in seq_along(ex_kurt)) {
  resAll$sim <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), resAll$sim)
  sim_models <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), sim_models)
}
resAll$sim <- sub("Cauchy", "CP", resAll$sim)
sim_models <- sub("Cauchy", "CP", sim_models)
models <- c("BM", "OUfixedRoot", "EB", "cauchy", "cauchyLambda")
models_pretty <- c("BM", "OU", "EB", "CP", "$\\lambda$ CP")

# Subsets Format
resAll <- subset(resAll, var_error == 0 & var_BM > 1)
resAll <- subset(resAll, name != "cauchyLambda")
sel <- subset(resAll, score == "AIC") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)

# colors
# colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), c(-0.8, -0.3, 0.2)),
#              lighten(rep(colorblindr:::palette_OkabeIto[1], 2), c(0.6, 0.8)))
colPlot <- c("#002418", "#026D4E", "#3AB589", "#FFD6AD", "#FFEBD9")

# labeller
# get_pretty_label <- function(string) paste0("$\\sigma^2_{BM} = ", string, "$")
# get_pretty_label <- as_labeller(get_pretty_label, default = label_parsed)
get_pretty_label <- function(string) paste0("sigma2_BM = ", string)
get_pretty_label <- as_labeller(get_pretty_label, default = label_value)

# Plot
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  facet_wrap(vars(var_BM), labeller = get_pretty_label) +
  theme_minimal() + 
  xlab("Simulation Model") + ylab("Count") + 
  scale_fill_manual(values = colPlot, name = "Selected Model:") + 
  theme(#legend.direction="vertical",
    # panel.grid = element_blank(),
    legend.position="top",
    legend.justification = c("left", "top"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.key.size = unit(12, 'pt'),
    text = element_text(size = 9))
# guides(fill = guide_legend(ncol = 4, byrow = TRUE))
colorblindr::cvd_grid()

################################################################################
## Plot : var_BM = 1, var_error changes
################################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

load(paste0(results_directory, "lizards_model_selection.rda")) 

# Merge error and no error for gaussian models
resAll$name <- sub("err", "", resAll$name)
# Names
for (i in seq_along(ex_kurt)) {
  resAll$sim <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), resAll$sim)
  sim_models <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), sim_models)
}
resAll$sim <- sub("Cauchy", "CP", resAll$sim)
sim_models <- sub("Cauchy", "CP", sim_models)
models <- c("BM", "OUfixedRoot", "EB", "cauchy", "cauchyLambda")
models_pretty <- c("BM", "OU", "EB", "CP", "$\\lambda$ CP")

# Subsets Format
resAll <- subset(resAll, var_error > 0 & var_BM == 1)
resAll <- subset(resAll, name != "cauchyLambda")
sel <- subset(resAll, score == "AIC") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)

# colors
# colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), c(-0.8, -0.3, 0.2)),
#              lighten(rep(colorblindr:::palette_OkabeIto[1], 2), c(0.6, 0.8)))
colPlot <- c("#002418", "#026D4E", "#3AB589", "#FFD6AD", "#FFEBD9")

# labeller
# get_pretty_label <- function(string) paste0("$\\sigma^2_{BM} = ", string, "$")
# get_pretty_label <- as_labeller(get_pretty_label, default = label_parsed)
get_pretty_label <- function(string) paste0("sigma2_error = ", string)
get_pretty_label <- as_labeller(get_pretty_label, default = label_value)

# Plot
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  facet_wrap(vars(var_error), labeller = get_pretty_label) +
  theme_minimal() + 
  xlab("Simulation Model") + ylab("Count") + 
  scale_fill_manual(values = colPlot, name = "Selected Model:") + 
  theme(#legend.direction="vertical",
    # panel.grid = element_blank(),
    legend.position="top",
    legend.justification = c("left", "top"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.key.size = unit(12, 'pt'),
    text = element_text(size = 9))
# guides(fill = guide_legend(ncol = 4, byrow = TRUE))
colorblindr::cvd_grid()

################################################################################
## Plot : var_BM = 1, var_error changes, with lambda CP
################################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

load(paste0(results_directory, "lizards_model_selection.rda")) 

# Merge error and no error for gaussian models
resAll$name <- sub("err", "", resAll$name)
# Names
for (i in seq_along(ex_kurt)) {
  resAll$sim <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), resAll$sim)
  sim_models <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), sim_models)
}
resAll$sim <- sub("Cauchy", "CP", resAll$sim)
sim_models <- sub("Cauchy", "CP", sim_models)
models <- c("BM", "OUfixedRoot", "EB", "cauchy", "cauchyLambda")
models_pretty <- c("BM", "OU", "EB", "CP", "$\\lambda$ CP")

# Subsets Format
resAll <- subset(resAll, var_error > 0 & var_BM == 1)
# resAll <- subset(resAll, name != "cauchyLambda")
sel <- subset(resAll, score == "AIC") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)

# colors
# colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), c(-0.8, -0.3, 0.2)),
#              lighten(rep(colorblindr:::palette_OkabeIto[1], 2), c(0.6, 0.8)))
colPlot <- c("#002418", "#026D4E", "#3AB589", "#FFD6AD", "#FFEBD9")

# labeller
# get_pretty_label <- function(string) paste0("$\\sigma^2_{BM} = ", string, "$")
# get_pretty_label <- as_labeller(get_pretty_label, default = label_parsed)
get_pretty_label <- function(string) paste0("sigma2_error = ", string)
get_pretty_label <- as_labeller(get_pretty_label, default = label_value)

# Plot
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  facet_wrap(vars(var_error), labeller = get_pretty_label) +
  theme_minimal() + 
  xlab("Simulation Model") + ylab("Count") + 
  scale_fill_manual(values = colPlot, name = "Selected Model:") + 
  theme(#legend.direction="vertical",
    # panel.grid = element_blank(),
    legend.position="top",
    legend.justification = c("left", "top"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.key.size = unit(12, 'pt'),
    text = element_text(size = 9))
# guides(fill = guide_legend(ncol = 4, byrow = TRUE))
colorblindr::cvd_grid()

################################################################################
## Plot : var_BM = 1, var_error = 0, BIC and AICc
################################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

load(paste0(results_directory, "lizards_model_selection.rda")) 

# Compute BIC
resAllwide <- tidyr::pivot_wider(resAll, names_from = "score", values_from = "value")
resAllwide$logLik <- do.call(c, resAllwide$logLik)
resAllwide$AIC <- do.call(c, resAllwide$AIC)
resAllwide$npars <- resAllwide$logLik + 0.5 * resAllwide$AIC
resAllwide$BIC <- -2 * resAllwide$logLik + log(length(phy$tip.label)) * resAllwide$npars
resAllwide$AICc <- resAllwide$AIC + (2 * resAllwide$npars^2 + 2 * resAllwide$npars) / (length(phy$tip.label) - resAllwide$npars - 1)
resAll <- tidyr::pivot_longer(resAllwide,
                              cols = all_of(c("BIC", "logLik", "AIC", "AICc")),
                              names_to = "score", values_to = "value")
  
# Merge error and no error for gaussian models
resAll$name <- sub("err", "", resAll$name)
# Names
for (i in seq_along(ex_kurt)) {
  resAll$sim <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), resAll$sim)
  sim_models <- sub(paste0("NIG", i), paste0("NIG ", ex_kurt[i]), sim_models)
}
resAll$sim <- sub("Cauchy", "CP", resAll$sim)
sim_models <- sub("Cauchy", "CP", sim_models)
models <- c("BM", "OUfixedRoot", "EB", "cauchy", "cauchyLambda")
models_pretty <- c("BM", "OU", "EB", "CP", "$\\lambda$ CP")

# Subsets Format
resAll <- subset(resAll, var_error == 0 & var_BM == 1)
resAll <- subset(resAll, name != "cauchyLambda")

## AIC
sel <- subset(resAll, score == "AIC") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)
selAIC <- sel
selAIC$score <- "AIC"

## BIC
sel <- subset(resAll, score == "BIC") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)
selBIC <- sel
selBIC$score <- "BIC"

## AICc
sel <- subset(resAll, score == "AICc") %>% group_by(rep, var_BM, var_error, sim) %>% summarise(model = name[which.min(value)])
sel$model <- factor(sel$model, levels = models)
levels(sel$model) <- models_pretty
sel$sim <- factor(sel$sim, levels = sim_models)
selAICc <- sel
selAICc$score <- "AICc"

sel <- rbind(selAIC, selAICc, selBIC)

# colors
# colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), c(-0.8, -0.3, 0.2)),
#              lighten(rep(colorblindr:::palette_OkabeIto[1], 2), c(0.6, 0.8)))
colPlot <- c("#002418", "#026D4E", "#3AB589", "#FFD6AD", "#FFEBD9")

# Plot
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  facet_wrap(vars(score), labeller = label_both) +
  theme_minimal() + 
  xlab("Simulation Model") + ylab("Count") + 
  scale_fill_manual(values = colPlot, name = "Selected Model:") + 
  theme(#legend.direction="vertical",
    # panel.grid = element_blank(),
    legend.position="top",
    legend.justification = c("left", "top"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.key.size = unit(12, 'pt'),
    text = element_text(size = 9))
# guides(fill = guide_legend(ncol = 4, byrow = TRUE))
# colorblindr::cvd_grid()

# colors
# colPlot <- c(lighten(rep(colorblindr:::palette_OkabeIto[3], 3), c(-0.8, -0.3, 0.2)),
#              lighten(rep(colorblindr:::palette_OkabeIto[1], 2), c(0.6, 0.8)))
colPlot <- c("#002418", "#026D4E", "#3AB589", "#FFD6AD", "#FFEBD9")

# Plot
ggplot(sel, aes(x = sim, fill = model)) + geom_bar() +
  theme_minimal() + 
  xlab("Simulation Model") + ylab("Count") + 
  scale_fill_manual(values = colPlot, name = "Selected Model:") + 
  theme(#legend.direction="vertical",
    # panel.grid = element_blank(),
    legend.position="top",
    legend.justification = c("left", "top"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.key.size = unit(12, 'pt'),
    text = element_text(size = 9))
# guides(fill = guide_legend(ncol = 4, byrow = TRUE))
colorblindr::cvd_grid()
