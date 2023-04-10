# cauchy_paper by PB, GD
#   
# To the extent possible under law, the person who associated CC0 with
# cauchy_paper has waived all copyright and related or neighboring rights
# to cauchy_paper.
# 
# You should have received a copy of the CC0 legal code along with this
# work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

# Program to conduct the simulations for the assessment of the parameters inference.

library(here)
library(ape)
library(cauphy)
library(phylolm)
library(foreach)
library(doParallel)
library(HDInterval)

################################################################################
## File management
################################################################################

datestamp <- format(Sys.time(), "%Y-%m-%d")
results_directory <- here("results", paste0(datestamp, "_"))

################################################################################
## Parameters
################################################################################
all_mu <- 0
values_root <- seq(-10, 10, 0.01)
all_disp_cau <- 1
all_lambda <- c(1)
all_n_tips <- c(10, 20, 30, 40, 50, 100, 150, 200, 250)
Nrep <- 500

################################################################################
## Data simulation
################################################################################

sim <- function(Nrep, n_tips, mu, disp, lambda) {
  dat <- NULL
  all_phy <- vector(mode = "list", length = Nrep)
  for (i in 1:Nrep) {
    # sim tree
    phy <- rphylo(n_tips, 0.1, 0)
    phy$edge.length <- phy$edge.length / max(vcv(phy))
    all_phy[[i]] <- phy
    # sim
    dat <- cbind(dat, rTraitCauchy(n = 1, phy = phy, model = "lambda",
                                   parameters = list(root.value = mu, disp = disp, lambda = lambda)))
  }
  return(list(dat = dat,
              phy = all_phy))
}

get_name <- function(Nrep, n_tips, mu, disp, lambda) {
  return(paste("Nrep", Nrep, "n_tips", n_tips, "mu", mu, "disp", disp, "lambda", lambda, sep = "_"))
}

all_sims <- NULL
set.seed(12891026)
for (n_tips in all_n_tips) {
  for (mu in all_mu) {
    for (disp in all_disp_cau) {
      for (lambda in all_lambda) {
        all_sims[[get_name(Nrep, n_tips, mu, disp, lambda)]] <- sim(Nrep, n_tips, mu, disp, lambda)
      }
    }
  }
}

################################################################################
## Estimation function
################################################################################

all_methods <- c("reml", "fixed.root")
all_models <- c("cauchy")
all_init <- c("Qn")

get_params <- function(phy, dat, model, method, init) {
  tt <- system.time(
    ff <- fitCauchy(phy, dat, model = model, method = method, optim = "local", method.init.disp = init)
  )
  ff <- compute_vcov(ff)
  suppressMessages(ci <- confint(ff))
  params <- list(x0 = NA, x0_min = NA, x0_max = NA,
                 disp = ff$disp, disp_min = ci["disp", 1], disp_max = ci["disp", 2],
                 lambda = NA, lambda_min = NA, lambda_max = NA)
  if (method == "fixed.root") {
    params$x0 <- ff$x0
    params$x0_min <- ci["x0", 1]
    params$x0_max <- ci["x0", 2]
  } else {
    anc_root <- ancestral(ff, node = length(phy$tip.label) + 1, values = values_root)
    dens_root <- list(x = as.numeric(colnames(anc_root)),
                      y = anc_root[1, ])
    class(dens_root) <- "density"
    hdi_root <- hdi(dens_root)
    estim_root <- as.numeric(names(which.max(anc_root[1, ])))
    params$x0 <- as.numeric(names(which.max(anc_root[1, ])))
    params$x0_min <- hdi_root[1]
    params$x0_max <- hdi_root[2]
  }
  if (model == "lambda") {
    params$lambda <- ff$lambda
    params$lambda_min <- ci["lambda", 1]
    params$lambda_max <- ci["lambda", 2]
  }
  params$loglik <- ff$logLik
  params$time <- unname(tt[3])
  return(params)
}

fit_params <- function(i, all_sims, Nrep, n_tips, mu, disp, lambda, model, method, init) {
  name_data <- get_name(Nrep, n_tips, mu, disp, lambda)
  all_dat <- all_sims[[name_data]]
  phy <- all_dat$phy[[i]]
  dat <- all_dat$dat[, i]
  names(dat) <- rownames(all_dat$dat)
  res <- get_params(phy, dat, model, method, init)
  res$Nrep <- Nrep
  res$rep <- i
  res$n_tips <- n_tips
  res$true_x0 <- mu
  res$true_disp <- disp
  res$true_lambda <- lambda
  res$model <- model
  res$method <- method
  res$init <- init
  return(res)
}

################################################################################
## Estimation in parallel
################################################################################

## Required packages to pass to all nodes
reqpckg <- c("cauphy", "ape", "HDInterval")
## Number of cores
Ncores <- 6
## Register nodes
cl <- makeCluster(Ncores, outfile = "")
registerDoParallel(cl)

res <- foreach(n_tips = all_n_tips, .packages = reqpckg, .combine = rbind) %:%
  foreach(mu = all_mu, .packages = reqpckg, .combine = rbind) %:%
  foreach(disp = all_disp_cau, .packages = reqpckg, .combine = rbind) %:%
  foreach(lambda = all_lambda, .packages = reqpckg, .combine = rbind) %:%
  foreach(method = all_methods, .packages = reqpckg, .combine = rbind) %:%
  foreach(model = all_models, .packages = reqpckg, .combine = rbind) %:%
  foreach(init = all_init, .packages = reqpckg, .combine = rbind) %:%
  foreach(i = seq_len(Nrep), .packages = reqpckg, .combine = rbind, .errorhandling = "remove") %dopar% {
    fit_params(i, all_sims, Nrep, n_tips, mu, disp, lambda, model, method, init)
  }

stopCluster(cl)

res <- as.data.frame(res)
res <- apply(res, 2, unlist, simplify = FALSE)
res <- as.data.frame(res)

res$x0_in <- (res$true_x0 >= res$x0_min) & (res$true_x0 <= res$x0_max)
res$disp_in <- (res$true_disp >= res$disp_min) & (res$true_disp <= res$disp_max)
res$lambda_in <- (res$true_lambda >= res$lambda_min) & (res$true_lambda <= res$lambda_max)

save.image(file = paste0(results_directory, "simulations_paramters.rda"))


################################################################################
## Plot the results
################################################################################
library(ggplot2)

load(file = paste0(results_directory, "simulations_paramters.rda"))
# load(file = "results/2023-02-03_simulations_paramters.rda")

################################################################################
## sigma estimation

dd <- subset(res, true_lambda == 1 & model == "cauchy")
nrow(dd) / 8 / 2 == Nrep
ggplot(dd, aes(y = disp, x = as.factor(n_tips), fill = method)) +
  # facet_grid(cols = vars(true_lambda), labeller = label_both) +
  geom_boxplot(aes(fill = method)) +
  geom_hline(aes(yintercept = true_disp), linetype = "dashed") +
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal()

coverrage_x0 <- tapply(dd$x0_in, list(dd$n_tips, dd$method), FUN = mean)
coverrage_disp <- tapply(dd$disp_in, list(dd$n_tips, dd$method), FUN = mean)

################################################################################
## root estimation

dd <- subset(res, true_lambda == 1 & model == "cauchy")
nrow(dd) / 2 / 8 == Nrep
ggplot(dd, aes(y = x0, x = as.factor(n_tips), fill = method)) +
  # facet_grid(cols = vars(true_lambda), labeller = label_both) +
  geom_boxplot(aes(fill = method)) +
  geom_hline(aes(yintercept = true_x0), linetype = "dashed") +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_minimal()
