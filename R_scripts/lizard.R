# cauchy_paper by PB, GD
#   
# To the extent possible under law, the person who associated CC0 with
# cauchy_paper has waived all copyright and related or neighboring rights
# to cauchy_paper.
# 
# You should have received a copy of the CC0 legal code along with this
# work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

# Program to re-analyse the data from:
# Mahler, D. Luke; Ingram, Travis; Revell, Liam J.; Losos, Jonathan B. (2013), Data from: Exceptional convergence on the macroevolutionary landscape in island lizard radiations, Dryad, Dataset, https://doi.org/10.5061/dryad.9g182

library(here)
library(ape)
library(phylolm)
library(Rphylopars)
library(cauphy)

################################################################################
## Data
################################################################################
## Phylogenetic Tree
phy <- read.tree(file = here("data", "Mahler_et_al_2013_Data", "GA_Anolis_MCC.tre"))
## Data
dat <- read.csv(file = here("data", "Mahler_et_al_2013_Data", "GA_Anolis_traits.csv"))
## ecomorphs
ecomorphs <- read.csv(here("data", "Mahler_et_al_2013_Data", "GA_Anolis_trad_ecomorph_class.csv"),
                      header = FALSE)
ecomorphs <- ecomorphs[match(phy$tip.label, ecomorphs[, 1]), ]
ecomorphs <- factor(ecomorphs[, 2])
levels(ecomorphs) <- c("Trunk-Ground",
                       "Trunk-Crown",
                       "Crown-Giant",
                       "Twig",
                       "Grass-Bush",
                       "Trunk",
                       "Unique")
ecomorphs <- factor(ecomorphs, level = sort(levels(ecomorphs)))
col_ecomorphs <- ecomorphs
levels(col_ecomorphs) <- hcl.colors(7, palette = "Dark 3")

## Keep only svl
svl <- dat[, "AVG.SVL"]
names(svl) <- dat$species
svl <- cauphy:::checkTraitTree(svl, phy)

## clades
cham_giant <- getMRCA(phy, 49:58)
cham_giant_nodes <- c(cham_giant, phytools::getDescendants(phy, cham_giant))
cham_giant_edges <- phy$edge[, 2] %in% cham_giant_nodes
cham_giant_edges <- (cham_giant_edges + 1) * 1.5

## Models
models <- c("BM", "lambda", "delta", "kappa", "OU", "EB", "cauchy", "cauchyLambda")
models_pretty <- c("BM", "$\\lambda$", "$\\delta$", "$\\kappa$", "OU", "EB", "CP", "Cau $\\lambda$")
rows_pretty <- c("logLik", "AIC", "$\\mu$", "$\\sigma$ ($\\times 10^{-3}$)", "$\\theta$", "time (s)")

################################################################################
## ML
################################################################################

fit_phylolm <- function(phy, svl, model) {
  if (model == "OU") model <- "OUfixedRoot"
  if (model == "cauchy") {
    tt <- system.time(
      res <- fitCauchy(phy, svl, model = "cauchy", method = "fixed.root", optim = "global")
    )
    return(c(logLik = res$logLik,
             AIC = res$aic,
             mu = res$x0,
             sigma = res$disp,
             param = NA,
             time = unname(tt[3])))
  } else if (model == "cauchyLambda") {
    tt <- system.time(
      res <- fitCauchy(phy, svl, model = "lambda", method = "fixed.root", optim = "global")
    )
    return(c(logLik = res$logLik,
             AIC = res$aic,
             mu = res$x0,
             sigma = res$disp,
             param = res$lambda,
             time = unname(tt[3])))
  } else {
    if (model == "kappa") {
      tt <- system.time(
        res <- phylolm(svl ~ 1, phy = phy, model = model, upper.bound = 3)
      )
    } else {
      tt <- system.time(
        res <- phylolm(svl ~ 1, phy = phy, model = model)
      )
    }
    return(c(logLik = res$logLik,
             AIC = res$aic,
             mu = unname(res$coefficients),
             sigma = sqrt(res$sigma2),
             param = ifelse(model == "BM", NA, res$optpar),
             time = unname(tt[3])))
  }
}

all_fits <- sapply(models, function(mm) fit_phylolm(phy, svl, mm))
all_fits

################################################################################
## REML
################################################################################

fit_reml <- function(phy, svl, model) {
  if (model == "cauchy") {
    tt <- system.time(
      res <- fitCauchy(phy, svl, model = "cauchy", method = "reml", optim = "global")
    )
    return(c(logLik = res$logLik,
             AIC = res$aic,
             mu = NA,
             sigma = res$disp,
             param = NA,
             time = unname(tt[3])))
  } else if (model == "cauchyLambda") {
    tt <- system.time(
      res <- fitCauchy(phy, svl, model = "lambda", method = "reml", optim = "global")
    )
    return(c(logLik = res$logLik,
             AIC = res$aic,
             mu = NA,
             sigma = res$disp,
             param = res$lambda,
             time = unname(tt[3])))
  } else {
    svl_pp <- data.frame(species = names(svl), svl = svl)
    if (model == "delta") {
      phy_norm <- phy
      h_tree <-  vcv(phy)[1, 1]
      phy_norm$edge.length <- phy_norm$edge.length / h_tree
      tt <- system.time(
        res <- phylopars(svl_pp, phy_norm, model = model, REML = TRUE)
      )
      res$pars$phylocov <- res$pars$phylocov / h_tree
    } else {
      tt <- system.time(
        res <- phylopars(svl_pp, phy, model = model, REML = TRUE)
      )
    }
    return(c(logLik = res$logLik,
             AIC = -2 * res$logLik + 2 * res$npars,
             mu = NA,
             sigma = sqrt(ifelse(model == "OU", 2 * res$model$alpha * res$pars$phylocov, res$pars$phylocov)),
             param = switch(model,
                            BM = NA,
                            OU = res$model$alpha,
                            lambda = res$model$lambda,
                            delta = res$model$delta,
                            kappa = res$model$kappa,
                            EB = res$model$rate),
             time = unname(tt[3])))
  }
}

all_fits_reml <- sapply(models, function(mm) fit_reml(phy, svl, mm))
all_fits_reml

################################################################################
## Ancestral Reconstruction
################################################################################
## fit
fit_Cau_reml <- fitCauchy(phy, svl, model = "cauchy", method = "reml", optim = "global")
## values
values <- seq(3, 5.5, length.out = 100)
values_inc <- seq(-1.5, 1.5, length.out = 100)
## edges long enough
min_length <- 3
nodes_inc <- phy$edge[phy$edge.length > min_length, 2]
## Reconstruction
time_anc <- system.time(
  anc_reml <- ancestral(fit_Cau_reml, values = values)
)
time_inc <- system.time(
  inc_reml <- increment(fit_Cau_reml, values = values_inc, node = nodes_inc, n_cores = 6)
)

## Plot
par(cex = 1, mar = c(0, 0, 0, 0), bg = "white")
plot_asr(fit_Cau_reml, anc = anc_reml, inc = inc_reml,
         show.tip.label = FALSE,
         width.node = 0.8, height.node = 1.8,
         width.edge = 1.5, height.edge = 0.8,
         x.lim = 65,
         x.legend = -2, y.legend = 10,
         edge.width = cham_giant_edges)
tiplabels(text = paste0("A. ", phy$tip.label), frame = "none",
          adj = c(0, 0.5), offset = 7,
          col = as.vector(col_ecomorphs), font = 4, cex = 0.85)
legend(x = -2, y = 103.5,
       legend = levels(ecomorphs),
       col = levels(col_ecomorphs),
       pch = 19, y.intersp = 0.8, cex = 0.9, ncol = 2)

