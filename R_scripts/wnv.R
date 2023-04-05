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

# Program to re-analyse the data from:
# PYBUS, Oliver G., SUCHARD, Marc A., LEMEY, Philippe, BERNARDIN, F. J., RAMBAUT, Andrew, CRAWFORD, F. W., GRAY, R. R., ARINAMINPATHY, N., STRAMER, S. L., BUSCH, M. P. y DELWART, E. L.. Unifying the spatial epidemiology and molecular evolution of emerging epidemics. 2012, 15066–15071. (109). ISSN: 0027-8424.ISBN: 1091-6490. 


################################################################################
## MCC tree
################################################################################
# These lines in a Terminal run the xml analysis, and then extract an MCC tree.
# Warning: the beast analysis takes a lot of time (~15 hours)
# Tree and data file resulting from this run are included in the data repository.

# beast -seed 1289 WNV_GTR_gamma.xml
# treeannotator -burninTrees 500 -heights ca WNV_GTR_gamma.trees WNV_GTR_gamma_MCC.tree

################################################################################
## R Setup
################################################################################

library(here)
library(ape)
library(cauphy)
library(phylolm)
library(HDInterval)

################################################################################
## Data from BEAST
################################################################################

## Read MCC tree
phy <- read.nexus(file = here("data/Pybus_et_al_2012_Data/WNV_GTR_gamma_MCC.tree"))
phy <- ladderize(phy, right = FALSE)
write.tree(phy, file = here("data/Pybus_et_al_2012_Data/WNV_GTR_gamma_MCC.newick"))

## Read tip data with noise
# We extract here the data with added jitter that was used in the tree building analysis.
dat <- read.table(file = here("data/Pybus_et_al_2012_Data/WNV_GTR_gamma_tips.log"), header = TRUE)
dat <- do.call(c, dat)
names(dat) <- sub("location.", "", names(dat))
dat <- dat[grepl("[A-Z]", names(dat))]
lat <- dat[grepl(".1$", names(dat))]
names(lat) <- sub(".1$", "", names(lat))
long <- dat[grepl(".2$", names(dat))]
names(long) <- sub(".2$", "", names(long))
lat <- cauphy:::checkTraitTree(lat, phy)
long <- cauphy:::checkTraitTree(long, phy)

write.table(data.frame(traits = names(lat), lat = lat, long = long),
            file = here("data/Pybus_et_al_2012_Data/WNV_GTR_gamma_data.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## plot
plot(phy, cex = 0.5); nodelabels(cex = 0.5)
plot(phy, show.tip.label = FALSE, x.lim = 25)
phydataplot(lat, phy, offset = 1, scaling = 0.1)
phydataplot(long, phy, offset = 15, scaling = 0.05)

## Accession numbers and states
accnum <- unlist(as.vector(read.table(here("data/Pybus_et_al_2012_Data/ordered_accession_numbers.txt"))))
states <- as.vector(read.table(here("data/Pybus_et_al_2012_Data/ordered_states.txt"), sep = "\t"))$V1
match_an <- unname(sapply(accnum, function(an) grep(paste0("^", an), phy$tip.label)))
tips_ordered <- phy$tip.label[match_an]
match_tips <- match(phy$tip.label, tips_ordered)
accnum <- unlist(accnum[match_tips])
states <- states[match_tips]

# cbind(phy$tip.label, accnum)
# cbind(phy$tip.label, states)

# https://www.ncbi.nlm.nih.gov/nuccore/AF202541
# https://www.ncbi.nlm.nih.gov/nuccore/FJ151394

################################################################################
## ML fits
################################################################################

models <- c("BM", "lambda", "delta", "kappa", "OU", "EB", "cauchy")

fit_phylolm <- function(phy, dat, model) {
  if (model == "OU") model <- "OUfixedRoot"
  if (model == "cauchy") {
    tt <- system.time(
      res <- fitCauchy(phy, dat, model = "cauchy", method = "fixed.root", optim = "global")
    )
    return(c(logLik = res$logLik,
             AIC = res$aic,
             mu = res$x0,
             sigma = res$disp,
             param = NA,
             time = unname(tt[3])))
  } else {
    if (model == "kappa") {
      tt <- system.time(
        res <- phylolm(dat ~ 1, phy = phy, model = model, upper.bound = 3)
      )
    } else {
      tt <- system.time(
        res <- phylolm(dat ~ 1, phy = phy, model = model)
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

## Lat
all_fits_lat <- sapply(models, function(mm) fit_phylolm(phy, lat, mm))
all_fits_lat

## Long
all_fits_long <- sapply(models, function(mm) fit_phylolm(phy, long, mm))
all_fits_long


################################################################################
## Comparison of independent beast runs with the CP on a square rooted tree
################################################################################
# Check that the "Cauchy RRW" gives the same results as a Cauchy process on a 
# tree with square-rooted branch lengths.

## beast fits
# beast -seed 1289 WNV_HMC_Cauchy.xml
# beast -seed 1289 WNV_HMC_Cauchy_lat.xml
# beast -seed 1289 WNV_HMC_Cauchy_long.xml

## beast logs
fit_beast <- read.table("data/Pybus_et_al_2012_Data/WNV_HMC_Cauchy.log", header = TRUE)
fit_beast_lat <- read.table("data/Pybus_et_al_2012_Data/WNV_HMC_Cauchy_lat.log", header = TRUE)
fit_beast_long <- read.table("data/Pybus_et_al_2012_Data/WNV_HMC_Cauchy_long.log", header = TRUE)
# 10% burnin
fit_beast <- fit_beast[floor(nrow(fit_beast)/10):nrow(fit_beast), ]
fit_beast_lat <- fit_beast_lat[floor(nrow(fit_beast_lat)/10):nrow(fit_beast_lat), ]
fit_beast_long <- fit_beast_long[floor(nrow(fit_beast_long)/10):nrow(fit_beast_long), ]

## square root tree
physqrt <- phy
physqrt$edge.length <- sqrt(physqrt$edge.length)

## CP fits
fit_Cau_sqrt_lat <- fitCauchy(physqrt, lat, model = "cauchy", method = "random.root", root.edge = 10000, hessian = TRUE)
fit_Cau_sqrt_long <- fitCauchy(physqrt, long, model = "cauchy", method = "random.root", root.edge = 10000, hessian = TRUE)

## Comparisons of parameter estimates
res_disp <- data.frame(CP_lat = c(fit_Cau_sqrt_lat$disp, confint(fit_Cau_sqrt_lat)),
                       RRW_lat = c(median(sqrt(fit_beast_lat$location.variance.diagonal)),
                                   hdi(sqrt(fit_beast_lat$location.variance.diagonal))),
                       CP_long = c(fit_Cau_sqrt_long$disp, confint(fit_Cau_sqrt_long)),
                       RRW_long = c(median(sqrt(fit_beast_long$location.variance.diagonal)),
                                   hdi(sqrt(fit_beast_long$location.variance.diagonal))))
res_disp

################################################################################
## Root Ancestral Reconstruction
################################################################################

## REML Root reconstructions
fit_Cau_lat <- fitCauchy(phy, lat, model = "cauchy", method = "reml", optim = "global")
values_lat <- seq(38, 45, length.out = 1000)
anc_root_lat <- ancestral(fit_Cau_lat, node = length(phy$tip.label) + 1, values = values_lat)
fit_Cau_long <- fitCauchy(phy, long, model = "cauchy", method = "reml", optim = "global")
values_long <- seq(-80, -70, length.out = 1000)
anc_root_long <- ancestral(fit_Cau_long, node = length(phy$tip.label) + 1, values = values_long)

df_CP <- data.frame(root = c(anc_root_lat, anc_root_long),
                    values = c(values_lat, values_long),
                    location = rep(c("latitude", "longitude"), each = length(anc_root_long)),
                    method = "Independent CP")

## BEAST root reconstruction
beast_dens_lat <- density(fit_beast$location.207.1)
beast_dens_long <- density(fit_beast$location.207.2)
df_beast <- rbind(data.frame(root = beast_dens_lat$y,
                             values = beast_dens_lat$x,
                             location = "latitude",
                             method = "Cauchy RRW"),
                  data.frame(root = beast_dens_long$y,
                             values = beast_dens_long$x,
                             location = "longitude",
                             method = "Cauchy RRW"))

df_plot <- rbind(df_CP, df_beast)

## Plot
library(ggplot2)
cols <- hcl.colors(10, palette = "plasma")[c(6,9)]
p <- ggplot(df_plot, aes(x = values, y = root, color = method)) +
  geom_line() + 
  facet_wrap(location ~ ., nrow = 2, scales = "free", strip.position = "bottom") + 
  theme_minimal() + 
  theme(strip.placement = "outside",) + xlab("") + 
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.99),
        legend.justification = c("left", "top"),
        legend.key.size = unit(12, 'pt'),
  ) + 
  scale_color_viridis_d(name = "Method", end = 0.7, option = "B")
p
colorblindr::cvd_grid()

## Peaks
beast_dens_lat$x[pracma::findpeaks(beast_dens_lat$y)[, 2]]
beast_dens_long$x[pracma::findpeaks(beast_dens_long$y)[, 2]]

################################################################################
## Full Ancestral Reconstruction
################################################################################

## values
values_lat <- seq(10, 50, length.out = 100)
values_inc_lat <- seq(-40, 40, length.out = 100)
values_long <- seq(-130, -60, length.out = 100)
values_inc_long <- seq(-70, 70, length.out = 100)

## Reconstruction
time_anc_lat <- system.time(
  anc_lat <- ancestral(fit_Cau_lat, values = values_lat)
)
time_anc_long <- system.time(
  anc_long <- ancestral(fit_Cau_long, values = values_long)
)
time_inc_lat <- system.time(
  inc_lat <- increment(fit_Cau_lat, values = values_inc_lat, n_cores = 6)
)
time_inc_long <- system.time(
  inc_long <- increment(fit_Cau_long, values = values_inc_long, n_cores = 6)
)

## Plot
plot_asr(fit_Cau_lat, anc = anc_lat, inc = inc_lat,
         show.tip.label = FALSE,
         width.node = 0.1, height.node = 2,
         width.edge = 0.2, height.edge = 0.8,
         x.legend = -2, y.legend = 10,
         scaling = 0.02)
axisPhylo(backward = FALSE, root.time = 1998.394)
tiplabels(states, adj = 0, frame = "none")

plot_asr(fit_Cau_long, anc = anc_long, inc = inc_long,
         show.tip.label = FALSE,
         width.node = 0.1, height.node = 2,
         width.edge = 0.2, height.edge = 0.8,
         x.legend = -2, y.legend = 10,
         scaling = 0.01)
axisPhylo(backward = FALSE, root.time = 1998.394)

################################################################################
## Ancestral Reconstruction on Given Nodes
################################################################################

nn <- c(194, 195, 196)

lim_lat <- c(20, 55)
lim_long <- c(-90, -70)
lim_lat_inc <- c(-20, 5)
lim_long_inc <- c(-15, 5)

anc_lat_nn <- ancestral(fit_Cau_lat, values = seq(lim_lat[1], lim_lat[2], length.out = 1000), node = nn)
anc_long_nn <- ancestral(fit_Cau_long, values = seq(lim_long[1], lim_long[2], length.out = 1000), node = nn)
inc_lat_nn <- increment(fit_Cau_lat, values = seq(lim_lat_inc[1], lim_lat_inc[2], length.out = 1000), node = nn, n_cores = 3)
inc_long_nn <- increment(fit_Cau_long, values = seq(lim_long_inc[1], lim_long_inc[2], length.out = 1000), node = nn, n_cores = 3)

plot(inc_lat_nn, node = nn, type = "l")
plot(anc_lat_nn, node = nn, type = "l")
plot(inc_long_nn, node = nn, type = "l")
plot(anc_long_nn, node = nn, type = "l")

## Matching with beast nodes
tree_RRW_0 <- treeio::read.beast(here("data/Pybus_et_al_2012_Data/WNV_HMC_Cauchy.trees"))[[1]]
datRRW_0 <- tidytree::get.data(tree_RRW_0)
fit_beast_0 <- read.table("data/Pybus_et_al_2012_Data/WNV_HMC_Cauchy.log", header = TRUE)

## Reconstructions
all_rec_nn <- NULL
for (i in seq_along(nn)) {
  node <- nn[i]
  par <- phy$edge[phy$edge[, 2] == node, 1]
  # latitude
  node_beast <- names(which(sapply(fit_beast_0, function(x) x[1]) == datRRW_0[datRRW_0[, "node"] == node, ]$location[[1]][1]))
  par_beast <- names(which(sapply(fit_beast_0, function(x) x[1]) == datRRW_0[datRRW_0[, "node"] == par, ]$location[[1]][1]))
  inc_dens <- density(fit_beast[[node_beast]] - fit_beast[[par_beast]])
  anc_dens <- density(fit_beast[[node_beast]])
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = inc_dens$y,
                                             values = inc_dens$x,
                                             location = "latitude",
                                             method = "Cauchy RRW",
                                             type = "increment",
                                             node = node))
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = anc_dens$y,
                                             values = anc_dens$x,
                                             location = "latitude",
                                             method = "Cauchy RRW",
                                             type = "ancestral",
                                             node = node))
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = anc_lat_nn[paste(node), ],
                                             values = as.numeric(colnames(anc_lat_nn)),
                                             location = "latitude",
                                             method = "Independent CP",
                                             type = "ancestral",
                                             node = node))
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = inc_lat_nn[paste(node), ],
                                             values = as.numeric(colnames(inc_lat_nn)),
                                             location = "latitude",
                                             method = "Independent CP",
                                             type = "increment",
                                             node = node))
  # longitude
  node_beast <- names(which(sapply(fit_beast_0, function(x) x[1]) == datRRW_0[datRRW_0[, "node"] == node, ]$location[[1]][2]))
  par_beast <- names(which(sapply(fit_beast_0, function(x) x[1]) == datRRW_0[datRRW_0[, "node"] == par, ]$location[[1]][2]))
  inc_dens <- density(fit_beast[[node_beast]] - fit_beast[[par_beast]])
  anc_dens <- density(fit_beast[[node_beast]])
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = inc_dens$y,
                                             values = inc_dens$x,
                                             location = "longitude",
                                             method = "Cauchy RRW",
                                             type = "increment",
                                             node = node))
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = anc_dens$y,
                                             values = anc_dens$x,
                                             location = "longitude",
                                             method = "Cauchy RRW",
                                             type = "ancestral",
                                             node = node))
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = anc_long_nn[paste(node), ],
                                             values = as.numeric(colnames(anc_long_nn)),
                                             location = "longitude",
                                             method = "Independent CP",
                                             type = "ancestral",
                                             node = node))
  all_rec_nn <- rbind(all_rec_nn, data.frame(density = inc_long_nn[paste(node), ],
                                             values = as.numeric(colnames(inc_long_nn)),
                                             location = "longitude",
                                             method = "Independent CP",
                                             type = "increment",
                                             node = node))
}

## Plot
library(ggplot2)
p <- ggplot(subset(all_rec_nn, location == "latitude"), aes(x = values, y = density, color = method)) +
  geom_line() + 
  facet_wrap(type ~ node, nrow = 2, scales = "free", strip.position = "bottom") +
  # facet_grid(location ~ node, scales = "free", space = "free") + 
  theme_minimal() + 
  theme(strip.placement = "outside",) + xlab("") + 
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.99),
        legend.justification = c("left", "top"),
        legend.key.size = unit(12, 'pt'),
  ) + 
  scale_color_viridis_d(name = "Method", end = 0.7, option = "B")
p
# colorblindr::cvd_grid()

p <- ggplot(subset(all_rec_nn, location == "longitude"), aes(x = values, y = density, color = method)) +
  geom_line() + 
  facet_wrap(type ~ node, nrow = 2, scales = "free", strip.position = "bottom") +
  # facet_grid(location ~ node, scales = "free", space = "free") + 
  theme_minimal() + 
  theme(strip.placement = "outside",) + xlab("") + 
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.99),
        legend.justification = c("left", "top"),
        legend.key.size = unit(12, 'pt'),
  ) + 
  scale_color_viridis_d(name = "Method", end = 0.7, option = "B")
p
