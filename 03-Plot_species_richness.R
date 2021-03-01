library(tidyverse)
library(R.utils)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod"
scripts.loc <- "forest-mapping-scripts/"
mod <- loadObject(mod.nam)
#______________________________________#

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing_community_occupancy.R"))

nsims <- dim(mod$mcmcOutput)[1] # number of posterior samples

# Compile values for plotting grid-level relationships #
vars.plt <- c("PACC10", "PACC40", "mnPerArRatio_Opn")
var.nams <- c("PACCGap", "PACCOpn", "PAROpn")

for(v in 1:length(vars.plt)) {
  vlab <- vars.plt[v]
  vind <- which(dimnames(X.beta)[[3]] == vlab)
  x <- seq(quantile(landscape_data[, vlab], prob = 0.01, type = 8),
           quantile(landscape_data[, vlab], prob = 0.99, type = 8),
           length.out = 20)
  z <- (x - mean(landscape_data[, vlab])) / sd(landscape_data[, vlab])
  if(str_detect(vlab, "PACC")) x <- x * 100
  dat.plt <- data.frame(x = x, z = z)
  z.arr <- z %>% array(dim = c(nrow(dat.plt), nsims, length(spp.list))) %>% aperm(c(2, 3, 1))
  omega <- mod$sims.list$omega %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  beta0 <- mod$sims.list$beta0 %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  beta1 <- mod$sims.list$beta1[,,vind] %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  psi <- QSLpersonal::expit(
    beta0  + beta1 * z.arr
  )
  y <- omega * psi
  SR <- apply(y, c(1, 3), sum)
  dat.plt$SR.pred <- apply(SR, 2, median)
  dat.plt$SR.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$SR.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))
  assign(str_c("dplt.", var.nams[v]), dat.plt)
}

# Compile values for plotting grid-level relationships #
x <- seq(quantile(Cov[, "CanCov"], prob = 0.01, type = 8),
         quantile(Cov[, "CanCov"], prob = 0.99, type = 8),
         length.out = 20)
z <- (x - mean(Cov[, "CanCov"])) / sd(Cov[, "CanCov"])
dat.plt <- data.frame(x = x, z = z)
z.arr <- z %>% array(dim = c(nrow(dat.plt), nsims, length(spp.list))) %>% aperm(c(2, 3, 1))
omega <- mod$sims.list$omega %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
psi <- QSLpersonal::expit(mod$sims.list$beta0) %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
alpha0 <- mod$sims.list$alpha0 %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
alpha1 <- mod$sims.list$alpha1[,,1] %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
alpha2 <- mod$sims.list$alpha1[,,2] %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
theta <- QSLpersonal::expit(
  alpha0  + alpha1 * z.arr + alpha2 * (z.arr^2)
)
y <- omega * psi * theta 
SR <- apply(y, c(1, 3), sum)
dat.plt$SR.pred <- apply(SR, 2, median)
dat.plt$SR.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
dat.plt$SR.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))
dplt.CanCov <- dat.plt

ymin <- min(dplt.PACCGap$SR.plo, dplt.PACCOpn$SR.plo, dplt.PAROpn$SR.plo)
ymax <- max(dplt.PACCGap$SR.phi, dplt.PACCOpn$SR.phi, dplt.PAROpn$SR.phi)

dat.plt <- dplt.PACCGap
p.PACCGap <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi), alpha = 0.3) +
  geom_line(aes(x = x, y = SR.pred), size = 1.5) +
  ylim(ymin, ymax) +
  labs(x = "Canopy gaps (%)", y = NULL)

dat.plt <- dplt.PACCOpn
p.PACCOpn <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi), alpha = 0.3) +
  geom_line(aes(x = x, y = SR.pred), size = 1.5) +
  ylim(ymin, ymax) +
  labs(x = "Open forest (%)", y = NULL)

dat.plt <- dplt.PAROpn
p.PAROpn <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi), alpha = 0.3) +
  geom_line(aes(x = x, y = SR.pred), size = 1.5) +
  ylim(ymin, ymax) +
  labs(x = str_wrap("Open forest perimeter-area ratio", 30), y = NULL)

dat.plt <- dplt.CanCov
p.CanCov <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi), alpha = 0.3) +
  geom_line(aes(x = x, y = SR.pred), size = 1.5) +
  labs(x = "Canopy cover (%)", y = NULL)

p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0.0500000, y = 0.525, width = 0.3166667, height = 0.475) +
  draw_plot(p.PACCOpn, x = 0.3666667, y = 0.525, width = 0.3166667, height = 0.475) +
  draw_plot(p.PAROpn,  x = 0.6833333, y = 0.525, width = 0.3166667, height = 0.475) +
  draw_plot(p.CanCov,  x = 0.0500000, y = 0.050, width = 1,         height = 0.475) +
  draw_plot_label(c("Point richness", "Grid cell richness"),
                  x = c(0, 0), y = c(0.3, 0.8), size = 15,
                  angle = c(90, 90), hjust = c(0.5, 0.5))

save_plot(str_c("Plot_species_richness.jpg"), p, ncol = 2, nrow = 2, dpi = 200)
