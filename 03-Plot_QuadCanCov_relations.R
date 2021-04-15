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
x.var <- "CanCov"
#______________________________________#

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing_community_occupancy.R"))

spp.ind <- which(apply(mod$sims.list$alpha1[, , 2], 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                   apply(mod$sims.list$alpha1[, , 2], 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plt <- spp.list[spp.ind]
spp.common.map <- spp.out$common_name
names(spp.common.map) <- spp.out$BirdCode
spp.common <- spp.common.map[spp.plt]

# Tabulate species estimates for plotting #
x <- seq(
  quantile(Cov[, x.var], probs = 0.01, type = 8),
  quantile(Cov[, x.var], probs = 0.99, type = 8),
  length.out = 20
)
z <- (x - mean(Cov[, x.var])) / sd(Cov[, x.var])
dat.template <- data.frame(x = x, z = z)
nsims <- dim(mod$mcmcOutput)[1]
z.arr <- dat.template$z %>% array(dim = c(nrow(dat.template), nsims)) %>% aperm(c(2, 1))

for(sp in 1:length(spp.plt)) {
  sp.plt <- spp.plt[sp]
  sp.ind <- spp.ind[sp]
  sp.com <- spp.common[sp]
  dat.sp <- dat.template
  psi <- QSLpersonal::expit(mod$sims.list$beta0[, sp.ind]) %>% array(dim = c(nsims, nrow(dat.sp)))
  alpha0 <- mod$sims.list$alpha0[, sp.ind] %>% array(dim = c(nsims, nrow(dat.sp)))
  alpha1 <- mod$sims.list$alpha1[, sp.ind, 1] %>% array(dim = c(nsims, nrow(dat.sp)))
  alpha2 <- mod$sims.list$alpha1[, sp.ind, 2] %>% array(dim = c(nsims, nrow(dat.sp)))
  theta <- QSLpersonal::expit(alpha0 + alpha1 * z.arr + alpha2 * (z.arr^2))
  dat.sp$y <- apply(psi * theta, 2, median)
  dat.sp$y.lo <- apply(psi * theta, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.sp$y.hi <- apply(psi * theta, 2, function(x) quantile(x, prob = 0.975, type = 8))
  assign(str_c("dat.", sp.plt), dat.sp)

  p <- ggplot(dat = dat.sp, aes(x = x, y = y)) +
    geom_ribbon(aes(ymin = y.lo, ymax = y.hi), alpha=0.3) +
    geom_line(size = 1.5) + 
    scale_y_continuous(lim = c(0, 1)) +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text.x=element_text(size=20)) +
    theme(axis.text.y=element_text(size=20)) +
    annotate(geom = "text", x = mean(x), y = 1, label = str_wrap(sp.com, 17), size = 6, vjust = 1)
  assign(str_c("p.", sp.plt), p)
}

p <- ggdraw() + 
  draw_plot(p.BTHU, x = 0,    y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(p.WEWP, x = 0.25, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(p.HAFL, x = 0.5,  y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(p.PLVI, x = 0.75, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(p.WAVI, x = 0,    y = 0.5,  width = 0.25, height = 0.25) +
  draw_plot(p.CLNU, x = 0.25, y = 0.5,  width = 0.25, height = 0.25) +
  draw_plot(p.VGSW, x = 0.5,  y = 0.5,  width = 0.25, height = 0.25) +
  draw_plot(p.PYNU, x = 0.75, y = 0.5,  width = 0.25, height = 0.25) +
  draw_plot(p.HOWR, x = 0,    y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(p.WEBL, x = 0.25, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(p.TOSO, x = 0.5,  y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(p.EVGR, x = 0.75, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(p.RECR, x = 0,    y = 0,    width = 0.25, height = 0.25) +
  draw_plot(p.CHSP, x = 0.25, y = 0,    width = 0.25, height = 0.25) +
  draw_plot(p.GTTO, x = 0.5,  y = 0,    width = 0.25, height = 0.25) +
  draw_plot(p.WETA, x = 0.75, y = 0,    width = 0.25, height = 0.25)
p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0.05, width = 0.95, height = 0.95) +
  draw_plot_label(c("Point occupancy", "Canopy cover (%)"),
                  x = c(0, 0.55), y = c(0.55, 0.05), size = c(28, 28),
                  angle = c(90, 0), hjust = c(0.5, 0.5))

save_plot(str_c("Plot_spp_", x.var, "_relations.jpg"), p, ncol = 3, nrow = 3, dpi = 200)
