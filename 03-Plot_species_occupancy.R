library(tidyverse)
library(R.utils)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

#__________ Script inputs _____________#
spp <- "WEBL" # Species for which occupancy plot is desired
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

spp.ind <- which(spp.list == spp)

for(v in 1:length(vars.plt)) {
  vlab <- vars.plt[v]
  vind <- which(dimnames(X.beta)[[3]] == vlab)
  x <- seq(quantile(landscape_data[, vlab], prob = 0.01, type = 8),
           quantile(landscape_data[, vlab], prob = 0.99, type = 8),
           length.out = 20)
  z <- (x - mean(landscape_data[, vlab])) / sd(landscape_data[, vlab])
  if(str_detect(vlab, "PACC")) x <- x * 100
  dat.plt <- data.frame(x = x, z = z)
  z.arr <- z %>% array(dim = c(nrow(dat.plt), nsims)) %>% aperm(c(2, 1))
  omega <- mod$sims.list$omega %>% array(dim = c(nsims, nrow(dat.plt)))
  beta0 <- mod$sims.list$beta0[,spp.ind] %>% array(dim = c(nsims, nrow(dat.plt)))
  beta1 <- mod$sims.list$beta1[,spp.ind,vind] %>% array(dim = c(nsims, nrow(dat.plt)))
  psi <- QSLpersonal::expit(
    beta0  + beta1 * z.arr
  )
  dat.plt$psi <- apply(psi, 2, median)
  dat.plt$psi.lo <- apply(psi, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$psi.hi <- apply(psi, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  assign(str_c("dplt.", var.nams[v]), dat.plt)
}