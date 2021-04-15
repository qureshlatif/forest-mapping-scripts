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
vars.plt <- c("PACC10", "PACC40")
var.nams <- c("PACCGap", "PACCOpn")

spp_cat <- read.csv("Spp_list_detected_&_categorized.csv", header = T, stringsAsFactors = F)
PIPO_assoc_ind <- which(spp.list %in% spp_cat$BirdCode[which(spp_cat$PIPO_associated)])
PIPO_spec_ind <- which(spp.list %in% spp_cat$BirdCode[which(spp_cat$PIPO_specialist)])
Riffraff_ind <- which(spp.list %in% spp_cat$BirdCode[which(!spp_cat$PIPO_associated)])

for(v in 1:length(vars.plt)) {
  vlab <- vars.plt[v]
  vind <- which(dimnames(X.beta)[[3]] == vlab)
  x <- rep(seq(quantile(landscape_data[, vlab], prob = 0.01, type = 8),
           quantile(landscape_data[, vlab], prob = 0.99, type = 8),
           length.out = 20), 2)
  z <- (x - mean(landscape_data[, vlab])) / sd(landscape_data[, vlab])
  x.LowMont <- c(rep(0, 20), rep(1, 20))
  z.LowMont <- (x.LowMont - mean(landscape_data$LowMont)) / sd(landscape_data$LowMont)
  
  if(str_detect(vlab, "PACC")) x <- x * 100
  dat.plt <- data.frame(x = x, z = z, Zone = c(rep("Upper", 20), rep("Lower", 20))) %>%
    mutate(Zone = factor(Zone, levels = c("Lower", "Upper")))
  z.arr <- z %>% array(dim = c(nrow(dat.plt), nsims, length(spp.list))) %>% aperm(c(2, 3, 1))
  zLM.arr <- z.LowMont %>% array(dim = c(nrow(dat.plt), nsims, length(spp.list))) %>% aperm(c(2, 3, 1))
  omega <- mod$sims.list$omega %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  beta0 <- mod$sims.list$beta0 %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  beta1 <- mod$sims.list$beta1[,,vind] %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  beta2 <- mod$sims.list$beta1[,,7] %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
  psi <- QSLpersonal::expit(
    beta0  + beta1 * z.arr + beta2 * zLM.arr
  )
  y <- omega * psi
  SR <- apply(y, c(1, 3), sum)
  dat.plt$SR.pred <- apply(SR, 2, median)
  dat.plt$SR.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$SR.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  y <- omega[, PIPO_assoc_ind, ] * psi[, PIPO_assoc_ind, ]
  SR <- apply(y, c(1, 3), sum)
  dat.plt$SR_assoc.pred <- apply(SR, 2, median)
  dat.plt$SR_assoc.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$SR_assoc.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))

  y <- omega[, PIPO_spec_ind, ] * psi[, PIPO_spec_ind, ]
  SR <- apply(y, c(1, 3), sum)
  dat.plt$SR_spec.pred <- apply(SR, 2, median)
  dat.plt$SR_spec.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$SR_spec.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  y_ass <- omega[, PIPO_assoc_ind, ] * psi[, PIPO_assoc_ind, ]
  y_riffraff <- omega[, Riffraff_ind, ] * psi[, Riffraff_ind, ]
  SR <- apply(y_ass, c(1, 3), sum) / apply(y_riffraff, c(1, 3), sum)
  dat.plt$Assoc_ratio.pred <- apply(SR, 2, median)
  dat.plt$Assoc_ratio.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$Assoc_ratio.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  y_spec <- omega[, PIPO_spec_ind, ] * psi[, PIPO_spec_ind, ]
  y_riffraff <- omega[, Riffraff_ind, ] * psi[, Riffraff_ind, ]
  SR <- apply(y_spec, c(1, 3), sum) / apply(y_riffraff, c(1, 3), sum)
  dat.plt$Special_ratio.pred <- apply(SR, 2, median)
  dat.plt$Special_ratio.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
  dat.plt$Special_ratio.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  assign(str_c("dplt.", var.nams[v]), dat.plt)
}

# Compile values for plotting point-level relationships #
x <- rep(seq(quantile(Cov[, "CanCov"], prob = 0.01, type = 8),
         quantile(Cov[, "CanCov"], prob = 0.99, type = 8),
         length.out = 20), 2)
z <- (x - mean(Cov[, "CanCov"])) / sd(Cov[, "CanCov"])
x.LowMont <- c(rep(0, 20), rep(1, 20))
z.LowMont <- (x.LowMont - mean(landscape_data$LowMont)) / sd(landscape_data$LowMont)

dat.plt <- data.frame(x = x, z = z, Zone = c(rep("Upper", 20), rep("Lower", 20))) %>%
  mutate(Zone = factor(Zone, levels = c("Lower", "Upper")))
z.arr <- z %>% array(dim = c(nrow(dat.plt), nsims, length(spp.list))) %>% aperm(c(2, 3, 1))
zLM.arr <- z.LowMont %>% array(dim = c(nrow(dat.plt), nsims, length(spp.list))) %>% aperm(c(2, 3, 1))
omega <- mod$sims.list$omega %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
beta0 <- mod$sims.list$beta0 %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
betaLM <- mod$sims.list$beta1[,,7] %>% array(dim = c(nsims, length(spp.list), nrow(dat.plt)))
psi <- QSLpersonal::expit(beta0 + betaLM * zLM.arr)
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

y <- omega[, PIPO_assoc_ind, ] * psi[, PIPO_assoc_ind, ] * theta[, PIPO_assoc_ind, ]
SR <- apply(y, c(1, 3), sum)
dat.plt$SR_assoc.pred <- apply(SR, 2, median)
dat.plt$SR_assoc.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
dat.plt$SR_assoc.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))

y <- omega[, PIPO_spec_ind, ] * psi[, PIPO_spec_ind, ] * theta[, PIPO_spec_ind, ]
SR <- apply(y, c(1, 3), sum)
dat.plt$SR_spec.pred <- apply(SR, 2, median)
dat.plt$SR_spec.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
dat.plt$SR_spec.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))

y_ass <- omega[, PIPO_assoc_ind, ] * psi[, PIPO_assoc_ind, ] * theta[, PIPO_assoc_ind, ]
y_riffraff <- omega[, Riffraff_ind, ] * psi[, Riffraff_ind, ] * theta[, Riffraff_ind, ]
SR <- apply(y_ass, c(1, 3), sum) / apply(y_riffraff, c(1, 3), sum)
dat.plt$Assoc_ratio.pred <- apply(SR, 2, median)
dat.plt$Assoc_ratio.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
dat.plt$Assoc_ratio.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))

y_spec <- omega[, PIPO_spec_ind, ] * psi[, PIPO_spec_ind, ] * theta[, PIPO_spec_ind, ]
y_riffraff <- omega[, Riffraff_ind, ] * psi[, Riffraff_ind, ] * theta[, Riffraff_ind, ]
SR <- apply(y_spec, c(1, 3), sum) / apply(y_riffraff, c(1, 3), sum)
dat.plt$Special_ratio.pred <- apply(SR, 2, median)
dat.plt$Special_ratio.plo <- apply(SR, 2, function(x) quantile(x, prob = 0.025, type = 8))
dat.plt$Special_ratio.phi <- apply(SR, 2, function(x) quantile(x, prob = 0.975, type = 8))

dplt.CanCov <- dat.plt

###################################
## Plot overall species richness ##
###################################

ymin <- min(dplt.PACCGap$SR.plo, dplt.PACCOpn$SR.plo)
ymax <- max(dplt.PACCGap$SR.phi, dplt.PACCOpn$SR.phi)

dat.plt <- dplt.PACCGap
p.PACCGap <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi, fill = Zone), alpha = 0.4) +
  geom_line(aes(x = x, y = SR.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  ylim(ymin, ymax) +
  labs(x = "Canopy gaps (%)", y = NULL) +
  guides(color = FALSE, fill = FALSE)

dat.plt <- dplt.PACCOpn
p.PACCOpn <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi, fill = Zone), alpha = 0.3) +
  geom_line(aes(x = x, y = SR.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  ylim(ymin, ymax) +
  labs(x = "Open forest (%)", y = NULL, color = "Life Zone", fill = "Life Zone")

dat.plt <- dplt.CanCov
p.CanCov <- ggplot(data = dat.plt, aes(x = x, y = SR.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR.plo, ymax = SR.phi, fill = Zone), alpha = 0.3) +
  geom_line(aes(x = x, y = SR.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  labs(x = "Canopy cover (%)", y = NULL) +
  guides(color = FALSE, fill = FALSE)

p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0,    y = 0.5, width = 0.47, height = 0.5) +
  draw_plot(p.PACCOpn, x = 0.47, y = 0.5, width = 0.53, height = 0.5) +
  draw_plot(p.CanCov,  x = 0,    y = 0,   width = 1,    height = 0.5)
p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0, width = 0.95, height = 0.95) +
  draw_plot_label(c("Point richness", "Grid cell richness", "All species"),
                  x = c(0, 0, 0.55), y = c(0.3, 0.75, 1), size = 15,
                  angle = c(90, 90, 0), hjust = c(0.5, 0.5, 0.5))

save_plot(str_c("Plot_species_richness.jpg"), p, ncol = 2, nrow = 2, dpi = 200)


###################################
## Plot PIPO specialist richness ##
###################################

ymin <- min(dplt.PACCGap$SR_spec.plo, dplt.PACCOpn$SR_spec.plo)
ymax <- max(dplt.PACCGap$SR_spec.phi, dplt.PACCOpn$SR_spec.phi)

dat.plt <- dplt.PACCGap
p.PACCGap <- ggplot(data = dat.plt, aes(x = x, y = SR_spec.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR_spec.plo, ymax = SR_spec.phi, fill = Zone), alpha = 0.4) +
  geom_line(aes(x = x, y = SR_spec.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  ylim(ymin, ymax) +
  labs(x = "Canopy gaps (%)", y = NULL) +
  guides(color = FALSE, fill = FALSE)

dat.plt <- dplt.PACCOpn
p.PACCOpn <- ggplot(data = dat.plt, aes(x = x, y = SR_spec.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR_spec.plo, ymax = SR_spec.phi, fill = Zone), alpha = 0.3) +
  geom_line(aes(x = x, y = SR_spec.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  ylim(ymin, ymax) +
  labs(x = "Open forest (%)", y = NULL, color = "Life Zone", fill = "Life Zone")

dat.plt <- dplt.CanCov
p.CanCov <- ggplot(data = dat.plt, aes(x = x, y = SR_spec.pred)) + 
  geom_ribbon(aes(x = x, ymin = SR_spec.plo, ymax = SR_spec.phi, fill = Zone), alpha = 0.3) +
  geom_line(aes(x = x, y = SR_spec.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  labs(x = "Canopy cover (%)", y = NULL) +
  guides(color = FALSE, fill = FALSE)

p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0,    y = 0.5, width = 0.47, height = 0.5) +
  draw_plot(p.PACCOpn, x = 0.47, y = 0.5, width = 0.53, height = 0.5) +
  draw_plot(p.CanCov,  x = 0,    y = 0,   width = 1,    height = 0.5)
p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0, width = 0.95, height = 0.95) +
  draw_plot_label(c("Point richness", "Grid cell richness", "Ponderosa pine forest specialists"),
                  x = c(0, 0, 0.55), y = c(0.3, 0.75, 1), size = 15,
                  angle = c(90, 90, 0), hjust = c(0.5, 0.5, 0.5))

save_plot(str_c("Plot_PIPO_specialist_richness.jpg"), p, ncol = 2, nrow = 2, dpi = 200)

####################################
## Plot specialist/riffraff ratio ##
####################################

ymin <- min(dplt.PACCGap$Special_ratio.plo, dplt.PACCOpn$Special_ratio.plo)
ymax <- max(dplt.PACCGap$Special_ratio.phi, dplt.PACCOpn$Special_ratio.phi)

dat.plt <- dplt.PACCGap
p.PACCGap <- ggplot(data = dat.plt, aes(x = x, y = Special_ratio.pred)) + 
  geom_ribbon(aes(x = x, ymin = Special_ratio.plo, ymax = Special_ratio.phi, fill = Zone), alpha = 0.4) +
  geom_line(aes(x = x, y = Special_ratio.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  ylim(ymin, ymax) +
  labs(x = "Canopy gaps (%)", y = NULL) +
  guides(color = FALSE, fill = FALSE)

dat.plt <- dplt.PACCOpn
p.PACCOpn <- ggplot(data = dat.plt, aes(x = x, y = Special_ratio.pred)) + 
  geom_ribbon(aes(x = x, ymin = Special_ratio.plo, ymax = Special_ratio.phi, fill = Zone), alpha = 0.3) +
  geom_line(aes(x = x, y = Special_ratio.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  ylim(ymin, ymax) +
  labs(x = "Open forest (%)", y = NULL, color = "Life Zone", fill = "Life Zone")

dat.plt <- dplt.CanCov
p.CanCov <- ggplot(data = dat.plt, aes(x = x, y = Special_ratio.pred)) + 
  geom_ribbon(aes(x = x, ymin = Special_ratio.plo, ymax = Special_ratio.phi, fill = Zone), alpha = 0.3) +
  geom_line(aes(x = x, y = Special_ratio.pred, color = Zone), size = 1.5) +
  scale_color_manual(values = c("saddlebrown", "#009E73")) +
  scale_fill_manual(values = c("saddlebrown", "#009E73")) +
  labs(x = "Canopy cover (%)", y = NULL) +
  guides(color = FALSE, fill = FALSE)

p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0,    y = 0.5, width = 0.47, height = 0.5) +
  draw_plot(p.PACCOpn, x = 0.47, y = 0.5, width = 0.53, height = 0.5) +
  draw_plot(p.CanCov,  x = 0,    y = 0,   width = 1,    height = 0.5)
p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0, width = 0.95, height = 0.95) +
  draw_plot_label(c("Point ratio", "Grid cell ratio", "PIPO specialist ratio"),
                  x = c(0, 0, 0.55), y = c(0.3, 0.75, 1), size = 15,
                  angle = c(90, 90, 0), hjust = c(0.5, 0.5, 0.5))

save_plot(str_c("Plot_PIPO_specialist_ratio.jpg"), p, ncol = 2, nrow = 2, dpi = 200)

# ##################################
# ## Plot PIPO associate richness ##
# ##################################
# 
# ymin <- min(dplt.PACCGap$SR_assoc.plo, dplt.PACCOpn$SR_assoc.plo, dplt.PAROpn$SR_assoc.plo)
# ymax <- max(dplt.PACCGap$SR_assoc.phi, dplt.PACCOpn$SR_assoc.phi, dplt.PAROpn$SR_assoc.phi)
# 
# dat.plt <- dplt.PACCGap
# p.PACCGap <- ggplot(data = dat.plt, aes(x = x, y = SR_assoc.pred)) + 
#   geom_ribbon(aes(x = x, ymin = SR_assoc.plo, ymax = SR_assoc.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = SR_assoc.pred), size = 1.5) +
#   ylim(ymin, ymax) +
#   labs(x = "Canopy gaps (%)", y = NULL)
# 
# dat.plt <- dplt.PACCOpn
# p.PACCOpn <- ggplot(data = dat.plt, aes(x = x, y = SR_assoc.pred)) + 
#   geom_ribbon(aes(x = x, ymin = SR_assoc.plo, ymax = SR_assoc.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = SR_assoc.pred), size = 1.5) +
#   ylim(ymin, ymax) +
#   labs(x = "Open forest (%)", y = NULL)
# 
# dat.plt <- dplt.PAROpn
# p.PAROpn <- ggplot(data = dat.plt, aes(x = x, y = SR_assoc.pred)) + 
#   geom_ribbon(aes(x = x, ymin = SR_assoc.plo, ymax = SR_assoc.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = SR_assoc.pred), size = 1.5) +
#   ylim(ymin, ymax) +
#   labs(x = str_wrap("Open forest perimeter-area ratio", 30), y = NULL)
# 
# dat.plt <- dplt.CanCov
# p.CanCov <- ggplot(data = dat.plt, aes(x = x, y = SR_assoc.pred)) + 
#   geom_ribbon(aes(x = x, ymin = SR_assoc.plo, ymax = SR_assoc.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = SR_assoc.pred), size = 1.5) +
#   labs(x = "Canopy cover (%)", y = NULL)
# 
# p <- ggdraw() + 
#   draw_plot(p.PACCGap, x = 0,      y = 0.5, width = 0.3333, height = 0.5) +
#   draw_plot(p.PACCOpn, x = 0.3333, y = 0.5, width = 0.3333, height = 0.5) +
#   draw_plot(p.PAROpn,  x = 0.6667, y = 0.5, width = 0.3333, height = 0.5) +
#   draw_plot(p.CanCov,  x = 0,      y = 0,   width = 1,      height = 0.5)
# p <- ggdraw() +
#   draw_plot(p, x = 0.05, y = 0, width = 0.95, height = 0.95) +
#   draw_plot_label(c("Point richness", "Grid cell richness", "Ponderosa pine forest associates"),
#                   x = c(0, 0, 0.55), y = c(0.3, 0.75, 1), size = 15,
#                   angle = c(90, 90, 0), hjust = c(0.5, 0.5, 0.5))
# 
# save_plot(str_c("Plot_PIPO_associate_richness.jpg"), p, ncol = 2, nrow = 2, dpi = 200)

# ###################################
# ## Plot associate/riffraff ratio ##
# ###################################
# 
# ymin <- min(dplt.PACCGap$Assoc_ratio.plo, dplt.PACCOpn$Assoc_ratio.plo, dplt.PAROpn$Assoc_ratio.plo)
# ymax <- max(dplt.PACCGap$Assoc_ratio.phi, dplt.PACCOpn$Assoc_ratio.phi, dplt.PAROpn$Assoc_ratio.phi)
# 
# dat.plt <- dplt.PACCGap
# p.PACCGap <- ggplot(data = dat.plt, aes(x = x, y = Assoc_ratio.pred)) + 
#   geom_ribbon(aes(x = x, ymin = Assoc_ratio.plo, ymax = Assoc_ratio.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = Assoc_ratio.pred), size = 1.5) +
#   ylim(ymin, ymax) +
#   labs(x = "Canopy gaps (%)", y = NULL)
# 
# dat.plt <- dplt.PACCOpn
# p.PACCOpn <- ggplot(data = dat.plt, aes(x = x, y = Assoc_ratio.pred)) + 
#   geom_ribbon(aes(x = x, ymin = Assoc_ratio.plo, ymax = Assoc_ratio.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = Assoc_ratio.pred), size = 1.5) +
#   ylim(ymin, ymax) +
#   labs(x = "Open forest (%)", y = NULL)
# 
# dat.plt <- dplt.PAROpn
# p.PAROpn <- ggplot(data = dat.plt, aes(x = x, y = Assoc_ratio.pred)) + 
#   geom_ribbon(aes(x = x, ymin = Assoc_ratio.plo, ymax = Assoc_ratio.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = Assoc_ratio.pred), size = 1.5) +
#   ylim(ymin, ymax) +
#   labs(x = str_wrap("Open forest perimeter-area ratio", 30), y = NULL)
# 
# dat.plt <- dplt.CanCov
# p.CanCov <- ggplot(data = dat.plt, aes(x = x, y = Assoc_ratio.pred)) + 
#   geom_ribbon(aes(x = x, ymin = Assoc_ratio.plo, ymax = Assoc_ratio.phi), alpha = 0.3) +
#   geom_line(aes(x = x, y = Assoc_ratio.pred), size = 1.5) +
#   labs(x = "Canopy cover (%)", y = NULL)
# 
# p <- ggdraw() + 
#   draw_plot(p.PACCGap, x = 0,      y = 0.5, width = 0.3333, height = 0.5) +
#   draw_plot(p.PACCOpn, x = 0.3333, y = 0.5, width = 0.3333, height = 0.5) +
#   draw_plot(p.PAROpn,  x = 0.6667, y = 0.5, width = 0.3333, height = 0.5) +
#   draw_plot(p.CanCov,  x = 0,      y = 0,   width = 1,      height = 0.5)
# p <- ggdraw() +
#   draw_plot(p, x = 0.05, y = 0, width = 0.95, height = 0.95) +
#   draw_plot_label(c("Point ratio", "Grid cell ratio", "PIPO associate ratio"),
#                   x = c(0, 0, 0.55), y = c(0.3, 0.75, 1), size = 15,
#                   angle = c(90, 90, 0), hjust = c(0.5, 0.5, 0.5))
# 
# save_plot(str_c("Plot_PIPO_associate_ratio.jpg"), p, ncol = 2, nrow = 2, dpi = 200)
