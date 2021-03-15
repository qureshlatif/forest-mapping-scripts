## Metrics desired for prediction & mapping ##
# 1.  Overall species richness for each scenario (SR_all).
# 2.  PIPO associate richness for each scenario (SR_assoc).
# 3.  PIPO specialist richness for each scenario (SR_spec).
# 4.  Ratio of PIPO associates vs non-associates for each scenario (PIPO_assoc_ratio).
# 5.  Individual species occupancy for PIPO associates with at least 30 detections.
# 6.  Gain/loss in overall richness for fuels reductin and restoration relative to do nothing.
# 7.  Gain/loss in PIPO associate richness for fuels reductin and restoration relative to do nothing.
# 8.  Gain/loss in PIPO specialist richness for fuels reductin and restoration relative to do nothing.
# 9.  Gain/loss in PIPO associate ratio for fuels reductin and restoration relative to do nothing.
# 10. Gain/loss in individual species occupancy.
#*** For all Gain/loss metrics, flag statistically supported differences.

library(tidyverse)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod"
scripts.loc <- "forest-mapping-scripts/"
mod <- loadObject(mod.nam)
#______________________________________#

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing_community_occupancy.R"))










# Tabulate parameter estimates
pars <- c(str_c("beta.", dimnames(X.beta)[[3]]),
          str_c("alpha.", dimnames(X.alpha)[[2]]),
          str_c("zeta.", dimnames(X.zeta)[[2]]))
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

pars.ind <- which(str_detect(pars, "beta"))
for(i in pars.ind) {
  parm <- mod$sims.list$beta1[,,i]
  tbl_pars[, pars[i]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  tbl_pars[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}

pars.ind <- which(str_detect(pars, "alpha"))
for(i in pars.ind) {
  parm <- mod$sims.list$alpha1[,,i-(min(pars.ind)-1)]
  tbl_pars[, pars[i]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  tbl_pars[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}

pars.ind <- which(str_detect(pars, "zeta"))
for(i in pars.ind) {
  parm <- mod$sims.list$zeta1[,,i-(min(pars.ind)-1)]
  tbl_pars[, pars[i]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  tbl_pars[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}

rm(parm)
spp.detected <- apply(Y.mat, 2, function(x) any(x > 0))
tbl_pars <- tbl_pars[which(spp.detected),]

#### Plot grid level vegetation relationships ####
pars.sub <- c("beta.PACC10", "beta.PACC40", "beta.mnPerArRatio_Opn", "alpha.CanCov", "alpha.CanCov2")
dat.plt <- tbl_pars %>% tbl_df() %>%
  select(beta.PACC10:beta.mnPerArRatio_Opn.hi, alpha.CanCov:alpha.CanCov2.hi) %>%
  mutate(Spp = spp.list[which(spp.detected)]) %>%
  mutate(index = row_number() %>% rev())

dat.plt.supp <- dat.plt %>%
  filter_at(vars(ends_with(".lo")), any_vars(. > 0)) %>%
  bind_rows(
    dat.plt %>%
      filter_at(vars(ends_with(".hi")), any_vars(. < 0))
  ) %>%
  distinct() %>%
  arrange(index %>% desc()) %>%
  mutate(index = row_number() %>% rev())

cols <- pars.sub %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.supp), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.supp[, which(str_detect(names(dat.plt.supp), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}
dat.plt.supp <- dat.plt.supp %>%
  bind_cols(
    dat.supp %>% data.frame(stringsAsFactors = F)
  )

min.y <- dat.plt.supp %>% select(ends_with(".lo")) %>% data.matrix() %>% min()
max.y <- dat.plt.supp %>% select(ends_with(".hi")) %>% data.matrix() %>% max()

p.PACC10 <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.PACC10)) +
  geom_errorbar(aes(ymin = beta.PACC10.lo, ymax = beta.PACC10.hi, color = beta.PACC10.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.PACC10.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp), expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["PACCGap"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.PACC40 <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.PACC40)) +
  geom_errorbar(aes(ymin = beta.PACC40.lo, ymax = beta.PACC40.hi, color = beta.PACC40.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.PACC40.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp), expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["PACCOpn"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.PAR40 <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.mnPerArRatio_Opn)) +
  geom_errorbar(aes(ymin = beta.mnPerArRatio_Opn.lo, ymax = beta.mnPerArRatio_Opn.hi, color = beta.mnPerArRatio_Opn.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.mnPerArRatio_Opn.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp), expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["PAROpn"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.CanCov <- ggplot(dat = dat.plt.supp, aes(x = index, y = alpha.CanCov)) +
  geom_errorbar(aes(ymin = alpha.CanCov.lo, ymax = alpha.CanCov.hi, color = alpha.CanCov.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = alpha.CanCov.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp), expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(alpha)["CanCov"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.CanCov2 <- ggplot(dat = dat.plt.supp, aes(x = index, y = alpha.CanCov2)) +
  geom_errorbar(aes(ymin = alpha.CanCov2.lo, ymax = alpha.CanCov2.hi, color = alpha.CanCov2.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = alpha.CanCov2.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp), expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(alpha)[CanCov^2])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.PACC10,  x = 0.05, y = 0, width = 0.19, height = 1) +
  draw_plot(p.PACC40,  x = 0.24, y = 0, width = 0.19, height = 1) +
  draw_plot(p.PAR40,   x = 0.43, y = 0, width = 0.19, height = 1) +
  draw_plot(p.CanCov,  x = 0.62, y = 0, width = 0.19, height = 1) +
  draw_plot(p.CanCov2, x = 0.81, y = 0, width = 0.19, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_effects.jpg", p, ncol = 5, nrow = 3.5, dpi = 200)
