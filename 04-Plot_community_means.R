# Metrics to summarize #
# 1. Area with supported gain/loss in species richness with restoration & fuels reduction
# 2. Area with supported gain/loss in specialist richness with restoration/fuels reduction
# 3. Area with supported gain/loss in specialist composition ratio with restoration/fuels reduction

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

dat_point <- read.table("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPpts_pred.txt", header = T, stringsAsFactors = F, sep = ",") %>%
  select(Id:DRt_RESp) %>%
  left_join(
    foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPptsPolyv2.dbf", as.is = T) %>%
      select(Id, CC_Nothing, CC_RES, FHR_CC, USNG_Code) %>%
      rename(CanCov_0 = CC_Nothing, CanCov_RES = CC_RES, CanCov_FHR = FHR_CC),
    by = "Id"
  )
dat_grid <- read.table("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPgrid_pred.txt", header = T, stringsAsFactors = F, sep = ",") %>%
  select(USNG_Code:DRt_RESp) %>%
  left_join(
    foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/Predict_RES_FHR_DoN_USNGv2.dbf", as.is = T) %>%
      select(USNG_Code, POpen2018:PGapRES, LowMont) %>%
      rename(percOpn_0 = POpen2018, percGap_0 = PGap2018, PAROpn_0 = Open2018PA,
             percOpn_FHR = POpenFHR, percGap_FHR = PGapFHR, PAROpn_FHR = OpenFHRPAR,
             percOpn_RES = POpenRES, percGap_RES = PGapRES, PAROpn_RES = OpenRESPAR) %>%
      mutate(PAROpn_0 = as.numeric(PAROpn_0),
             PAROpn_FHR = as.numeric(PAROpn_FHR),
             PAROpn_RES = as.numeric(PAROpn_RES)),
    by = "USNG_Code"
  ) %>%
  mutate(PAROpn_0 = ifelse(percOpn_0 == 0, mean(landscape_data$mnPerArRatio_Opn), PAROpn_0), # Treat PAR as undefined when there is no open forest.
         PAROpn_FHR = ifelse(percOpn_FHR == 0, mean(landscape_data$mnPerArRatio_Opn), PAROpn_FHR),
         PAROpn_RES = ifelse(percOpn_RES == 0, mean(landscape_data$mnPerArRatio_Opn), PAROpn_RES))
dat_point <- dat_point %>% left_join(
  dat_grid %>% select(USNG_Code, LowMont),
  by = "USNG_Code"
)

## Summarize values for plotting ##
dat_grid_long <- dat_grid %>%
  select(LowMont, SR0, Spec0, Rat0) %>%
  rename(SR = SR0, Spec = Spec0, Rat = Rat0) %>%
  mutate(Scenario = "Baseline") %>%
  bind_rows(
    dat_grid %>%
      select(LowMont, SR_FHR, Spec_FHR, Rat_FHR) %>%
      rename(SR = SR_FHR, Spec = Spec_FHR, Rat = Rat_FHR) %>%
      mutate(Scenario = "Fuels reduction")
  ) %>%
  bind_rows(
    dat_grid %>%
      select(LowMont, SR_RES, Spec_RES, Rat_RES) %>%
      rename(SR = SR_RES, Spec = Spec_RES, Rat = Rat_RES) %>%
      mutate(Scenario = "Restoration")
    ) %>%
  mutate(LifeZone = ifelse(LowMont == 1, "Lower montane", "Upper montane")) %>%
  mutate(LifeZone = factor(LifeZone, levels = c("Lower montane", "Upper montane")),
         Scenario = factor(Scenario, levels = c("Baseline", "Fuels reduction", "Restoration")))

dat_point_long <- dat_point %>%
  select(LowMont, SR0, Spec0, Rat0) %>%
  rename(SR = SR0, Spec = Spec0, Rat = Rat0) %>%
  mutate(Scenario = "Baseline") %>%
  bind_rows(
    dat_point %>%
      select(LowMont, SR_FHR, Spec_FHR, Rat_FHR) %>%
      rename(SR = SR_FHR, Spec = Spec_FHR, Rat = Rat_FHR) %>%
      mutate(Scenario = "Fuels reduction")
  ) %>%
  bind_rows(
    dat_point %>%
      select(LowMont, SR_RES, Spec_RES, Rat_RES) %>%
      rename(SR = SR_RES, Spec = Spec_RES, Rat = Rat_RES) %>%
      mutate(Scenario = "Restoration")
  ) %>%
  mutate(LifeZone = ifelse(LowMont == 1, "Lower montane", "Upper montane")) %>%
  mutate(LifeZone = factor(LifeZone, levels = c("Lower montane", "Upper montane")),
         Scenario = factor(Scenario, levels = c("Baseline", "Fuels reduction", "Restoration")))

## Define functions for customizing boxplot ##
f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), type = 8)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

o <- function(x) {
  subset(x, x < quantile(x, prob = 0.05, type = 8) | quantile(x, prob = 0.95, type = 8) < x)
}

## Generate box plots ##
p.SR.grid <- ggplot(dat_grid_long, aes(x = Scenario, y = SR)) +
  facet_grid(. ~ LifeZone) +
  stat_summary(fun.data = f, geom="boxplot") +
  stat_summary(fun = o, geom="point", alpha = 0.3, size = 0.2, positiion = "jitter") +
  ylab("Species richness") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 10))

p.Spec.grid <- ggplot(dat_grid_long, aes(x = Scenario, y = Spec)) +
  facet_grid(. ~ LifeZone) +
  stat_summary(fun.data = f, geom="boxplot") +
  stat_summary(fun = o, geom="point", alpha = 0.3, size = 0.2) +
  ylab("Specialist richness") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 10))

p.Rat.grid <- ggplot(dat_grid_long, aes(x = Scenario, y = Rat)) +
  facet_grid(. ~ LifeZone) +
  stat_summary(fun.data = f, geom="boxplot") +
  stat_summary(fun = o, geom="point", alpha = 0.3, size = 0.2) +
  ylab("Specialist-generalist ratio") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 10))


p.SR.point <- ggplot(dat_point_long, aes(x = Scenario, y = SR)) +
  facet_grid(. ~ LifeZone) +
  stat_summary(fun.data = f, geom="boxplot") +
  stat_summary(fun = o, geom="point", alpha = 0.3, size = 0.2) +
  ylab("Species richness") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 10))

p.Spec.point <- ggplot(dat_point_long, aes(x = Scenario, y = Spec)) +
  facet_grid(. ~ LifeZone) +
  stat_summary(fun.data = f, geom="boxplot") +
  stat_summary(fun = o, geom="point", alpha = 0.3, size = 0.2) +
  ylab("Specialist richness") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 10))

p.Rat.point <- ggplot(dat_point_long, aes(x = Scenario, y = Rat)) +
  facet_grid(. ~ LifeZone) +
  stat_summary(fun.data = f, geom="boxplot") +
  stat_summary(fun = o, geom="point", alpha = 0.3, size = 0.2) +
  ylab("Specialist-generalist ratio") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 10))


# Assemble panels #
p <- ggdraw() + 
  draw_plot(p.SR.grid,    x = 0,   y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(p.Spec.grid,  x = 0,   y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(p.Rat.grid,   x = 0,   y = 0,      width = 0.5, height = 0.3333) +
  draw_plot(p.SR.point,   x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(p.Spec.point, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(p.Rat.point,  x = 0.5, y = 0,      width = 0.5, height = 0.3333)
p <- ggdraw() +
  draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
  draw_plot_label(c("Grid scale (100 ha)", "Point scale (6.25 ha)"),
                  x = c(0.3, 0.78),
                  y = c(1, 1), size = 20,
                  angle = c(0, 0),
                  hjust = c(0.5, 0.5))

#save_plot("Plot_community_predictions.jpg", p, ncol = 1.75, nrow = 3, dpi = 200)
save_plot("manuscript/Figure7.jpg", p, ncol = 1.75, nrow = 3, dpi = 600)
