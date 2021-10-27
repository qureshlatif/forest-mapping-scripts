# Correlates ecosystem services gains/losses at the catchment level #
#***Requries GIS work to summarize and consolidate values at the catchment level***#

library(tidyverse)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
dat.catch <- read.csv("Catchments_metrics.csv", header=TRUE, stringsAsFactors = F) %>%
  mutate(D_H_FHR = H_FHR - H_no_action,
         D_H_RES = H_RES - H_no_action,
         D_RC2_FHR = RC2_FHR - RC2_no_action,
         D_RC2_RES = RC2_RES - RC2_no_action) %>%
  select(FEATURE, AreSqKM, life_zn, D_H_FHR:D_RC2_RES, DSR_FHR_pt:DRat_Res_grid, D__3_FH, D__3_RE, SF_FHR:IFH_RES) %>%
  rename(D_Sed3_FHR = D__3_FH,
         D_Sed3_RES = D__3_RE,
         LHSPI_FHR = LHSPI_F,
         LHSPI_RES = LHSPI_R) %>%
  filter(across(D_H_FHR:IFH_RES, function(x) !is.na(x)))

# Complete correlation matrices #
cor.FHR <- dat.catch %>%
  select(contains("FHR")) %>%
  data.matrix() %>%
  cor

cor.FHR.LM <- dat.catch %>%
  filter(life_zn == "Lower montane") %>%
  select(contains("FHR")) %>%
  data.matrix() %>%
  cor

cor.FHR.UM <- dat.catch %>%
  filter(life_zn == "Upper montane") %>%
  select(contains("FHR")) %>%
  data.matrix() %>%
  cor

cor.RES <- dat.catch %>%
  select(contains("RES")) %>%
  data.matrix() %>%
  cor

cor.RES.LM <- dat.catch %>%
  filter(life_zn == "Lower montane") %>%
  select(contains("RES")) %>%
  data.matrix() %>%
  cor

cor.RES.UM <- dat.catch %>%
  filter(life_zn == "Upper montane") %>%
  select(contains("RES")) %>%
  data.matrix() %>%
  cor

# Correlations relating avian to non-avian metrics #
  # Fuels reduction
dat.avian <- dat.catch %>%
  select(DSR_FHR_grid, DSR_FHR_pt, DSpec_FHR_grid, DSpec_FHR_pt, DRat_FHR_grid, DRat_FHR_pt) %>%
  data.matrix()
dat.nonav <- dat.catch %>%
  select(D_H_FHR, D_RC2_FHR, D_Sed3_FHR, SF_FHR, LHSPI_FHR, IFH_FHR) %>%
  data.matrix()
cor.FHR <- cor(dat.avian, dat.nonav)
write.csv(cor.FHR, "Catchment_Av-NonAv_correlations_FHR.csv", row.names = T)

dat.avian <- dat.catch %>%
  select(DSR_Res_grid, DSR_RES_pt, DSpec_Res_grid, DSpec_RES_pt, DRat_Res_grid, DRat_RES_pt) %>%
  data.matrix()
dat.nonav <- dat.catch %>%
  select(D_H_RES, D_RC2_RES, D_Sed3_RES, SF_RES, LHSPI_RES, IFH_RES) %>%
  data.matrix()
cor.RES <- cor(dat.avian, dat.nonav)
write.csv(cor.RES, "Catchment_Av-NonAv_correlations_RES.csv", row.names = T)
