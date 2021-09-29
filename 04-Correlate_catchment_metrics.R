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

