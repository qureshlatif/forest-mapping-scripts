# Metrics to summarize #
# 1. Area with supported gain/loss in species richness with restoration & fuels reduction
# 2. Area with supported gain/loss in specialist richness with restoration/fuels reduction
# 3. Area with supported gain/loss in specialist composition ratio with restoration/fuels reduction

library(tidyverse)

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

sum.fn <- function(x, dig = 1) {
  median(x) %>%
    round(digits = dig) %>%
    str_c(
      " (",
      quantile(x, prob = 0.025, type = 8) %>% round(digits = dig),
      ",",
      quantile(x, prob = 0.975, type = 8) %>% round(digits = dig),
      ")"
    )
}

## Summaries of the landscape ##
cols <- c("Ref_Low", "FHR_Low", "RES_Low", "Ref_Upp", "FHR_Upp", "RES_Upp")
rows <- c("SRGrid", "SRPoint", "SpecGrid", "SpecPoint", "RatGrid", "RatPoint", "percGap", "percOpn", "PAROpn", "CanCov")
out <- matrix("", nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

out["SRGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(SR0, SR_FHR, SR_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(SR0, SR_FHR, SR_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["SpecGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(Spec0, Spec_FHR, Spec_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(Spec0, Spec_FHR, Spec_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["RatGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(Rat0, Rat_FHR, Rat_RES) %>%
  summarise_all(function(x) sum.fn(x, dig = 2)) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(Rat0, Rat_FHR, Rat_RES) %>%
      summarise_all(function(x) sum.fn(x, dig = 2)) %>% as.matrix() %>% as.character())

out["SRPoint", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(SR0, SR_FHR, SR_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(SR0, SR_FHR, SR_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["SpecPoint", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(Spec0, Spec_FHR, Spec_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(Spec0, Spec_FHR, Spec_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["RatPoint", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(Rat0, Rat_FHR, Rat_RES) %>%
  summarise_all(function(x) sum.fn(x, dig = 2)) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(Rat0, Rat_FHR, Rat_RES) %>%
      summarise_all(function(x) sum.fn(x, dig = 2)) %>% as.matrix() %>% as.character())

out["percGap", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(percGap_0, percGap_FHR, percGap_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(percGap_0, percGap_FHR, percGap_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["percOpn", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(percOpn_0, percOpn_FHR, percOpn_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(percOpn_0, percOpn_FHR, percOpn_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["PAROpn", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(PAROpn_0, PAROpn_FHR, PAROpn_RES) %>%
  summarise_all(function(x) sum.fn(x, dig = 2)) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(PAROpn_0, PAROpn_FHR, PAROpn_RES) %>%
      summarise_all(function(x) sum.fn(x, dig = 2)) %>% as.matrix() %>% as.character())
out["CanCov", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(CanCov_0, CanCov_FHR, CanCov_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(CanCov_0, CanCov_FHR, CanCov_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())

write.csv(out, "Map_summaries.csv")

## Summaries of the landscape ##
cols <- c("DFHR_low", "DRES_low", "DFHR_upp", "DRES_upp")
rows <- c("SRGrid", "SRPoint", "SpecGrid", "SpecPoint", "RatGrid", "RatPoint")
out <- matrix("", nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

out["SRGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(DSR_FHR, DSR_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(DSR_FHR, DSR_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["SpecGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(DSpc_FHR, DSpc_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(DSpc_FHR, DSpc_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["RatGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(DRat_FHR, DRat_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(DRat_FHR, DRat_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())

out["SRPoint", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(DSR_FHR, DSR_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(DSR_FHR, DSR_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["SpecPoint", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(DSpc_FHR, DSpc_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(DSpc_FHR, DSpc_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())
out["RatPoint", ] <- dat_point %>% filter(LowMont == 1) %>%
  select(DRat_FHR, DRat_RES) %>%
  summarise_all(sum.fn) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(DRat_FHR, DRat_RES) %>%
      summarise_all(sum.fn) %>% as.matrix() %>% as.character())

write.csv(out, "Map_summaries_diff.csv")

## Area with supported gains or losses ##
cols <- c("FHR_Gain", "FHR_Loss", "RES_Gain", "RES_Loss")
rows <- c("SRGrid", "SRPoint", "SpecGrid", "SpecPoint", "RatGrid", "RatPoint")
out <- matrix("", nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

out[c(1, 3, 5), "FHR_Gain"] <- dat_grid %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric()
out[c(1, 3, 5), "FHR_Loss"] <- dat_grid %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric()
out[c(1, 3, 5), "RES_Gain"] <- dat_grid %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric()
out[c(1, 3, 5), "RES_Loss"] <- dat_grid %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric()

out[c(2, 4, 6), "FHR_Gain"] <- dat_point %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) round(sum(x >= 0.975) / 16)) %>% data.matrix() %>% as.numeric()
out[c(2, 4, 6), "FHR_Loss"] <- dat_point %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) round(sum(x <= 0.025) / 16)) %>% data.matrix() %>% as.numeric()
out[c(2, 4, 6), "RES_Gain"] <- dat_point %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) round(sum(x >= 0.975) / 16)) %>% data.matrix() %>% as.numeric()
out[c(2, 4, 6), "RES_Loss"] <- dat_point %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) round(sum(x <= 0.025) / 16)) %>% data.matrix() %>% as.numeric()

write.csv(out, "Map_summaries_area.csv")
