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

## Summaries of differences with reference scenario ##
cols <- c("DFHR_low", "DRES_low", "DFHR_upp", "DRES_upp")
rows <- c("SRGrid", "SRPoint", "SpecGrid", "SpecPoint", "RatGrid", "RatPoint")
out <- matrix("", nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

out["SRGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(DSR_FHR, DSR_RES) %>%
  summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(DSR_FHR, DSR_RES) %>%
      summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character())
out["SpecGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(DSpc_FHR, DSpc_RES) %>%
  summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(DSpc_FHR, DSpc_RES) %>%
      summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character())
out["RatGrid", ] <- dat_grid %>% filter(LowMont == 1) %>%
  select(DRat_FHR, DRat_RES) %>%
  summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character() %>%
  c(dat_grid %>% filter(LowMont == 0) %>%
      select(DRat_FHR, DRat_RES) %>%
      summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character())

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
  summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character() %>%
  c(dat_point %>% filter(LowMont == 0) %>%
      select(DRat_FHR, DRat_RES) %>%
      summarise_all(sum.fn, dig = 3) %>% as.matrix() %>% as.character())

write.csv(out, "Map_summaries_diff.csv")

## Area with supported gains or losses ##
cols <- c("FHR_Gain", "FHR_Loss", "RES_Gain", "RES_Loss")
rows <- c("SRGrid_LowMont", "SRPoint_LowMont", "SRGrid_UppMont", "SRPoint_UppMont",
          "SpecGrid_LowMont", "SpecPoint_LowMont", "SpecGrid_UppMont", "SpecPoint_UppMont",
          "RatGrid_LowMont", "RatPoint_LowMont", "RatGrid_UppMont", "RatPoint_UppMont")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "FHR_Gain"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric()
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "FHR_Loss"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric()
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "RES_Gain"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric()
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "RES_Loss"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric()

out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "FHR_Gain"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) round(sum(x >= 0.975) / 16)) %>% data.matrix() %>% as.numeric()
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "FHR_Loss"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) round(sum(x <= 0.025) / 16)) %>% data.matrix() %>% as.numeric()
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "RES_Gain"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) round(sum(x >= 0.975) / 16)) %>% data.matrix() %>% as.numeric()
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "RES_Loss"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) round(sum(x <= 0.025) / 16)) %>% data.matrix() %>% as.numeric()

out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "FHR_Gain"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric()
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "FHR_Loss"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric()
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "RES_Gain"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric()
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "RES_Loss"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric()

out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "FHR_Gain"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) round(sum(x >= 0.975) / 16)) %>% data.matrix() %>% as.numeric()
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "FHR_Loss"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) round(sum(x <= 0.025) / 16)) %>% data.matrix() %>% as.numeric()
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "RES_Gain"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) round(sum(x >= 0.975) / 16)) %>% data.matrix() %>% as.numeric()
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "RES_Loss"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) round(sum(x <= 0.025) / 16)) %>% data.matrix() %>% as.numeric()

write.csv(out, "Map_summaries_area.csv")

## Proportion area with supported gains or losses ##
cols <- c("FHR_Gain", "FHR_Loss", "RES_Gain", "RES_Loss")
rows <- c("SRGrid_LowMont", "SRPoint_LowMont", "SRGrid_UppMont", "SRPoint_UppMont",
          "SpecGrid_LowMont", "SpecPoint_LowMont", "SpecGrid_UppMont", "SpecPoint_UppMont",
          "RatGrid_LowMont", "RatPoint_LowMont", "RatGrid_UppMont", "RatPoint_UppMont")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

max.extent <- sum(dat_grid$LowMont == 1)
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "FHR_Gain"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "FHR_Loss"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "RES_Gain"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRGrid_LowMont", "SpecGrid_LowMont", "RatGrid_LowMont"), "RES_Loss"] <- dat_grid %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)

max.extent <- sum(dat_point$LowMont == 1)
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "FHR_Gain"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "FHR_Loss"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "RES_Gain"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRPoint_LowMont", "SpecPoint_LowMont", "RatPoint_LowMont"), "RES_Loss"] <- dat_point %>%
  filter(LowMont == 1) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)

max.extent <- sum(dat_grid$LowMont == 0)
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "FHR_Gain"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "FHR_Loss"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "RES_Gain"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRGrid_UppMont", "SpecGrid_UppMont", "RatGrid_UppMont"), "RES_Loss"] <- dat_grid %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)

max.extent <- sum(dat_point$LowMont == 0)
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "FHR_Gain"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "FHR_Loss"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_FHRp, DSp_FHRp, DRt_FHRp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "RES_Gain"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x >= 0.975)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)
out[c("SRPoint_UppMont", "SpecPoint_UppMont", "RatPoint_UppMont"), "RES_Loss"] <- dat_point %>%
  filter(LowMont == 0) %>% select(DSR_RESp, DSp_RESp, DRt_RESp) %>%
  summarise_all(function(x) sum(x <= 0.025)) %>% data.matrix() %>% as.numeric() %>%
  (function(x) x / max.extent)

write.csv(out, "Map_summaries_prparea.csv")

## Relate with heterogeneity ##
dat_catch <- read.csv("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/heterogeneity_metrics.csv", header = T, stringsAsFactors = F) %>%
  select(FEATURE, H, RC2, TRT)
dat_catch <- dat_catch %>% filter(TRT == "UNT") %>% select(-TRT) %>%
  rename(H_UNT = H, RC2_UNT = RC2) %>%
  left_join(
    dat_catch %>% filter(TRT == "FHR") %>% select(-TRT) %>%
              rename(H_FHR = H, RC2_FHR = RC2),
            by = "FEATURE"
    ) %>%
  left_join(
    dat_catch %>% filter(TRT == "RES") %>% select(-TRT) %>%
      rename(H_RES = H, RC2_RES = RC2),
    by = "FEATURE"
  ) %>%
  mutate(H_FHRd = H_FHR - H_UNT,
         RC2_FHRd = RC2_FHR - RC2_UNT,
         H_RESd = H_RES - H_UNT,
         RC2_RESd = RC2_RES - RC2_UNT)
#hist(dat_catch$H_FHRd)
#hist(dat_catch$H_RESd)
#hist(dat_catch$RC2_FHRd)
#hist(dat_catch$RC2_RESd)


##*** Need to rejoin catchment IDs with point files containing predictions if we want this again ***##

# dat_point_sum <- dat_point %>% left_join(
#   foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPpts_pred.dbf", as.is = T) %>%
#     select(Id, CatchID),
#   by = "Id"
#   ) %>%
#   left_join(dat_catch, by = c("CatchID" = "FEATURE")) %>%
#   dplyr::group_by(CatchID) %>%
#   summarise(SR_FHR = mean(SR_FHR),
#             SR_RES = mean(SR_RES),
#             DSR_FHR = mean(SR_FHR - SR0),
#             DSR_RES = mean(SR_RES - SR0),
#             H_FHR = mean(H_FHR),
#             H_RES = mean(H_RES),
#             H_FHRd = mean(H_FHRd),
#             H_RESd = mean(H_RESd),
#             RC2_FHR = mean(RC2_FHR),
#             RC2_RES = mean(RC2_RES),
#             RC2_FHRd = mean(RC2_FHRd),
#             RC2_RESd = mean(RC2_RESd),
#             LowMont = first(LowMont))
# 
# cor(dat_point_sum$H_FHR[which(dat_point_sum$LowMont == 1)], dat_point_sum$SR_FHR[which(dat_point_sum$LowMont == 1)])
# cor(dat_point_sum$H_FHR[which(dat_point_sum$LowMont == 0)], dat_point_sum$SR_FHR[which(dat_point_sum$LowMont == 0)])
# 
# cor(dat_point_sum$H_RES[which(dat_point_sum$LowMont == 1)], dat_point_sum$SR_RES[which(dat_point_sum$LowMont == 1)])
# cor(dat_point_sum$H_RES[which(dat_point_sum$LowMont == 0)], dat_point_sum$SR_RES[which(dat_point_sum$LowMont == 0)])
# 
# cor(dat_point_sum$H_FHRd[which(dat_point_sum$LowMont == 1)], dat_point_sum$DSR_FHR[which(dat_point_sum$LowMont == 1)])
# cor(dat_point_sum$H_FHRd[which(dat_point_sum$LowMont == 0)], dat_point_sum$DSR_FHR[which(dat_point_sum$LowMont == 0)])
# 
# cor(dat_point_sum$H_RESd[which(dat_point_sum$LowMont == 1)], dat_point_sum$DSR_RES[which(dat_point_sum$LowMont == 1)])
# cor(dat_point_sum$H_RESd[which(dat_point_sum$LowMont == 0)], dat_point_sum$DSR_RES[which(dat_point_sum$LowMont == 0)])
# 
# dat_grid_sum <- dat_grid %>% left_join(
#   foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPgrid_pred.dbf", as.is = T) %>%
#     select(USNG_Code, CatchID),
#   by = "USNG_Code"
# ) %>%
#   left_join(dat_catch, by = c("CatchID" = "FEATURE")) %>%
#   dplyr::group_by(CatchID) %>%
#   summarise(SR_FHR = mean(SR_FHR),
#             SR_RES = mean(SR_RES),
#             DSR_FHR = mean(SR_FHR - SR0),
#             DSR_RES = mean(SR_RES - SR0),
#             H_FHR = mean(H_FHR),
#             H_RES = mean(H_RES),
#             H_FHRd = mean(H_FHRd),
#             H_RESd = mean(H_RESd),
#             RC2_FHR = mean(RC2_FHR),
#             RC2_RES = mean(RC2_RES),
#             RC2_FHRd = mean(RC2_FHRd),
#             RC2_RESd = mean(RC2_RESd),
#             LowMont = first(LowMont))
# 
# cor(dat_grid_sum$H_FHR[which(dat_grid_sum$LowMont == 1)], dat_grid_sum$SR_FHR[which(dat_grid_sum$LowMont == 1)], use = "complete")
# cor(dat_grid_sum$H_FHR[which(dat_grid_sum$LowMont == 0)], dat_grid_sum$SR_FHR[which(dat_grid_sum$LowMont == 0)], use = "complete")
# 
# cor(dat_grid_sum$H_RES[which(dat_grid_sum$LowMont == 1)], dat_grid_sum$SR_RES[which(dat_grid_sum$LowMont == 1)], use = "complete")
# cor(dat_grid_sum$H_RES[which(dat_grid_sum$LowMont == 0)], dat_grid_sum$SR_RES[which(dat_grid_sum$LowMont == 0)], use = "complete")
# 
# cor(dat_grid_sum$H_FHRd[which(dat_grid_sum$LowMont == 1)], dat_grid_sum$DSR_FHR[which(dat_grid_sum$LowMont == 1)], use = "complete")
# cor(dat_grid_sum$H_FHRd[which(dat_grid_sum$LowMont == 0)], dat_grid_sum$DSR_FHR[which(dat_grid_sum$LowMont == 0)], use = "complete")
# 
# # At grid cell level, positive relationship between change in heterogeneity and change in species richness under restoration scenario
# cor(dat_grid_sum$H_RESd[which(dat_grid_sum$LowMont == 1)], dat_grid_sum$DSR_RES[which(dat_grid_sum$LowMont == 1)], use = "complete")
# plot(dat_grid_sum$H_RESd[which(dat_grid_sum$LowMont == 1)], dat_grid_sum$DSR_RES[which(dat_grid_sum$LowMont == 1)])
# cor(dat_grid_sum$H_RESd[which(dat_grid_sum$LowMont == 0)], dat_grid_sum$DSR_RES[which(dat_grid_sum$LowMont == 0)], use = "complete")
# plot(dat_grid_sum$H_RESd[which(dat_grid_sum$LowMont == 0)], dat_grid_sum$DSR_RES[which(dat_grid_sum$LowMont == 0)])
# 
# # I suspect with restoration, heterogeneity is positively related with treatment intensity, whereas with fuels reduction the opposite is true.
# # This is probably why you have opposite relationships between change in hetergeneity and change in species richness with FHR (negative) vs RES (positive).
# 
# ## Summarize and tabulate heterogeneity ##
# dat_catch <- dat_catch %>%
#   left_join(
#     dat_point_sum %>% select(CatchID, LowMont) %>%
#       distinct(),
#     by = c("FEATURE" = "CatchID")
#   )
# 
# cols <- c("Ref", "FHR", "FHR_Diff", "RES", "RES_Diff")
# rows <- c("H_LowM", "RC2_LowM", "H_UppM", "RC2_UppM")
# out <- matrix("", nrow = length(rows), ncol = length(cols),
#               dimnames = list(rows, cols))
# 
# out["H_LowM", ] <- dat_catch %>% filter(LowMont == 1) %>%
#   select(H_UNT, H_FHR, H_FHRd, H_RES, H_RESd) %>%
#   summarise_all(sum.fn, dig = 2) %>% as.matrix() %>% as.character()
# out["H_UppM", ] <- dat_catch %>% filter(LowMont == 0) %>%
#   select(H_UNT, H_FHR, H_FHRd, H_RES, H_RESd) %>%
#   summarise_all(sum.fn, dig = 2) %>% as.matrix() %>% as.character()
# out["RC2_LowM", ] <- dat_catch %>% filter(LowMont == 1) %>%
#   select(RC2_UNT, RC2_FHR, RC2_FHRd, RC2_RES, RC2_RESd) %>%
#   summarise_all(sum.fn, dig = 2) %>% as.matrix() %>% as.character()
# out["RC2_UppM", ] <- dat_catch %>% filter(LowMont == 0) %>%
#   select(RC2_UNT, RC2_FHR, RC2_FHRd, RC2_RES, RC2_RESd) %>%
#   summarise_all(sum.fn, dig = 2) %>% as.matrix() %>% as.character()
# 
# write.csv(out, "Heterogeneity_summaries.csv")
