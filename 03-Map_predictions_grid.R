## Metrics desired for prediction & mapping ##
# 1. Overall species richness for each scenario (SR_all; m = 3).
# 2. PIPO specialist richness for each scenario (SR_spec; m = 3).
# 3. Ratio of PIPO specialists vs non-associates for each scenario (PIPO_special_ratio; m = 3).
# 4. Individual species occupancy for PIPO specialists with at least 30 detections (m = 54).
# 5. Gain/loss in overall richness for fuels reductin and restoration relative to do nothing (m = 2).
# 6. Gain/loss in PIPO specialist richness for fuels reductin and restoration relative to do nothing (m = 2).
# 7. Gain/loss in PIPO specialist ratio for fuels reductin and restoration relative to do nothing (m = 2).
# 8. Gain/loss in individual species occupancy (m = 36).
#*** For all Gain/loss metrics, flag statistically supported differences.
#*** Produce maps for all metrics at grid cell scale.

library(tidyverse)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod"
scripts.loc <- "forest-mapping-scripts/"
chunks.loc <- "predChunks/grid/"
mod <- loadObject(mod.nam)
#______________________________________#

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing_community_occupancy.R"))

# Get covariates for spatial prediction #
dat_grid <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/Predict_RES_FHR_DoN_USNGv2.dbf", as.is = T) %>%
  rename(Latitude = CP_Y,
         percGap_2018 = PGap2018, percGap_RES = PGapRES, percGap_FHR = PGapFHR,
         percOpn_2018 = POpen2018, percOpn_RES = POpenRES, percOpn_FHR = POpenFHR,
         PAROpn_2018 = Open2018PA, PAROpn_RES = OpenRESPAR, PAROpn_FHR = OpenFHRPAR) %>%
  mutate(TWI = ifelse(TWI == -9999, NA, TWI),
         PAROpn_2018 = PAROpn_2018 %>% as.numeric(),
         PAROpn_FHR = PAROpn_FHR %>% as.numeric(),
         PAROpn_RES = PAROpn_RES %>% as.numeric()) %>%
  mutate(PAROpn_2018 = ifelse(percOpn_2018 == 0, mean(landscape_data$mnPerArRatio_Opn), PAROpn_2018), # Treat PAR as undefined when there is no open forest.
         PAROpn_FHR = ifelse(percOpn_FHR == 0, mean(landscape_data$mnPerArRatio_Opn), PAROpn_FHR),
         PAROpn_RES = ifelse(percOpn_RES == 0, mean(landscape_data$mnPerArRatio_Opn), PAROpn_RES))

# Filter by study area #
dat_grid <- dat_grid %>%
  select(-StudyArea) %>%
  left_join(
    foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPptsPolyv2.dbf", as.is = T) %>%
      select(USNG_Code, StudyArea) %>%
      dplyr::group_by(USNG_Code) %>%
      summarise(StudyArea = max(StudyArea)),
    by = "USNG_Code"
  ) %>% filter(StudyArea == 1)

# Index PIPO specialist species #
spp_cat <- read.csv("Spp_list_detected_&_categorized.csv", header = T, stringsAsFactors = F)
PIPO_spec_ind <- which(spp.list %in% spp_cat$BirdCode[which(spp_cat$PIPO_specialist)])
Riffraff_ind <- which(spp.list %in% spp_cat$BirdCode[which(!spp_cat$PIPO_associated)])
spp_pred <- spp_cat %>% filter(Detections >= 30 & BirdCode %in% spp.list[PIPO_spec_ind]) %>% pull(BirdCode)
spp_pred_ind <- 1:length(spp.list)
names(spp_pred_ind) <- spp.list
spp_pred_ind <- spp_pred_ind[spp_pred]

# Calculate predictions #
nsamp <- length(mod$sims.list$alpha0.mu)
dat_grid$SR0 <- dat_grid$Spec0 <- dat_grid$Rat0 <-
  dat_grid$SR_FHR <- dat_grid$Spec_FHR <- dat_grid$Rat_FHR <-
  dat_grid$SR_RES <- dat_grid$Spec_RES <- dat_grid$Rat_RES <-
  dat_grid$DSR_FHR <- dat_grid$DSR_FHRp <- dat_grid$DSpc_FHR <- dat_grid$DSp_FHRp <-
  dat_grid$DRat_FHR <- dat_grid$DRt_FHRp <- dat_grid$DSR_RES <- dat_grid$DSR_RESp <-
  dat_grid$DSpc_RES <- dat_grid$DSp_RESp <- dat_grid$DRat_RES <- dat_grid$DRt_RESp <- NA
dat_grid <- dat_grid %>%
  select(USNG_Code:LowMont, SR0:DRt_RESp)
Psi_spp_pred0 <- Psi_spp_pred_FHR <- Psi_spp_pred_RES <-
  Diff_spp_FHR <- Diff_spp_RES <- Diff_spp_FHRp <- Diff_spp_RESp <-
  matrix(NA, nrow = nrow(dat_grid), ncol = length(spp_pred),
         dimnames = list(NULL, spp_pred))

# For transfer to analysis server #
#save.image("Mapping_workspace_grid.RData")
#load("Mapping_workspace_grid.RData")

# Break job up into chunks #
chunk.size <- 100

omega <-         mod$sims.list$omega %>%       array(dim = c(nsamp, length(spp.list), chunk.size))
beta0 <-         mod$sims.list$beta0 %>%       array(dim = c(nsamp, length(spp.list), chunk.size))
beta.percGap <-  mod$sims.list$beta1[,,1] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.percOpn <-  mod$sims.list$beta1[,,2] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.PAROpn <-   mod$sims.list$beta1[,,3] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.Latitude <- mod$sims.list$beta1[,,4] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.heatload <- mod$sims.list$beta1[,,5] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.TWI <-      mod$sims.list$beta1[,,6] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.LowMont <-  mod$sims.list$beta1[,,7] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))

chunks <- 1:ceiling(nrow(dat_grid) / chunk.size)
chnk.st <- (chunks * chunk.size) - (chunk.size - 1)
chnk.end <- chunks * chunk.size
chnk.end[length(chunks)] <- nrow(dat_grid)
for(chnk in chunks) {
  dat_grid_chunk <- dat_grid %>% slice(chnk.st[chnk]:chnk.end[chnk])
  
  n <- length(chnk.st[chnk]:chnk.end[chnk])
  
  Latitude <- ((dat_grid_chunk$Latitude - mean(landscape_data$Latitude)) / sd(landscape_data$Latitude)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  heatload <- ((dat_grid_chunk$heatload - mean(landscape_data$heatload)) / sd(landscape_data$heatload)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  TWI <- ((dat_grid_chunk$TWI - mean(landscape_data$TWI)) / sd(landscape_data$TWI)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  LowMont <- ((dat_grid_chunk$LowMont - mean(landscape_data$LowMont)) / sd(landscape_data$LowMont)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  
  # 2018 #
  percGap <- ((dat_grid_chunk$percGap_2018 - mean(landscape_data$PACC10 * 100)) / sd(landscape_data$PACC10 * 100)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  percOpn <- ((dat_grid_chunk$percOpn_2018 - mean(landscape_data$PACC40 * 100)) / sd(landscape_data$PACC40 * 100)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  PAROpn <- (dat_grid_chunk$PAROpn_2018 %>%
               (function(x) ifelse(is.na(x), 0, (x - mean(landscape_data$mnPerArRatio_Opn)) / sd(landscape_data$mnPerArRatio_Opn)))) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))

  psi <- (beta0[,,1:n] + beta.percGap[,,1:n] * percGap + beta.percOpn[,,1:n] * percOpn + beta.PAROpn[,,1:n] * PAROpn +
            beta.Latitude[,,1:n] * Latitude + beta.heatload[,,1:n] * heatload + beta.TWI[,,1:n] * TWI +
            beta.LowMont[,,1:n] * LowMont) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  occ0 <- omega[,,1:n] * psi
  SR_all0 <- apply(occ0, c(1, 3), sum)
  dat_grid_chunk$SR0 <- apply(SR_all0, 2, median)
  SR_spec0 <- apply(occ0[,PIPO_spec_ind,], c(1, 3), sum)
  dat_grid_chunk$Spec0 <- apply(SR_spec0, 2, median)
  SR_rifraf0 <- apply(occ0[,Riffraff_ind,], c(1, 3), sum)
  Rat0 <- SR_spec0 / SR_rifraf0
  dat_grid_chunk$Rat0 <- apply(Rat0, 2, median)
  occ_cond0 <- psi
  apply(occ_cond0[,spp_pred_ind,], c(3, 2), median) %>% saveObject(str_c(chunks.loc, "Psi_spp_pred0_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  # FHR #
  percGap <- ((dat_grid_chunk$percGap_FHR - mean(landscape_data$PACC10 * 100)) / sd(landscape_data$PACC10 * 100)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  percOpn <- ((dat_grid_chunk$percOpn_FHR - mean(landscape_data$PACC40 * 100)) / sd(landscape_data$PACC40 * 100)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  PAROpn <- (dat_grid_chunk$PAROpn_FHR %>%
               (function(x) ifelse(is.na(x), 0, (x - mean(landscape_data$mnPerArRatio_Opn)) / sd(landscape_data$mnPerArRatio_Opn)))) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))

  psi <- (beta0[,,1:n] + beta.percGap[,,1:n] * percGap + beta.percOpn[,,1:n] * percOpn + beta.PAROpn[,,1:n] * PAROpn +
            beta.Latitude[,,1:n] * Latitude + beta.heatload[,,1:n] * heatload + beta.TWI[,,1:n] * TWI +
            beta.LowMont[,,1:n] * LowMont) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  occFHR <- omega[,,1:n] * psi
  SR_all_FHR <- apply(occFHR, c(1, 3), sum)
  dat_grid_chunk$SR_FHR <- apply(SR_all_FHR, 2, median)
  SR_spec_FHR <- apply(occFHR[,PIPO_spec_ind,], c(1, 3), sum)
  dat_grid_chunk$Spec_FHR <- apply(SR_spec_FHR, 2, median)
  SR_rifraf_FHR <- apply(occFHR[,Riffraff_ind,], c(1, 3), sum)
  Rat_FHR <- SR_spec_FHR / SR_rifraf_FHR
  dat_grid_chunk$Rat_FHR <- apply(Rat_FHR, 2, median)
  occ_cond_FHR <- psi
  apply(occ_cond_FHR[,spp_pred_ind,], c(3, 2), median) %>% saveObject(str_c(chunks.loc, "Psi_spp_pred_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  dat_grid_chunk$DSR_FHR <- apply(SR_all_FHR - SR_all0, 2, median)
  dat_grid_chunk$DSR_FHRp <- apply(SR_all_FHR - SR_all0, 2, function(x) sum(x > 0) / length(x))
  dat_grid_chunk$DSpc_FHR <- apply(SR_spec_FHR - SR_spec0, 2, median)
  dat_grid_chunk$DSp_FHRp <- apply(SR_spec_FHR - SR_spec0, 2, function(x) sum(x > 0) / length(x))
  dat_grid_chunk$DRat_FHR <- apply(Rat_FHR - Rat0, 2, median)
  dat_grid_chunk$DRt_FHRp <- apply(Rat_FHR - Rat0, 2, function(x) sum(x > 0) / length(x))
  apply(occ_cond_FHR[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), median) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  apply(occ_cond_FHR[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), function(x) sum(x > 0) / length(x)) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_FHRp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  # RES #
  percGap <- ((dat_grid_chunk$percGap_RES - mean(landscape_data$PACC10 * 100)) / sd(landscape_data$PACC10 * 100)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  percOpn <- ((dat_grid_chunk$percOpn_RES - mean(landscape_data$PACC40 * 100)) / sd(landscape_data$PACC40 * 100)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  PAROpn <- (dat_grid_chunk$PAROpn_RES %>%
               (function(x) ifelse(is.na(x), 0, (x - mean(landscape_data$mnPerArRatio_Opn)) / sd(landscape_data$mnPerArRatio_Opn)))) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))

  psi <- (beta0[,,1:n] + beta.percGap[,,1:n] * percGap + beta.percOpn[,,1:n] * percOpn + beta.PAROpn[,,1:n] * PAROpn +
            beta.Latitude[,,1:n] * Latitude + beta.heatload[,,1:n] * heatload + beta.TWI[,,1:n] * TWI +
            beta.LowMont[,,1:n] * LowMont) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  occRES <- omega[,,1:n] * psi
  SR_all_RES <- apply(occRES, c(1, 3), sum)
  dat_grid_chunk$SR_RES <- apply(SR_all_RES, 2, median)
  SR_spec_RES <- apply(occRES[,PIPO_spec_ind,], c(1, 3), sum)
  dat_grid_chunk$Spec_RES <- apply(SR_spec_RES, 2, median)
  SR_rifraf_RES <- apply(occRES[,Riffraff_ind,], c(1, 3), sum)
  Rat_RES <- SR_spec_RES / SR_rifraf_RES
  dat_grid_chunk$Rat_RES <- apply(Rat_RES, 2, median)
  occ_cond_RES <- psi
  apply(occ_cond_RES[,spp_pred_ind,], c(3, 2), median) %>% saveObject(str_c(chunks.loc, "Psi_spp_pred_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  dat_grid_chunk$DSR_RES <- apply(SR_all_RES - SR_all0, 2, median)
  dat_grid_chunk$DSR_RESp <- apply(SR_all_RES - SR_all0, 2, function(x) sum(x > 0) / length(x))
  dat_grid_chunk$DSpc_RES <- apply(SR_spec_RES - SR_spec0, 2, median)
  dat_grid_chunk$DSp_RESp <- apply(SR_spec_RES - SR_spec0, 2, function(x) sum(x > 0) / length(x))
  dat_grid_chunk$DRat_RES <- apply(Rat_RES - Rat0, 2, median)
  dat_grid_chunk$DRt_RESp <- apply(Rat_RES - Rat0, 2, function(x) sum(x > 0) / length(x))
  apply(occ_cond_RES[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), median) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  apply(occ_cond_RES[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), function(x) sum(x > 0) / length(x)) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_RESp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  write.csv(dat_grid_chunk, str_c(chunks.loc, "dat_grid_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"), ".csv"), row.names = F)
  gc(verbose = F)
}

# Pull together the chunks and save for joining with shapefile in ArcGIS #
load("Mapping_workspace_grid.RData")
chunks.loc <- "predChunks/grid/"
chunk.size <- 100
chunks <- 1:ceiling(nrow(dat_grid) / chunk.size)
dat_grid <- read.csv(str_c(chunks.loc, "dat_grid_chnk0001.csv"), header = T, stringsAsFactors = F)
for(chnk in chunks[-1]) dat_grid <- dat_grid %>%
  bind_rows(read.csv(str_c(chunks.loc, "dat_grid_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"),".csv"), header = T, stringsAsFactors = F))

Psi_spp <- loadObject(str_c(chunks.loc, "Psi_spp_pred0_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Psi_spp_pred0_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "0")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Psi_spp_pred_FHR_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Psi_spp_pred_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "_FHR")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Psi_spp_pred_RES_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Psi_spp_pred_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "_RES")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_FHR_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "FHRd")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_FHRp_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_FHRp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "FHRp")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_RES_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "RESd")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_RESp_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_RESp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "RESp")
dat_grid <- dat_grid %>% bind_cols(Psi_spp)

dat_grid <- dat_grid %>% select(USNG_Code, SR0:WETARESp)
dat_grid <- dat_grid %>%
  left_join(foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/Predict_RES_FHR_DoN_USNGv2.dbf", as.is = T) %>%
              select(USNG_Code, CP_X, CP_Y),
            by = "USNG_Code")
dat_grid <- dat_grid %>%
  select(USNG_Code, CP_X, CP_Y, SR0:WETARESp)

write.table(dat_grid, "C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPgrid_pred.txt", row.names = F, sep = ",")
