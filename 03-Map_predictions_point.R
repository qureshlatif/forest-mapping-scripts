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
#*** Produce maps for all metrics at point scale integrating grid cell and point level relationships.

library(tidyverse)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod"
scripts.loc <- "forest-mapping-scripts/"
chunks.loc <- "predChunks/point/"
mod <- loadObject(mod.nam)
#______________________________________#

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing_community_occupancy.R"))

# Get covariates for spatial prediction #
dat_grid <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/Predict_RES_FHR_DoN_USNGv2.dbf", as.is = T) %>%
  rename(Latitude = CP_X,
         percGap_2018 = PGap2018, percGap_RES = PGapRES, percGap_FHR = PGapFHR,
         percOpn_2018 = POpen2018, percOpn_RES = POpenRES, percOpn_FHR = POpenFHR,
         PAROpn_2018 = Open2018PA, PAROpn_RES = OpenRESPAR, PAROpn_FHR = OpenFHRPAR) %>%
  mutate(TWI = ifelse(TWI == -9999, NA, TWI),
         PAROpn_2018 = PAROpn_2018 %>% as.numeric(),
         PAROpn_FHR = PAROpn_FHR %>% as.numeric(),
         PAROpn_RES = PAROpn_RES %>% as.numeric())
dat_point <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPptsPolyv2.dbf", as.is = T) %>%
  rename(CanCov_2018 = CC_Nothing, CanCov_RES = CC_RES, CanCov_FHR = FHR_CC) %>%
  left_join(
    dat_grid %>%
      select(USNG_Code, Latitude, percOpn_2018:TWI),
    by = "USNG_Code"
  )

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
dat_point$SR0 <- dat_point$Spec0 <- dat_point$Rat0 <-
  dat_point$SR_FHR <- dat_point$Spec_FHR <- dat_point$Rat_FHR <-
  dat_point$SR_RES <- dat_point$Spec_RES <- dat_point$Rat_RES <-
  dat_point$DSR_FHR <- dat_point$DSR_FHRp <- dat_point$DSpc_FHR <- dat_point$DSp_FHRp <-
  dat_point$DRat_FHR <- dat_point$DRt_FHRp <- dat_point$DSR_RES <- dat_point$DSR_RESp <-
  dat_point$DSpc_RES <- dat_point$DSp_RESp <- dat_point$DRat_RES <- dat_point$DRt_RESp <- NA
dat_point <- dat_point %>%
  select(Id:TWI, SR0:DRt_RESp)
Psi_spp_pred0 <- Psi_spp_pred_FHR <- Psi_spp_pred_RES <-
  Diff_spp_FHR <- Diff_spp_RES <- Diff_spp_FHRp <- Diff_spp_RESp <-
  matrix(NA, nrow = nrow(dat_point), ncol = length(spp_pred),
         dimnames = list(NULL, spp_pred))

# For transfer to analysis server #
#save.image("Mapping_workspace_point.RData")
#load("Mapping_workspace_point.RData")

# Break job up into chunks #
chunk.size <- 250

omega <-         mod$sims.list$omega %>%       array(dim = c(nsamp, length(spp.list), chunk.size))
beta0 <-         mod$sims.list$beta0 %>%       array(dim = c(nsamp, length(spp.list), chunk.size))
beta.percGap <-  mod$sims.list$beta1[,,1] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.percOpn <-  mod$sims.list$beta1[,,2] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.PAROpn <-   mod$sims.list$beta1[,,3] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.Latitude <- mod$sims.list$beta1[,,4] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.heatload <- mod$sims.list$beta1[,,5] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
beta.TWI <-      mod$sims.list$beta1[,,6] %>%  array(dim = c(nsamp, length(spp.list), chunk.size))
alpha0 <-        mod$sims.list$alpha0 %>%      array(dim = c(nsamp, length(spp.list), chunk.size))
alpha.CanCov <-  mod$sims.list$alpha1[,,1] %>% array(dim = c(nsamp, length(spp.list), chunk.size))
alpha.CanCov2 <- mod$sims.list$alpha1[,,1] %>% array(dim = c(nsamp, length(spp.list), chunk.size))

chunks <- 1:ceiling(nrow(dat_point) / chunk.size)
chnk.st <- (chunks * chunk.size) - (chunk.size - 1)
chnk.end <- chunks * chunk.size
chnk.end[length(chunks)] <- nrow(dat_point)
for(chnk in chunks) {
  dat_point_chunk <- dat_point %>% slice(chnk.st[chnk]:chnk.end[chnk])
  
  n <- length(chnk.st[chnk]:chnk.end[chnk])
  
  Latitude <- ((dat_point_chunk$Latitude - mean(landscape_data$Latitude)) / sd(landscape_data$Latitude)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  heatload <- ((dat_point_chunk$heatload - mean(landscape_data$heatload)) / sd(landscape_data$heatload)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  TWI <- ((dat_point_chunk$TWI - mean(landscape_data$TWI)) / sd(landscape_data$TWI)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  
  # 2018 #
  percGap <- ((dat_point_chunk$percGap_2018 - mean(landscape_data$PACC10)) / sd(landscape_data$PACC10)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  percOpn <- ((dat_point_chunk$percOpn_2018 - mean(landscape_data$PACC40)) / sd(landscape_data$PACC40)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  PAROpn <- (dat_point_chunk$PAROpn_2018 %>%
               (function(x) ifelse(is.na(x), 0, (x - mean(landscape_data$mnPerArRatio_Opn)) / sd(landscape_data$mnPerArRatio_Opn)))) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  CanCov <- ((dat_point_chunk$CanCov_2018 - mean(Cov[, "CanCov"])) / sd(Cov[, "CanCov"])) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  
  psi <- (beta0[,,1:n] + beta.percGap[,,1:n] * percGap + beta.percOpn[,,1:n] * percOpn + beta.PAROpn[,,1:n] * PAROpn +
            beta.Latitude[,,1:n] * Latitude + beta.heatload[,,1:n] * heatload + beta.TWI[,,1:n] * TWI) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  theta <- (alpha0[,,1:n] + alpha.CanCov[,,1:n] * CanCov + alpha.CanCov2[,,1:n] * CanCov^2) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  occ0 <- omega[,,1:n] * psi * theta
  SR_all0 <- apply(occ0, c(1, 3), sum)
  dat_point_chunk$SR0 <- apply(SR_all0, 2, median)
  SR_spec0 <- apply(occ0[,PIPO_spec_ind,], c(1, 3), sum)
  dat_point_chunk$Spec0 <- apply(SR_spec0, 2, median)
  SR_rifraf0 <- apply(occ0[,Riffraff_ind,], c(1, 3), sum)
  Rat0 <- SR_spec0 / SR_rifraf0
  dat_point_chunk$Rat0 <- apply(Rat0, 2, median)
  occ_cond0 <- psi * theta
  apply(occ_cond0[,spp_pred_ind,], c(3, 2), median) %>% saveObject(str_c(chunks.loc, "Psi_spp_pred0_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  # FHR #
  percGap <- ((dat_point_chunk$percGap_FHR - mean(landscape_data$PACC10)) / sd(landscape_data$PACC10)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  percOpn <- ((dat_point_chunk$percOpn_FHR - mean(landscape_data$PACC40)) / sd(landscape_data$PACC40)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  PAROpn <- (dat_point_chunk$PAROpn_FHR %>%
               (function(x) ifelse(is.na(x), 0, (x - mean(landscape_data$mnPerArRatio_Opn)) / sd(landscape_data$mnPerArRatio_Opn)))) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  CanCov <- ((dat_point_chunk$CanCov_FHR - mean(Cov[, "CanCov"])) / sd(Cov[, "CanCov"])) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  
  psi <- (beta0[,,1:n] + beta.percGap[,,1:n] * percGap + beta.percOpn[,,1:n] * percOpn + beta.PAROpn[,,1:n] * PAROpn +
            beta.Latitude[,,1:n] * Latitude + beta.heatload[,,1:n] * heatload + beta.TWI[,,1:n] * TWI) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  theta <- (alpha0[,,1:n] + alpha.CanCov[,,1:n] * CanCov + alpha.CanCov2[,,1:n] * CanCov^2) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  occFHR <- omega[,,1:n] * psi * theta
  SR_all_FHR <- apply(occFHR, c(1, 3), sum)
  dat_point_chunk$SR_FHR <- apply(SR_all_FHR, 2, median)
  SR_spec_FHR <- apply(occFHR[,PIPO_spec_ind,], c(1, 3), sum)
  dat_point_chunk$Spec_FHR <- apply(SR_spec_FHR, 2, median)
  SR_rifraf_FHR <- apply(occFHR[,Riffraff_ind,], c(1, 3), sum)
  Rat_FHR <- SR_spec_FHR / SR_rifraf_FHR
  dat_point_chunk$Rat_FHR <- apply(Rat_FHR, 2, median)
  occ_cond_FHR <- psi * theta
  apply(occ_cond_FHR[,spp_pred_ind,], c(3, 2), median) %>% saveObject(str_c(chunks.loc, "Psi_spp_pred_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  dat_point_chunk$DSR_FHR <- apply(SR_all_FHR - SR_all0, 2, median)
  dat_point_chunk$DSR_FHRp <- apply(SR_all_FHR - SR_all0, 2, function(x) sum(x > 0) / length(x))
  dat_point_chunk$DSpc_FHR <- apply(SR_spec_FHR - SR_spec0, 2, median)
  dat_point_chunk$DSp_FHRp <- apply(SR_spec_FHR - SR_spec0, 2, function(x) sum(x > 0) / length(x))
  dat_point_chunk$DRat_FHR <- apply(Rat_FHR - Rat0, 2, median)
  dat_point_chunk$DRt_FHRp <- apply(Rat_FHR - Rat0, 2, function(x) sum(x > 0) / length(x))
  apply(occ_cond_FHR[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), median) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  apply(occ_cond_FHR[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), function(x) sum(x > 0) / length(x)) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_FHRp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  # RES #
  percGap <- ((dat_point_chunk$percGap_RES - mean(landscape_data$PACC10)) / sd(landscape_data$PACC10)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  percOpn <- ((dat_point_chunk$percOpn_RES - mean(landscape_data$PACC40)) / sd(landscape_data$PACC40)) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  PAROpn <- (dat_point_chunk$PAROpn_RES %>%
               (function(x) ifelse(is.na(x), 0, (x - mean(landscape_data$mnPerArRatio_Opn)) / sd(landscape_data$mnPerArRatio_Opn)))) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  CanCov <- ((dat_point_chunk$CanCov_RES - mean(Cov[, "CanCov"])) / sd(Cov[, "CanCov"])) %>%
    array(dim = c(n, nsamp, length(spp.list))) %>% aperm(c(2, 3, 1))
  
  psi <- (beta0[,,1:n] + beta.percGap[,,1:n] * percGap + beta.percOpn[,,1:n] * percOpn + beta.PAROpn[,,1:n] * PAROpn +
            beta.Latitude[,,1:n] * Latitude + beta.heatload[,,1:n] * heatload + beta.TWI[,,1:n] * TWI) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  theta <- (alpha0[,,1:n] + alpha.CanCov[,,1:n] * CanCov + alpha.CanCov2[,,1:n] * CanCov^2) %>%
    (function(x) ifelse(x > 709, 709, x)) %>% # Truncate to avoid NAs
    QSLpersonal::expit()
  occRES <- omega[,,1:n] * psi * theta
  SR_all_RES <- apply(occRES, c(1, 3), sum)
  dat_point_chunk$SR_RES <- apply(SR_all_RES, 2, median)
  SR_spec_RES <- apply(occRES[,PIPO_spec_ind,], c(1, 3), sum)
  dat_point_chunk$Spec_RES <- apply(SR_spec_RES, 2, median)
  SR_rifraf_RES <- apply(occRES[,Riffraff_ind,], c(1, 3), sum)
  Rat_RES <- SR_spec_RES / SR_rifraf_RES
  dat_point_chunk$Rat_RES <- apply(Rat_RES, 2, median)
  occ_cond_RES <- psi * theta
  apply(occ_cond_RES[,spp_pred_ind,], c(3, 2), median) %>% saveObject(str_c(chunks.loc, "Psi_spp_pred_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  dat_point_chunk$DSR_RES <- apply(SR_all_RES - SR_all0, 2, median)
  dat_point_chunk$DSR_RESp <- apply(SR_all_RES - SR_all0, 2, function(x) sum(x > 0) / length(x))
  dat_point_chunk$DSpc_RES <- apply(SR_spec_RES - SR_spec0, 2, median)
  dat_point_chunk$DSp_RESp <- apply(SR_spec_RES - SR_spec0, 2, function(x) sum(x > 0) / length(x))
  dat_point_chunk$DRat_RES <- apply(Rat_RES - Rat0, 2, median)
  dat_point_chunk$DRt_RESp <- apply(Rat_RES - Rat0, 2, function(x) sum(x > 0) / length(x))
  apply(occ_cond_RES[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), median) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  apply(occ_cond_RES[,spp_pred_ind,] - occ_cond0[,spp_pred_ind,], c(3, 2), function(x) sum(x > 0) / length(x)) %>%
    saveObject(str_c(chunks.loc, "Diff_spp_RESp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0")))
  
  write.csv(dat_point_chunk, str_c(chunks.loc, "dat_point_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"), ".csv"), row.names = F)
  gc(verbose = F)
}

# Pull together the chunks and save for joining with shapefile in ArcGIS #
load("Mapping_workspace_point.RData")
chunks.loc <- "predChunks/point/"
chunk.size <- 250
chunks <- 1:ceiling(nrow(dat_point) / chunk.size)
dat_point <- read.csv(str_c(chunks.loc, "dat_point_chnk0001.csv"), header = T, stringsAsFactors = F)
for(chnk in chunks[-1]) dat_point <- dat_point %>%
  bind_rows(read.csv(str_c(chunks.loc, "dat_point_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"),".csv"), header = T, stringsAsFactors = F))

Psi_spp <- loadObject(str_c(chunks.loc, "Psi_spp_pred0_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Psi_spp_pred0_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "0")
dat_point <- dat_point %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Psi_spp_pred_FHR_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Psi_spp_pred_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "_FHR")
dat_point <- dat_point %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Psi_spp_pred_RES_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Psi_spp_pred_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "_RES")
dat_point <- dat_point %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_FHR_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_FHR_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "FHRd")
dat_point <- dat_point %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_FHRp_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_FHRp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "FHRp")
dat_point <- dat_point %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_RES_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_RES_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "RESd")
dat_point <- dat_point %>% bind_cols(Psi_spp)

Psi_spp <- loadObject(str_c(chunks.loc, "Diff_spp_RESp_chnk0001"))
for(chnk in chunks[-1]) Psi_spp <- Psi_spp %>%
  rbind(loadObject(str_c(chunks.loc, "Diff_spp_RESp_chnk", str_pad(chnk, width = 4, side = "left", pad = "0"))))
Psi_spp <- data.frame(Psi_spp)
names(Psi_spp) <- str_c(spp_pred, "RESp")
dat_point <- dat_point %>% bind_cols(Psi_spp)

dat_point <- dat_point %>% select(Id, SR0:WETARESp)
dat_point <- dat_point %>%
  left_join(foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPptsPolyv2.dbf", as.is = T) %>%
              select(Id, UTMX, UTMY, StudyArea),
            by = "Id") %>%
  filter(StudyArea == 1) %>%
  select(-StudyArea)
dat_point <- dat_point %>%
  select(Id, UTMX, UTMY, SR0:WETARESp)

write.table(dat_point, "C:/Users/Quresh.Latif/files/GIS/CEAP/PredictGrid/CEAPpts_pred.txt", row.names = F, sep = ",")
