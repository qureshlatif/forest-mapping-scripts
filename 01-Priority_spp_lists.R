library(stringr)
library(BCRDataAPI)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")

#__________ Script inputs _____________#
Spp_priority <- read.csv("C:/Users/Quresh.Latif/files/data/Priority Species 2020v5.csv", stringsAsFactors = F) %>%
  as_tibble() %>%
  filter(Agency == "USFS Region 2")
Spp_list <- read.csv("Spp_list.csv", stringsAsFactors = F) %>% as_tibble() %>%
  rename(Detections = Detections..max...4184.)
BBS_declining <- read.csv("C:/Users/Quresh.Latif/files/data/BBS_Trend_Estimates_2017.csv", stringsAsFactors = F) %>%
  as_tibble() %>%
  filter(Region.Code == "SUR") %>%
  rename(Species = Species.Name,
         Trend = X1966.2015.Trend.Estimates,
         Trend_CI = X1966.2015.Credible.Interval.for.Trend.Estimate) %>%
  mutate(Trend_UCL = str_split(Trend_CI, ",", simplify = T)[, 2] %>%
           str_sub(1, -2) %>% as.numeric()) %>%
  filter(Trend_UCL < 0) %>%
  select(Species, Trend, Trend_CI)
#______________________________________#

# Mark priority species #
Spp_list <- Spp_list %>%
  mutate(Priority = common_name %in% Spp_priority$Species)

# Mark declining species #
Spp_list <- Spp_list %>%
  mutate(Declining = common_name %in% BBS_declining$Species)

# Compile specialist index based on raw counts (uncomment and run to update cache file from database) #
#' BCRDataAPI::set_api_server('analysis.api.bcr.eco')
#' 
#' SampDesign <- c("IMBCR", "GRTS")
#' ss <- BCRDataAPI::subspecies()
#' BCRDataAPI::reset_api()
#' BCRDataAPI::set_api_server('analysis.api.bcr.eco')
#' BCRDataAPI::add_columns(c('TransectNum|str',
#'                           'Point|int',
#'                           'Year|int',
#'                           'CL_count|int',
#'                           'BirdCode|str',
#'                           'Species|str',
#'                           'How|str',
#'                           'Sex|str'
#' ))
#' 
#' BCRDataAPI::filter_on(c(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
#'                         #'BCR = 16',
#'                         'ninetynine = 0',
#'                         'eightyeight = 0',
#'                         'How <> F',
#'                         'Sex <> J',
#'                         'Migrant = 0',
#'                         'TimePeriod > -1',
#'                         'radialDistance < 125'))
#' grab <- BCRDataAPI::get_data(interpolate_effort = T) %>%
#'   mutate(BirdCode = ss[BirdCode] %>% as.character) %>%
#'   filter(BirdCode %in% Spp_list$BirdCode)
#' 
#' BCRDataAPI::reset_api()
#' BCRDataAPI::set_api_server('analysis.api.bcr.eco')
#' BCRDataAPI::add_columns(c('TransectNum|str',
#'                           'Point|int',
#'                           'Year|int',
#'                           'primaryHabitat|str'
#' ))
#' 
#' BCRDataAPI::filter_on(c(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
#'                         'BCR = 16'))
#' BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'primaryHabitat'))
#' grab_phab <- BCRDataAPI::get_data()
#' 
#' # Drop CO-BLM-BR and CO-BLM-GU because they were entered twice and I don't want to have to deal with removing duplicates #
#' grab <- grab %>%
#'   filter(!str_sub(TransectNum, 1, 9) %in% c("CO-BLM-GR", "CO-BLM-GU"))
#' grab_phab <- grab_phab %>%
#'   filter(!str_sub(TransectNum, 1, 9) %in% c("CO-BLM-GR", "CO-BLM-GU"))
#' 
#' # Identify species with at least 10 detections
#' spp_d10 <- grab %>%
#'   dplyr::group_by(BirdCode) %>%
#'   summarise(sumCount = sum(CL_count)) %>%
#'   filter(sumCount >= 10) %>%
#'   pull(BirdCode)
#' 
#' dat <- grab %>% filter(BirdCode %in% spp_d10) # Drop species with <10 detections.
#' dat <- dat %>% # Join primary habitat to bird detections & drop surveys with no primary habitat recorded.
#'   left_join(grab_phab, by = c("TransectNum", "Point", "Year")) %>%
#'   filter(!is.na(primaryHabitat))
#' 
#' nPIPO <- dat %>%
#'   filter(primaryHabitat == "PP") %>%
#'   select(TransectNum, Point, Year) %>%
#'   distinct() %>% nrow
#' nNonPIPO <- dat %>%
#'   filter(primaryHabitat != "PP") %>%
#'   select(TransectNum, Point, Year) %>%
#'   distinct() %>% nrow
#' 
#' dat_sum <- dat %>% # Calculate relative abundance by species within PIPO forest, everything except PIPO forest, and overall.
#'   filter(primaryHabitat == "PP") %>%
#'   dplyr::group_by(BirdCode) %>%
#'   summarise(RA_pp = sum(CL_count) / nPIPO) %>% # PIPO forest
#'   left_join(dat %>%
#'               filter(primaryHabitat != "PP") %>%
#'               dplyr::group_by(BirdCode) %>%
#'               summarise(RA_nonpp = sum(CL_count) / nNonPIPO), # non-PIPO forest
#'             by = "BirdCode") %>% # non-PIPO forest
#'   mutate(PSI = RA_pp / (RA_nonpp + RA_pp))
#' 
#' write.csv(dat_sum, "PIPO_associate_spp_cache.csv", row.names = F)


#sum(dat_sum$PSI > 0.5) # Number of species with specialization index > 0.5
#dat_sum$BirdCode[which(dat_sum$PSI > 0.5)] # Specialized species

# Mark PIPO forest associates #
dat_sum <- read.csv("PIPO_associate_spp_cache.csv", stringsAsFactors = F) %>% as_tibble()
Spp_list <- Spp_list %>%
  mutate(PIPO_associated = BirdCode %in% dat_sum$BirdCode[which(dat_sum$PSI > 0.5)],
         PIPO_specialist = BirdCode %in% dat_sum$BirdCode[which(dat_sum$PSI > 0.66)]) %>%
  filter(Detections > 0)

write.csv(Spp_list, "Spp_list_detected_&_categorized.csv", row.names = F)
