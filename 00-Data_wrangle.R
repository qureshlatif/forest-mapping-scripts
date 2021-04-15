library(stringr)
library(BCRDataAPI)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")

################## Inputs ####################
# Get latest AOU checklist with tax names and order #
aou.checklist <- read.csv("C:/Users/Quresh.Latif/files/data/NACC_list_bird_species_downloaded_20200416.csv",
                          header = T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(tax_ord = row_number())

spp.exclude <- c("Squirrel, Red", "Ruffed Grouse", "Turkey Vulture", "Wild Turkey",
                 "Sandhill Crane", "Bald Eagle", "American Kestrel", "Red-tailed Hawk",
                 "Great Blue Heron", "Swainson's Hawk", "Canada Goose", "Squirrel, Abert's",
                 "Northern Pygmy-Owl", "Northern Goshawk", "Sharp-shinned Hawk", "Green-winged Teal",
                 "Cooper's Hawk", "Great Horned Owl", "Pika", "Gambel's Quail", "Osprey",
                 "Common Merganser", "White-tailed Ptarmigan", "Peregrine Falcon",
                 "Boreal Owl", "Spotted Owl", "Black-crowned Night-Heron", "Ring-necked Duck",
                 "California Gull", "Northern Saw-whet Owl", "Long-eared Owl", "Flammulated Owl",
                 "Prairie Falcon", "Northern Harrier", "American White Pelican", "Western Screech-Owl",
                 "Double-crested Cormorant", "Bufflehead", "Thicket Tinamou", "Dusky Grouse", "Mallard",
                 "Golden Eagle", "Gadwall", "Virginia Rail", "Chukar")
strata <- c("CO-CFLRP-CF", "CO-BCR16-RC", "CO-BCR16-PC")
SampDesign <- c("IMBCR", "GRTS")
##############################################

#### Compile species list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('BirdCode|str',
                          'Species|str')
)
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('BirdCode <> NOBI')
BCRDataAPI::filter_on('BCR = 16')
BCRDataAPI::filter_on('primaryHabitat in LP,MC,II,PP')
BCRDataAPI::filter_on(str_c('Year in ', str_c(2008:2020, collapse = ",")))
BCRDataAPI::group_by(c('BirdCode', 'Species'))
grab <- BCRDataAPI::get_data()

spp.out <- grab %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(!Species %in% spp.exclude)

# Collapsing sub-species and renamed species #
ss <- BCRDataAPI::subspecies()
spp.out <- spp.out %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)%>%
  dplyr::group_by(BirdCode) %>%
  mutate(min_length = min(nchar(Species))) %>%
  mutate(Species = str_sub(Species, 1, min_length)) %>%
  select(BirdCode, Species) %>%
  # Additional tweaks #
  mutate(Species = replace(Species, which(Species %in% c("Western Scrub-Jay", "Woodhouse's Scrub")), "Woodhouse's Scrub-Jay")) %>%
  ungroup %>%
  unique

#sum(!spp.out$Species %in% aou.checklist$common_name) # check - should be zero
spp.out <- spp.out %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

# Remove additional implausible members of the metacommunity (based on review of BNA range maps and habitat accounts) #
spp.out <- spp.out %>%
  filter(!BirdCode %in% c("RUHU", "PAWR", "OLWA", "AMPI", "WWCR", "SABS"))

spp.excluded <- grab %>%
  select(BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(Species %in% spp.exclude) %>%
  select(BirdCode, Species) %>%
  unique %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

#### Detection data ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'PointLatitude|num',
                          'PointLongitude|num',
                          'zone|int',
                          'Year|int',
                          'Date|str',
                          'PointVisitStartTime|str',
                          'easting|int',
                          'northing|int',
                          'zone|int',
                          'Stratum|str',
                          'radialDistance|int',
                          'CL_count|int',
                          'BirdCode|str',
                          'Species|str',
                          'How|str',
                          'Sex|str',
                          'TimePeriod|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'Year in 2014,2016,2018',
                        str_c('BirdCode in ', str_c(spp.out$BirdCode, collapse = ",")),
                        'ninetynine = 0',
                        'eightyeight = 0',
                        'How <> F',
                        'Sex <> J',
                        'Migrant = 0',
                        'TimePeriod > -1',
                        'radialDistance < 125'))
grab <- BCRDataAPI::get_data(interpolate_effort = T) %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)

point.coords <- grab %>%
  select(TransectNum, Point, easting, northing, zone) %>%
  unique
point.list <- unique(str_c(point.coords$TransectNum,
                           str_pad(point.coords$Point, width = 2, side = "left", pad = "0"), sep = "-")) %>%
  sort
grid.list <- unique(point.coords$TransectNum) %>% sort

## Point X years surveyed ##
pointXyears.list <- unique(str_c(grab$TransectNum,
                                 str_pad(grab$Point, width = 2,
                                         side = "left", pad = "0"),
                                 grab$Year, sep = "-")) %>% sort

## Add number of detections and count summaries to spp.out by stratum ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>%
  unique %>% dplyr::group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, CL_count) %>%
  dplyr::group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

spp.out <- spp.out %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount), (function(x) replace(x, is.na(x), 0)))

maxDetPossible <- length(pointXyears.list) # max possible by stratum
names(spp.out)[which(names(spp.out) == "Detections")] <-
  str_c("Detections (max = ", maxDetPossible, ")")

write.csv(spp.out, "Spp_list.csv", row.names = F)
rm(smry)

## Add number of detections and count summaries to excluded species ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>% unique %>%
  dplyr::group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, TransectNum, Point, CL_count) %>%
  dplyr::group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

spp.excluded <- spp.excluded %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount),
            (function(x) replace(x, is.na(x), 0)))

write.csv(spp.excluded, "Spp_excluded.csv", row.names = F)
rm(smry)

bird_data <- grab %>%  # Store bird survey data for later use.
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year))

## Get GIS covariates ##
dat.gis <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/FS/CFLRP/Bird_survey_point_coords.dbf", as.is = T) %>%
  tbl_df() %>%
  rename(Grid = TransectNu) %>%
  select(Grid, Point, Northing, heatload, TWI, LifeZone)

## Get grid-level landscape structure (LANDFIRE) covariates ##
landscape_data <- data.frame(Grid = pointXyears.list %>% str_sub(1, -9),
                             Year = pointXyears.list %>% str_sub(-4, -1), stringsAsFactors = F) %>%
  mutate(gridIndex = Grid %>% as.factor %>% as.integer,
         YearInd = Year %>% as.factor %>% as.integer) %>%
  unique

PA_data <- read.csv("data/CFLRP Percent Areas 2014.csv", header = T, stringsAsFactors = F) %>%
  tbl_df %>%
  mutate(Year = "2014") %>%
  bind_rows(
    read.csv("data/CFLRP Percent Areas 2016.csv", header = T, stringsAsFactors = F) %>%
      tbl_df %>%
      mutate(Year = "2016")
  ) %>%
  bind_rows(
    read.csv("data/Gaps_3km_2018.csv", header = T, stringsAsFactors = F) %>%
      tbl_df() %>%
      rename(Gap_3km_square = PercentArea,
             TransectNum = TransNum) %>%
      select(TransectNum, Gap_3km_square) %>%
      left_join(read.csv("data/Open_3km_2018.csv", header = T, stringsAsFactors = F) %>%
                  tbl_df() %>%
                  rename(Open_3km_square = PercentArea,
                         TransectNum = TransNum) %>%
                  select(TransectNum, Open_3km_square),
                by = "TransectNum") %>%
      mutate(Year = "2018")
  ) %>%
  rename(PACC10 = Gap_3km_square,
         PACC40 = Open_3km_square) %>%
  select(TransectNum, Year, PACC10, PACC40)

landscape_data <- landscape_data %>%
  left_join(
    PA_data %>%
      select(TransectNum, Year, PACC10, PACC40),
    by = c("Grid" = "TransectNum", "Year" = "Year")
  ) %>%
  left_join(read.csv("data/Patch_struct_&_config.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              mutate(Year = as.character(Year),
                     # Rescale area and distance covariates to ha and km
                     mnPtchAr_Gap = mnPtchAr_Gap3km / 10000,
                     mnPtchAr_Opn = mnPtchAr_Opn3km / 10000,
                     mnPerArRatio_Gap = mnPerArRatio_Gap3km,
                     mnPerArRatio_Opn = mnPerArRatio_Opn3km,
                     NNdist_Gap = NNdist_Gap3km / 1000,
                     NNdist_Opn = NNdist_Opn3km / 1000) %>%
              select(TransNum, Year, mnPtchAr_Gap:NNdist_Opn) %>%
              bind_rows(
                read.csv("data/Gaps_3km_2018.csv", header = T, stringsAsFactors = F) %>%
                  tbl_df() %>%
                  rename(mnPtchAr_Gap = MeanPatch,
                         mnPerArRatio_Gap = para_mn,
                         NNdist_Gap = enn_mn) %>%
                  mutate(NNdist_Gap = NNdist_Gap / 1000) %>%
                  select(TransNum, mnPtchAr_Gap, NNdist_Gap, mnPerArRatio_Gap) %>%
                  left_join(
                    read.csv("data/Open_3km_2018.csv", header = T, stringsAsFactors = F) %>%
                      tbl_df() %>%
                      rename(mnPtchAr_Opn = MeanPatch,
                             mnPerArRatio_Opn = para_mn,
                             NNdist_Opn = enn_mn) %>%
                      mutate(NNdist_Opn = NNdist_Opn / 1000) %>%
                      select(TransNum, mnPtchAr_Opn, NNdist_Opn, mnPerArRatio_Opn),
                    by = "TransNum"
                  ) %>%
                  mutate(Year = "2018")
              ),
            by = c("Grid" = "TransNum", "Year" = "Year")) %>%
  select(Grid:NNdist_Opn)

Mode <- function(x) { # Define function for getting the most frequent level
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
dat.gis <- dat.gis %>%
  dplyr::group_by(Grid) %>%
  summarize(Latitude = mean(Northing),
            heatload = mean(heatload),
            TWI = mean(TWI),
            LifeZone = Mode(LifeZone)) %>%
  mutate(LowMont = as.integer(LifeZone == "Lower montane"))

landscape_data <- landscape_data %>%
  left_join(dat.gis, by = "Grid")

rm(PA_data)
rm(dat.gis)

#cor(landscape_data %>% select(PACC10:TWI), use = "complete")
# Remove mnPtchAr_Gap & mnPtchAr_Opn from planned analysis model

## Get point level canopy cover values ##
veg_data <- data.frame(Point_year = pointXyears.list, stringsAsFactors = F) %>%
  tbl_df() %>%
  mutate(Point = str_sub(Point_year, 1, -6),
         Year = str_sub(Point_year, -4, -1)) %>%
  left_join(
    read.csv("data/PointLevel_ContinuousCC_2014.csv", header = T, stringsAsFactors = F) %>%
      tbl_df() %>%
      select(TransectNum, Point, Year, cc2014ext) %>%
      mutate(Point = str_c(TransectNum, "-", str_pad(Point, side = "left", width = 2, pad = "0")),
             Year = as.character(Year)) %>%
      rename(CanCov = cc2014ext) %>%
      bind_rows(
        read.csv("data/PointLevel_ContinuousCC_2016.csv", header = T, stringsAsFactors = F) %>%
          tbl_df() %>%
          select(TransectNum, Point, Year, cc2016ext) %>%
          mutate(Point = str_c(TransectNum, "-", str_pad(Point, side = "left", width = 2, pad = "0")),
                 Year = as.character(Year)) %>%
          rename(CanCov = cc2016ext)
      ) %>%
      bind_rows(
        read.csv("data/PointLevel_ContinuousCC_2018.csv", header = T, stringsAsFactors = F) %>%
          tbl_df() %>%
          select(TransectNum, Point, Year, cc2018ext) %>%
          mutate(Point = str_c(TransectNum, "-", str_pad(Point, side = "left", width = 2, pad = "0")),
                 Year = as.character(Year)) %>%
          rename(CanCov = cc2018ext)
      ) %>%
      select(-TransectNum),
    by = c("Point", "Year")
  ) %>%
  mutate(CanCov = ifelse(is.na(CanCov), 0, CanCov))

## Trim dates, compile day of year & start time in minutes ##
library(lubridate)
bird_data <- bird_data %>%
  mutate(Date = str_sub(Date, 6, -14) %>% dmy) %>%
  mutate(DOY = yday(Date)) %>%
  mutate(PointVisitStartTime = PointVisitStartTime %>%
           replace(which(PointVisitStartTime == "0"), NA)) %>%
  mutate(HR = PointVisitStartTime %>% str_sub(1, -3) %>% as.integer) %>%
  mutate(MIN = PointVisitStartTime %>% str_sub(-2, -1) %>% as.integer) %>%
  mutate(Time = str_c(HR, MIN, "00", sep = ":")) %>%
  mutate(dateTime = str_c(Date, Time, sep = " ")) %>%
  mutate(Time_ssr = QSLpersonal::tssr(PointLatitude, PointLongitude, dateTime)) %>%
  select(TransectNum:Date, DOY, Time_ssr, PointVisitStartTime:TimePeriod)

## Compile multidimensional detection data array ##
spp.list <- spp.out$BirdCode
cov.names <- c("gridIndex", "YearInd", "DayOfYear", "Time_ssr", names(veg_data)[-c(1:3)])

bird_data <- bird_data %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year))
Y.mat <- matrix(NA, nrow = length(pointXyears.list), ncol = length(spp.list),
               dimnames = list(pointXyears.list, spp.list))
TR.mat <- matrix(6, nrow = length(pointXyears.list), ncol = length(spp.list),
               dimnames = list(pointXyears.list, spp.list))
for(sp in 1:length(spp.list)) {
  obs <- bird_data %>% filter(BirdCode == spp.list[sp])
  if(nrow(obs) > 0) {
    Y.mat[, sp] <- (pointXyears.list %in% obs$Point_year) %>% as.integer
    tvec <- tapply(obs$TimePeriod, obs$Point_year, min)
    tvec <- tvec[order(names(tvec))]
    TR.mat[which(pointXyears.list %in% obs$Point_year), sp] <- tvec
  } else {
    Y.mat[, sp] <- 0
  }
}

Cov <- matrix(NA, nrow = length(pointXyears.list), ncol = length(cov.names),
              dimnames = list(pointXyears.list, cov.names))
Cov[, "gridIndex"] <- pointXyears.list %>% str_sub(1, -9) %>% as.factor %>% as.integer
Cov[, "YearInd"] <- pointXyears.list %>% str_sub(-4, -1) %>% as.factor %>% as.integer
Cov[, "DayOfYear"] <- bird_data %>%
  select(Point_year, DOY) %>% distinct() %>% arrange(Point_year) %>% pull(DOY)
Cov[, "Time_ssr"] <- bird_data %>%
  select(Point_year, Time_ssr) %>% distinct() %>% arrange(Point_year) %>%
  pull(Time_ssr)
Cov[, -c(1:4)] <- veg_data %>%
  arrange(Point_year) %>%
  select(-Point_year) %>%
  select(-Point) %>%
  select(-Year) %>%
  data.matrix()

rm(obs, maxDetPossible, sp, ss, tvec, grab)
save.image("Data_compiled.RData")

