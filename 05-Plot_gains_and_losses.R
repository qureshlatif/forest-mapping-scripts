library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

dat <- read.csv("Map_summaries_area.csv", header = T, stringsAsFactors = F) %>%
  rename(Class = X) %>%
  mutate(Zone = str_sub(Class, -7, -1)) %>%
  mutate(MetricScale = str_sub(Class, 1, -9)) %>%
  mutate(Metric = ifelse(str_detect(Class, "SR"), "All species",
                         ifelse(str_detect(Class, "Spec"), "Specialists", "S-G ratio"))) %>%
  mutate(Scale = ifelse(str_detect(Class, "Grid"), "1km cells", "250m cells")) %>%
  select(Metric, Scale, Zone, FHR_Gain:RES_Loss)

# Upper montane, fuels reduction
dat.plot <- dat %>% filter(Zone == "UppMont") %>%
  select(Metric, Scale, FHR_Gain) %>%
  mutate(Direction = "Gains") %>%
  rename(Area = FHR_Gain) %>%
  bind_rows(
    dat %>% filter(Zone == "UppMont") %>%
      select(Metric, Scale, FHR_Loss) %>%
      mutate(Direction = "Losses") %>%
      rename(Area = FHR_Loss)
  ) %>%
  select(Metric, Scale, Direction, Area) %>%
  mutate(Metric = factor(Metric, levels = c("All species", "Specialists", "S-G ratio")),
         Scale = factor(Scale, levels = c("1km cells", "250m cells")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Upper_FHR <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=expression("Area ("*km^2*")"), fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_y_continuous(limits = c(0, 3200)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1),
        axis.text.x = element_text(angle = 10))

# Upper montane, restoration
dat.plot <- dat %>% filter(Zone == "UppMont") %>%
  select(Metric, Scale, RES_Gain) %>%
  mutate(Direction = "Gains") %>%
  rename(Area = RES_Gain) %>%
  bind_rows(
    dat %>% filter(Zone == "UppMont") %>%
      select(Metric, Scale, RES_Loss) %>%
      mutate(Direction = "Losses") %>%
      rename(Area = RES_Loss)
  ) %>%
  select(Metric, Scale, Direction, Area) %>%
  mutate(Metric = factor(Metric, levels = c("All species", "Specialists", "S-G ratio")),
         Scale = factor(Scale, levels = c("1km cells", "250m cells")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Upper_RES <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=expression("Area ("*km^2*")"), fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  guides(fill = F) +
  scale_y_continuous(limits = c(0, 3200)) +
  theme(axis.text.x = element_text(angle = 10))

# Lower montane, fuels reduction
dat.plot <- dat %>% filter(Zone == "LowMont") %>%
  select(Metric, Scale, FHR_Gain) %>%
  mutate(Direction = "Gains") %>%
  rename(Area = FHR_Gain) %>%
  bind_rows(
    dat %>% filter(Zone == "UppMont") %>%
      select(Metric, Scale, FHR_Loss) %>%
      mutate(Direction = "Losses") %>%
      rename(Area = FHR_Loss)
  ) %>%
  select(Metric, Scale, Direction, Area) %>%
  mutate(Metric = factor(Metric, levels = c("All species", "Specialists", "S-G ratio")),
         Scale = factor(Scale, levels = c("1km cells", "250m cells")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Lower_FHR <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=expression("Area ("*km^2*")"), fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  guides(fill = F) +
  scale_y_continuous(limits = c(0, 1900)) +
  theme(axis.text.x = element_text(angle = 10))

# Lower montane, restoration
dat.plot <- dat %>% filter(Zone == "LowMont") %>%
  select(Metric, Scale, RES_Gain) %>%
  mutate(Direction = "Gains") %>%
  rename(Area = RES_Gain) %>%
  bind_rows(
    dat %>% filter(Zone == "UppMont") %>%
      select(Metric, Scale, RES_Loss) %>%
      mutate(Direction = "Losses") %>%
      rename(Area = RES_Loss)
  ) %>%
  select(Metric, Scale, Direction, Area) %>%
  mutate(Metric = factor(Metric, levels = c("All species", "Specialists", "S-G ratio")),
         Scale = factor(Scale, levels = c("1km cells", "250m cells")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Lower_RES <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=expression("Area ("*km^2*")"), fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  guides(fill = F) +
  scale_y_continuous(limits = c(0, 1900)) +
  theme(axis.text.x = element_text(angle = 10))

# Assemble panels #
p <- ggdraw() + 
  draw_plot(p_Upper_FHR, x = 0,   y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(p_Upper_RES, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(p_Lower_FHR, x = 0,   y = 0,   width = 0.5, height = 0.5) +
  draw_plot(p_Lower_RES, x = 0.5, y = 0,   width = 0.5, height = 0.5)
p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0, width = 0.95, height = 0.95) +
  draw_plot_label(c("Fuels reduction", "Restoration",
                    "Lower Montane", "Upper Montane"),
                  x = c(0.35, 0.8, 0, 0),
                  y = c(1, 1, 0.25, 0.72), size = 20,
                  angle = c(0, 0, 90, 90),
                  hjust = c(0.5, 0.5, 0.5, 0.5))

save_plot("Plot_gains_losses.jpg", p, ncol = 2, nrow = 3, dpi = 200)
