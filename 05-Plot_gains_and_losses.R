library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

dat <- read.csv("Map_summaries_prparea.csv", header = T, stringsAsFactors = F) %>%
  rename(Class = X) %>%
  mutate(Zone = str_sub(Class, -7, -1)) %>%
  mutate(MetricScale = str_sub(Class, 1, -9)) %>%
  mutate(Metric = ifelse(str_detect(Class, "SR"), "All species",
                         ifelse(str_detect(Class, "Spec"), "Specialists", "S-G ratio"))) %>%
  mutate(Scale = ifelse(str_detect(Class, "Grid"), "Grid scale (100 ha)", "Point scale (6.25 ha)")) %>%
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
         Scale = factor(Scale, levels = c("Grid scale (100 ha)", "Point scale (6.25 ha)")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Upper_FHR <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=NULL, fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(fill=guide_legend(reverse=TRUE)) +
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
         Scale = factor(Scale, levels = c("Grid scale (100 ha)", "Point scale (6.25 ha)")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Upper_RES <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=NULL, fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  guides(fill = F) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 10))

# Lower montane, fuels reduction
dat.plot <- dat %>% filter(Zone == "LowMont") %>%
  select(Metric, Scale, FHR_Gain) %>%
  mutate(Direction = "Gains") %>%
  rename(Area = FHR_Gain) %>%
  bind_rows(
    dat %>% filter(Zone == "LowMont") %>%
      select(Metric, Scale, FHR_Loss) %>%
      mutate(Direction = "Losses") %>%
      rename(Area = FHR_Loss)
  ) %>%
  select(Metric, Scale, Direction, Area) %>%
  mutate(Metric = factor(Metric, levels = c("All species", "Specialists", "S-G ratio")),
         Scale = factor(Scale, levels = c("Grid scale (100 ha)", "Point scale (6.25 ha)")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Lower_FHR <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=NULL, fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  guides(fill = F) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 10))

# Lower montane, restoration
dat.plot <- dat %>% filter(Zone == "LowMont") %>%
  select(Metric, Scale, RES_Gain) %>%
  mutate(Direction = "Gains") %>%
  rename(Area = RES_Gain) %>%
  bind_rows(
    dat %>% filter(Zone == "LowMont") %>%
      select(Metric, Scale, RES_Loss) %>%
      mutate(Direction = "Losses") %>%
      rename(Area = RES_Loss)
  ) %>%
  select(Metric, Scale, Direction, Area) %>%
  mutate(Metric = factor(Metric, levels = c("All species", "Specialists", "S-G ratio")),
         Scale = factor(Scale, levels = c("Grid scale (100 ha)", "Point scale (6.25 ha)")),
         Direction = factor(Direction, levels = c("Losses", "Gains")))
p_Lower_RES <- ggplot(dat.plot, aes(x = Metric, y = Area, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Scale) +
  labs(x=NULL, y=NULL, fill=NULL) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  guides(fill = F) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 10))

# Assemble panels #
p <- ggdraw() + 
  draw_plot(p_Upper_FHR, x = 0,   y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(p_Upper_RES, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(p_Lower_FHR, x = 0,   y = 0,   width = 0.5, height = 0.5) +
  draw_plot(p_Lower_RES, x = 0.5, y = 0,   width = 0.5, height = 0.5)
p <- ggdraw() +
  draw_plot(p, x = 0.07, y = 0, width = 0.93, height = 0.95) +
  draw_plot_label(c("Fuels reduction", "Restoration",
                    "Proportion landscape",
                    "Lower Montane", "Upper Montane"),
                  x = c(0.32, 0.79, 0.02, 0, 0),
                  y = c(1, 1, 0.5, 0.1, 0.9),
                  size = c(20, 20, 23, 20, 20),
                  angle = c(0, 0, 90, 90, 90),
                  hjust = c(0.5, 0.5, 0.5, 0, 1),
                  fontface = c("bold", "bold", "italic", "bold", "bold"))

save_plot("manuscript/Figure8.jpg", p, ncol = 2, nrow = 3, dpi = 600)
