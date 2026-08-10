##################################################################################
## Get select graphs for 20 and 100%
## 
## Author Daniel Anstett
## 
## 
## Written July 2026
###################################################################################
#Import libraries
library(tidyverse)

#Import Data
PS2_all <- read_csv("data/PS2_all.csv")
PS2_all <- PS2_all %>% mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
                              Timepoint = factor(Timepoint, levels = c("Early", "Late")))

#Colors: solid line/point colors by Year, and lighter fill colors for geom_smooth ribbons
year_colors <- c("2010" = "skyblue3", "2014" = "#FF7700")
year_fill_colors <- c("2010" = "#8FC1E3", "2014" = "#FFA858")

############################################################################################################
###################################### Fv.Fm: 20% vs 100%, Both Sites ######################################
############################################################################################################
p_fvfm <- PS2_all %>%
  dplyr::filter(Treatment %in% c("20%", "100%")) %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "100%")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
  ggplot(aes(x = TOE, y = Fv.Fm, color = Year, fill = Year)) +
  geom_point() +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values = year_fill_colors, name = "Year") +
  geom_smooth() +
  facet_grid(Treatment ~ Site,
             labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California"))) +
  coord_cartesian(ylim = c(0.75, NA)) +
  labs(x = "Time of Experiment", y = "Quantum Yield of PSII (Fv/Fm)") +
  theme_bw() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(color = "black", size = 26, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 26, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 26, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_text(size = 24, face = "bold"),
    strip.background = element_blank()
  )

ggsave("graphs/multi_join/select/fvfm_20_100.pdf",
       plot = p_fvfm, width = 10, height = 10, units = "in")

############################################################################################################
###################################### ChlIdx: 20% vs 100%, Both Sites #####################################
############################################################################################################

p_chlidx <- PS2_all %>%
  dplyr::filter(Treatment %in% c("20%", "100%")) %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "100%")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
  ggplot(aes(x = TOE, y = ChlIdx, color = Year, fill = Year)) +
  geom_point() +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values = year_fill_colors, name = "Year") +
  geom_smooth() +
  facet_grid(Treatment ~ Site,
             labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California"))) +
  labs(x = "Time of Experiment", y = "Chlorophyll Index") +
  theme_bw() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(color = "black", size = 26, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 26, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_text(size = 24, face = "bold"),
    strip.background = element_blank()
  )

ggsave("graphs/multi_join/select/chlidx_20_100.pdf",
       plot = p_chlidx, width = 10, height = 10, units = "in")

library(cowplot)

############################################################################################################
###################################### Combined: Side-by-side (wide) #######################################
############################################################################################################

# Drop the legend and the Treatment facet strip from the left-hand plot so
# only one shared legend and one set of Treatment labels (on the right plot) appear
p_fvfm_nolegend <- p_fvfm + theme(
  legend.position = "none",
  strip.background.y = element_blank(),
  strip.text.y = element_blank()
)

# Move the legend into the empty space in the upper part of the lower-right
# facet panel (Treatment = 100%, Site = SCalifornia), instead of it sitting
# off to the side. This is only applied for the combined figure -- the
# standalone chlidx_20_100.pdf keeps its normal side legend.
# >>> ADJUST THESE TWO NUMBERS (x, y) TO REPOSITION THE LEGEND <<<
p_chlidx_combined <- p_chlidx + theme(
  legend.position = "inside",
  legend.position.inside = c(0.99, 0.49),   # <-- tweak x/y (0-1 range) here
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = NA),
  legend.text = element_text(size = 22, face = "bold"),
  legend.key.size = unit(3, "lines"),
  legend.key.height = unit(2.4, "lines")
)

p_combined_wide <- plot_grid(p_fvfm_nolegend, p_chlidx_combined, nrow = 1)

ggsave("graphs/multi_join/select/fvfm_chlidx_2x4.pdf",
       plot = p_combined_wide, width = 20, height = 10, units = "in")

############################################################################################################
###################################### Combined: Stacked (tall) #############################################
############################################################################################################

p_combined_tall <- plot_grid(p_fvfm, p_chlidx, ncol = 1)

ggsave("graphs/multi_join/select/fvfm_chlidx_4x2.pdf",
       plot = p_combined_tall, width = 10, height = 20, units = "in")



