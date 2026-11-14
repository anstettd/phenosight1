library(tidyverse)
library(car)
library(patchwork)

# ---- Data -----------------------------------------------------------------
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California")))

# ---- Lag-2 predictors (24 h = 2 steps of 12 TOE) --------------------------
lag_hours <- 24

lagged <- PS2_all %>%
  select(Population, Treatment, Site, Year, TOE, ChlIdx, AriIdx) %>%
  mutate(TOE = TOE + lag_hours) %>%
  rename(ChlIdx_lag2 = ChlIdx, AriIdx_lag2 = AriIdx)

PS2_lag <- PS2_all %>%
  left_join(lagged, by = c("Population", "Treatment", "Site", "Year", "TOE"))

PS2_early_lagged <- PS2_lag %>% filter(Timepoint == "Early")
PS2_late_lagged  <- PS2_lag %>% filter(Timepoint == "Late")

# ---- Models ---------------------------------------------------------------
lm_chl_early <- lm(Fv.Fm ~ ChlIdx_lag2 + Site * Year + Treatment, data = PS2_early_lagged)
lm_ari_early <- lm(Fv.Fm ~ AriIdx_lag2 + Site * Year + Treatment, data = PS2_early_lagged)
lm_chl_late  <- lm(Fv.Fm ~ ChlIdx_lag2 + Site * Year + Treatment, data = PS2_late_lagged)
lm_ari_late  <- lm(Fv.Fm ~ AriIdx_lag2 + Site * Year + Treatment, data = PS2_late_lagged)

Anova(lm_chl_early)
Anova(lm_ari_early)
Anova(lm_chl_late)
Anova(lm_ari_late)

# ---- Theme ----------------------------------------------------------------
fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 11, face = "bold"),
  axis.text.y = element_text(size = 11, face = "bold", margin = margin(r = 1)),
  axis.title = element_text(color = "black", size = 14, face = "bold"),
  axis.title.y = element_text(color = "black", size = 14, face = "bold", margin = margin(r = 2)),
  axis.ticks.length = unit(2, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 8, face = "bold"),
  legend.key.size = unit(0.8, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  strip.text = element_text(size = 13, face = "bold"),
  strip.background = element_blank(),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
)

# ---- Plot function --------------------------------------------------------
# legend = FALSE -> no legend (Early)
# legend = TRUE  -> vertical legend inside bottom-right of California 2010 (Late)
make_xy <- function(data, xvar, xlab, legend = FALSE) {
  d <- data %>% filter(!is.na(.data[[xvar]]), !is.na(Fv.Fm))
  
  p <- ggplot(d, aes(x = .data[[xvar]], y = Fv.Fm,
                     color = Treatment, fill = Treatment, group = Treatment)) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.9) +
    facet_grid(Year ~ Site) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(x = xlab, y = "Fv/Fm") +
    fig_theme
  
  if (legend) {
    # x = 0.93: right edge of the California column; y = 0.55: bottom of the 2010 (top) row
    p <- p + theme(legend.position = "inside",
                   legend.position.inside = c(0.93, 0.55),
                   legend.justification.inside = c(1, 0),
                   legend.direction = "vertical")
  } else {
    p <- p + theme(legend.position = "none")
  }
  p
}

# ---- Individual plots -----------------------------------------------------
p_chl_early <- make_xy(PS2_early_lagged, "ChlIdx_lag2", "Chlorophyll Index (24 h prior)")
p_chl_late  <- make_xy(PS2_late_lagged,  "ChlIdx_lag2", "Chlorophyll Index (24 h prior)", legend = TRUE)
p_ari_early <- make_xy(PS2_early_lagged, "AriIdx_lag2", "Anthocyanin Index (24 h prior)")
p_ari_late  <- make_xy(PS2_late_lagged,  "AriIdx_lag2", "Anthocyanin Index (24 h prior)", legend = TRUE)

# ---- Combine early | late side by side -------------------------------------
tag_theme <- theme(plot.tag = element_text(size = 20, face = "bold"))

Fig_chl <- (p_chl_early | p_chl_late) +
  plot_annotation(tag_levels = "A", theme = tag_theme)

Fig_ari <- (p_ari_early | p_ari_late) +
  plot_annotation(tag_levels = "A", theme = tag_theme)

Fig_chl
Fig_ari

# ---- Save -----------------------------------------------------------------
dir.create("paper_graphs/z_Extra", showWarnings = FALSE, recursive = TRUE)

ggsave("paper_graphs/fvfm_chlidx_early_late.pdf",        plot = Fig_chl, width = 12, height = 6, units = "in")
ggsave("paper_graphs/z_Extra/fvfm_aridx_early_late.pdf", plot = Fig_ari, width = 12, height = 6, units = "in")