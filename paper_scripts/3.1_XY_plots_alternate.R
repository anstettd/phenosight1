library(tidyverse)
library(car)

# ---- Data -----------------------------------------------------------------
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California")))

year_colors <- c("2010" = "skyblue3", "2014" = "#FF7700")

# ---- Lag-2 predictors (24 h = 2 steps of 12 TOE) --------------------------
# Treatment stays in the join keys only so each row picks up its lag from the
# correct series; it is not used in any model or plot.
lag_hours <- 24

lagged <- PS2_all %>%
  select(Population, Treatment, Site, Year, TOE, ChlIdx, AriIdx) %>%
  mutate(TOE = TOE + lag_hours) %>%
  rename(ChlIdx_lag2 = ChlIdx, AriIdx_lag2 = AriIdx)

PS2_lag <- PS2_all %>%
  left_join(lagged, by = c("Population", "Treatment", "Site", "Year", "TOE"))

PS2_early_lagged <- PS2_lag %>% filter(Timepoint == "Early")
PS2_late_lagged  <- PS2_lag %>% filter(Timepoint == "Late")

# ---- Models (no Treatment) ------------------------------------------------
lm_chl_early <- lm(Fv.Fm ~ ChlIdx_lag2 + Site * Year, data = PS2_early_lagged)
lm_ari_early <- lm(Fv.Fm ~ AriIdx_lag2 + Site * Year, data = PS2_early_lagged)
lm_chl_late  <- lm(Fv.Fm ~ ChlIdx_lag2 + Site * Year, data = PS2_late_lagged)
lm_ari_late  <- lm(Fv.Fm ~ AriIdx_lag2 + Site * Year, data = PS2_late_lagged)

Anova(lm_chl_early)
Anova(lm_ari_early)
Anova(lm_chl_late)
Anova(lm_ari_late)

# ---- Theme (same as Fig 1) ------------------------------------------------
legend_space <- 5

fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 11, face = "bold"),
  axis.text.y = element_text(size = 11, face = "bold", margin = margin(r = 1)),
  axis.title = element_text(color = "black", size = 12, face = "bold"),
  axis.title.y = element_text(color = "black", size = 12, face = "bold", margin = margin(r = 2)),
  axis.ticks.length = unit(2, "pt"),
  legend.position = "none",
  strip.text = element_text(size = 13, face = "bold"),
  strip.background = element_blank(),
  plot.margin = margin(t = 0, r = legend_space, b = 0, l = 5)
)

# ---- Plot function --------------------------------------------------------
# Same structure as the original XY plots: Year fill, black-outlined points,
# black lm line with the default shaded band.
make_xy <- function(data, xvar, xlab) {
  d <- data %>% filter(!is.na(.data[[xvar]]), !is.na(Fv.Fm))
  
  ggplot(d, aes(x = .data[[xvar]], y = Fv.Fm, fill = Year)) +
    geom_point(shape = 21, color = "black", size = 2) +
    geom_smooth(method = "lm", color = "black") +
    scale_fill_manual(values = year_colors, name = "Year") +
    facet_grid(Year ~ Site) +
    labs(x = xlab, y = "Fv/Fm") +
    fig_theme
}

# ---- Plots ----------------------------------------------------------------
p_chl_early <- make_xy(PS2_early_lagged, "ChlIdx_lag2", "Chlorophyll Index (24 h prior)")
p_ari_early <- make_xy(PS2_early_lagged, "AriIdx_lag2", "Anthocyanin Index (24 h prior)")
p_chl_late  <- make_xy(PS2_late_lagged,  "ChlIdx_lag2", "Chlorophyll Index (24 h prior)")
p_ari_late  <- make_xy(PS2_late_lagged,  "AriIdx_lag2", "Anthocyanin Index (24 h prior)")

p_chl_early
p_ari_early
p_chl_late
p_ari_late

ggsave("paper_graphs/z_Extra/fvfm_chlidx_lag2_early_noTrt.pdf", plot = p_chl_early, width = 7, height = 6, units = "in")
ggsave("paper_graphs/z_Extra/fvfm_aridx_lag2_early_noTrt.pdf",  plot = p_ari_early, width = 7, height = 6, units = "in")
ggsave("paper_graphs/z_Extra/fvfm_chlidx_lag2_late_noTrt.pdf",  plot = p_chl_late,  width = 7, height = 6, units = "in")
ggsave("paper_graphs/z_Extra/fvfm_aridx_lag2_late_noTrt.pdf",   plot = p_ari_late,  width = 7, height = 6, units = "in")