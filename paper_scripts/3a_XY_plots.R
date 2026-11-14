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
tag_y <- 0.955   # >>> vertical position of the A / B tags (1 = very top); lower number = lower <<<
tag_x <- 0.005   # >>> horizontal position of the tags <<<

fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(size = 11, face = "bold", margin = margin(r = 1)),
  axis.title = element_text(color = "black", size = 14, face = "bold"),
  axis.title.y = element_text(color = "black", size = 14, face = "bold", margin = margin(r = 2)),
  axis.ticks.length = unit(2, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 11, face = "bold"),
  legend.key.size = unit(1, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  strip.text = element_text(size = 12, face = "bold"),
  strip.background = element_blank(),
  # A / B tags sit inside each plot, level with the 20%, 40%... strip labels
  plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 0.5),
  plot.tag.location = "plot",
  plot.tag.position = c(tag_x, tag_y),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
)

# Same Site x Year groups and styling as Fig 1C: California in reds, Oregon in blues
grp_levels <- c("California 2014", "California 2010", "Oregon 2014", "Oregon 2010")
grp_cols   <- c("California 2014" = "#B2182B", "California 2010" = "#EF6B4F",
                "Oregon 2014"     = "#2166AC", "Oregon 2010"     = "#6BAED6")
grp_lines  <- c("California 2014" = "dashed", "California 2010" = "solid",
                "Oregon 2014"     = "dashed", "Oregon 2010"     = "solid")

# ---- Plot function --------------------------------------------------------
# Columns = % WHC treatment; Oregon and California share each panel (4 lines: Site x Year)
# early = TRUE -> no x-axis title, no legend, tag A (top row)
# early = FALSE -> x-axis title, vertical legend inside bottom right, tag B (bottom row)
make_xy <- function(data, xvar, xlab, early = FALSE) {
  d <- data %>%
    filter(!is.na(.data[[xvar]]), !is.na(Fv.Fm)) %>%
    mutate(grp = factor(paste(Site, Year), levels = grp_levels))
  
  p <- ggplot(d, aes(x = .data[[xvar]], y = Fv.Fm,
                     color = grp, fill = grp, linetype = grp, group = grp)) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.9) +
    facet_grid(. ~ Treatment) +
    scale_color_manual(values = grp_cols,  breaks = grp_levels, name = NULL) +
    scale_fill_manual(values = grp_cols,   breaks = grp_levels, name = NULL) +
    scale_linetype_manual(values = grp_lines, breaks = grp_levels, name = NULL) +
    scale_x_continuous(n.breaks = 3) +
    fig_theme
  
  if (early) {
    p <- p + labs(x = NULL, y = "Fv/Fm", tag = "A") +
      theme(legend.position = "none")
  } else {
    p <- p + labs(x = xlab, y = "Fv/Fm", tag = "B") +
      theme(legend.position = "inside",
            legend.position.inside = c(0.99, 0.02),   # bottom right; adjust x/y here
            legend.justification.inside = c(1, 0),
            legend.direction = "vertical",
            legend.background = element_blank(),
            legend.key.width = unit(2, "lines"))
  }
  p
}

# ---- Individual plots -----------------------------------------------------
p_chl_early <- make_xy(PS2_early_lagged, "ChlIdx_lag2", "Chlorophyll Index (24 h prior)", early = TRUE)
p_chl_late  <- make_xy(PS2_late_lagged,  "ChlIdx_lag2", "Chlorophyll Index (24 h prior)")
p_ari_early <- make_xy(PS2_early_lagged, "AriIdx_lag2", "Anthocyanin Index (24 h prior)", early = TRUE)
p_ari_late  <- make_xy(PS2_late_lagged,  "AriIdx_lag2", "Anthocyanin Index (24 h prior)")

# ---- Combine: A (Early) on top, B (Late) below, tightly stacked ----------------
Fig_chl <- p_chl_early / p_chl_late
Fig_ari <- p_ari_early / p_ari_late

Fig_chl
Fig_ari

# ---- Save -----------------------------------------------------------------
dir.create("paper_graphs/z_Extra", showWarnings = FALSE, recursive = TRUE)

ggsave("paper_graphs/fvfm_chlidx_early_late.pdf",        plot = Fig_chl, width = 14, height = 6.5, units = "in")
ggsave("paper_graphs/z_Extra/fvfm_aridx_early_late.pdf", plot = Fig_ari, width = 14, height = 6.5, units = "in")