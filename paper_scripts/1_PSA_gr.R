##################################################################################
## Fig 1: Projected Shoot Area (A), Growth Rate Boxplots (B), Growth Rate vs SWHC (C)
## In format for paper
##
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(patchwork)   # needs ggplot2 >= 3.5.0 for legend.position.inside
library(multcompView) # multcompLetters() for the Tukey letters in B

#Import Data
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year      = factor(Year, levels = c(2010, 2014)),
         Site      = factor(Site, levels = c("Oregon", "California")))

GR_all <- read_csv("data/GR_all.csv") %>%
  mutate(SWHC      = parse_number(Treatment),                        # 20, 40, ... 100 for regression in C
         Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year      = factor(Year, levels = c(2010, 2014)),
         Site      = factor(Site, levels = c("Oregon", "California")))

#Overall model
mod_PSA <- lm(PSA ~ Site * Year + Treatment, data = PS2_all)
Anova(mod_PSA)

############################################################################################################
#Shared theme (used by A, B and C)

fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  panel.spacing = unit(0.8, "lines"),
  axis.text = element_text(size = 12, face = "bold"),
  axis.text.y = element_text(size = 12, face = "bold", margin = margin(r = 1)),
  axis.title = element_text(color = "black", size = 14, face = "bold"),
  axis.title.y = element_text(color = "black", size = 14, face = "bold", margin = margin(r = 2)),
  axis.ticks.length = unit(3, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 11, face = "bold"),
  legend.key.size = unit(1, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  legend.justification = c(0, 0.5),
  strip.text = element_text(size = 14, face = "bold"),
  strip.background = element_blank(),
  plot.tag = element_text(size = 24, face = "bold", hjust = 0, vjust = 1),
  plot.tag.location = "plot",
  plot.tag.position = c(0.005, 1),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
)

# Larger x and y tick labels, plus a little more space between the y-axis title and the axis (B and C only)
tick_size_BC <- 15   # >>> tick label size for B and C; adjust here <<<
ytitle_gap_BC <- 7   # >>> gap (pt) between the y-axis title and the tick labels in B and C; adjust here <<<
tick_theme_BC <- theme(
  axis.text.x = element_text(size = tick_size_BC, face = "bold"),
  axis.text.y = element_text(size = tick_size_BC, face = "bold", margin = margin(r = 1)),
  axis.title.y = element_text(color = "black", size = 14, face = "bold", margin = margin(r = ytitle_gap_BC))
)

# Larger x tick labels for the A panels
tick_size_A <- 15    # >>> x tick label size for A; adjust here <<<
tick_theme_A <- theme(
  axis.text.x = element_text(size = tick_size_A, face = "bold")
)

# Larger y tick labels for all three A graphs
tick_size_A_y <- 13  # >>> y tick label size for all of A; adjust here <<<
tick_theme_A_y <- theme(
  axis.text.y = element_text(size = tick_size_A_y, face = "bold", margin = margin(r = 1))
)

# Extra gap between the y-axis title and the tick labels on the upper A graph
ytitle_gap_A <- 8    # >>> gap (pt) for the upper A graph's y title; adjust here <<<

# Pull the y titles of the middle (-log(P)) and lower (R2) A graphs in toward the axis
ytitle_gap_A_lower <- -4   # >>> pt; negative moves the title inward. Try -2 for less, -6 for more <<<

############################################################################################################
#Fig 1A - Projected Shoot Area over time, with P-values and adjusted R2 across TOE

PSA_summ <- PS2_all %>%
  group_by(Site, Year, Treatment, TOE) %>%
  summarise(mean_PSA = mean(PSA, na.rm = TRUE),
            se_PSA   = sd(PSA, na.rm = TRUE) / sqrt(sum(!is.na(PSA))),
            .groups = "drop")

PSA_plot <- ggplot(PSA_summ, aes(x = TOE, y = mean_PSA, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_ribbon(aes(ymin = mean_PSA - se_PSA, ymax = mean_PSA + se_PSA), alpha = 0.2, color = NA) +
  geom_vline(xintercept = 410, color = "black", linetype = "11",
             linewidth = 0.5, alpha = 0.6, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  facet_grid(Year ~ Site) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = NULL, y = "Projected Shoot Area", tag = "A") +
  fig_theme + tick_theme_A + tick_theme_A_y +
  theme(axis.title.y = element_text(color = "black", size = 14, face = "bold",
                                    margin = margin(r = ytitle_gap_A)),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.01),   # bottom right of California 2014; adjust x/y here
        legend.justification = c(1, 0),
        legend.background = element_blank())

fit_one <- function(d) {
  m <- lm(PSA ~ Site * Year + Treatment, data = d,
          contrasts = list(Site = "contr.sum", Year = "contr.sum", Treatment = "contr.sum"))
  a <- as.data.frame(Anova(m, type = 3)) %>%
    rownames_to_column("Term") %>%
    filter(Term %in% c("Site", "Year", "Treatment", "Site:Year")) %>%
    select(Term, P = `Pr(>F)`)
  a$adjR2 <- summary(m)$adj.r.squared
  a
}

toe_stats <- PS2_all %>%
  group_by(TOE) %>%
  group_modify(~ tryCatch(fit_one(.x),
                          error = function(e) tibble(Term = NA_character_, P = NA_real_, adjR2 = NA_real_))) %>%
  ungroup() %>%
  filter(!is.na(Term)) %>%
  mutate(Term = factor(Term, levels = c("Site:Year", "Treatment", "Site", "Year")),
         negLogP = -log10(P))

term_cols   <- c("Site:Year" = "royalblue3", "Treatment" = "#D95F02",
                 "Site" = "#7570B3", "Year" = "#E7298A")
term_labels <- c("Site:Year" = "Site*Year", "Treatment" = "Treatment",
                 "Site" = "Site", "Year" = "Year")

#Middle: -log10(P); its legend is drawn inside the R2 graph below, so none here
P_plot <- ggplot(toe_stats, aes(x = TOE, y = negLogP, color = Term)) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "11", linewidth = 0.8) +
  geom_vline(xintercept = 410, color = "black", linetype = "11",
             linewidth = 0.5, alpha = 0.6, show.legend = FALSE) +
  geom_line(linewidth = 0.9, linetype = "solid") +
  scale_color_manual(values = term_cols, labels = term_labels) +
  labs(x = NULL, y = "-log(P)") +
  fig_theme + tick_theme_A_y +
  theme(legend.position = "none",
        axis.title.y = element_text(color = "black", size = 14, face = "bold",
                                    margin = margin(r = ytitle_gap_A_lower)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Bottom: adjusted R2, carrying the legend for the P-value lines above it.
#The invisible layer (all-NA y, so nothing is drawn) only exists to create the legend entries.
R2_plot <- toe_stats %>% distinct(TOE, adjR2) %>%
  ggplot(aes(x = TOE, y = adjR2)) +
  geom_vline(xintercept = 410, color = "black", linetype = "11",
             linewidth = 0.5, alpha = 0.6, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_line(data = toe_stats %>% mutate(adjR2 = NA_real_),
            aes(color = Term), linewidth = 0.9, na.rm = TRUE) +
  scale_color_manual(values = term_cols, labels = term_labels) +
  labs(x = "Time of Experiment (h)", y = expression(bold(R^2))) +
  fig_theme + tick_theme_A + tick_theme_A_y +
  theme(axis.title.y = element_text(color = "black", size = 14, face = "bold",
                                    margin = margin(r = ytitle_gap_A_lower)),
        legend.position = "inside",
        legend.position.inside = c(0.995, 0.93),   # upper right of R2; adjust x/y here
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.85, "lines"),
        legend.key.width = unit(1.2, "lines"))

col_A <- (PSA_plot / P_plot / R2_plot) + plot_layout(heights = c(4, 1, 1))

############################################################################################################
#Fig 1B - Growth rate box plots with Tukey letters
#Tukey: one ANOVA per Site x Timepoint panel across all 10 Treatment x Year groups

gr_letters <- GR_all %>%
  mutate(tgrp = paste(Treatment, Year, sep = "_")) %>%
  group_by(Site, Timepoint) %>%
  group_modify(function(d, key) {
    fit <- aov(GR ~ tgrp, data = d)
    tk  <- TukeyHSD(fit)$tgrp[, "p adj"]
    L   <- multcompLetters(tk)$Letters
    tibble(tgrp = names(L), Letters = unname(L))
  }) %>%
  ungroup()

# Letters sit just above the top of each box's data (including outliers)
y_rng <- diff(range(GR_all$GR, na.rm = TRUE))

GR_letter_pos <- GR_all %>%
  mutate(tgrp = paste(Treatment, Year, sep = "_")) %>%
  group_by(Site, Timepoint, Treatment, Year, tgrp) %>%
  summarise(top = max(GR, na.rm = TRUE), .groups = "drop") %>%
  left_join(gr_letters, by = c("Site", "Timepoint", "tgrp")) %>%
  mutate(y = top + 0.03 * y_rng)

# Shared y range with headroom at the top, where the legend sits (California, Late panel)
y_lo <- min(GR_all$GR, na.rm = TRUE) - 0.02 * y_rng
y_hi <- max(GR_all$GR, na.rm = TRUE) + 0.30 * y_rng

# Hand-drawn legend, stacked vertically (2014 above 2010), top right of the California / Late panel only
lb_site <- factor("California", levels = levels(GR_all$Site))
lb_time <- factor("Late",       levels = levels(GR_all$Timepoint))
key_h   <- 0.045 * y_rng                                   # half-height of each legend key
lb_mid14 <- y_hi - 0.07 * y_rng                            # >>> move legend up/down: edit these two rows <<<
lb_mid10 <- y_hi - 0.22 * y_rng
key_x0 <- 4.40; key_x1 <- 4.75; txt_x <- 4.82              # >>> move legend left/right: edit these <<<

lb_keys <- tibble(Site = lb_site, Timepoint = lb_time,
                  Year = c("2014", "2010"),
                  xmin = key_x0, xmax = key_x1,
                  mid  = c(lb_mid14, lb_mid10),
                  ymin = mid - key_h, ymax = mid + key_h)

Fig1B <- ggplot(GR_all, aes(x = Treatment, y = GR, fill = Year)) +
  geom_boxplot(position = position_dodge(width = 0.8), linewidth = 0.6, outlier.size = 1.2,
               fatten = 1) +                                   # median line same thickness as the box outline
  geom_text(data = GR_letter_pos,
            aes(x = Treatment, y = y, label = Letters, group = Year),
            position = position_dodge(width = 0.8),
            size = 5, fontface = "bold", vjust = 0, inherit.aes = FALSE) +
  # legend layers: colored key, horizontal median line through its middle, then the label
  geom_rect(data = lb_keys, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Year),
            color = "black", linewidth = 0.6, inherit.aes = FALSE, show.legend = FALSE) +
  geom_segment(data = lb_keys, aes(x = xmin, xend = xmax, y = mid, yend = mid),
               color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  geom_text(data = lb_keys, aes(x = txt_x, y = mid, label = Year),
            size = 5.5, fontface = "bold", hjust = 0, vjust = 0.5, inherit.aes = FALSE) +
  facet_grid(Timepoint ~ Site) +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20)) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  labs(x = "% Soil Water-Holding Capacity", y = "Growth Rate", tag = "B") +
  fig_theme + tick_theme_BC +
  theme(legend.position = "none")

############################################################################################################
#Fig 1C - Growth rate regressed on % Soil Water-Holding Capacity within each Site x Year group, Early and Late

# One group per Site x Year; California in reds, Oregon in blues
grp_levels <- c("California 2014", "California 2010", "Oregon 2014", "Oregon 2010")
grp_cols   <- c("California 2014" = "#B2182B", "California 2010" = "#EF6B4F",
                "Oregon 2014"     = "#2166AC", "Oregon 2010"     = "#6BAED6")
grp_shapes <- c("California 2014" = 17, "California 2010" = 16,
                "Oregon 2014"     = 17, "Oregon 2010"     = 16)
grp_lines  <- c("California 2014" = "dashed", "California 2010" = "solid",
                "Oregon 2014"     = "dashed", "Oregon 2010"     = "solid")

GR_C <- GR_all %>%
  mutate(grp = factor(paste(Site, Year), levels = grp_levels))

gr_head <- max(GR_C$GR, na.rm = TRUE) * 1.45

# Regression of growth rate on % SWHC for each Site x Year group within each Timepoint
gr_fits <- GR_C %>%
  group_by(Site, Year, Timepoint, grp) %>%
  group_modify(function(d, key) {
    m <- lm(GR ~ SWHC, data = d)
    s <- summary(m)
    tibble(slope = coef(m)[2], R2 = s$r.squared, P = s$coefficients[2, 4], n = nrow(d))
  }) %>%
  ungroup()

# One coloured text line per group, stacked in legend order, top left of each panel
fit_lab <- gr_fits %>%
  mutate(lab = paste0("R2 = ", formatC(R2, format = "f", digits = 2),
                      "   P = ", formatC(P, format = "g", digits = 2)),
         y   = gr_head * (0.975 - 0.06 * (as.integer(grp) - 1)))

Fig1C <- ggplot(GR_C, aes(x = SWHC, y = GR, color = grp, shape = grp, linetype = grp)) +
  geom_point(size = 2.8, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 0.8) +
  geom_text(data = fit_lab, aes(x = 18, y = y, label = lab, color = grp),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4.2, fontface = "bold",
            show.legend = FALSE) +
  facet_wrap(~ Timepoint) +
  scale_color_manual(values = grp_cols,   breaks = grp_levels, name = NULL) +
  scale_shape_manual(values = grp_shapes, breaks = grp_levels, name = NULL) +
  scale_linetype_manual(values = grp_lines, breaks = grp_levels, name = NULL) +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4.2, linewidth = 1.1))) +
  coord_cartesian(ylim = c(0, gr_head)) +
  labs(x = "% Soil Water-Holding Capacity", y = "Growth Rate", tag = "C") +
  fig_theme + tick_theme_BC +
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.97),   # top right of the Late panel; adjust x/y here
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.size = unit(1.5, "lines"),
        legend.key.width = unit(2.4, "lines"),
        legend.spacing.y = unit(2, "pt"))

col_BC <- (Fig1B / Fig1C) + plot_layout(heights = c(1.7, 1.3))

############################################################################################################
#Combine: A on the left, B over C on the right

Fig1 <- (col_A | col_BC) + plot_layout(widths = c(1, 1.15))
Fig1

dir.create("paper_graphs", showWarnings = FALSE)
ggsave("paper_graphs/Fig1.pdf", Fig1, width = 18, height = 11, units = "in", device = "pdf")