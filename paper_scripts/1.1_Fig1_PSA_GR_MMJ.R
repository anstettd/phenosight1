##################################################################################
## Fig 1: Projected Shoot Area + Growth Rate, with per-treatment statistics
## (MMJ, Nov 2026) — follows the style of 1_PSA_gr.R
##
##  A  PSA over time, maternal-line means, one line per SWHC regime,
##     with a -log10(P) panel beneath: each SWHC regime against 100% SWHC
##  B  Growth rate per interval, with the ANOVA P and Tukey HSD letters
##  C  Growth rate regressed on % SWHC, per population group
##
## Everything runs off the line-mean tables written by 0_data_setup.R:
##   data/PS2_all.csv, data/GR_all.csv
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(ggh4x)
library(patchwork)
library(multcompView)

#The dashed vertical line in every panel: the setup changed here because the
#plants outgrew the trays, so the step at this point is technical, not biological.
RUN_BREAK <- (408 + 411) / 2

#Import Data
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year = factor(Year, levels = c(2010, 2014)),      # 2010 on top row
         Site = factor(Site, levels = c("Oregon", "California")))   # Oregon on left

GR_all <- read_csv("data/GR_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California")),
         SWHC = as.numeric(gsub("%", "", Treatment)))

############################################################################################################
#Shared theme, copied from 1_PSA_gr.R so the panels match. If we keep both scripts
#it would be worth pulling this into its own file and sourcing it from each.

#Same look as 1_PSA_gr.R but with the type scaled down. His figure has 4 panels at
#6 x 6 in; this one has ~18 panels, so his absolute point sizes swamp the plots.
#Change TEXT_SCALE alone to move all of it together.
TEXT_SCALE <- 1.00

fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 11 * TEXT_SCALE, face = "bold"),
  axis.text.x = element_text(size = 11 * TEXT_SCALE, face = "bold", margin = margin(t = 1)),
  axis.text.y = element_text(size = 11 * TEXT_SCALE, face = "bold", margin = margin(r = 1)),
  axis.title = element_text(color = "black", size = 12 * TEXT_SCALE, face = "bold"),
  axis.title.x = element_text(color = "black", size = 12 * TEXT_SCALE, face = "bold",
                              margin = margin(t = 1)),
  axis.title.y = element_text(color = "black", size = 12 * TEXT_SCALE, face = "bold",
                              margin = margin(r = 1)),
  axis.ticks.length = unit(2, "pt"),
  panel.spacing = unit(3, "pt"),          # default is 5.5pt
  legend.title = element_blank(),
  legend.text = element_text(size = 8 * TEXT_SCALE, face = "bold"),
  legend.key.size = unit(0.8, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  legend.justification = c(0, 0.5),
  strip.text = element_text(size = 13 * TEXT_SCALE, face = "bold",
                            margin = margin(2, 2, 2, 2)),
  strip.text.y.left = element_text(angle = 90, size = 13 * TEXT_SCALE, face = "bold",
                                   margin = margin(1, 1, 1, 1)),
  strip.placement = "outside",            # row strips sit left of the axis labels
  strip.background = element_blank(),
  plot.margin = margin(t = 2, r = 4, b = 2, l = 4),  # legends sit in the guide area now
  strip.switch.pad.grid = unit(1, "pt")
)

############################################################################################################
# Fig 1A — PSA over time + P-values of each SWHC regime against 100%

#Top: mean +/- SE across maternal lines
PSA_summ <- PS2_all %>%
  group_by(Site, Year, Treatment, TOE) %>%
  summarise(mean_PSA = mean(PSA, na.rm = TRUE),
            se_PSA = sd(PSA, na.rm = TRUE) / sqrt(sum(!is.na(PSA))),
            n_lines = sum(!is.na(PSA)),
            .groups = "drop")

#(the curves and their P panel are drawn together further down, so that each
# P panel sits directly beneath the curves it belongs to)

#One-way ANOVA across the five regimes at each TOE within each Site x Year,
#then Tukey HSD, keeping only the contrasts against 100% SWHC.
tukey_vs_ref <- function(d, response, ref = "100%") {
  d <- d[is.finite(d[[response]]), ]
  keep <- names(which(table(droplevels(d$Treatment)) >= 2))
  d <- d[d$Treatment %in% keep, ]
  if (!(ref %in% keep) || length(keep) < 3) return(tibble())
  d$Treatment <- droplevels(d$Treatment)

  fit <- aov(d[[response]] ~ d$Treatment)
  tk  <- TukeyHSD(fit)[[1]]

  parts <- strsplit(rownames(tk), "-", fixed = TRUE)
  keep_rows <- vapply(parts, function(x) ref %in% x, logical(1))
  tibble(Treatment = vapply(parts[keep_rows], function(x) setdiff(x, ref)[1], character(1)),
         P = unname(tk[keep_rows, "p adj"]),
         P_anova = summary(fit)[[1]][["Pr(>F)"]][1])
}

toe_stats <- PS2_all %>%
  group_by(Site, Year, TOE) %>%
  group_modify(~ tukey_vs_ref(.x, "PSA")) %>%
  ungroup() %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         negLogP = -log10(P))

write_csv(toe_stats, "data/Fig1A_PSA_vs100_pvalues.csv")

#Curves and P values in ONE facet grid, rows interleaved as
#   2010 PSA / 2010 -log(P) / 2014 PSA / 2014 -log(P)
#so each P panel sits immediately under the curves it describes. independent = "y"
#lets the PSA and P rows carry their own scales; force_panelsizes keeps the P rows
#thin, the way they looked when they were a separate strip.
PSA_long <- PSA_summ %>%
  transmute(Site, Year, TOE, Treatment,
            value = mean_PSA, lo = mean_PSA - se_PSA, hi = mean_PSA + se_PSA,
            Panel = "PSA")

P_long <- toe_stats %>%
  transmute(Site, Year, TOE, Treatment, value = negLogP,
            lo = NA_real_, hi = NA_real_, Panel = "-log(P)")

row_levels <- c("2010  PSA", "2010  -log(P)", "2014  PSA", "2014  -log(P)")

#The P rows are thin, and the right-hand strip text is rotated, so the full
#"2010  -log(P)" gets clipped. The year is already on the PSA row directly above.
row_labels <- c("2010  PSA" = "2010  PSA (a.u.)", "2010  -log(P)" = "-log P",
                "2014  PSA" = "2014  PSA (a.u.)", "2014  -log(P)" = "-log P")

A_dat <- bind_rows(PSA_long, P_long) %>%
  mutate(RowFacet = factor(paste0(Year, "  ", Panel), levels = row_levels))

A_ribbon <- A_dat %>% filter(Panel == "PSA")

#Force one shared y range across all four PSA panels, and another across both
#-log(P) panels: an invisible geom_blank layer carrying the global min/max of each
#panel type. With scales = "free_y" (and NOT independent) the two PSA rows then
#resolve to the same limits, and so do the two P rows.
A_range <- bind_rows(
  A_ribbon %>% summarise(Panel = "PSA", lo = 0, hi = max(hi, na.rm = TRUE) * 1.02),
  A_dat %>% filter(Panel == "-log(P)") %>%
    summarise(Panel = "-log(P)", lo = 0, hi = max(value, na.rm = TRUE) * 1.05))

A_blank <- A_dat %>%
  distinct(Site, RowFacet, Panel) %>%
  left_join(A_range, by = "Panel") %>%
  pivot_longer(c(lo, hi), names_to = "end", values_to = "value") %>%
  mutate(TOE = min(A_dat$TOE, na.rm = TRUE))

#P = 0.05 reference, only in the two P rows
A_hline <- A_dat %>%
  filter(Panel == "-log(P)") %>%
  distinct(Site, RowFacet) %>%
  mutate(y = -log10(0.05))

Fig1A <- ggplot(A_dat, aes(x = TOE, y = value,
                           color = Treatment, fill = Treatment, group = Treatment)) +
  geom_blank(data = A_blank, aes(x = TOE, y = value), inherit.aes = FALSE) +
  geom_ribbon(data = A_ribbon, aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_hline(data = A_hline, aes(yintercept = y), inherit.aes = FALSE,
             color = "grey50", linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = RUN_BREAK, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  facet_grid2(RowFacet ~ Site, scales = "free_y",
              labeller = labeller(RowFacet = row_labels), switch = "y") +
  force_panelsizes(rows = c(4, 1.3, 4, 1.3)) +
  scale_color_viridis_d(drop = FALSE) +
  scale_fill_viridis_d(drop = FALSE) +
  labs(x = "Time of Experiment", y = NULL) +
  fig_theme +
  theme(strip.text.y.left = element_text(angle = 90, size = 11 * TEXT_SCALE, face = "bold",
                                         margin = margin(1, 1, 1, 1)),
        legend.position = "right")

############################################################################################################
# Fig 1B — Growth rate per interval, with ANOVA P and Tukey letters

#Compact letter display, using the same multcompLetters route as the other notebooks
gr_letters <- GR_all %>%
  group_by(Site, Year, Timepoint) %>%
  group_modify(function(d, key) {
    d$Treatment <- droplevels(d$Treatment)
    if (nlevels(d$Treatment) < 3 || any(table(d$Treatment) < 2))
      return(tibble(Treatment = levels(d$Treatment), Letters = "", P_anova = NA_real_))
    fit <- aov(GR ~ Treatment, data = d)
    P   <- summary(fit)[[1]][["Pr(>F)"]][1]
    tk  <- TukeyHSD(fit)$Treatment[, "p adj"]
    L   <- multcompLetters(tk)$Letters
    tibble(Treatment = names(L), Letters = unname(L), P_anova = P)
  }) %>%
  ungroup() %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")))

write_csv(gr_letters, "data/Fig1B_GR_ANOVA_Tukey_letters.csv")

#Position the letters just above the top of each treatment's spread
GR_summ <- GR_all %>%
  group_by(Site, Year, Timepoint, Treatment) %>%
  summarise(mean_GR = mean(GR, na.rm = TRUE),
            se_GR = sd(GR, na.rm = TRUE) / sqrt(sum(!is.na(GR))),
            top = max(GR, na.rm = TRUE), .groups = "drop") %>%
  left_join(gr_letters, by = c("Site", "Year", "Timepoint", "Treatment"))

gr_head <- max(GR_all$GR, na.rm = TRUE) * 1.32   # headroom for the letters and the ANOVA P

#One ANOVA P per panel, printed top left
anova_lab <- gr_letters %>%
  distinct(Site, Year, Timepoint, P_anova) %>%
  mutate(lab = ifelse(is.na(P_anova), "",
                      paste0("ANOVA P = ", formatC(P_anova, format = "g", digits = 2))))

Fig1B <- ggplot(GR_summ, aes(x = Treatment, y = mean_GR, color = Treatment)) +
  geom_jitter(data = GR_all, aes(x = Treatment, y = GR, color = Treatment),
              width = 0.12, size = 1, alpha = 0.45, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_GR - se_GR, ymax = mean_GR + se_GR), width = 0.18, linewidth = 0.7) +
  geom_point(size = 2.4) +
  #Tukey letters: black, and all on one line well clear of the points rather than
  #each riding on its own treatment's maximum
  geom_text(aes(label = Letters), y = gr_head * 0.84, color = "black", size = 4.2,
            fontface = "bold", show.legend = FALSE) +
  geom_text(data = anova_lab, aes(x = 0.6, y = gr_head * 0.99, label = lab),
            inherit.aes = FALSE, hjust = 0, size = 3.4, fontface = "bold", color = "grey25") +
  facet_nested(Timepoint ~ Site + Year, switch = "y") +
  scale_color_viridis_d() +
  coord_cartesian(ylim = c(0, gr_head)) +
  labs(x = "% SWHC", y = expression(bold(Growth~rate~(a.u.~h^-1)))) +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8 * TEXT_SCALE, face = "bold", angle = 45, hjust = 1))

############################################################################################################
# Fig 1C — Growth rate regressed on % SWHC, per population group

gr_fits <- GR_all %>%
  group_by(Site, Year, Timepoint) %>%
  group_modify(function(d, key) {
    m <- lm(GR ~ SWHC, data = d)
    s <- summary(m)
    tibble(slope = coef(m)[2], R2 = s$r.squared, P = s$coefficients[2, 4], n = nrow(d))
  }) %>%
  ungroup()

write_csv(gr_fits, "data/Fig1C_GR_vs_SWHC_regression.csv")

#Solid where the regression is significant, dotted where it is not
gr_plotdat <- GR_all %>%
  left_join(gr_fits %>% select(Site, Year, Timepoint, P, R2), by = c("Site", "Year", "Timepoint")) %>%
  mutate(Sig = ifelse(P < 0.05, "P < 0.05", "n.s."))

fit_lab <- gr_fits %>%
  mutate(lab = paste0(Site, " ", Year, "   R2 = ", formatC(R2, format = "f", digits = 2),
                      "   P = ", formatC(P, format = "g", digits = 2))) %>%
  group_by(Timepoint) %>%
  summarise(lab = paste(lab, collapse = "\n"), .groups = "drop")

Fig1C <- ggplot(gr_plotdat, aes(x = SWHC, y = GR, color = Site, linetype = Year)) +
  geom_point(aes(shape = Year), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 0.55) +
  geom_text(data = fit_lab, aes(x = 18, y = gr_head * 0.98, label = lab),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2, color = "grey25") +
  facet_wrap(~ Timepoint) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  scale_shape_manual(values = c("2010" = 16, "2014" = 17)) +
  guides(linetype = guide_legend(override.aes = list(color = "grey30")),
         shape    = guide_legend(override.aes = list(color = "grey30"))) +
  coord_cartesian(ylim = c(0, gr_head)) +
  labs(x = "% SWHC", y = expression(bold(Growth~rate~(a.u.~h^-1)))) +
  fig_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 8 * TEXT_SCALE, face = "bold"))

############################################################################################################
#Assemble

#Fig1A is now a single ggplot (curves and P rows interleaved in one facet grid),
#so this is a plain three-panel stack.
Fig1 <- Fig1A / Fig1B / Fig1C +
  plot_layout(heights = c(5.8, 3, 2.6)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

#NB: do not print Fig1 after this. Under Rscript a bare `Fig1` opens the default
#device and dumps a second, 7x7 in copy as Rplots.pdf in the project root, which
#is cramped and easy to mistake for the real figure.
ggsave("paper_graphs/Figure 1.pdf", Fig1, width = 13.5, height = 18)
