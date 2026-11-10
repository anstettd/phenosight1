##################################################################################
## Fig 2: Fv/Fm and Chlorophyll Index over time, with per-treatment
## statistics (MMJ, Nov 2026) — same construction as panel A of
## 1.1_Fig1_PSA_GR_MMJ.R, and the same style as 1_PSA_gr.R
##
##  A  Fv/Fm over time, maternal-line means, one line per SWHC regime, with a
##     -log10(P) row directly beneath each row of curves: each regime vs 100% SWHC
##  B  the same for Chlorophyll Index
##
## Runs off the line-mean table written by 0_data_setup.R: data/PS2_all.csv
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(ggh4x)
library(patchwork)

#The dashed vertical line in every panel: the setup changed here because the
#plants outgrew the trays, so the step at this point is technical, not biological.
RUN_BREAK <- (408 + 411) / 2
REF <- "100%"

#Import Data
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year = factor(Year, levels = c(2010, 2014)),      # 2010 on top row
         Site = factor(Site, levels = c("Oregon", "California")))   # Oregon on left

############################################################################################################
#Shared theme, kept identical to 1.1_Fig1_PSA_GR_MMJ.R so the two figures match.
#If we keep both scripts it is worth pulling this into its own file and sourcing it.

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
  panel.spacing = unit(3, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 8 * TEXT_SCALE, face = "bold"),
  legend.key.size = unit(0.8, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  legend.justification = c(0, 0.5),
  strip.text = element_text(size = 13 * TEXT_SCALE, face = "bold",
                            margin = margin(2, 2, 2, 2)),
  strip.text.y.left = element_text(angle = 90, size = 13 * TEXT_SCALE, face = "bold",
                                   margin = margin(2, 2, 2, 2)),
  strip.placement = "outside",            # row strips sit left of the axis labels
  strip.background = element_blank(),
  plot.margin = margin(t = 2, r = 4, b = 2, l = 4),
  strip.switch.pad.grid = unit(1, "pt")
)

############################################################################################################
#Builds one Figure-1A style block for a given trait: curves per watering regime in
#a Year x Site grid, with a thin -log10(P) row directly beneath each row of curves,
#all in ONE facet grid so each P panel sits under the curves it belongs to.
#
#  trait  column name in PS2_all, e.g. "Fv.Fm"
#  label  strip text, e.g. "Fv/Fm"
#
#Returns a list: $plot (the ggplot) and $stats (per-timepoint Tukey table).

trait_block <- function(dat, trait, label) {

  summ <- dat %>%
    group_by(Site, Year, Treatment, TOE) %>%
    summarise(mean_v = mean(.data[[trait]], na.rm = TRUE),
              se_v = sd(.data[[trait]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[trait]]))),
              .groups = "drop")

  #One-way ANOVA across the five regimes at each TOE within each Site x Year,
  #then Tukey HSD, keeping only the contrasts against 100% SWHC.
  tk_one <- function(d) {
    d <- d[is.finite(d[[trait]]), ]
    keep <- names(which(table(droplevels(d$Treatment)) >= 2))
    d <- d[d$Treatment %in% keep, ]
    if (!(REF %in% keep) || length(keep) < 3) return(tibble())
    d$Treatment <- droplevels(d$Treatment)
    fit <- aov(d[[trait]] ~ d$Treatment)
    tk  <- TukeyHSD(fit)[[1]]
    parts <- strsplit(rownames(tk), "-", fixed = TRUE)
    kr <- vapply(parts, function(x) REF %in% x, logical(1))
    tibble(Treatment = vapply(parts[kr], function(x) setdiff(x, REF)[1], character(1)),
           P = unname(tk[kr, "p adj"]),
           P_anova = summary(fit)[[1]][["Pr(>F)"]][1])
  }

  stats_tbl <- dat %>%
    group_by(Site, Year, TOE) %>%
    group_modify(~ tk_one(.x)) %>%
    ungroup() %>%
    mutate(Treatment = factor(Treatment, levels = levels(dat$Treatment)),
           negLogP = -log10(P))

  row_levels <- c("2010  trait", "2010  -log(P)", "2014  trait", "2014  -log(P)")
  row_labels <- setNames(c(paste0("2010  ", label), "-log P",
                           paste0("2014  ", label), "-log P"), row_levels)

  curves <- summ %>%
    transmute(Site, Year, TOE, Treatment, value = mean_v,
              lo = mean_v - se_v, hi = mean_v + se_v, Panel = "trait")
  pvals <- stats_tbl %>%
    transmute(Site, Year, TOE, Treatment, value = negLogP,
              lo = NA_real_, hi = NA_real_, Panel = "-log(P)")

  d_all <- bind_rows(curves, pvals) %>%
    mutate(RowFacet = factor(paste0(Year, "  ", Panel), levels = row_levels))
  d_rib <- d_all %>% filter(Panel == "trait")

  #One shared y range across both trait panels, and another across both -log(P)
  #panels, so every panel of the same kind is directly comparable.
  span <- diff(range(c(d_rib$lo, d_rib$hi), na.rm = TRUE))
  rng <- bind_rows(
    tibble(Panel = "trait",
           lo = min(d_rib$lo, na.rm = TRUE) - 0.03 * span,
           hi = max(d_rib$hi, na.rm = TRUE) + 0.03 * span),
    tibble(Panel = "-log(P)", lo = 0,
           hi = max(d_all$value[d_all$Panel == "-log(P)"], na.rm = TRUE) * 1.05))

  d_blank <- d_all %>%
    distinct(Site, RowFacet, Panel) %>%
    left_join(rng, by = "Panel") %>%
    pivot_longer(c(lo, hi), names_to = "end", values_to = "value") %>%
    mutate(TOE = min(d_all$TOE, na.rm = TRUE))

  d_hline <- d_all %>%
    filter(Panel == "-log(P)") %>%
    distinct(Site, RowFacet) %>%
    mutate(y = -log10(0.05))

  p <- ggplot(d_all, aes(x = TOE, y = value,
                         color = Treatment, fill = Treatment, group = Treatment)) +
    geom_blank(data = d_blank, aes(x = TOE, y = value), inherit.aes = FALSE) +
    geom_ribbon(data = d_rib, aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
    geom_hline(data = d_hline, aes(yintercept = y), inherit.aes = FALSE,
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
    theme(strip.text.y.left = element_text(angle = 90, size = 11 * TEXT_SCALE,
                                           face = "bold", margin = margin(2, 2, 2, 2)),
          legend.position = "right")

  list(plot = p, stats = stats_tbl)
}

############################################################################################################
# Fig 2A — Fv/Fm      Fig 2B — Chlorophyll Index

A <- trait_block(PS2_all, "Fv.Fm",  "Fv/Fm (a.u.)")
B <- trait_block(PS2_all, "ChlIdx", "Chlorophyll Index (a.u.)")

write_csv(A$stats, "data/Fig2A_FvFm_vs100_pvalues.csv")
write_csv(B$stats, "data/Fig2B_ChlIdx_vs100_pvalues.csv")

#How many timepoints reach P < 0.05 for each regime, as a quick sanity check
cat("\nFv/Fm: timepoints with Tukey P < 0.05 vs 100% SWHC\n")
print(A$stats %>% filter(P < 0.05) %>% count(Site, Year, Treatment))
cat("\nChlIdx: timepoints with Tukey P < 0.05 vs 100% SWHC\n")
print(B$stats %>% filter(P < 0.05) %>% count(Site, Year, Treatment))

Fig2 <- A$plot / B$plot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

#NB: do not print Fig2 after this. Under Rscript a bare `Fig2` opens the default
#device and dumps a second, 7x7 in copy as Rplots.pdf in the project root.
ggsave("paper_graphs/Figure 2.pdf", Fig2, width = 13.5, height = 17)
