##################################################################################
## Supplemental Fig S1: Projected Shoot Area + Growth Rate, with per-treatment
## statistics (MMJ, Nov 2026) — follows the style of 1_PSA_gr.R
##
##  A  PSA over time, one line per SWHC regime, with a -log10(P) panel beneath:
##     each SWHC regime against 100% SWHC
##  B  Growth rate per interval, with the ANOVA P and Tukey HSD letters
##  C  Growth rate regressed on % SWHC, per population group
##
## The script writes the figure TWICE, once per replication unit:
##   _lineMeans  the replicate is the maternal line  (data/PS2_all.csv, GR_all.csv)
##   _plants     the replicate is the individual plant
##                                 (data/PS2_all_plants.csv, GR_all_plants.csv)
## Same figure, same statistics, different n. The line-level version is the
## conservative one; the plant-level version has 4-5x the replication and will
## show far more significance, but the plants within a maternal line are not
## independent of one another, so it overstates the evidence. Both are in the
## supplement so the reader can see which conclusions survive the conservative
## treatment.
##
## All four input tables are written by 0_data_setup.R.
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

OUT_DIR <- "paper_graphs/Supplemental/magda"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#The two replication units. Add or drop entries here to change what gets built.
UNITS <- list(
  list(tag   = "lineMeans",
       ps2   = "data/PS2_all.csv",
       gr    = "data/GR_all.csv",
       label = "Replicate = maternal line (3 lines per group; 2 for Oregon 2014)",
       pt    = 1.5, alpha = 0.70, jit = 1.0),
  list(tag   = "plants",
       ps2   = "data/PS2_all_plants.csv",
       gr    = "data/GR_all_plants.csv",
       label = "Replicate = individual plant (3-5 plants per line per regime)",
       pt    = 0.8, alpha = 0.40, jit = 0.6)
)

read_unit <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
           Timepoint = factor(Timepoint, levels = c("Early", "Late")),
           Year = factor(Year, levels = c(2010, 2014)),      # 2010 on top row
           Site = factor(Site, levels = c("Oregon", "California")))   # Oregon on left
}

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

#One-way ANOVA across the five regimes at each TOE within each Site x Year,
#then Tukey HSD, keeping only the contrasts against 100% SWHC. Works unchanged
#whether the rows are line means or plants - only n changes.
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

############################################################################################################
#Builds the whole of Figure S1 for one replication unit and writes it out.

build_fig1 <- function(u) {

  PS2 <- read_unit(u$ps2)
  GR  <- read_unit(u$gr) %>% mutate(SWHC = as.numeric(gsub("%", "", Treatment)))

  cat("\n================ ", u$tag, " ================\n", sep = "")
  cat("PS2 rows:", nrow(PS2), " GR rows:", nrow(GR), "\n")
  cat("replicates per Site x Year x Treatment x TOE (PS2):\n")
  print(PS2 %>% count(Site, Year, Treatment, TOE) %>% count(Site, Year, n, name = "cells"))

  ##########################################################################################
  # S1A — PSA over time + P-values of each SWHC regime against 100%

  PSA_summ <- PS2 %>%
    group_by(Site, Year, Treatment, TOE) %>%
    summarise(mean_PSA = mean(PSA, na.rm = TRUE),
              se_PSA = sd(PSA, na.rm = TRUE) / sqrt(sum(!is.na(PSA))),
              n_reps = sum(!is.na(PSA)),
              .groups = "drop")

  toe_stats <- PS2 %>%
    group_by(Site, Year, TOE) %>%
    group_modify(~ tukey_vs_ref(.x, "PSA")) %>%
    ungroup() %>%
    mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
           negLogP = -log10(P))

  write_csv(toe_stats, paste0("data/FigS1A_PSA_vs100_pvalues_", u$tag, ".csv"))

  #Curves and P values in ONE facet grid, rows interleaved as
  #   2010 PSA / 2010 -log(P) / 2014 PSA / 2014 -log(P)
  #so each P panel sits immediately under the curves it describes. force_panelsizes
  #keeps the P rows thin, the way they looked when they were a separate strip.
  PSA_long <- PSA_summ %>%
    transmute(Site, Year, TOE, Treatment,
              value = mean_PSA, lo = mean_PSA - se_PSA, hi = mean_PSA + se_PSA,
              Panel = "PSA")

  P_long <- toe_stats %>%
    transmute(Site, Year, TOE, Treatment, value = negLogP,
              lo = NA_real_, hi = NA_real_, Panel = "-log(P)")

  row_levels <- c("2010  PSA", "2010  -log(P)", "2014  PSA", "2014  -log(P)")

  #The P rows are thin and the strip text is rotated, so the full "2010  -log(P)"
  #gets clipped. The year is already on the PSA row directly above.
  row_labels <- c("2010  PSA" = "2010  PSA (a.u.)", "2010  -log(P)" = "-log P",
                  "2014  PSA" = "2014  PSA (a.u.)", "2014  -log(P)" = "-log P")

  A_dat <- bind_rows(PSA_long, P_long) %>%
    mutate(RowFacet = factor(paste0(Year, "  ", Panel), levels = row_levels))

  A_ribbon <- A_dat %>% filter(Panel == "PSA")

  #Force one shared y range across all four PSA panels, and another across both
  #-log(P) panels: an invisible geom_blank layer carrying the global min/max of
  #each panel type. With scales = "free_y" (and NOT independent) the two PSA rows
  #then resolve to the same limits, and so do the two P rows.
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

  ##########################################################################################
  # S1B — Growth rate per interval, with ANOVA P and Tukey letters

  gr_letters <- GR %>%
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

  write_csv(gr_letters, paste0("data/FigS1B_GR_ANOVA_Tukey_letters_", u$tag, ".csv"))

  GR_summ <- GR %>%
    group_by(Site, Year, Timepoint, Treatment) %>%
    summarise(mean_GR = mean(GR, na.rm = TRUE),
              se_GR = sd(GR, na.rm = TRUE) / sqrt(sum(!is.na(GR))),
              top = max(GR, na.rm = TRUE), .groups = "drop") %>%
    left_join(gr_letters, by = c("Site", "Year", "Timepoint", "Treatment"))

  gr_head <- max(GR$GR, na.rm = TRUE) * 1.32   # headroom for the letters and the ANOVA P

  #One ANOVA P per panel, printed top left
  anova_lab <- gr_letters %>%
    distinct(Site, Year, Timepoint, P_anova) %>%
    mutate(lab = ifelse(is.na(P_anova), "",
                        paste0("ANOVA P = ", formatC(P_anova, format = "g", digits = 2))))

  Fig1B <- ggplot(GR_summ, aes(x = Treatment, y = mean_GR, color = Treatment)) +
    geom_jitter(data = GR, aes(x = Treatment, y = GR, color = Treatment),
                width = 0.12, size = u$pt * 0.7, alpha = u$alpha * 0.65, show.legend = FALSE) +
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

  ##########################################################################################
  # S1C — Growth rate regressed on % SWHC, per population group

  gr_fits <- GR %>%
    group_by(Site, Year, Timepoint) %>%
    group_modify(function(d, key) {
      m <- lm(GR ~ SWHC, data = d)
      s <- summary(m)
      tibble(slope = coef(m)[2], R2 = s$r.squared, P = s$coefficients[2, 4], n = nrow(d))
    }) %>%
    ungroup()

  write_csv(gr_fits, paste0("data/FigS1C_GR_vs_SWHC_regression_", u$tag, ".csv"))

  gr_plotdat <- GR %>%
    left_join(gr_fits %>% select(Site, Year, Timepoint, P, R2), by = c("Site", "Year", "Timepoint"))

  fit_lab <- gr_fits %>%
    mutate(lab = paste0(Site, " ", Year, "   R2 = ", formatC(R2, format = "f", digits = 2),
                        "   P = ", formatC(P, format = "g", digits = 2))) %>%
    group_by(Timepoint) %>%
    summarise(lab = paste(lab, collapse = "\n"), .groups = "drop")

  Fig1C <- ggplot(gr_plotdat, aes(x = SWHC, y = GR, color = Site, linetype = Year)) +
    geom_point(aes(shape = Year), size = u$pt, alpha = u$alpha) +
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

  ##########################################################################################
  #Assemble. The title names the replication unit, so the two versions of this
  #figure cannot be confused for one another.

  Fig1 <- Fig1A / Fig1B / Fig1C +
    plot_layout(heights = c(5.8, 3, 2.6)) +
    plot_annotation(title = u$label,
                    tag_levels = list(c("A", "B", "C"))) &
    theme(plot.tag = element_text(size = 18, face = "bold"))

  out <- file.path(OUT_DIR, paste0("FigS1_PSA_GR_", u$tag, ".pdf"))

  #NB: do not print Fig1 after this. Under Rscript a bare `Fig1` opens the default
  #device and dumps a second, 7x7 in copy as Rplots.pdf in the project root, which
  #is cramped and easy to mistake for the real figure.
  ggsave(out, Fig1, width = 13.5, height = 18)
  cat("wrote", out, "\n")

  cat("timepoints with Tukey P < 0.05 vs 100% SWHC (PSA):\n")
  print(toe_stats %>% filter(P < 0.05) %>% count(Site, Year, Treatment))

  invisible(NULL)
}

for (u in UNITS) build_fig1(u)
