# All 14 traits over time, at MATERNAL LINE level
# ------------------------------------------------------------------------------
# Same construction as Mimulus_AllTraits_TimeSeries_byPopulation.R, but one panel
# per maternal line (Population A-K) instead of per site x year group.
#
#   rows    = population group (Oregon above S. California, 2010 above 2014)
#   columns = the maternal lines within that group
#   panel   = all five watering regimes, mean +/- SE, Early and Late on one axis
#   strip   = one-way ANOVA across the five regimes at that timepoint, then Tukey
#             HSD; one colored line per regime against 100% SWHC,
#             labelled by that regime
#
# NOTE ON POWER: there are only 3-5 plants per watering regime per maternal line
# per run, so these Tukey tests are far weaker than the pooled-group version and
# most timepoints will not reach p 0.05. Read the curves, not the strip.
#
# Run from .../collaborations/Anstett_Debbas/20251209_Experiment01
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggsci)      # scale_color_d3
library(cowplot)

REF          <- "100%"
treat_levels <- c("20%", "40%", "60%", "80%", "100%")
run_break    <- (408 + 411) / 2
toe_min      <- 140          # keeps 144 h onwards, the projected-shoot-area trough
outdir       <- "paper_graphs/magda/TimeSeries_byMaternalLine"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

traits <- c("Roi.Mask" = "Projected shoot area (a.u.)", "Fv.Fm" = "Fv/Fm",
            "Fq.Fm" = "Fq'/Fm'", "NPQ" = "NPQ", "qP" = "qP", "qN" = "qN",
            "qL" = "qL", "qI" = "qI", "phiNO" = "phiNO", "phiNPQ" = "phiNPQ",
            "npq.t" = "NPQ(t)", "ChlIdx" = "Chlorophyll index",
            "AriIdx" = "Anthocyanin reflectance index", "NDVI" = "NDVI")
file_tag <- c("Roi.Mask" = "Area", "Fv.Fm" = "FvFm", "Fq.Fm" = "FqFm", "npq.t" = "npqt")

# rows of the figure: group -> its maternal lines (NA = blank cell)
blocks <- list("Oregon 2010"        = c("G","H","I"),
               "Oregon 2014"        = c("J","K", NA),
               "California 2010" = c("A","B","C"),
               "California 2014" = c("D","E","F"))

# --- data ---------------------------------------------------------------------

ps2 <- rbind(read.csv("data/Experiment01_EarlyPS2_Mimulus.csv"),
             read.csv("data/Experiment01_LatePS2_Mimulus.csv"))
ps2 <- ps2[!duplicated(ps2[, c("Plant.ID", "TOE")]), ]
ps2 <- subset(ps2, TOE >= toe_min)
ps2$Treatment <- factor(ps2$Treatment, levels = treat_levels)
ps2$Run <- ifelse(ps2$TOE < run_break, "Exp20243", "Exp20244")

# --- ANOVA + Tukey vs REF, per trait x maternal line x timepoint --------------

tukey_vs_ref <- function(d, trait) {
  d <- d[is.finite(d[[trait]]), ]
  keep <- names(which(table(droplevels(d$Treatment)) >= 3))
  d <- d[d$Treatment %in% keep, ]
  if (!(REF %in% keep) || length(keep) < 3) return(data.frame())
  d$Treatment <- droplevels(factor(d$Treatment, levels = treat_levels))

  fit     <- aov(d[[trait]] ~ d$Treatment)
  p_anova <- summary(fit)[[1]][["Pr(>F)"]][1]
  tk      <- TukeyHSD(fit)[[1]]

  parts     <- strsplit(rownames(tk), "-", fixed = TRUE)
  keep_rows <- vapply(parts, function(x) REF %in% x, logical(1))
  other     <- vapply(parts[keep_rows], function(x) setdiff(x, REF)[1], character(1))

  data.frame(Treatment = other,
             n_ref     = sum(d$Treatment == REF),
             p_tukey   = unname(tk[keep_rows, "p adj"]),
             p_anova   = p_anova)
}

pvals <- do.call(rbind, lapply(names(traits), function(trait) {
  out <- ps2 %>%
    group_by(Population, TOE) %>%
    group_modify(function(d, key) tukey_vs_ref(d, trait)) %>%
    ungroup()
  out$Trait <- trait
  out
}))

pvals$LOD       <- -log10(pmax(pvals$p_tukey, 1e-16))
pvals$Treatment <- factor(pvals$Treatment, levels = treat_levels)
pvals$Run       <- ifelse(pvals$TOE < run_break, "Exp20243", "Exp20244")

write.csv(pvals, "data/Mimulus_AllTraits_vs100pct_pvalues_byPopulation.csv",
          row.names = FALSE)

# --- one time course + p-value strip per maternal line ------------------------

make_pair <- function(pop, gname, trait, ylim_trait, lodmax) {
  d <- subset(ps2,   Population == pop)
  q <- subset(pvals, Population == pop & Trait == trait)

  lgraph <- ggplot(d, aes(x = TOE, y = .data[[trait]], color = Treatment, fill = Treatment)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", aes(group = interaction(Treatment, Run)),
                 linetype = 0, alpha = 0.20) +
    stat_summary(fun = mean, geom = "line", aes(group = interaction(Treatment, Run)),
                 linewidth = 0.6) +
    geom_vline(xintercept = run_break, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_d3("category10") + scale_fill_d3("category10") +
    coord_cartesian(ylim = ylim_trait) +
    labs(x = "", y = unname(traits[trait]),
         title = paste0(gname, " — line ", pop,
                        "  (n = ", length(unique(d$Plant.ID)), ")")) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_blank(),
          plot.title = element_text(size = 10, hjust = 0))

  pplot <- ggplot(q, aes(x = TOE, y = LOD, color = Treatment,
                         group = interaction(Treatment, Run))) +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = run_break, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_d3("category10") +
    coord_cartesian(ylim = c(0, lodmax)) +
    labs(x = "Hours Of Imaging", y = "-log10(p-value)") +
    theme_minimal() + theme_bw() +
    theme(legend.position = "none",
          strip.text = element_blank(), strip.background = element_blank())

  plot_grid(lgraph, pplot, rel_heights = c(4, 1), ncol = 1, align = "v", axis = "l")
}

legend_treat <- get_legend(
  ggplot(ps2, aes(x = TOE, y = Roi.Mask, color = Treatment)) + geom_line() +
    scale_color_d3("category10") + labs(color = "% SWHC") + theme_minimal())

for (trait in names(traits)) {
  v   <- ps2[[trait]][is.finite(ps2[[trait]])]
  rng <- quantile(v, c(0.001, 0.999)); pad <- diff(rng) * 0.06
  ylim_trait <- if (trait == "Roi.Mask") c(0, rng[2] + pad) else c(rng[1] - pad, rng[2] + pad)
  lodmax <- max(pvals$LOD[pvals$Trait == trait & is.finite(pvals$LOD)], na.rm = TRUE) * 1.10

  # a blank placeholder keeps Oregon 2014's third column empty
  cells <- unlist(lapply(names(blocks), function(g)
    lapply(blocks[[g]], function(p)
      if (is.na(p)) ggplot() + theme_void()
      else make_pair(p, g, trait, ylim_trait, lodmax))),
    recursive = FALSE)

  figure <- plot_grid(plot_grid(plotlist = cells, ncol = 3),
                      legend_treat, ncol = 2, rel_widths = c(1, 0.09))

  tag <- if (trait %in% names(file_tag)) file_tag[[trait]] else trait
  pdf(file.path(outdir, paste0("Mimulus_", tag, "_TimeSeries_byMaternalLine.pdf")),
      height = 21, width = 15)
  print(figure)
  dev.off()
}
