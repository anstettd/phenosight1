# All 14 traits over time, Early + Late combined, split by population group
# ------------------------------------------------------------------------------
# One figure per trait. In each: rows = collection site (Oregon on top,
# S. California below), columns = collection year (2010, 2014); all five watering
# regimes as mean +/- SE. Under each panel the usual -log10(p) strip, from a
# one-way ANOVA across the five regimes at that timepoint followed by Tukey HSD,
# with one colored line per regime compared against 100% SWHC.
# Each line is labelled by the regime it represents, not by "100% vs ...".
#
# Writes <trait> PDFs plus one combined p-value table.
# Run from .../collaborations/Anstett_Debbas/20251209_Experiment01
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggsci)      # scale_color_d3
library(cowplot)

REF       <- "100%"                            # reference watering regime
treat_levels <- c("20%", "40%", "60%", "80%", "100%")
run_break <- (408 + 411) / 2                   # Exp 20243 ends 408 h, Exp 20244 starts 411 h
toe_min   <- 140     # keeps 144 h onwards: the projected-shoot-area trough in all
                     # four groups, after which the splines are well behaved
outdir    <- "paper_graphs/magda/TimeSeries_byGroup"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

traits <- c("Roi.Mask" = "Projected shoot area (a.u.)",
            "Fv.Fm"    = "Fv/Fm",
            "Fq.Fm"    = "Fq'/Fm'",
            "NPQ"      = "NPQ",
            "qP"       = "qP",
            "qN"       = "qN",
            "qL"       = "qL",
            "qI"       = "qI",
            "phiNO"    = "phiNO",
            "phiNPQ"   = "phiNPQ",
            "npq.t"    = "NPQ(t)",
            "ChlIdx"   = "Chlorophyll index",
            "AriIdx"   = "Anthocyanin reflectance index",
            "NDVI"     = "NDVI")
file_tag <- c("Roi.Mask" = "Area", "Fv.Fm" = "FvFm", "Fq.Fm" = "FqFm", "npq.t" = "npqt")

# --- data ---------------------------------------------------------------------

decode <- read.csv("data/decode.csv", fileEncoding = "UTF-8-BOM")

ps2 <- rbind(read.csv("data/Experiment01_EarlyPS2_Mimulus.csv"),
             read.csv("data/Experiment01_LatePS2_Mimulus.csv"))
ps2 <- ps2[!duplicated(ps2[, c("Plant.ID", "TOE")]), ]   # `comp` repeats the control rows
ps2 <- subset(ps2, TOE >= toe_min)

ps2$Site <- decode$Site[match(ps2$Population, decode$Population)]
ps2$Year <- decode$Year[match(ps2$Population, decode$Population)]
ps2$Treatment <- factor(ps2$Treatment, levels = treat_levels)
ps2$Run   <- ifelse(ps2$TOE < run_break, "Exp20243", "Exp20244")
ps2$Group <- paste(ps2$Site, ps2$Year)

group_order <- c("Oregon 2010", "Oregon 2014", "SCalifornia 2010", "SCalifornia 2014")
group_title <- c("Oregon 2010" = "Oregon 2010", "Oregon 2014" = "Oregon 2014",
                 "SCalifornia 2010" = "California 2010",
                 "SCalifornia 2014" = "California 2014")

# --- ANOVA + Tukey vs REF, per trait x group x timepoint ----------------------

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
    group_by(Group, TOE) %>%
    group_modify(function(d, key) tukey_vs_ref(d, trait)) %>%
    ungroup()
  out$Trait <- trait
  out
}))

pvals$LOD       <- -log10(pmax(pvals$p_tukey, 1e-16))
pvals$Treatment <- factor(pvals$Treatment, levels = treat_levels)
pvals$Run       <- ifelse(pvals$TOE < run_break, "Exp20243", "Exp20244")

write.csv(pvals, "data/Mimulus_AllTraits_vs100pct_pvalues.csv",
          row.names = FALSE)

# --- one time course + p-value strip per group, one figure per trait ----------

make_pair <- function(g, trait, ylim_trait, lodmax) {
  d <- subset(ps2,   Group == g)
  q <- subset(pvals, Group == g & Trait == trait)

  lgraph <- ggplot(d, aes(x = TOE, y = .data[[trait]], color = Treatment, fill = Treatment)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", aes(group = interaction(Treatment, Run)),
                 linetype = 0, alpha = 0.22) +
    stat_summary(fun = mean, geom = "line", aes(group = interaction(Treatment, Run)),
                 linewidth = 0.7) +
    geom_vline(xintercept = run_break, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_d3("category10") + scale_fill_d3("category10") +
    coord_cartesian(ylim = ylim_trait) +
    labs(x = "", y = unname(traits[trait]), title = unname(group_title[g])) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_blank(),
          plot.title = element_text(size = 11, hjust = 0))

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
  rng <- quantile(v, c(0.001, 0.999))
  pad <- diff(rng) * 0.06
  ylim_trait <- if (trait == "Roi.Mask") c(0, rng[2] + pad) else c(rng[1] - pad, rng[2] + pad)
  lodmax <- max(pvals$LOD[pvals$Trait == trait & is.finite(pvals$LOD)], na.rm = TRUE) * 1.10

  pairs  <- lapply(group_order, make_pair, trait = trait,
                   ylim_trait = ylim_trait, lodmax = lodmax)
  figure <- plot_grid(plot_grid(plotlist = pairs, ncol = 2),
                      legend_treat, ncol = 2, rel_widths = c(1, 0.13))

  tag <- if (trait %in% names(file_tag)) file_tag[[trait]] else trait
  pdf(file.path(outdir, paste0("Mimulus_", tag, "_TimeSeries_byPopulation.pdf")),
      height = 11, width = 15)
  print(figure)
  dev.off()
}
