# Mimulus tree-ring dose-response figure, at MATERNAL LINE level
# ------------------------------------------------------------------------------
# Same construction as Mimulus_TreeRing_DoseResponse.R, but the dose regression is
# fitted within each maternal line (Population A-K) instead of pooling the lines
# of a site x year group.
#
#   radius = TOE (center = 0 h, edge = 747 h)
#   angle  = trait, 14 in the same order inside every line's wedge
#   fill   = -log10(p) of the dose slope for that line x trait x timepoint
#   stars  = FDR < 0.05, BH within line x trait across timepoints
#
# Layout: the four site x year quadrants are kept in the same corners as the
# pooled figure (2014 up, 2010 down, S. California right, Oregon left). Each
# quadrant is three maternal-line wedges wide; Oregon 2014 has only two lines, so
# its third wedge is left blank and every quadrant still spans exactly 90 degrees.
# Trait names are drawn only on the first wedge of each quadrant - at 180 slots
# there is no room to label all 154.
#
# NOTE ON POWER: each line x timepoint has only about 20 plants across the five
# watering regimes, so these fits are far noisier than the pooled-group version.
#
# Run from .../collaborations/Anstett_Debbas/20251209_Experiment01
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

trait_order <- c("Roi.Mask","Fv.Fm","Fq.Fm","NPQ","qP","qN","qL","qI",
                 "phiNO","phiNPQ","npq.t","ChlIdx","AriIdx","NDVI")
trait_label <- setNames(trait_order, trait_order)
trait_label["Roi.Mask"] <- "Digital biomass"

# quadrant -> its maternal lines, clockwise from the top; NA is the blank wedge
quadrants <- list("California 2014" = c("D","E","F"),
                  "California 2010" = c("A","B","C"),
                  "Oregon 2010"        = c("G","H","I"),
                  "Oregon 2014"        = c("J","K", NA))

# --- data ---------------------------------------------------------------------

ps2 <- rbind(read.csv("data/Experiment01_EarlyPS2_Mimulus.csv"),
             read.csv("data/Experiment01_LatePS2_Mimulus.csv"))
ps2 <- ps2[!duplicated(ps2[, c("Plant.ID", "TOE")]), ]   # `comp` repeats the control rows
ps2$StressDose <- 100 - as.numeric(sub("%", "", ps2$Treatment))

# --- dose-response at every timepoint, within each maternal line ---------------

dose_response_by_timepoint <- function(df, traits, line_label) {
  toes <- sort(unique(df$TOE))
  out  <- vector("list", length(toes) * length(traits)); i <- 0
  for (t in toes) {
    d_t <- df[df$TOE == t, ]
    for (tr in traits) {
      i <- i + 1
      d_use <- d_t[is.finite(d_t[[tr]]) & is.finite(d_t$StressDose), c(tr, "StressDose")]
      if (nrow(d_use) < 6 || length(unique(d_use$StressDose)) < 3) {
        out[[i]] <- data.frame(Population = line_label, Trait = tr, TOE = t,
                               slope = NA_real_, p_value = NA_real_, r2 = NA_real_,
                               n = nrow(d_use))
        next
      }
      fit <- lm(d_use[[tr]] ~ d_use$StressDose); s <- summary(fit)
      out[[i]] <- data.frame(Population = line_label, Trait = tr, TOE = t,
                             slope = coef(fit)[2], p_value = s$coefficients[2, 4],
                             r2 = s$r.squared, n = nrow(d_use))
    }
  }
  do.call(rbind, out)
}

pop_order <- unlist(quadrants, use.names = FALSE)
pop_order <- pop_order[!is.na(pop_order)]

dose_response_all <- do.call(rbind, lapply(pop_order, function(p)
  dose_response_by_timepoint(subset(ps2, Population == p), trait_order, p)))

dose_response_all$LOD <- -log10(dose_response_all$p_value)
dose_response_all <- dose_response_all %>%
  group_by(Population, Trait) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(dose_response_all, "data/Mimulus_dose_response_byPopulation.csv", row.names = FALSE)

# --- tree-ring layout ---------------------------------------------------------

n_traits  <- length(trait_order)
n_gap     <- 3
pops_per_q <- 3
block     <- pops_per_q * n_traits + n_gap        # 45 slots per quadrant
gap_id    <- paste0("__gap", seq_len(n_gap))

slot_levels <- unlist(lapply(names(quadrants), function(g) {
  pops <- quadrants[[g]]
  c(unlist(lapply(seq_along(pops), function(i) {
      p <- pops[i]
      if (is.na(p)) paste0(g, "__blank", i, ":", trait_order)   # blank wedge
      else paste(p, trait_order, sep = ": ")
    })),
    paste0(g, gap_id))
}))

dose_response_all$Slot <- factor(
  paste(dose_response_all$Population, dose_response_all$Trait, sep = ": "),
  levels = slot_levels)

max_toe   <- max(dose_response_all$TOE, na.rm = TRUE)
n_slots   <- length(slot_levels)
run_break <- (408 + 411) / 2

q_boundaries <- 0.5 + block * (0:length(quadrants))
divider_q    <- data.frame(x = q_boundaries, y = 0, yend = max_toe)
pop_start    <- 0.5 + sort(unique(unlist(lapply(seq_along(quadrants) - 1, function(qi)
                  qi * block + (0:(pops_per_q - 1)) * n_traits))))
divider_pop  <- data.frame(x = pop_start, y = 0, yend = max_toe)

# quadrant names, horizontal, in the corner each quadrant points into
q_center <- q_boundaries[-length(q_boundaries)] + pops_per_q * n_traits / 2
q_lab <- data.frame(x = q_center, label = names(quadrants), y = max_toe * 1.235)
deg <- 90 - 360 * (q_lab$x - 0.5) / n_slots
q_lab$hjust <- ifelse(cos(deg * pi / 180) > 0, 0, 1)
q_lab$vjust <- ifelse(sin(deg * pi / 180) > 0, 0, 1)

# maternal-line letter on every wedge
line_lab <- do.call(rbind, lapply(seq_along(quadrants), function(qi) {
  pops <- quadrants[[qi]]
  keep <- !is.na(pops)
  data.frame(x = (qi - 1) * block + (which(keep) - 1) * n_traits + n_traits / 2 + 0.5,
             label = pops[keep])
}))
line_lab$y     <- max_toe * 1.165
line_lab$angle <- 90 - 360 * (line_lab$x - 0.5) / n_slots
line_lab$hjust <- ifelse(line_lab$angle < -90, 1, 0)
line_lab$angle <- ifelse(line_lab$angle < -90, line_lab$angle + 180, line_lab$angle)

# trait names, only on the first wedge of each quadrant
trait_lab <- do.call(rbind, lapply(seq_along(quadrants), function(qi)
  data.frame(x = (qi - 1) * block + seq_len(n_traits) - 0.5 + 0.5,
             label = unname(trait_label[trait_order]))))
trait_lab$y     <- max_toe * 1.045
trait_lab$angle <- 90 - 360 * (trait_lab$x - 0.5) / n_slots
trait_lab$hjust <- ifelse(trait_lab$angle < -90, 1, 0)
trait_lab$angle <- ifelse(trait_lab$angle < -90, trait_lab$angle + 180, trait_lab$angle)

toe_breaks   <- c(100, 200, 300, 400, 500, 600, 700)
toe_breaks   <- toe_breaks[toe_breaks <= max_toe]
toe_label_x  <- pops_per_q * n_traits + ceiling(n_gap / 2)
toe_label_df <- data.frame(x = toe_label_x, y = toe_breaks, label = paste0(toe_breaks, "h"))

p_treering <- ggplot(dose_response_all, aes(x = Slot, y = TOE, fill = LOD)) +
  geom_tile() +
  geom_hline(yintercept = toe_breaks, color = "grey55", linewidth = 0.3, alpha = 0.45) +
  geom_hline(yintercept = run_break, color = "grey25", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = subset(dose_response_all, FDR < 0.05), shape = 8, size = 0.35,
             color = "red", stroke = 0.25) +
  geom_segment(data = divider_pop, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, color = "white", linewidth = 0.4) +
  geom_segment(data = divider_q, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, color = "white", linewidth = 1.1) +
  geom_text(data = q_lab, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, fontface = "bold", size = 4.6, angle = 0) +
  geom_text(data = line_lab, aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
            inherit.aes = FALSE, fontface = "bold", size = 3.8) +
  geom_text(data = trait_lab, aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
            inherit.aes = FALSE, size = 1.7, color = "grey15") +
  geom_text(data = toe_label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 2.2, hjust = 0.5, color = "grey30") +
  scale_fill_distiller(palette = "GnBu", direction = 1, name = "-log10(p)\ndose slope",
                       na.value = "grey90") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, max_toe * 1.32), expand = c(0, 0)) +
  coord_polar(theta = "x") +
  labs(title = "Dose-response across 14 traits and 11 Mimulus maternal lines",
       subtitle = paste0("radius = Time Of Experiment (center = TOE 0, edge = ", max_toe,
          " h); angle = trait. Each quadrant is one population group and is divided into\n",
          "its maternal lines; the 14 traits run in the same order in every wedge, named on the first wedge of each quadrant.\n",
          "Red stars: FDR < 0.05 within that trait x line. Oregon 2014 has two lines, so its third wedge is blank.\n",
          "Only about 20 plants per line per timepoint.")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), axis.title = element_blank(),
        plot.subtitle = element_text(size = 8))

pdf("paper_graphs/magda/Dose_Response_TreeRing_Mimulus_byPopulation.pdf", height = 13, width = 13)
p_treering
dev.off()
