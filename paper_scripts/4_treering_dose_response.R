# Mimulus drought-gradient tree-ring figure
# ------------------------------------------------------------------------------
# Same construction as Dose_Response_TreeRing_AllTraits.pdf in the comparative
# drought physiology notebook (20260707_DroughtCompPhysiology_Magda.Rmd), applied
# to Dan Anstett's Mimulus drought-gradient experiment (Exp 20243 + 20244).
#
#   radius = TOE (center = start of experiment, edge = 747 h)
#   angle  = trait, in four 90-degree quadrants, one per population group
#            (collection site x collection year), same 14-trait order in each
#   fill   = -log10(p) of the dose slope at that trait x timepoint x group
#   stars  = FDR < 0.05, BH within group x trait across timepoints
#
# Stress dose follows the SWHC convention of the original: larger = more severe,
# so StressDose = 100 - % SWHC (100% -> 0, 20% -> 80).
#
# Run from the experiment folder:
#   .../collaborations/Anstett_Debbas/20251209_Experiment01
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

trait_order <- c("Roi.Mask","Fv.Fm","Fq.Fm","NPQ","qP","qN","qL","qI",
                 "phiNO","phiNPQ","npq.t","ChlIdx","AriIdx","NDVI")

# display names: the column stays Roi.Mask, the figure says what it measures
trait_label <- setNames(trait_order, trait_order)
trait_label["Roi.Mask"] <- "Digital biomass"
# clockwise from the top: 2014 in the upper half, 2010 in the lower half,
# S. California on the right, Oregon on the left
group_order <- c("SCalifornia 2014","SCalifornia 2010","Oregon 2010","Oregon 2014")
group_label <- c("SCalifornia 2010" = "California 2010", "SCalifornia 2014" = "California 2014",
                 "Oregon 2010" = "Oregon 2010", "Oregon 2014" = "Oregon 2014")

# --- data ---------------------------------------------------------------------

decode <- read.csv("data/decode.csv", fileEncoding = "UTF-8-BOM")
decode$Group <- paste(decode$Site, decode$Year)

early <- read.csv("data/Experiment01_EarlyPS2_Mimulus.csv")   # Exp 20243, 0-408 h
late  <- read.csv("data/Experiment01_LatePS2_Mimulus.csv")    # Exp 20244, 411-747 h
ps2   <- rbind(early, late)

# the `comp` column repeats the 100% control rows once per pairwise comparison
ps2 <- ps2[!duplicated(ps2[, c("Plant.ID","TOE")]), ]

ps2$Group      <- decode$Group[match(ps2$Population, decode$Population)]
ps2$StressDose <- 100 - as.numeric(sub("%", "", ps2$Treatment))

# --- dose-response at every timepoint, within each population group ------------
# identical to dose_response_by_timepoint() in the comparative physiology
# notebook: >= 3 distinct dose levels and >= 6 plants required to fit.

dose_response_by_timepoint <- function(df, traits, group_label) {
  toes <- sort(unique(df$TOE))
  out  <- vector("list", length(toes) * length(traits))
  i <- 0
  for (t in toes) {
    d_t <- df[df$TOE == t, ]
    for (tr in traits) {
      i <- i + 1
      d_use <- d_t[is.finite(d_t[[tr]]) & is.finite(d_t$StressDose), c(tr, "StressDose")]
      if (nrow(d_use) < 6 || length(unique(d_use$StressDose)) < 3) {
        out[[i]] <- data.frame(Group = group_label, Trait = tr, TOE = t,
                               slope = NA_real_, p_value = NA_real_, r2 = NA_real_,
                               n = nrow(d_use))
        next
      }
      fit <- lm(d_use[[tr]] ~ d_use$StressDose)
      s   <- summary(fit)
      out[[i]] <- data.frame(Group = group_label, Trait = tr, TOE = t,
                             slope = coef(fit)[2], p_value = s$coefficients[2, 4],
                             r2 = s$r.squared, n = nrow(d_use))
    }
  }
  do.call(rbind, out)
}

dose_response_all <- do.call(rbind, lapply(group_order, function(g)
  dose_response_by_timepoint(subset(ps2, Group == g), trait_order, g)))

dose_response_all$LOD <- -log10(dose_response_all$p_value)
dose_response_all <- dose_response_all %>%
  group_by(Group, Trait) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(dose_response_all, "data/Mimulus_dose_response_by_timepoint.csv", row.names = FALSE)

# --- tree-ring layout ---------------------------------------------------------

n_traits <- length(trait_order)
n_gap    <- 3
block    <- n_traits + n_gap
gap_id   <- paste0("__gap", seq_len(n_gap))

slot_levels <- unlist(lapply(group_order, function(g)
  c(paste(g, trait_order, sep = ": "), paste0(g, gap_id))))
dose_response_all$Slot <- factor(paste(dose_response_all$Group, dose_response_all$Trait, sep = ": "),
                                 levels = slot_levels)

max_toe   <- max(dose_response_all$TOE, na.rm = TRUE)
n_slots   <- length(slot_levels)
run_break <- (408 + 411) / 2     # Exp 20243 ends 408 h, Exp 20244 starts 411 h

boundaries     <- 0.5 + block * (0:length(group_order))
divider_df     <- data.frame(x = boundaries, xend = boundaries, y = 0, yend = max_toe)
label_center_x <- boundaries[-length(boundaries)] + n_traits / 2

# quadrant names set horizontally in the corner each quadrant points into,
# rather than curved around the rim - easier to read
group_label_df <- data.frame(x = label_center_x,
                             label = unname(group_label[group_order]),
                             y = max_toe * 1.17)
deg <- 90 - 360 * (group_label_df$x - 0.5) / n_slots
group_label_df$hjust <- ifelse(cos(deg * pi / 180) > 0, 0, 1)
group_label_df$vjust <- ifelse(sin(deg * pi / 180) > 0, 0, 1)

toe_breaks   <- c(100, 200, 300, 400, 500, 600, 700)
toe_breaks   <- toe_breaks[toe_breaks <= max_toe]
toe_label_x  <- match(paste0(group_order[1], gap_id[ceiling(n_gap / 2)]), slot_levels)
toe_label_df <- data.frame(x = toe_label_x, y = toe_breaks, label = paste0(toe_breaks, "h"))
brk_label_df <- data.frame(x = match(paste0(group_order[2], gap_id[ceiling(n_gap / 2)]), slot_levels),
                           y = run_break, label = "run break")

trait_label_df <- data.frame(
  Slot = factor(paste(rep(group_order, each = n_traits),
                      rep(trait_order, times = length(group_order)), sep = ": "),
                levels = slot_levels),
  trait_label = rep(unname(trait_label[trait_order]), times = length(group_order)))
trait_label_df$slot_index <- as.integer(trait_label_df$Slot)
trait_label_df$angle <- 90 - 360 * (trait_label_df$slot_index - 0.5) / n_slots
trait_label_df$hjust <- ifelse(trait_label_df$angle < -90, 1, 0)
trait_label_df$angle <- ifelse(trait_label_df$angle < -90, trait_label_df$angle + 180, trait_label_df$angle)
trait_label_df$y     <- max_toe * 1.05

p_treering <- ggplot(dose_response_all, aes(x = Slot, y = TOE, fill = LOD)) +
  geom_tile() +
  geom_hline(yintercept = toe_breaks, color = "grey55", linewidth = 0.3, alpha = 0.45) +
  geom_hline(yintercept = run_break, color = "grey25", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = subset(dose_response_all, FDR < 0.05), shape = 8, size = 0.5,
             color = "red", stroke = 0.3) +
  geom_segment(data = divider_df, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, color = "white", linewidth = 0.9) +
  geom_text(data = group_label_df, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, fontface = "bold", size = 5, angle = 0) +
  geom_text(data = toe_label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 2.5, hjust = 0.5, color = "grey30") +
  geom_text(data = brk_label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 2.3, hjust = 0.5, color = "grey20", fontface = "italic") +
  geom_text(data = trait_label_df, aes(x = Slot, y = y, label = trait_label, angle = angle, hjust = hjust),
            inherit.aes = FALSE, size = 2.4, color = "grey15") +
  scale_fill_distiller(palette = "GnBu", direction = 1, name = "-log10(p)\ndose slope",
                       na.value = "grey90") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, max_toe * 1.28), expand = c(0, 0)) +
  coord_polar(theta = "x") +
  labs(title = "Dose-response across 14 traits and 4 Mimulus population groups",
       subtitle = paste0("radius = Time Of Experiment (center = TOE 0, edge = ", max_toe,
                         " h); angle = trait, grouped into one 90-degree quadrant per population group\n",
                         "(collection site x collection year), trait names labeled at each wedge's outer edge. ",
                         "Red stars: FDR < 0.05 within that trait x group.\n",
                         "Stress dose = 100 - % SWHC, so a larger dose is a drier pot. ",
                         "Dashed ring: Exp 20243 (inner) and Exp 20244 (outer) are separate plant cohorts.")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), axis.title = element_blank(),
        plot.subtitle = element_text(size = 9))

pdf("paper_graphs/magda/Dose_Response_TreeRing_Mimulus.pdf", height = 12, width = 12)
p_treering
dev.off()
