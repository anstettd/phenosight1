##################################################################################
## Fig 4: trait ~ SWHC regression, at one timepoint and across the whole
## experiment (MMJ, Nov 2026)
##
##  A  NDVI regressed on % SWHC at the first timepoint at which the regression
##     reaches P < 0.05 in at least one population group, all four groups shown
##  B  the tree-ring: the same regression repeated at every timepoint and for
##     every trait at once. Radius = time, angle = trait, one 90-degree quadrant
##     per population group, fill = the SIGNED strength of the SWHC slope, so
##     direction is visible: green = the trait rises as the pots get wetter,
##     purple = it falls (i.e. the trait is higher under drought).
##
## Every regression is fitted on MATERNAL-LINE MEANS, not on individual plants:
## 15 points per group per timepoint (5 SWHC regimes x 3 lines), 10 for Oregon
## 2014, which has only 2 lines.
##
## Runs off data/PS2_all.csv, written by 0_data_setup.R.
###################################################################################
#Import libraries
library(tidyverse)
library(patchwork)

RUN_BREAK <- (408 + 411) / 2
PANEL_A_TRAIT <- "NDVI"
PANEL_A_LABEL <- "NDVI (a.u.)"
#How many groups must reach P < 0.05 for a timepoint to qualify for panel A.
#1 gives the first timepoint the regression appears anywhere; 2 gives the first
#timepoint at which it is corroborated by a second group.
N_GROUPS_REQUIRED <- 1

#What the ring is coloured by. All three are diverging except "logp".
#  "signed_logp"  sign(slope) x -log10(P)  - direction + strength of evidence
#  "signed_r"     sign(slope) x sqrt(R2)   - direction + effect size (the
#                 correlation between trait and SWHC; comparable across traits,
#                 which a raw slope is not, since the traits have different units)
#  "logp"         -log10(P) alone, the earlier version of this figure: strength
#                 of evidence only, direction invisible
FILL_MODE <- "signed_logp"

#NB: 0_data_setup.R renames Roi.Mask to PSA, so the column here is PSA.
trait_order <- c("PSA","Fv.Fm","Fq.Fm","NPQ","qP","qN","qL","qI",
                 "phiNO","phiNPQ","npq.t","ChlIdx","AriIdx","NDVI")
trait_label <- setNames(trait_order, trait_order)
#PSA keeps its abbreviation here, as in Fig 1; the traits are spelled out in the legend.

#Quadrants clockwise from the top: 2014 upper half, 2010 lower half,
#California on the right, Oregon on the left - as in Figure 1.
group_order <- c("California 2014", "California 2010", "Oregon 2010", "Oregon 2014")

#Import Data (already maternal-line means per Population x Treatment x TOE)
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Group = paste(Site, Year),
         SWHC = as.numeric(gsub("%", "", Treatment)))

cat("points per group per timepoint:\n")
print(PS2_all %>% count(Group, TOE) %>% count(Group, n, name = "timepoints"))

############################################################################################################
#One regression per Group x Trait x TOE, on the line means

fit_one <- function(d, trait) {
  if (!trait %in% names(d)) stop("column not found in PS2_all: ", trait)
  u <- d[is.finite(d[[trait]]) & is.finite(d$SWHC), ]
  if (nrow(u) < 6 || length(unique(u$SWHC)) < 3)
    return(tibble(slope = NA_real_, P = NA_real_, R2 = NA_real_, n = nrow(u)))
  m <- lm(u[[trait]] ~ u$SWHC); s <- summary(m)
  tibble(slope = coef(m)[2], P = s$coefficients[2, 4], R2 = s$r.squared, n = nrow(u))
}

reg <- map_dfr(trait_order, function(tr) {
  PS2_all %>%
    group_by(Group, TOE) %>%
    group_modify(~ fit_one(.x, tr)) %>%
    ungroup() %>%
    mutate(Trait = tr)
})

reg <- reg %>%
  mutate(LOD = -log10(P),
         Group = factor(Group, levels = group_order),
         Trait = factor(Trait, levels = trait_order)) %>%
  group_by(Group, Trait) %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>%      # across timepoints, within trait x group
  ungroup()

#The quantity the ring is filled with. A raw slope cannot be put on one colour
#scale across the 14 traits - PSA is an area in a.u., Fv/Fm is a ratio between 0
#and 1 - so the sign of the slope is carried by the sign of the fill and the
#magnitude by a quantity that IS comparable across traits.
reg <- reg %>%
  mutate(Fill = switch(FILL_MODE,
                       signed_logp = sign(slope) * LOD,
                       signed_r    = sign(slope) * sqrt(R2),
                       logp        = LOD,
                       stop("unknown FILL_MODE: ", FILL_MODE)))

FILL_NAME <- switch(FILL_MODE,
                    signed_logp = "signed -log10(P)\nSWHC slope",
                    signed_r    = "signed sqrt(R2)\nSWHC slope",
                    logp        = "-log10(P)\nSWHC slope")

write_csv(reg, "data/Fig4_trait_vs_SWHC_lineMeans.csv")

cat("\nfill =", FILL_MODE, "| negative (purple) slopes:",
    sum(reg$slope < 0, na.rm = TRUE), "of", sum(!is.na(reg$slope)), "cells\n")

############################################################################################################
# Panel A — one timepoint

qualifying <- reg %>%
  filter(Trait == PANEL_A_TRAIT, P < 0.05) %>%
  count(TOE, name = "n_groups") %>%
  filter(n_groups >= N_GROUPS_REQUIRED) %>%
  arrange(TOE)

TOE_A <- qualifying$TOE[1]
cat("\nPanel A timepoint: TOE", TOE_A, "h  (first with >=", N_GROUPS_REQUIRED,
    "group(s) at P < 0.05 for", PANEL_A_TRAIT, ")\n")
cat("   next few qualifying timepoints:", head(qualifying$TOE, 6), "\n\n")

A_dat <- PS2_all %>%
  filter(TOE == TOE_A) %>%
  mutate(Site = factor(Site, levels = c("Oregon", "California")),
         Year = factor(Year, levels = c(2010, 2014)))

A_fit <- reg %>%
  filter(Trait == PANEL_A_TRAIT, TOE == TOE_A) %>%
  mutate(lab = paste0(Group, "   R2 = ", formatC(R2, format = "f", digits = 2),
                      "   P = ", formatC(P, format = "g", digits = 2)))
print(A_fit %>% select(Group, n, slope, R2, P))

A_lab <- paste(A_fit$lab[match(group_order, A_fit$Group)], collapse = "\n")

Fig4A <- ggplot(A_dat, aes(x = SWHC, y = .data[[PANEL_A_TRAIT]],
                           color = Site, linetype = Year)) +
  geom_point(aes(shape = Year), size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 0.7) +
  annotate("text", x = 18, y = min(A_dat[[PANEL_A_TRAIT]], na.rm = TRUE),
           label = A_lab, hjust = 0, vjust = 0, size = 4.0, color = "grey25") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("2010" = 16, "2014" = 17)) +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  guides(linetype = guide_legend(override.aes = list(color = "grey30")),
         shape    = guide_legend(override.aes = list(color = "grey30"))) +
  labs(x = "% SWHC", y = PANEL_A_LABEL,
       title = paste0(PANEL_A_LABEL, " vs % SWHC at ", TOE_A, " h")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 13.5, face = "bold"),
        plot.margin = margin(2, 2, 0, 2),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"))

############################################################################################################
# Panel B — the tree-ring

n_traits <- length(trait_order)
n_gap    <- 3
block    <- n_traits + n_gap
gap_id   <- paste0("__gap", seq_len(n_gap))

slot_levels <- unlist(lapply(group_order, function(g)
  c(paste(g, trait_order, sep = ": "), paste0(g, gap_id))))
reg$Slot <- factor(paste(reg$Group, reg$Trait, sep = ": "), levels = slot_levels)

max_toe <- max(reg$TOE, na.rm = TRUE)
n_slots <- length(slot_levels)

#geom_tile takes its height from the SMALLEST gap between TOE values. The grid is
#12 h everywhere except across the run break (408 -> 411 = 3 h), so left to itself
#it draws every tile 3 h tall and the ring comes out as thin slivers separated by
#white. Set the height to the usual spacing instead.
toe_steps <- diff(sort(unique(reg$TOE)))
tile_h <- as.numeric(names(sort(table(toe_steps), decreasing = TRUE))[1])
cat("tile height set to", tile_h, "h (modal TOE spacing; min gap is",
    min(toe_steps), "h)\n")

boundaries     <- 0.5 + block * (0:length(group_order))
divider_df     <- data.frame(x = boundaries, y = 0, yend = max_toe)
label_center_x <- boundaries[-length(boundaries)] + n_traits / 2

#Quadrant names horizontal, in the corner each quadrant points into
group_label_df <- data.frame(x = label_center_x, label = group_order, y = max_toe * 1.17)
deg <- 90 - 360 * (group_label_df$x - 0.5) / n_slots
group_label_df$hjust <- ifelse(cos(deg * pi / 180) > 0, 0, 1)
group_label_df$vjust <- ifelse(sin(deg * pi / 180) > 0, 0, 1)

toe_breaks   <- c(200, 300, 400, 500, 600, 700)
toe_breaks   <- toe_breaks[toe_breaks <= max_toe]
toe_label_x  <- match(paste0(group_order[1], gap_id[ceiling(n_gap / 2)]), slot_levels)
toe_label_df <- data.frame(x = toe_label_x, y = toe_breaks, label = paste0(toe_breaks, "h"))
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

#Mark the timepoint shown in panel A
a_ring <- data.frame(y = TOE_A)

#Diverging scale centred on zero, so a sign change is a colour change and a
#non-responding trait fades to near-white rather than to one end of the ramp.
#Limits are symmetric, otherwise the white would not sit at slope = 0.
if (FILL_MODE == "logp") {
  fill_scale <- scale_fill_distiller(palette = "GnBu", direction = 1,
                                     name = FILL_NAME, na.value = "grey90")
} else {
  fill_max <- max(abs(reg$Fill), na.rm = TRUE)
  fill_scale <- scale_fill_gradient2(low = "#762A83", mid = "grey96", high = "#1B7837",
                                     midpoint = 0, limits = c(-fill_max, fill_max),
                                     name = FILL_NAME, na.value = "grey90")
}

Fig4B <- ggplot(reg, aes(x = Slot, y = TOE, fill = Fill)) +
  geom_tile(height = tile_h) +
  geom_hline(yintercept = toe_breaks, color = "grey55", linewidth = 0.3, alpha = 0.45) +
  geom_hline(yintercept = RUN_BREAK, color = "grey25", linewidth = 0.4, linetype = "dashed") +
  geom_hline(data = a_ring, aes(yintercept = y), color = "red", linewidth = 0.3,
             linetype = "dotdash", alpha = 0.7) +
  geom_point(data = subset(reg, FDR < 0.05), shape = 8, size = 1.05,
             color = "white", stroke = 0.65) +
  geom_point(data = subset(reg, FDR < 0.05), shape = 8, size = 0.7,
             color = "red", stroke = 0.35) +
  geom_segment(data = divider_df, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, color = "white", linewidth = 0.9) +
  geom_text(data = group_label_df, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, fontface = "bold", size = 6.2, angle = 0) +
  geom_text(data = toe_label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.4, hjust = 0.5, color = "grey30") +
  geom_text(data = trait_label_df, aes(x = Slot, y = y, label = trait_label,
                                       angle = angle, hjust = hjust),
            inherit.aes = FALSE, size = 3.6, color = "grey15") +
  fill_scale +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, max_toe * 1.22), expand = c(0, 0)) +
  coord_polar(theta = "x") +
  theme_minimal() +
  guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.9,
                                barwidth = unit(15, "lines"), barheight = unit(0.8, "lines"))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "bottom",       # frees the full page width for the ring
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0))

############################################################################################################
#Assemble

Fig4 <- Fig4A / Fig4B +
  plot_layout(heights = c(1, 5.2)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

#NB: do not print Fig4 after this - a bare `Fig4` under Rscript writes Rplots.pdf.
ggsave("paper_graphs/Figure 4.pdf", Fig4, width = 13, height = 16)

cat("\nFDR < 0.05 cells:", sum(reg$FDR < 0.05, na.rm = TRUE), "of", sum(!is.na(reg$P)), "\n")
cat("sign of the slope in the FDR-significant cells, per trait",
    "(neg = higher under drought):\n")
print(reg %>% filter(FDR < 0.05) %>%
        group_by(Trait) %>%
        summarise(neg = sum(slope < 0), pos = sum(slope > 0), .groups = "drop"))
cat("significant timepoints per trait x group:\n")
print(reg %>% mutate(s = FDR < 0.05) %>%
        group_by(Trait, Group) %>% summarise(n = sum(s, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Group, values_from = n))
