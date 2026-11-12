##################################################################################
## Fig 4 (MAIN): trait ~ SWHC regression across the whole experiment, as a
## tree-ring, with the single-timepoint regression inset in each corner
## (MMJ, Nov 2026)
##
## One figure, not two panels. The ring is the result: the trait ~ SWHC
## regression repeated at every timepoint and for every trait. Radius = time,
## angle = trait, one 90-degree quadrant per population group, fill = the SIGNED
## strength of the SWHC slope, so direction is visible: green = the trait rises
## as the pots get wetter, purple = it falls (i.e. the trait is higher under
## drought).
##
## In each corner, in the corner its own quadrant points into, a small
## regression of PANEL_A_TRAIT on % SWHC at one timepoint for that group - the
## worked example of the single cell the ring summarises. These replace the old
## full-width panel A, which took a third of the page to illustrate a concept.
##
## The cell each inset comes from is marked three ways, all of them exact in
## polar coordinates: outlined in red, sitting on the red dot-dash ring that
## carries its timepoint right round the figure, and at the end of a short red
## radial pointer running out from the centre. The ring and the pointer cross on
## the cell.
##
## The quadrant names run along the arc, centred on each quadrant's own span,
## just outside the trait labels. Tangential rather than horizontal because a
## horizontal name centred on the quadrant bisector either dips into the ring at
## its ends or runs out into the corner where the inset is; along the arc it does
## neither.
##
## Written twice, once per replication unit:
##   Figure4_TreeRing_lineMeans.pdf  replicate = maternal line
##                                              (data/PS2_all.csv)
##   Figure4_TreeRing_plants.pdf     replicate = individual plant
##                                              (data/PS2_all_plants.csv)
## Same regressions, different n: 15 points per group per timepoint at line level
## (5 SWHC regimes x 3 lines; 10 for Oregon 2014), roughly 50-60 at plant level.
## The plant-level version will look far more significant, but plants within a
## maternal line are not independent, so it overstates the evidence.
##
## Both tables are written by 0_data_setup.R.
###################################################################################
#Import libraries
library(tidyverse)
library(patchwork)

RUN_BREAK <- (408 + 411) / 2
PANEL_A_TRAIT <- "NDVI"
PANEL_A_LABEL <- "NDVI (a.u.)"
#How many groups must reach P < 0.05 for a timepoint to qualify for the insets.
#1 gives the first timepoint the regression appears anywhere; 2 gives the first
#timepoint at which it is corroborated by a second group. The timepoint is chosen
#ONCE, from the line means, and used for both versions of the figure, so that the
#only thing that differs between them is the replication unit.
N_GROUPS_REQUIRED <- 1

#What the ring is coloured by. All three are diverging except "logp".
#  "signed_logp"  sign(slope) x -log10(P)  - direction + strength of evidence
#  "signed_r"     sign(slope) x sqrt(R2)   - direction + effect size (the
#                 correlation between trait and SWHC; comparable across traits,
#                 which a raw slope is not, since the traits have different units)
#  "logp"         -log10(P) alone, the earlier version of this figure: strength
#                 of evidence only, direction invisible
FILL_MODE <- "signed_logp"

NEG_COL <- "#762A83"   # purple: trait falls as the pots get wetter
POS_COL <- "#1B7837"   # green:  trait rises as the pots get wetter

#Outer padding of the radial axis: the plotted radius runs to YLIM_PAD * max_toe
#while the tiles stop at max_toe, which is what leaves room for the trait labels
#and the quadrant names. The leader lines below depend on this number, so change
#it in one place only.
YLIM_PAD <- 1.22

#Side of each corner inset, as a fraction of the (square) ring panel. The ring's
#drawn content, trait labels included, reaches 1.05/1.22 = 0.86 of the panel's
#half-width, i.e. a radius of 0.43 panel sides. A square inset of side s in a
#corner first touches that circle when sqrt(2) * (0.5 - s) = 0.43, i.e. at
#s = 0.196 - so keep this below ~0.19. Raise it for legibility, lower it if the
#trait labels start to collide.
INSET_SIZE <- 0.18

#Radius at which the quadrant names sit, as a multiple of max_toe: centred on
#each quadrant's own span and just outside the trait labels, which start at 1.05
#and run out to about 1.12 for the longest of them.
GROUP_LABEL_R <- 1.16

#NB: 0_data_setup.R renames Roi.Mask to PSA, so the column here is PSA.
trait_order <- c("PSA","Fv.Fm","Fq.Fm","NPQ","qP","qN","qL","qI",
                 "phiNO","phiNPQ","npq.t","ChlIdx","AriIdx","NDVI")
trait_label <- setNames(trait_order, trait_order)
#PSA keeps its abbreviation here, as in Fig 1; the traits are spelled out in the legend.

#Quadrants clockwise from the top: 2014 upper half, 2010 lower half,
#California on the right, Oregon on the left - as in Figure 1.
group_order <- c("California 2014", "California 2010", "Oregon 2010", "Oregon 2014")

#Which corner each quadrant points into, as c(left, bottom, right, top) in
#fractions of the ring panel. Follows group_order above: the first quadrant runs
#clockwise from 12 o'clock, so it occupies the upper right, and so on round.
S <- INSET_SIZE
CORNERS <- list(
  "California 2014" = c(1 - S, 1 - S, 1,     1    ),   # upper right
  "California 2010" = c(1 - S, 0,     1,     S    ),   # lower right
  "Oregon 2010"     = c(0,     0,     S,     S    ),   # lower left
  "Oregon 2014"     = c(0,     1 - S, S,     1    )    # upper left
)

#Main figure, so this one goes in paper_graphs/ itself, not the supplement.
OUT_DIR <- "paper_graphs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

UNITS <- list(
  list(tag = "lineMeans", ps2 = "data/PS2_all.csv",
       label = "Replicate = maternal line (3 lines per group; 2 for Oregon 2014)",
       pt = 1.8, alpha = 0.75),
  list(tag = "plants",    ps2 = "data/PS2_all_plants.csv",
       label = "Replicate = individual plant (3-5 plants per line per regime)",
       pt = 0.9, alpha = 0.40)
)

read_unit <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(Group = paste(Site, Year),
           SWHC = as.numeric(gsub("%", "", Treatment)))
}

############################################################################################################
#One regression per Group x Trait x TOE. Works unchanged whether the rows are
#line means or plants - only n changes.

fit_one <- function(d, trait) {
  if (!trait %in% names(d)) stop("column not found in the input table: ", trait)
  u <- d[is.finite(d[[trait]]) & is.finite(d$SWHC), ]
  if (nrow(u) < 6 || length(unique(u$SWHC)) < 3)
    return(tibble(slope = NA_real_, P = NA_real_, R2 = NA_real_, n = nrow(u)))
  m <- lm(u[[trait]] ~ u$SWHC); s <- summary(m)
  tibble(slope = coef(m)[2], P = s$coefficients[2, 4], R2 = s$r.squared, n = nrow(u))
}

regress_all <- function(dat) {
  reg <- map_dfr(trait_order, function(tr) {
    dat %>%
      group_by(Group, TOE) %>%
      group_modify(~ fit_one(.x, tr)) %>%
      ungroup() %>%
      mutate(Trait = tr)
  })

  reg %>%
    mutate(LOD = -log10(P),
           Group = factor(Group, levels = group_order),
           Trait = factor(Trait, levels = trait_order)) %>%
    group_by(Group, Trait) %>%
    mutate(FDR = p.adjust(P, method = "BH")) %>%   # across timepoints, within trait x group
    ungroup() %>%
    #The quantity the ring is filled with. A raw slope cannot be put on one colour
    #scale across the 14 traits - PSA is an area in a.u., Fv/Fm is a ratio between
    #0 and 1 - so the sign of the slope is carried by the sign of the fill and the
    #magnitude by a quantity that IS comparable across traits.
    mutate(Fill = switch(FILL_MODE,
                         signed_logp = sign(slope) * LOD,
                         signed_r    = sign(slope) * sqrt(R2),
                         logp        = LOD,
                         stop("unknown FILL_MODE: ", FILL_MODE)))
}

FILL_NAME <- switch(FILL_MODE,
                    signed_logp = "signed -log10(P)\nSWHC slope",
                    signed_r    = "signed sqrt(R2)\nSWHC slope",
                    logp        = "-log10(P)\nSWHC slope")

############################################################################################################
#The timepoint shown in the insets: chosen once, from the line means, so both
#versions of the figure illustrate the same cell of the ring.

reg_lines <- regress_all(read_unit("data/PS2_all.csv"))

qualifying <- reg_lines %>%
  filter(Trait == PANEL_A_TRAIT, P < 0.05) %>%
  count(TOE, name = "n_groups") %>%
  filter(n_groups >= N_GROUPS_REQUIRED) %>%
  arrange(TOE)

TOE_A <- qualifying$TOE[1]
cat("\nInset timepoint: TOE", TOE_A, "h  (first with >=", N_GROUPS_REQUIRED,
    "group(s) at P < 0.05 for", PANEL_A_TRAIT, ", on line means)\n")
cat("   next few qualifying timepoints:", head(qualifying$TOE, 6), "\n\n")

############################################################################################################
#One corner inset: that group's PANEL_A_TRAIT ~ SWHC regression at TOE_A.
#The fitted line takes the ring's colour semantics - purple for a negative
#slope, green for a positive one - so the inset and the wedge it sits beside
#say the same thing.

mini_reg <- function(dat, fits, g, u) {
  d <- dat %>% filter(Group == g)
  f <- fits %>% filter(Group == g)
  line_col <- if (nrow(f) == 1 && is.finite(f$slope) && f$slope < 0) NEG_COL else POS_COL
  sub <- if (nrow(f) == 1)
    paste0("R2 = ", formatC(f$R2, format = "f", digits = 2),
           "   P = ", formatC(f$P, format = "g", digits = 2),
           "   n = ", f$n) else ""

  ggplot(d, aes(x = SWHC, y = .data[[PANEL_A_TRAIT]])) +
    geom_point(size = u$pt, alpha = u$alpha, color = "grey35") +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
                color = line_col, linewidth = 0.9) +
    scale_x_continuous(breaks = c(20, 60, 100), limits = c(12, 108)) +
    labs(title = g, subtitle = sub, x = "% SWHC", y = PANEL_A_LABEL) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8.5, color = "grey25"),
          axis.title = element_text(size = 9, face = "bold"),
          axis.text = element_text(size = 8),
          plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
          plot.margin = margin(4, 6, 4, 4))
}

############################################################################################################
#The ring, and the whole figure, for one replication unit.

build_fig4 <- function(u) {

  dat <- read_unit(u$ps2)
  reg <- if (u$tag == "lineMeans") reg_lines else regress_all(dat)

  cat("\n================ ", u$tag, " ================\n", sep = "")
  cat("points per group per timepoint:\n")
  print(dat %>% count(Group, TOE) %>% count(Group, n, name = "timepoints"))
  cat("fill =", FILL_MODE, "| negative (purple) slopes:",
      sum(reg$slope < 0, na.rm = TRUE), "of", sum(!is.na(reg$slope)), "cells\n")

  write_csv(reg, paste0("data/Fig4_trait_vs_SWHC_", u$tag, ".csv"))

  ##########################################################################################
  #Geometry of the ring

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

  boundaries <- 0.5 + block * (0:length(group_order))
  divider_df <- data.frame(x = boundaries, y = 0, yend = max_toe)

  #Quadrant names, centred on each quadrant's own span - the middle of its block
  #of traits, not the gap that closes it. That gap sits against the NEXT
  #quadrant's first trait, so a name placed there reads as belonging to that one.
  #
  #Set along the arc. A tangential label keeps a constant radius, so it neither
  #dips into the ring at its ends (as a horizontal label centred on the bisector
  #does) nor runs out into the corner where the inset is. The trait labels read
  #radially at 90 - 360 * (slot - 0.5) / n_slots; tangential is that minus 90,
  #then folded into (-90, 90] so no name comes out upside down.
  group_label_df <- data.frame(
    x = 0.5 + block * (seq_along(group_order) - 1) + n_traits / 2,
    label = group_order,
    y = max_toe * GROUP_LABEL_R)
  group_label_df$angle <- -360 * (group_label_df$x - 0.5) / n_slots
  group_label_df$angle <- ((group_label_df$angle + 90) %% 180) - 90

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
  trait_label_df$angle <- ifelse(trait_label_df$angle < -90,
                                 trait_label_df$angle + 180, trait_label_df$angle)
  trait_label_df$y     <- max_toe * 1.05

  #Mark the timepoint shown in the insets, and the four cells themselves.
  #The pointer runs from the centre out to the cell at the same slot, so it is
  #drawn in the ring's own polar coordinates and lands on the cell by
  #construction rather than by arithmetic. A line at constant angle is straight
  #in polar coordinates, which a line between two arbitrary points is not. Most
  #of it falls in the empty hole at the centre, since no tile is drawn inside the
  #first timepoint of the series.
  a_ring <- data.frame(y = TOE_A)
  hl <- reg %>% filter(Trait == PANEL_A_TRAIT, TOE == TOE_A)
  spoke <- data.frame(x = match(paste(group_order, PANEL_A_TRAIT, sep = ": "),
                                slot_levels),
                      y = 0, yend = TOE_A - tile_h / 2)

  #Diverging scale centred on zero, so a sign change is a colour change and a
  #non-responding trait fades to near-white rather than to one end of the ramp.
  #Limits are symmetric, otherwise the white would not sit at slope = 0.
  if (FILL_MODE == "logp") {
    fill_scale <- scale_fill_distiller(palette = "GnBu", direction = 1,
                                       name = FILL_NAME, na.value = "grey90")
  } else {
    fill_max <- max(abs(reg$Fill), na.rm = TRUE)
    fill_scale <- scale_fill_gradient2(low = NEG_COL, mid = "grey96", high = POS_COL,
                                       midpoint = 0, limits = c(-fill_max, fill_max),
                                       name = FILL_NAME, na.value = "grey90")
  }

  ring <- ggplot(reg, aes(x = Slot, y = TOE, fill = Fill)) +
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
    #The one cell each inset is drawn from: a radial pointer out to it, then its
    #outline. Both get a white underlay so they read against any fill.
    geom_segment(data = spoke, aes(x = x, xend = x, y = y, yend = yend),
                 inherit.aes = FALSE, color = "white", linewidth = 1.6) +
    geom_segment(data = spoke, aes(x = x, xend = x, y = y, yend = yend),
                 inherit.aes = FALSE, color = "red", linewidth = 0.7) +
    geom_tile(data = hl, aes(x = Slot, y = TOE), inherit.aes = FALSE,
              height = tile_h, fill = NA, color = "red", linewidth = 1.0) +
    geom_text(data = group_label_df, aes(x = x, y = y, label = label,
                                         angle = angle),
              inherit.aes = FALSE, fontface = "bold", size = 5.6) +
    geom_text(data = toe_label_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 3.4, hjust = 0.5, color = "grey30") +
    geom_text(data = trait_label_df, aes(x = Slot, y = y, label = trait_label,
                                         angle = angle, hjust = hjust),
              inherit.aes = FALSE, size = 3.6, color = "grey15") +
    fill_scale +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, max_toe * YLIM_PAD), expand = c(0, 0)) +
    coord_polar(theta = "x") +
    labs(title = u$label) +
    theme_minimal() +
    guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.9,
                                  barwidth = unit(15, "lines"), barheight = unit(0.8, "lines"))) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "bottom",       # frees the full page width for the ring
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 10),
          plot.margin = margin(0, 0, 0, 0))

  ##########################################################################################
  #The four corner insets

  A_dat <- dat %>% filter(TOE == TOE_A)
  A_fit <- reg %>% filter(Trait == PANEL_A_TRAIT, TOE == TOE_A)
  print(A_fit %>% select(Group, n, slope, R2, P))

  #No leader line from the inset to the cell. It was drawn on a transparent
  #overlay laid over the panel, which assumes that overlay's 0-1 box and the
  #ring's drawn circle share a centre AND a radius. coord_polar forces the panel
  #square, so as soon as the page is not exactly square the drawn circle is inset
  #within the panel box and the two disagree: the line lands at the wrong radius,
  #which reads as pointing at a different timepoint from the cell it belongs to.
  #The radial pointer above does the same job inside the polar coordinates, where
  #it cannot drift.
  Fig4 <- ring
  for (g in group_order) {
    cr <- CORNERS[[g]]
    Fig4 <- Fig4 + inset_element(mini_reg(A_dat, A_fit, g, u),
                                 left = cr[1], bottom = cr[2], right = cr[3], top = cr[4],
                                 align_to = "panel")
  }

  out <- file.path(OUT_DIR, paste0("Figure4_TreeRing_", u$tag, ".pdf"))

  #coord_polar forces a square panel, so the page has to be close to square or the
  #ring shrinks and leaves bands of white. The extra height is the colourbar.
  #NB: do not print Fig4 after this - a bare `Fig4` under Rscript writes Rplots.pdf.
  ggsave(out, Fig4, width = 15, height = 15.8)
  cat("wrote", out, "\n")

  cat("FDR < 0.05 cells:", sum(reg$FDR < 0.05, na.rm = TRUE), "of", sum(!is.na(reg$P)), "\n")
  cat("significant timepoints per trait x group:\n")
  print(reg %>% mutate(s = FDR < 0.05) %>%
          group_by(Trait, Group) %>% summarise(n = sum(s, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = Group, values_from = n))
  cat("sign of the slope in the FDR-significant cells, per trait",
      "(neg = higher under drought):\n")
  print(reg %>% filter(FDR < 0.05) %>%
          group_by(Trait) %>%
          summarise(neg = sum(slope < 0), pos = sum(slope > 0), .groups = "drop"))

  invisible(NULL)
}

for (u in UNITS) build_fig4(u)
