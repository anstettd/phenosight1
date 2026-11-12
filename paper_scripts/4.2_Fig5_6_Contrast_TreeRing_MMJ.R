##################################################################################
## Figs 5 and 6: contrasts between collections, as half-rings (MMJ, Nov 2026)
##
## Derivatives of Figure 4 - same geometry, same reading - but each cell carries a
## CONTRAST BETWEEN TWO COLLECTIONS rather than a trait. Two contrasts, two
## quantities, two replication units: eight figures from one run.
##
## THE CONTRAST (what is divided or subtracted, and what the quadrants are)
##
##   "year"  Figure 5. 2014 against 2010, within each site.
##           One quadrant per SITE: California upper right, Oregon upper left,
##           keeping Figure 4's left/right convention.
##
##   "site"  Figure 6. California against Oregon, within each collection year.
##           One quadrant per YEAR: 2010 upper left, 2014 upper right, so the
##           older collection reads on the left.
##
## THE QUANTITY (what fills a cell)
##
##   "dslope"  CHANGE IN PLASTICITY. Within a quadrant, trait and timepoint,
##             trait ~ SWHC is fitted separately in each level of the contrast -
##             these are Figure 4's regressions - and the cell carries
##                   (slope_comparison - slope_reference) / residual SD
##             i.e. how much more steeply the comparison collection responds to
##             watering, in residual standard deviations per % SWHC. Dividing by
##             the residual SD of the same model makes it an effect size and so
##             comparable across traits carrying different units.
##             The asterisks are the contrast x SWHC interaction from
##             lm(trait ~ C * SWHC): the formal test of exactly what the colour
##             shows, so fill and asterisks test the same thing.
##
##   "ratio"   LEVEL CONTRAST, the "evolutionary index" for the year version.
##             EI = mean(comparison) / mean(reference) per quadrant, trait,
##             timepoint and watering regime; EI is then regressed on % SWHC and
##             the cell carries that slope.
##             The asterisks are the P of the C term in lm(trait ~ C + SWHC) -
##             "do the two collections differ here at all", holding watering
##             constant. Note this is NOT a test of what the colour shows.
##
## The two quantities measure different things and both are worth a look:
## "ratio" is mostly a LEVEL shift, "dslope" a PLASTICITY shift. A population can
## evolve either without the other.
##
## Reading plasticity acquired over the experiment: radius is time, so a wedge
## that starts near-white at the centre and deepens outward is a difference that
## BUILDS as the experiment runs, rather than one present from the first image.
##
## REPLICATION. Each figure is written twice:
##   *_lineMeans.pdf   replicate = maternal line    (data/PS2_all.csv)
##   *_plants.pdf      replicate = individual plant (data/PS2_all_plants.csv)
## Both tables come from 0_data_setup.R. The line-level version is the
## conservative one; plants within a line are not independent replicates of the
## line's response, so the plant-level P values are optimistic.
##
## CAVEATS worth keeping in the text
##   - Figure 5: "year" is confounded with "these particular maternal lines".
##     Oregon 2010 is G, H, I and Oregon 2014 is J, K, so with three and two
##     lines the contrast is between-line as much as between-year.
##   - Figure 6: "site" is likewise confounded with line sampling, and the two
##     sites were not collected in balanced numbers in both years.
##   - "ratio" divides by a trait mean, so it is unstable wherever that mean
##     approaches zero; see MIN_DENOM_FRAC. "dslope" never divides by a trait
##     value and needs no such guard.
###################################################################################
#Import libraries
library(tidyverse)
library(patchwork)

#Everything below is looped over. Drop entries to build fewer figures.
CONTRASTS <- c("year", "site")
MODES     <- c("dslope", "ratio")

RUN_BREAK <- (408 + 411) / 2
YLIM_PAD  <- 1.22

#"ratio" mode only: "ratio" = comparison/reference as specified; "log2" = log2 of
#it, symmetric in the two directions and better behaved as a regression
#response, but it needs both means positive.
INDEX_SCALE <- "ratio"

#"ratio" mode only. A ratio of means blows up when the denominator approaches
#zero, and several of these traits legitimately sit near zero (phiNO, qI, NPQ
#early on). Cells whose reference mean is under this fraction of the trait's own
#typical magnitude are dropped rather than drawn as a spike.
MIN_DENOM_FRAC <- 0.05

#Colour limits sit at this quantile of |fill|; anything past them is squished to
#the end of the ramp, so one unstable cell cannot flatten the rest.
FILL_QUANTILE <- 0.98

#The trait and timepoint shown in the two lower insets.
INSET_TRAIT <- "NDVI"
DRAW_INSETS <- TRUE

#Where the left inset sits in the ring's panel, c(left, bottom, right, top); the
#right one is its mirror image. Keep top at or below 0.5 - the tiles occupy the
#two UPPER quadrants, so anything above the half way line lands on them.
INSET_BOX <- c(0.04, 0.23, 0.36, 0.48)

#NB: 0_data_setup.R renames Roi.Mask to PSA, so the column here is PSA.
trait_order <- c("PSA","Fv.Fm","Fq.Fm","NPQ","qP","qN","qL","qI",
                 "phiNO","phiNPQ","npq.t","ChlIdx","AriIdx","NDVI")
trait_label <- setNames(trait_order, trait_order)

NEG_COL <- "#762A83"   # purple
POS_COL <- "#1B7837"   # green
LVL_COL <- c("#2166AC", "#D95F02")   # reference, comparison

OUT_DIR <- "paper_graphs"      # main figures, beside Figure 4
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#The two contrasts. cvar is contrasted, fvar becomes the quadrants, and the
#FIRST facet listed lands in the upper right, since the slot order runs
#clockwise from 12 o'clock.
SPECS <- list(
  year = list(fig = "Figure5", cvar = "Year", fvar = "Site",
              ref = "2010", cmp = "2014",
              facets = c("California", "Oregon"),
              #the level contrast between years is the "evolutionary index"
              ratio_stem = "EvolIndex",
              what = "2014 relative to 2010, within each site"),
  site = list(fig = "Figure6", cvar = "Site", fvar = "Year",
              ref = "Oregon", cmp = "California",
              facets = c("2014", "2010"),
              #between sites it is a divergence, not an evolutionary change, so
              #it is not called an evolutionary index
              ratio_stem = "LevelRatio",
              what = "California relative to Oregon, within each collection year")
)

UNITS <- list(
  list(tag = "lineMeans", ps2 = "data/PS2_all.csv",
       label = "Replicate = maternal line (3 lines per group; 2 for Oregon 2014)"),
  list(tag = "plants",    ps2 = "data/PS2_all_plants.csv",
       label = "Replicate = individual plant (3-5 plants per line per regime)")
)

read_unit <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(SWHC = as.numeric(gsub("%", "", Treatment)),
           Year = as.character(Year),
           Site = as.character(Site))
}

#Con is the contrasted variable, reference level first; Fac is the quadrant.
prep <- function(dat, sp) {
  dat %>%
    mutate(Con = factor(.data[[sp$cvar]], levels = c(sp$ref, sp$cmp)),
           Fac = factor(.data[[sp$fvar]], levels = sp$facets)) %>%
    filter(!is.na(Con), !is.na(Fac))
}

long_form <- function(d) {
  d %>%
    select(Fac, Con, SWHC, TOE, all_of(trait_order)) %>%
    pivot_longer(all_of(trait_order), names_to = "Trait", values_to = "value")
}

############################################################################################################
#Both tests and both per-level slopes, from the replicate-level rows of one
#Fac x Trait x TOE.
#
#  P_con        C term of lm(value ~ C + SWHC)  - do the two levels differ at all
#  P_con_swhc   C:SWHC term of lm(value ~ C * SWHC) - do they differ in HOW they
#               respond to watering. This is the test of "dslope".
#  sigma        residual SD of the interaction model, standardises dslope

test_contrast <- function(d) {
  u <- d[is.finite(d$value), ]
  blank <- tibble(beta_con = NA_real_, P_con = NA_real_, P_con_swhc = NA_real_,
                  sigma = NA_real_, slope_ref = NA_real_, slope_cmp = NA_real_,
                  n_con = nrow(u))
  if (nrow(u) < 8) return(blank)
  u$Con <- droplevels(u$Con)
  if (nlevels(u$Con) < 2) return(blank)
  #each level needs at least three watering regimes for a slope of its own
  if (any(tapply(u$SWHC, u$Con, function(x) length(unique(x))) < 3)) return(blank)

  lv  <- levels(u$Con)
  co  <- summary(lm(value ~ Con + SWHC, data = u))$coefficients
  m2  <- lm(value ~ Con * SWHC, data = u)
  co2 <- summary(m2)$coefficients
  ic  <- grep("^Con", rownames(co))[1]
  ii  <- grep("^Con.*:SWHC$|^SWHC:Con", rownames(co2))[1]

  s_ref <- coef(lm(value ~ SWHC, data = u[u$Con == lv[1], ]))[2]
  s_cmp <- coef(lm(value ~ SWHC, data = u[u$Con == lv[2], ]))[2]

  tibble(beta_con = co[ic, 1], P_con = co[ic, 4],
         P_con_swhc = if (!is.na(ii)) co2[ii, 4] else NA_real_,
         sigma = summary(m2)$sigma,
         slope_ref = unname(s_ref), slope_cmp = unname(s_cmp),
         n_con = nrow(u))
}

############################################################################################################
#The level contrast: mean(comparison) / mean(reference), with the near-zero
#denominator guard, then regressed on SWHC.

group_means <- function(lf) {
  lf %>%
    group_by(Fac, Con, Trait, SWHC, TOE) %>%
    summarise(m = mean(value, na.rm = TRUE),
              n_reps = sum(is.finite(value)), .groups = "drop")
}

level_index <- function(gm, sp) {
  typ <- gm %>% group_by(Trait) %>%
    summarise(typ = median(abs(m), na.rm = TRUE), .groups = "drop")

  wide <- gm %>%
    select(Fac, Con, Trait, SWHC, TOE, m) %>%
    pivot_wider(names_from = Con, values_from = m, names_prefix = "lvl_") %>%
    left_join(typ, by = "Trait") %>%
    mutate(m_ref = .data[[paste0("lvl_", sp$ref)]],
           m_cmp = .data[[paste0("lvl_", sp$cmp)]])

  ok <- is.finite(wide$m_ref) & is.finite(wide$m_cmp) &
        abs(wide$m_ref) > MIN_DENOM_FRAC * wide$typ
  if (INDEX_SCALE == "log2") ok <- ok & wide$m_ref > 0 & wide$m_cmp > 0

  wide$EI <- ifelse(ok, wide$m_cmp / wide$m_ref, NA_real_)
  if (INDEX_SCALE == "log2") wide$EI <- log2(wide$EI)

  cat("  index cells kept:", sum(ok), "of", nrow(wide), "\n")
  if (any(!ok)) print(wide %>% filter(!ok) %>% count(Trait, name = "dropped"))
  wide
}

fit_index <- function(d) {
  u <- d[is.finite(d$EI) & is.finite(d$SWHC), ]
  if (nrow(u) < 3 || length(unique(u$SWHC)) < 3)
    return(tibble(ei_slope = NA_real_, ei_P = NA_real_, ei_R2 = NA_real_, n_swhc = nrow(u)))
  m <- lm(EI ~ SWHC, data = u); s <- summary(m)
  tibble(ei_slope = coef(m)[2], ei_P = s$coefficients[2, 4],
         ei_R2 = s$r.squared, n_swhc = nrow(u))
}

#Everything for one contrast and one replication unit, both quantities at once -
#the regressions are the same for either, only the column that becomes the fill
#differs, so they are computed once and shared.
build_reg <- function(dat, sp) {
  d  <- prep(dat, sp)
  lf <- long_form(d)

  con <- lf %>% group_by(Fac, Trait, TOE) %>%
    group_modify(~ test_contrast(.x)) %>% ungroup()

  ei  <- level_index(group_means(lf), sp)
  idx <- ei %>% group_by(Fac, Trait, TOE) %>%
    group_modify(~ fit_index(.x)) %>% ungroup()

  reg <- con %>%
    left_join(idx, by = c("Fac", "Trait", "TOE")) %>%
    mutate(dslope = slope_cmp - slope_ref,
           dslope_std = ifelse(is.finite(sigma) & sigma > 0, dslope / sigma, NA_real_),
           Fac = factor(Fac, levels = sp$facets),
           Trait = factor(Trait, levels = trait_order))

  list(reg = reg, ei = ei, lf = lf)
}

#Pick the quantity and the P the asterisks come from.
with_mode <- function(reg, mode) {
  reg$Fill   <- if (mode == "dslope") reg$dslope_std else reg$ei_slope
  reg$P_mark <- if (mode == "dslope") reg$P_con_swhc else reg$P_con
  reg %>%
    group_by(Fac, Trait) %>%
    mutate(FDR = p.adjust(P_mark, method = "BH")) %>%    # across timepoints
    ungroup()
}

############################################################################################################
#Geometry of the half ring.
#
#coord_polar always runs a discrete x right round the circle, so the half is made
#by giving the scale twice the slots it needs and leaving the second half empty.
#Clockwise from 12 o'clock: the first quadrant's block (0-90 degrees), the empty
#half, then the second quadrant's block (270-360, i.e. upper left). Each block
#carries a gap wedge at both ends, so the quadrants are separated at 12 o'clock
#and meet the empty half the same way on both sides.

n_traits <- length(trait_order)
n_gap    <- 3
block    <- n_traits + 2 * n_gap          # 20
n_slots  <- 4 * block                     # 80

fac_slots <- function(s) c(paste0(s, "__gapA", seq_len(n_gap)),
                           paste(s, trait_order, sep = ": "),
                           paste0(s, "__gapB", seq_len(n_gap)))

make_slots <- function(facets) {
  sl <- c(fac_slots(facets[1]),
          paste0("__empty", seq_len(2 * block)),
          fac_slots(facets[2]))
  stopifnot(length(sl) == n_slots)
  sl
}

############################################################################################################

run_contrast <- function(cname) {

  sp <- SPECS[[cname]]
  slot_levels <- make_slots(sp$facets)
  fac_first <- setNames(match(paste(sp$facets, trait_order[1], sep = ": "), slot_levels),
                        sp$facets)
  fac_mid   <- fac_first + (n_traits - 1) / 2
  lvl_col   <- setNames(LVL_COL, c(sp$ref, sp$cmp))
  lvl_shape <- setNames(c(16, 17), c(sp$ref, sp$cmp))
  ratio_lab <- paste(sp$cmp, "/", sp$ref)

  cat("\n############ ", sp$fig, " - contrast: ", sp$what, "\n", sep = "")

  #One build per replication unit, shared by both quantities
  BUILDS <- lapply(UNITS, function(u) {
    cat("\n-- building", cname, u$tag, "\n")
    build_reg(read_unit(u$ps2), sp)
  })
  names(BUILDS) <- vapply(UNITS, function(u) u$tag, character(1))

  for (u_i in seq_along(UNITS)) {
    for (mode in MODES) {

      u   <- UNITS[[u_i]]
      bld <- BUILDS[[u_i]]
      reg <- with_mode(bld$reg, mode)
      ei  <- bld$ei
      lf  <- bld$lf

      STEM <- paste0(sp$fig, "_",
                     if (mode == "dslope") "DeltaPlasticity" else sp$ratio_stem)
      FILL_NAME <- if (mode == "dslope")
        paste0("change in plasticity\nslope(", sp$cmp, ") - slope(", sp$ref,
               ")\nresidual SD per % SWHC") else
        paste0("slope of\n", if (INDEX_SCALE == "log2")
               paste0("log2(", ratio_lab, ")") else paste0("(", ratio_lab, ")"),
               "\nper % SWHC")

      #The timepoint the insets show: first at which the marked P drops below
      #0.05 in at least one quadrant, taken from the LINE MEANS so that both
      #replication versions of the same figure illustrate the same cell;
      #failing that, the smallest P.
      ref_reg <- with_mode(BUILDS[["lineMeans"]]$reg, mode)
      cand <- ref_reg %>% filter(Trait == INSET_TRAIT, is.finite(P_mark)) %>% arrange(TOE)
      sig  <- cand %>% filter(P_mark < 0.05)
      TOE_I <- if (nrow(sig)) sig$TOE[1] else cand$TOE[which.min(cand$P_mark)]

      cat("\n== ", STEM, " ", u$tag, " | inset timepoint TOE ", TOE_I, " h for ",
          INSET_TRAIT, if (nrow(sig)) "" else " (no P < 0.05; smallest instead)",
          "\n", sep = "")

      #One CSV per contrast and replication unit: it carries BOTH quantities and
      #both P values, so it does not need writing once per mode.
      if (mode == MODES[1])
        write_csv(bld$reg, paste0("data/", sp$fig, "_contrast_regressions_",
                                  u$tag, ".csv"))

      reg$Slot <- factor(paste(reg$Fac, reg$Trait, sep = ": "), levels = slot_levels)
      max_toe <- max(reg$TOE, na.rm = TRUE)

      #geom_tile takes its height from the SMALLEST gap between TOE values, and
      #the grid is 12 h everywhere except across the run break (408 -> 411 = 3 h),
      #so left alone every tile comes out a quarter of its proper size.
      toe_steps <- diff(sort(unique(reg$TOE)))
      tile_h <- as.numeric(names(sort(table(toe_steps), decreasing = TRUE))[1])

      divider_df <- data.frame(x = c(0.5, block + 0.5, 3 * block + 0.5),
                               y = 0, yend = max_toe)

      #Quadrant names along the arc, centred on each quadrant's own traits
      fac_label_df <- data.frame(x = unname(fac_mid), label = sp$facets,
                                 y = max_toe * 1.16)
      fac_label_df$angle <- -360 * (fac_label_df$x - 0.5) / n_slots
      fac_label_df$angle <- ((fac_label_df$angle + 90) %% 180) - 90

      toe_breaks  <- c(200, 300, 400, 500, 600, 700)
      toe_breaks  <- toe_breaks[toe_breaks <= max_toe]
      toe_label_x <- match(paste0(sp$facets[1], "__gapB", ceiling(n_gap / 2)), slot_levels)
      toe_label_df <- data.frame(x = toe_label_x, y = toe_breaks,
                                 label = paste0(toe_breaks, "h"))

      trait_label_df <- data.frame(
        Slot = factor(paste(rep(sp$facets, each = n_traits),
                            rep(trait_order, times = 2), sep = ": "),
                      levels = slot_levels),
        trait_label = rep(unname(trait_label[trait_order]), times = 2))
      trait_label_df$slot_index <- as.integer(trait_label_df$Slot)
      trait_label_df$angle <- 90 - 360 * (trait_label_df$slot_index - 0.5) / n_slots
      trait_label_df$hjust <- ifelse(trait_label_df$angle < -90, 1, 0)
      trait_label_df$angle <- ifelse(trait_label_df$angle < -90,
                                     trait_label_df$angle + 180, trait_label_df$angle)
      trait_label_df$y <- max_toe * 1.05

      #Rings drawn as ARCS over the two occupied quadrants only. geom_hline in a
      #polar coord sweeps the whole 360 degrees, which puts gridlines and the
      #reference rings across the empty lower half and under the insets.
      #geom_segment is no use either: in polar coordinates it draws the straight
      #chord between its ends, not the arc. So each arc is a geom_path over a
      #fine sequence of slots at constant radius.
      arc_df <- function(ys) {
        rngs <- list(c(0.5, block + 0.5), c(3 * block + 0.5, 4 * block + 0.5))
        do.call(rbind, lapply(seq_along(rngs), function(i)
          do.call(rbind, lapply(ys, function(yy)
            data.frame(x = seq(rngs[[i]][1], rngs[[i]][2], length.out = 80),
                       y = yy, grp = paste0(i, "_", yy))))))
      }

      hl <- reg %>% filter(Trait == INSET_TRAIT, TOE == TOE_I)
      spoke <- data.frame(x = match(paste(sp$facets, INSET_TRAIT, sep = ": "), slot_levels),
                          y = 0, yend = TOE_I - tile_h / 2)

      lim <- quantile(abs(reg$Fill), FILL_QUANTILE, na.rm = TRUE)
      cat("   colour limits +/-", signif(lim, 3), "| squished:",
          sum(abs(reg$Fill) > lim, na.rm = TRUE), "cells\n")

      ring <- ggplot(reg, aes(x = Slot, y = TOE, fill = Fill)) +
        geom_tile(height = tile_h) +
        geom_path(data = arc_df(toe_breaks), aes(x = x, y = y, group = grp),
                  inherit.aes = FALSE, color = "grey55", linewidth = 0.3, alpha = 0.45) +
        geom_path(data = arc_df(RUN_BREAK), aes(x = x, y = y, group = grp),
                  inherit.aes = FALSE, color = "grey25", linewidth = 0.4,
                  linetype = "dashed") +
        geom_path(data = arc_df(TOE_I), aes(x = x, y = y, group = grp),
                  inherit.aes = FALSE, color = "red", linewidth = 0.4,
                  linetype = "dotdash", alpha = 0.8) +
        geom_point(data = subset(reg, FDR < 0.05), shape = 8, size = 1.05,
                   color = "white", stroke = 0.65) +
        geom_point(data = subset(reg, FDR < 0.05), shape = 8, size = 0.7,
                   color = "red", stroke = 0.35) +
        geom_segment(data = divider_df, aes(x = x, xend = x, y = y, yend = yend),
                     inherit.aes = FALSE, color = "white", linewidth = 0.9) +
        geom_segment(data = spoke, aes(x = x, xend = x, y = y, yend = yend),
                     inherit.aes = FALSE, color = "white", linewidth = 1.6) +
        geom_segment(data = spoke, aes(x = x, xend = x, y = y, yend = yend),
                     inherit.aes = FALSE, color = "red", linewidth = 0.7) +
        geom_tile(data = hl, aes(x = Slot, y = TOE), inherit.aes = FALSE,
                  height = tile_h, fill = NA, color = "red", linewidth = 1.0) +
        geom_text(data = fac_label_df, aes(x = x, y = y, label = label, angle = angle),
                  inherit.aes = FALSE, fontface = "bold", size = 6.2) +
        geom_text(data = toe_label_df, aes(x = x, y = y, label = label),
                  inherit.aes = FALSE, size = 3.4, hjust = 0.5, color = "grey30") +
        geom_text(data = trait_label_df, aes(x = Slot, y = y, label = trait_label,
                                             angle = angle, hjust = hjust),
                  inherit.aes = FALSE, size = 3.6, color = "grey15") +
        scale_fill_gradient2(low = NEG_COL, mid = "grey96", high = POS_COL,
                             midpoint = 0, limits = c(-lim, lim),
                             oob = scales::squish, name = FILL_NAME,
                             na.value = "grey90") +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(limits = c(0, max_toe * YLIM_PAD), expand = c(0, 0)) +
        coord_polar(theta = "x") +
        labs(title = paste0(sp$what, "\n", u$label)) +
        theme_minimal() +
        guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                      barwidth = unit(13, "lines"),
                                      barheight = unit(0.8, "lines"))) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.title = element_blank(),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                        lineheight = 1.15),
              legend.position = "inside",        # the lower half is empty
              legend.position.inside = c(0.5, 0.10),
              legend.direction = "horizontal",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 9),
              plot.margin = margin(0, 0, 0, 0))

      #Insets. In "dslope" mode the two REACTION NORMS whose slopes the cell
      #differences; in "ratio" mode, the index against SWHC.
      mini_norm <- function(f) {
        d <- lf  %>% filter(Fac == f, Trait == INSET_TRAIT, TOE == TOE_I, is.finite(value))
        r <- reg %>% filter(Fac == f, Trait == INSET_TRAIT, TOE == TOE_I)
        sub <- if (nrow(r) == 1)
          paste0("slope ", sp$ref, " = ", formatC(r$slope_ref, format = "g", digits = 2),
                 "    ", sp$cmp, " = ", formatC(r$slope_cmp, format = "g", digits = 2),
                 "\ndifference = ", formatC(r$dslope_std, format = "f", digits = 2),
                 " SD per % SWHC,  P = ", formatC(r$P_con_swhc, format = "g", digits = 2)) else ""
        ggplot(d, aes(x = SWHC, y = value, color = Con, shape = Con)) +
          geom_point(size = 1.9, alpha = 0.8) +
          geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 1.0) +
          scale_color_manual(values = lvl_col, drop = FALSE) +
          scale_shape_manual(values = lvl_shape, drop = FALSE) +
          scale_x_continuous(breaks = c(20, 60, 100), limits = c(12, 108)) +
          labs(title = paste0(f, " - ", INSET_TRAIT, " at ", TOE_I, " h"),
               subtitle = sub, x = "% SWHC", y = paste0(INSET_TRAIT, " (a.u.)")) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                plot.title = element_text(size = 11, face = "bold"),
                plot.subtitle = element_text(size = 8, color = "grey25", lineheight = 1.05),
                axis.title = element_text(size = 9, face = "bold"),
                axis.text = element_text(size = 8),
                legend.position = "inside", legend.position.inside = c(0.99, 0.99),
                legend.justification = c(1, 1), legend.title = element_blank(),
                legend.text = element_text(size = 8),
                legend.background = element_rect(fill = "white", color = NA),
                legend.key.size = unit(0.7, "lines"),
                plot.background = element_rect(fill = "white", color = "grey70",
                                               linewidth = 0.3),
                plot.margin = margin(4, 6, 4, 4))
      }

      mini_index <- function(f) {
        d <- ei  %>% filter(Fac == f, Trait == INSET_TRAIT, TOE == TOE_I, is.finite(EI))
        r <- reg %>% filter(Fac == f, Trait == INSET_TRAIT, TOE == TOE_I)
        lc <- if (nrow(r) == 1 && is.finite(r$ei_slope) && r$ei_slope < 0) NEG_COL else POS_COL
        ref0 <- if (INDEX_SCALE == "log2") 0 else 1
        sub <- if (nrow(r) == 1)
          paste0("slope = ", formatC(r$ei_slope, format = "g", digits = 2),
                 "   R2 = ", formatC(r$ei_R2, format = "f", digits = 2),
                 "\nP(slope) = ", formatC(r$ei_P, format = "g", digits = 2),
                 "   P(", sp$cvar, ") = ", formatC(r$P_con, format = "g", digits = 2)) else ""
        ggplot(d, aes(x = SWHC, y = EI)) +
          geom_hline(yintercept = ref0, linetype = "dotted", color = "grey45") +
          geom_point(size = 2.2, alpha = 0.85, color = "grey25") +
          geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
                      color = lc, linewidth = 0.9) +
          scale_x_continuous(breaks = c(20, 60, 100), limits = c(12, 108)) +
          labs(title = paste0(f, " - ", INSET_TRAIT, " at ", TOE_I, " h"),
               subtitle = sub, x = "% SWHC",
               y = if (INDEX_SCALE == "log2") paste0("log2(", ratio_lab, ")") else ratio_lab) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                plot.title = element_text(size = 11, face = "bold"),
                plot.subtitle = element_text(size = 8.5, color = "grey25", lineheight = 1.05),
                axis.title = element_text(size = 9, face = "bold"),
                axis.text = element_text(size = 8),
                plot.background = element_rect(fill = "white", color = "grey70",
                                               linewidth = 0.3),
                plot.margin = margin(4, 6, 4, 4))
      }

      inset_for <- function(f) if (mode == "dslope") mini_norm(f) else mini_index(f)

      Fig <- ring
      if (DRAW_INSETS) {
        #The second quadrant is the upper LEFT one, so its inset goes lower left;
        #the first is upper right, so its inset goes lower right. Both sit
        #strictly below the half way line, clear of every tile.
        Fig <- Fig +
          inset_element(inset_for(sp$facets[2]),
                        left = INSET_BOX[1], bottom = INSET_BOX[2],
                        right = INSET_BOX[3], top = INSET_BOX[4], align_to = "panel") +
          inset_element(inset_for(sp$facets[1]),
                        left = 1 - INSET_BOX[3], bottom = INSET_BOX[2],
                        right = 1 - INSET_BOX[1], top = INSET_BOX[4], align_to = "panel")
      }

      f_out <- file.path(OUT_DIR, paste0(STEM, "_", u$tag, ".pdf"))
      #coord_polar forces a square panel, so the page stays square even though
      #only the upper half carries tiles; the lower half holds the key and insets.
      #NB: no bare `Fig` after this - under Rscript that writes Rplots.pdf.
      ggsave(f_out, Fig, width = 14, height = 14)
      cat("   wrote", f_out, "\n")

      cat("   asterisks (", if (mode == "dslope") paste0(sp$cvar, " x SWHC") else sp$cvar,
          "), FDR < 0.05: ", sum(reg$FDR < 0.05, na.rm = TRUE), " of ",
          sum(!is.na(reg$P_mark)), "\n", sep = "")
      print(reg %>% mutate(s = FDR < 0.05) %>%
              group_by(Trait, Fac) %>%
              summarise(n = sum(s, na.rm = TRUE), .groups = "drop") %>%
              pivot_wider(names_from = Fac, values_from = n))
      cat("   direction where it survives FDR (neg =", sp$cmp,
          "lower / less responsive):\n")
      print(reg %>% filter(FDR < 0.05) %>%
              group_by(Fac) %>%
              summarise(neg = sum(Fill < 0, na.rm = TRUE),
                        pos = sum(Fill > 0, na.rm = TRUE), .groups = "drop"))
    }
  }
  invisible(NULL)
}

for (cn in CONTRASTS) run_contrast(cn)
