##################################################################################
## Fig 2 : Fv/FM and Chlorophyl graphs with P-values
## Get select graphs for 20 and 100%
## 
###################################################################################
#Import libraries
library(tidyverse)

#Import Data
PS2_all <- read_csv("data/PS2_all.csv")
PS2_all <- PS2_all %>% mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
                              Timepoint = factor(Timepoint, levels = c("Early", "Late")))

#Colors: solid line/point colors by Year, and lighter fill colors for geom_smooth ribbons
year_colors <- c("2010" = "skyblue3", "2014" = "#FF7700")
year_fill_colors <- c("2010" = "#8FC1E3", "2014" = "#FFA858")


###################################################################################
#Linear Models

library(car)
library(patchwork)

# Data for the models: 20% and 100% only (same subset as the top graphs)
PS2_sub <- PS2_all %>%
  dplyr::filter(Treatment %in% c("20%", "100%")) %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "100%")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California")))

#Overall models
mod_FvFm <- lm(Fv.Fm ~ Site * Year + Treatment, data = PS2_sub)
Anova(mod_FvFm)

mod_ChlIdx <- lm(ChlIdx ~ Site * Year + Treatment, data = PS2_sub)
Anova(mod_ChlIdx)



###################################### Fv.Fm: 20% vs 100%, Both Sites ######################################
############################################################################################################
p_fvfm <- PS2_all %>%
  dplyr::filter(Treatment %in% c("20%", "100%")) %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "100%")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California"))) %>%
  ggplot(aes(x = TOE, y = Fv.Fm, color = Year, fill = Year)) +
  geom_point() +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values = year_fill_colors, name = "Year") +
  geom_smooth() +
  facet_grid(Treatment ~ Site) +
  coord_cartesian(ylim = c(0.75, NA)) +
  labs(x = "Time of Experiment", y = "Quantum Yield of PSII (Fv/Fm)") +
  theme_bw() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(color = "black", size = 26, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 26, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 26, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_text(size = 24, face = "bold"),
    strip.background = element_blank()
  )

#ggsave("paper_graphs/multi_join/select/fvfm_20_100.pdf",
#       plot = p_fvfm, width = 10, height = 10, units = "in")


###################################### ChlIdx: 20% vs 100%, Both Sites #####################################
############################################################################################################

p_chlidx <- PS2_all %>%
  dplyr::filter(Treatment %in% c("20%", "100%")) %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "100%")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California"))) %>%
  ggplot(aes(x = TOE, y = ChlIdx, color = Year, fill = Year)) +
  geom_point() +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values = year_fill_colors, name = "Year") +
  geom_smooth() +
  facet_grid(Treatment ~ Site) +
  labs(x = "Time of Experiment", y = "Chlorophyll Index") +
  theme_bw() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(color = "black", size = 26, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 26, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_text(size = 24, face = "bold"),
    strip.background = element_blank()
  )

#ggsave("paper_graphs/chlidx_20_100.pdf",
#       plot = p_chlidx, width = 10, height = 10, units = "in")

library(cowplot)


###################################### Combined: Side-by-side (wide) #######################################
############################################################################################################

# Drop the legend and the Treatment facet strip from the left-hand plot so
# only one shared legend and one set of Treatment labels (on the right plot) appear
p_fvfm_nolegend <- p_fvfm + theme(
  legend.position = "none",
  strip.background.y = element_blank(),
  strip.text.y = element_blank()
)

# Move the legend into the empty space in the upper part of the lower-right
# facet panel (Treatment = 100%, Site = SCalifornia), instead of it sitting
# off to the side. This is only applied for the combined figure -- the
# standalone chlidx_20_100.pdf keeps its normal side legend.
# >>> ADJUST THESE TWO NUMBERS (x, y) TO REPOSITION THE LEGEND <<<
p_chlidx_combined <- p_chlidx + theme(
  legend.position = "inside",
  legend.position.inside = c(0.99, 0.49),   # <-- tweak x/y (0-1 range) here
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = NA),
  legend.text = element_text(size = 22, face = "bold"),
  legend.key.size = unit(3, "lines"),
  legend.key.height = unit(2.4, "lines")
)

p_combined_wide <- plot_grid(p_fvfm_nolegend, p_chlidx_combined, nrow = 1)

ggsave("paper_graphs/fvfm_chlidx_2x4.pdf",
       plot = p_combined_wide, width = 20, height = 10, units = "in")

############################################################################################################
#P-value and R2 graphs
# Fit the same model at each TOE; keep P-value per term and adjusted R2 of the full model
get_toe_stats <- function(data, response) {
  data %>%
    group_by(TOE) %>%
    group_modify(~ {
      out <- tryCatch({
        m <- lm(as.formula(paste0("`", response, "` ~ Site * Year + Treatment")), data = .x)
        a <- Anova(m)
        tibble(Term = rownames(a),
               p = a$`Pr(>F)`,
               adjR2 = summary(m)$adj.r.squared)
      }, error = function(e) {
        tibble(Term = c("Site", "Year", "Treatment", "Site:Year", "Residuals"),
               p = NA_real_, adjR2 = NA_real_)
      })
      out
    }) %>%
    ungroup() %>%
    dplyr::filter(Term != "Residuals") %>%
    mutate(negLogP = -log10(p),
           Term = dplyr::recode(Term, "Site:Year" = "Site*Year"),
           Term = factor(Term, levels = c("Site*Year", "Site", "Year", "Treatment")))
}

toe_stats_fvfm   <- get_toe_stats(PS2_sub, "Fv.Fm")
toe_stats_chlidx <- get_toe_stats(PS2_sub, "ChlIdx")

############################################################################################################
###################################### Formatting (middle + bottom) ########################################
############################################################################################################
legend_space <- 5.5   # right margin; 5.5 is the ggplot default so the right edges match the top graphs

# >>> Gap (pt) between the y-axis title and the tick text on the P-value and R2 graphs <<<
# Smaller = title sits closer to the tick text (was 12); raise it to move the title back left
ytitle_gap <- 2

#Shared theme: no grid lines, bold text, text sizes matched to the top figures, zero vertical margins
fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
  axis.text.y = element_text(size = 22, face = "bold", margin = margin(r = 4)),      # <-- gap between tick text and plot
  axis.title = element_text(color = "black", size = 26, face = "bold"),
  axis.title.y = element_text(color = "black", size = 26, face = "bold", margin = margin(r = ytitle_gap)),  # <-- gap between title and tick text
  axis.ticks.length = unit(2, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 8, face = "bold"),
  legend.key.size = unit(0.8, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  legend.justification = c(0, 0.5),
  strip.text = element_text(size = 24, face = "bold"),
  strip.background = element_blank(),
  plot.margin = margin(t = 0, r = legend_space, b = 0, l = 5)
)

# Thin dotted black vertical line at TOE = 410 (used on all six graphs)
# >>> linewidth = thickness, alpha = how dark (1 = solid black, lower = fainter) <<<
vline_410 <- function() {
  geom_vline(xintercept = 410, color = "black", linetype = "dotted",
             linewidth = 0.5, alpha = 0.6, show.legend = FALSE)
}

#Term colors, matched to Fig 1a (Site*Year = royalblue3; others keep their Dark2 colors from Fig 1a)
term_colors <- c("Site*Year" = "royalblue3",
                 "Treatment" = "#D95F02",
                 "Site"      = "#7570B3",
                 "Year"      = "#E7298A")

#Middle: -log10(P-value) across TOE (black dashed line = P of 0.05, same as Fig 1a)
# show_legend = TRUE puts the legend inside the upper right corner; FALSE removes it
make_P_plot <- function(toe_stats, show_legend = TRUE) {
  p <- ggplot(toe_stats, aes(x = TOE, y = negLogP, color = Term)) +
    vline_410() +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "11", linewidth = 0.8) +
    geom_line(linewidth = 0.9, linetype = "solid") +
    scale_color_manual(values = term_colors) +
    scale_y_continuous(breaks = seq(0, 10, 2.5)) +   # <-- ticks at 0, 2.5, 5, 7.5, 10
    coord_cartesian(ylim = c(0, 10)) +               # <-- upper limit of the axis
    labs(x = NULL, y = "\u2010log(P)") +     # \u2010 = short hyphen instead of the wide minus sign
    fig_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  if (show_legend) {
    p + theme(
      legend.position = "inside",
      legend.position.inside = c(0.99, 0.99),    # <-- tweak x/y (0-1 range) to reposition
      legend.justification = c(1, 1),            # legend's upper right corner sits at the position above
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 14, face = "bold"),   # <-- legend text size
      legend.key.size = unit(1.2, "lines"),
      legend.key.width = unit(2, "lines")
    )
  } else {
    p + theme(legend.position = "none")
  }
}

#Bottom: adjusted R2 of the full model across TOE (only graph with the x-axis title)
R2_title_hjust <- 0.56   # <-- 0.5 = centered; raise to move the R2 title up, lower to move it down

make_R2_plot <- function(toe_stats, breaks = ggplot2::waiver()) {
  toe_stats %>% distinct(TOE, adjR2) %>%
    ggplot(aes(x = TOE, y = adjR2)) +
    vline_410() +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(breaks = breaks) +
    expand_limits(y = 0.3) +     # forces the axis down to 0.3 so the 0.3 tick shows even with no data there
    labs(x = "Time of Experiment (h)", y = expression(bold(R^2))) +
    fig_theme +
    theme(axis.title.y = element_text(hjust = R2_title_hjust))
}

P_fvfm    <- make_P_plot(toe_stats_fvfm,   show_legend = FALSE)
P_chlidx  <- make_P_plot(toe_stats_chlidx, show_legend = TRUE)
R2_fvfm   <- make_R2_plot(toe_stats_fvfm,   breaks = c(0.3, 0.5, 0.7, 0.9))   # fewer ticks for Fv/Fm R2
R2_chlidx <- make_R2_plot(toe_stats_chlidx, breaks = c(0.3, 0.5, 0.7, 0.9))   # same ticks as Fv/Fm R2

############################################################################################################
###################################### Top row: A / B labels, no x-axis title ##############################
############################################################################################################
# New copies of the top graphs so the originals above stay untouched

# >>> A / B horizontal position (0-1 across each plot); only x changes, y stays at the top <<<
tag_x_A <- 0.025   # <-- raise to move A right, lower to move it left
tag_x_B <- 0.020   # <-- same for B

tag_theme <- function(x) {
  theme(
    plot.tag = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 1),  # same size/bolding as the y-axis titles
    plot.tag.location = "plot",
    plot.tag.position = c(x, 1),       # x = horizontal placement, 1 = top
    axis.title.x = element_blank()     # x-axis title removed from the top graphs
  )
}

# Put the dotted line at TOE = 410 behind the points and curves of a top graph
add_vline_back <- function(p) {
  p$layers <- append(list(vline_410()), p$layers)
  p
}

p_A <- add_vline_back(p_fvfm_nolegend) + labs(tag = "A") + tag_theme(tag_x_A)

# vjust on the Chlorophyll y-title raised from 1.7 to 1.9 to add a little space from the tick labels
p_B <- add_vline_back(p_chlidx_combined) + labs(tag = "B") + tag_theme(tag_x_B) +
  theme(axis.title.y = element_text(vjust = 1.9))

############################################################################################################
###################################### Combine: flat 4 rows x 2 columns grid ###############################
############################################################################################################
# One flat grid (not rows nested inside rows) so patchwork aligns the panel edges across every row
# Rows: top graphs, P-value, blank spacer, R2
# >>> Heights: top = 4, P-value = 1, spacer = 0.075 (white space between P-value and R2), R2 = 0.5 (half the P-value panel) <<<
# Raise the spacer number for more white space, lower it for less
Fig2_full <- wrap_plots(p_A, p_B,
                        P_fvfm, P_chlidx,
                        plot_spacer(), plot_spacer(),
                        R2_fvfm, R2_chlidx,
                        ncol = 2) +
  plot_layout(heights = c(4, 1, 0.075, 0.5))
Fig2_full

ggsave("paper_graphs/Fig2_fvfm_chlidx_stats.pdf",
       plot = Fig2_full, width = 20, height = 16, units = "in")