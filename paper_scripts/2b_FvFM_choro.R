##################################################################################
## Fig 2 : Fv/FM and Chlorophyl graphs with P-values
## Get select graphs for 20 and 100%
## Four-color version: California 2014/2010 in reds, Oregon 2014/2010 in blues
###################################################################################
#Import libraries
library(tidyverse)

#Import Data
PS2_all <- read_csv("data/PS2_all.csv")
PS2_all <- PS2_all %>% mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
                              Timepoint = factor(Timepoint, levels = c("Early", "Late")))

#Colors: same four Site x Year groups and colors as Fig 1C (legend order too)
grp_levels <- c("California 2014", "California 2010", "Oregon 2014", "Oregon 2010")
grp_cols   <- c("California 2014" = "#B2182B", "California 2010" = "#EF6B4F",
                "Oregon 2014"     = "#2166AC", "Oregon 2010"     = "#6BAED6")

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

# Plot data: 20% and 100% only, with the Site x Year group used for color
plot_data <- PS2_all %>%
  dplyr::filter(Treatment %in% c("20%", "100%")) %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "100%")),
         Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "California")),
         grp  = factor(paste(Site, Year), levels = grp_levels))

###################################### Fv.Fm: 20% vs 100%, Both Sites ######################################
############################################################################################################
p_fvfm <- plot_data %>%
  ggplot(aes(x = TOE, y = Fv.Fm, color = grp, fill = grp)) +
  geom_point() +
  scale_color_manual(values = grp_cols, breaks = grp_levels, name = NULL) +
  scale_fill_manual(values = grp_cols,  breaks = grp_levels, name = NULL) +
  geom_smooth(alpha = 0.25) +
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

###################################### ChlIdx: 20% vs 100%, Both Sites #####################################
############################################################################################################

p_chlidx <- plot_data %>%
  ggplot(aes(x = TOE, y = ChlIdx, color = grp, fill = grp)) +
  geom_point() +
  scale_color_manual(values = grp_cols, breaks = grp_levels, name = NULL) +
  scale_fill_manual(values = grp_cols,  breaks = grp_levels, name = NULL) +
  geom_smooth(alpha = 0.25) +
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

# Legend inside the right-hand plot; four entries now, so text and keys are smaller than the 2-color version
# >>> ADJUST THESE NUMBERS (x, y) TO REPOSITION THE LEGEND, AND THE SIZES TO RESIZE IT <<<
p_chlidx_combined <- p_chlidx + theme(
  legend.position = "inside",
  legend.position.inside = c(0.99, 0.49),   # <-- tweak x/y (0-1 range) here
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = NA),
  legend.text = element_text(size = 18, face = "bold"),
  legend.key.size = unit(2.2, "lines"),
  legend.key.height = unit(1.5, "lines")
)

############################################################################################################
#P-value and R2 graphs (unchanged: these show model terms, not Site x Year groups)
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
ytitle_gap <- 2

#Shared theme: no grid lines, bold text, text sizes matched to the top figures, zero vertical margins
fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
  axis.text.y = element_text(size = 22, face = "bold", margin = margin(r = 4)),
  axis.title = element_text(color = "black", size = 26, face = "bold"),
  axis.title.y = element_text(color = "black", size = 26, face = "bold", margin = margin(r = ytitle_gap)),
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
vline_410 <- function() {
  geom_vline(xintercept = 410, color = "black", linetype = "dotted",
             linewidth = 0.5, alpha = 0.6, show.legend = FALSE)
}

#Term colors, matched to Fig 1a
term_colors <- c("Site*Year" = "royalblue3",
                 "Treatment" = "#D95F02",
                 "Site"      = "#7570B3",
                 "Year"      = "#E7298A")

#Middle: -log10(P-value) across TOE (black dashed line = P of 0.05, same as Fig 1a)
make_P_plot <- function(toe_stats, show_legend = TRUE) {
  p <- ggplot(toe_stats, aes(x = TOE, y = negLogP, color = Term)) +
    vline_410() +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "11", linewidth = 0.8) +
    geom_line(linewidth = 0.9, linetype = "solid") +
    scale_color_manual(values = term_colors) +
    scale_y_continuous(breaks = seq(0, 10, 2.5)) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(x = NULL, y = "\u2010log(P)") +
    fig_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  if (show_legend) {
    p + theme(
      legend.position = "inside",
      legend.position.inside = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 14, face = "bold"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.width = unit(2, "lines")
    )
  } else {
    p + theme(legend.position = "none")
  }
}

#Bottom: adjusted R2 of the full model across TOE (only graph with the x-axis title)
R2_title_hjust <- 0.56

make_R2_plot <- function(toe_stats, breaks = ggplot2::waiver()) {
  toe_stats %>% distinct(TOE, adjR2) %>%
    ggplot(aes(x = TOE, y = adjR2)) +
    vline_410() +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(breaks = breaks) +
    expand_limits(y = 0.3) +
    labs(x = "Time of Experiment (h)", y = expression(bold(R^2))) +
    fig_theme +
    theme(axis.title.y = element_text(hjust = R2_title_hjust))
}

P_fvfm    <- make_P_plot(toe_stats_fvfm,   show_legend = FALSE)
P_chlidx  <- make_P_plot(toe_stats_chlidx, show_legend = TRUE)
R2_fvfm   <- make_R2_plot(toe_stats_fvfm,   breaks = c(0.3, 0.5, 0.7, 0.9))
R2_chlidx <- make_R2_plot(toe_stats_chlidx, breaks = c(0.3, 0.5, 0.7, 0.9))

############################################################################################################
###################################### Top row: A / B labels, no x-axis title ##############################
############################################################################################################
tag_x_A <- 0.025
tag_x_B <- 0.020

tag_theme <- function(x) {
  theme(
    plot.tag = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 1),
    plot.tag.location = "plot",
    plot.tag.position = c(x, 1),
    axis.title.x = element_blank()
  )
}

# Put the dotted line at TOE = 410 behind the points and curves of a top graph
add_vline_back <- function(p) {
  p$layers <- append(list(vline_410()), p$layers)
  p
}

p_A <- add_vline_back(p_fvfm_nolegend) + labs(tag = "A") + tag_theme(tag_x_A)

p_B <- add_vline_back(p_chlidx_combined) + labs(tag = "B") + tag_theme(tag_x_B) +
  theme(axis.title.y = element_text(vjust = 1.9))

############################################################################################################
###################################### Combine: flat 4 rows x 2 columns grid ###############################
############################################################################################################
Fig2_full <- wrap_plots(p_A, p_B,
                        P_fvfm, P_chlidx,
                        plot_spacer(), plot_spacer(),
                        R2_fvfm, R2_chlidx,
                        ncol = 2) +
  plot_layout(heights = c(4, 1, 0.075, 0.5))
Fig2_full

ggsave("paper_graphs/Fig2b_fvfm_chlidx_stats.pdf",
       plot = Fig2_full, width = 20, height = 16, units = "in")