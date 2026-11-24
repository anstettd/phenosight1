##################################################################################
## Fig 2 : Fv/FM and Chlorophyll graphs with P-values
## For each response variable: all-treatments plot (top) + -log10(P) (middle)
## + adjusted R2 (bottom), stacked like Fig 1A. Dashed vertical line at TOE = 410.
###################################################################################
# Import libraries
library(tidyverse)
library(car)
library(patchwork)   # needs ggplot2 >= 3.5.0 for legend.position.inside

# Import Data
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year      = factor(Year, levels = c(2010, 2014)),
         Site      = factor(Site, levels = c("Oregon", "California")))

# Colors by Year
year_colors <- c("2010" = "skyblue3", "2014" = "#FF7700")

# Vertical dashed line position (Time of Experiment)
toe_line <- 410   # >>> adjust here <<<

# Output folders
dir.create("paper_graphs/Supplemental/all_treatments", recursive = TRUE, showWarnings = FALSE)

############################################################################################################
# Shared theme (used by top, middle and lower panels)
fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  panel.spacing = unit(0.8, "lines"),
  axis.text = element_text(size = 13, face = "bold"),
  axis.title = element_text(color = "black", size = 16, face = "bold"),
  axis.ticks.length = unit(3, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 14, face = "bold"),
  strip.text = element_text(size = 14, face = "bold"),
  strip.background = element_blank(),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
)

# Pull y titles of the middle and lower panels in toward the axis (negative = inward)
ytitle_gap_lower <- -4

# Term colors / labels for the p-value lines
term_cols   <- c("Site:Year" = "royalblue3", "Treatment" = "#D95F02",
                 "Site" = "#7570B3", "Year" = "#E7298A")
term_labels <- c("Site:Year" = "Site*Year", "Treatment" = "Treatment",
                 "Site" = "Site", "Year" = "Year")

# Dashed vertical line used on every panel
vline_410 <- geom_vline(xintercept = toe_line, color = "black", linetype = "11",
                        linewidth = 0.5, alpha = 0.6, show.legend = FALSE)

# Fit the model at one TOE; returns P for each term and adjusted R2
fit_one <- function(d, var) {
  m <- lm(as.formula(paste(var, "~ Site*Year + Treatment")), data = d,
          contrasts = list(Site = "contr.sum", Year = "contr.sum", Treatment = "contr.sum"))
  a <- as.data.frame(Anova(m, type = 3)) %>%
    rownames_to_column("Term") %>%
    filter(Term %in% c("Site", "Year", "Treatment", "Site:Year")) %>%
    select(Term, P = `Pr(>F)`)
  a$adjR2 <- summary(m)$adj.r.squared
  a
}

############################################################################################################
###################################### Multi-Spectral Data  ################################################
############################################################################################################

response_vars <- c("Fv.Fm", "Fq.Fm", "NPQ", "qP", "qN", "qL", "qI",
                   "phiNO", "phiNPQ", "npq.t", "ChlIdx", "AriIdx", "NDVI")

for (var in response_vars) {
  
  # ---- stats at each TOE ----
  toe_stats <- PS2_all %>%
    group_by(TOE) %>%
    group_modify(~ tryCatch(fit_one(.x, var),
                            error = function(e) tibble(Term = NA_character_, P = NA_real_, adjR2 = NA_real_))) %>%
    ungroup() %>%
    filter(!is.na(Term)) %>%
    mutate(Term = factor(Term, levels = c("Site:Year", "Treatment", "Site", "Year")),
           negLogP = -log10(P))
  
  # ---- top: all treatments ----
  p_top <- ggplot(PS2_all, aes(x = TOE, y = .data[[var]], color = Year)) +
    geom_point() +
    scale_color_manual(values = year_colors, name = "Year") +
    geom_smooth() +
    vline_410 +
    facet_grid(Treatment ~ Site,
               labeller = labeller(Site = c("Oregon" = "Oregon", "California" = "Southern California"))) +
    labs(x = NULL, y = var) +
    fig_theme +
    theme(legend.position = "none")
  
  # ---- middle: -log10(P); legend lives in the R2 panel ----
  p_mid <- ggplot(toe_stats, aes(x = TOE, y = negLogP, color = Term)) +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "11", linewidth = 0.8) +
    vline_410 +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = term_cols, labels = term_labels) +
    labs(x = NULL, y = "-log(P)") +
    fig_theme +
    theme(legend.position = "none",
          axis.title.y = element_text(color = "black", size = 16, face = "bold",
                                      margin = margin(r = ytitle_gap_lower)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # ---- bottom: adjusted R2 (carries the legend for the P-value lines) ----
  p_r2 <- toe_stats %>% distinct(TOE, adjR2) %>%
    ggplot(aes(x = TOE, y = adjR2)) +
    vline_410 +
    geom_line(linewidth = 0.9) +
    geom_line(data = toe_stats %>% mutate(adjR2 = NA_real_),   # invisible; only creates legend entries
              aes(color = Term), linewidth = 0.9, na.rm = TRUE) +
    scale_color_manual(values = term_cols, labels = term_labels) +
    labs(x = "Time of Experiment (h)", y = expression(bold(R^2))) +
    fig_theme +
    theme(axis.title.y = element_text(color = "black", size = 16, face = "bold",
                                      margin = margin(r = ytitle_gap_lower)),
          legend.position = "inside",
          legend.position.inside = c(0.995, 0.93),   # upper right of R2; adjust here
          legend.justification = c(1, 1),
          legend.background = element_blank(),
          legend.text = element_text(size = 10, face = "bold"),
          legend.key.size = unit(0.85, "lines"),
          legend.key.width = unit(1.2, "lines"))
  
  fig <- (p_top / p_mid / p_r2) + plot_layout(heights = c(5, 1, 1))
  
  ggsave(paste0("paper_graphs/Supplemental/all_treatments/", tolower(var), "_all.pdf"),
         plot = fig, width = 12, height = 18, units = "in")
  
  message("Done: ", var)
}