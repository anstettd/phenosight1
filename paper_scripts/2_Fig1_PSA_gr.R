##################################################################################
## Fig 1: Make Projected Shoot Area and Growth Rate Graphs
## In format for paper
##
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(ggh4x)
library(patchwork)   # needs ggplot2 >= 3.5.0 for legend.position.inside

#Import Data
PS2_all <- read_csv("data/PS2_all.csv") %>%
  mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
         Timepoint = factor(Timepoint, levels = c("Early", "Late")),
         Year = factor(Year, levels = c(2010, 2014)),            # 2010 on top row
         Site = factor(Site, levels = c("Oregon", "California"))) # Oregon on left
#GR_early3 <- read_csv(data/GR_early3)
#GR_late3 <- read_csv(data/GR_late3)

#Overall model
mod_PSA <- lm(PSA ~ Site * Year + Treatment, data = PS2_all)
Anova(mod_PSA)

############################################################################################################
#Graphs

legend_space <- 55   # white space (pt) on the right that holds the legends; same on all three graphs so panels align

#Shared theme: no grid lines, bold text, zero vertical margins so stacked graphs sit flush
fig_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 11, face = "bold"),
  axis.text.y = element_text(size = 11, face = "bold", margin = margin(r = 1)),
  axis.title = element_text(color = "black", size = 12, face = "bold"),
  axis.title.y = element_text(color = "black", size = 12, face = "bold", margin = margin(r = 2)),
  axis.ticks.length = unit(2, "pt"),
  legend.title = element_blank(),
  legend.text = element_text(size = 8, face = "bold"),
  legend.key.size = unit(0.8, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.y = unit(0, "pt"),
  legend.justification = c(0, 0.5),          # legend's left edge, vertical middle, sits at the position below
  strip.text = element_text(size = 13, face = "bold"),
  strip.background = element_blank(),
  plot.margin = margin(t = 0, r = legend_space, b = 0, l = 5)
)

#Fig 1a - Projected Shoot Area (PSA) over time, with P-values and adjusted R2 across TOE

#Top: PSA line graphs (mean +/- SE per Treatment)
PSA_summ <- PS2_all %>%
  group_by(Site, Year, Treatment, TOE) %>%
  summarise(mean_PSA = mean(PSA, na.rm = TRUE),
            se_PSA   = sd(PSA, na.rm = TRUE) / sqrt(sum(!is.na(PSA))),
            .groups = "drop")

PSA_plot <- ggplot(PSA_summ, aes(x = TOE, y = mean_PSA, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_ribbon(aes(ymin = mean_PSA - se_PSA, ymax = mean_PSA + se_PSA), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(Year ~ Site) +     # columns = Site, rows = Year (2010 above 2014)
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = NULL, y = "Projected Shoot Area") +
  fig_theme +
  theme(axis.text.y = element_text(size = 8, face = "bold", margin = margin(r = 1)),
        legend.position = "inside",
        legend.position.inside = c(1.07, 0.5))   # just right of the 2010/2014 labels, centered between them

#Fit the full model at each TOE; extract P-values (Type III) and adjusted R2
fit_one <- function(d) {
  m <- lm(PSA ~ Site * Year + Treatment, data = d,
          contrasts = list(Site = "contr.sum", Year = "contr.sum", Treatment = "contr.sum"))
  a <- as.data.frame(Anova(m, type = 3)) %>%
    rownames_to_column("Term") %>%
    filter(Term %in% c("Site", "Year", "Treatment", "Site:Year")) %>%
    select(Term, P = `Pr(>F)`)
  a$adjR2 <- summary(m)$adj.r.squared
  a
}

toe_stats <- PS2_all %>%
  group_by(TOE) %>%
  group_modify(~ tryCatch(fit_one(.x),
                          error = function(e) tibble(Term = NA_character_, P = NA_real_, adjR2 = NA_real_))) %>%
  ungroup() %>%
  filter(!is.na(Term)) %>%
  mutate(Term = factor(Term, levels = c("Site:Year", "Treatment", "Site", "Year")),
         negLogP = -log10(P))

#Middle: -log10(P-value) across TOE (dotted line = P of 0.05); shares the x-axis with the graph below
P_plot <- ggplot(toe_stats, aes(x = TOE, y = negLogP, color = Term, linetype = Term)) +
  geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "dotted", linewidth = 0.8) +
  geom_line(linewidth = 0.9) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = NULL, y = "-log(P)") +
  fig_theme +
  theme(legend.key.width = unit(1.2, "lines"),
        legend.position = "inside",
        legend.position.inside = c(1.01, 0.5),   # directly beside B
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Bottom: adjusted R2 of the full model across TOE (only graph with the x-axis title)
R2_plot <- toe_stats %>% distinct(TOE, adjR2) %>%
  ggplot(aes(x = TOE, y = adjR2)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Time of Experiment", y = expression(bold(R^2))) +
  fig_theme

#Combine: PSA graph much taller than the other two, which sit flush below it
Fig1a <- PSA_plot / P_plot / R2_plot +
  plot_layout(heights = c(4, 1, 1))
Fig1a

ggsave("paper_graphs/Fig1a.pdf", Fig1a, width = 6, height = 8)