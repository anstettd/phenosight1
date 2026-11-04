##################################################################################
## Fig 2 : Fv/FM and Chlorophyl graphs with P-values
## Get select graphs for 20 and 100%
## 
###################################################################################
# Import libraries
library(tidyverse)  
library(car)        
library(ggh4x)      
library(pdftools) 

# Import Data
PS2_all <- read_csv("data/PS2_all.csv")
PS2_all <- PS2_all %>% mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
                              Timepoint = factor(Timepoint, levels = c("Early", "Late")))

# Colors: solid line/point colors by Year, and lighter fill colors for geom_smooth ribbons
year_colors <- c("2010" = "skyblue3", "2014" = "#FF7700")
year_fill_colors <- c("2010" = "#8FC1E3", "2014" = "#FFA858")



############################################################################################################
###################################### Multi-Spectral Data  ################################################
############################################################################################################

response_vars <- c("Fv.Fm", "Fq.Fm", "NPQ", "qP", "qN", "qL", "qI",
                   "phiNO", "phiNPQ", "npq.t", "ChlIdx", "AriIdx", "NDVI")

###############################################################################
## All Multi_spectral across all TOE
###############################################################################
for(var in response_vars) {
  
  time_PS2 <- unique(PS2_all$TOE)
  results <- data.frame()
  
  for(i in 1:length(time_PS2)){
    toe_select <- time_PS2[i]
    PS2_select <- PS2_all %>% dplyr::filter(TOE == toe_select)
    
    formula <- as.formula(paste(var, "~ Site*Year + Treatment"))
    model <- lm(formula, data = PS2_select)
    anova_out <- Anova(model, type = 3)
    
    results[i, 1] <- toe_select
    results[i, 2] <- anova_out[["Pr(>F)"]][2]
    results[i, 3] <- anova_out[["Pr(>F)"]][3]
    results[i, 4] <- anova_out[["Pr(>F)"]][4]
    results[i, 5] <- anova_out[["Pr(>F)"]][5]
    results[i, 6] <- summary(model)$r.squared
  }
  
  colnames(results) <- c("TOE", "Site", "Year", "Treatment", "Site*Year", "R2")
  
  # --- p-value plot ---
  p1 <- results %>%
    pivot_longer(cols = c(Treatment, Site, Year, `Site*Year`, R2),
                 names_to = "Variable",
                 values_to = "Value") %>%
    mutate(Variable = factor(Variable, levels = c("Treatment", "Site*Year", "Site", "Year", "R2")),
           Value = ifelse(Variable %in% c("Treatment", "Site", "Year", "Site*Year"), -log10(Value), Value)) %>%
    ggplot(aes(x = TOE, y = Value)) +
    geom_point() +
    geom_line() +
    geom_hline(data = data.frame(Variable = factor(c("Treatment", "Site", "Year", "Site*Year"),
                                                   levels = c("Treatment", "Site*Year", "Site", "Year", "R2")),
                                 yintercept = -log10(0.05)),
               aes(yintercept = yintercept),
               linetype = "dashed", color = "red") +
    geom_vline(xintercept = 400, color = "black", linetype = "dotted",
               linewidth = 0.3, alpha = 0.3, show.legend = FALSE) +
    facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
    facetted_pos_scales(y = list(
      Variable == "Treatment" ~ scale_y_continuous(name = "-log10(p-value)"),
      Variable == "Site*Year" ~ scale_y_continuous(name = "-log10(p-value)"),
      Variable == "Site"      ~ scale_y_continuous(name = "-log10(p-value)"),
      Variable == "Year"      ~ scale_y_continuous(name = "-log10(p-value)"),
      Variable == "R2"        ~ scale_y_continuous(name = "R²")
    )) +
    labs(x = "Hours of Imaging", y = "", title = var) +
    theme_bw() + theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 16, face = "bold", hjust = 0.4),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(color = "black", size = 20, vjust = 0.5, face = "bold"),
      axis.title.y = element_text(color = "black", size = 20, vjust = 1.7, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 14, face = "bold")
    )
  ggsave(paste0("graphs/multi_join/p-val/", tolower(var), "_p_val.pdf"),
         plot = p1, width = 6, height = 8, units = "in")
  
  # --- all treatments plot ---
  p2 <- PS2_all %>%
    mutate(Year = factor(Year, levels = c(2010, 2014)),
           Site = factor(Site, levels = c("Oregon", "California"))) %>%
    ggplot(aes(x = TOE, y = .data[[var]], color = factor(Year))) +
    geom_point() +
    scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"), name = "Year") +
    geom_smooth() +
    geom_vline(xintercept = 400, color = "black", linetype = "dotted",
               linewidth = 0.3, alpha = 0.3, show.legend = FALSE) +
    facet_grid(Treatment ~ Site,
               labeller = labeller(Site = c("Oregon" = "Oregon", "California" = "Southern California"))) +
    labs(x = "Time of Experiment", y = var) +
    theme_bw() + theme(
      axis.text.x = element_text(size = 16, face = "bold", hjust = 0.4),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(color = "black", size = 20, vjust = 0.5, face = "bold"),
      axis.title.y = element_text(color = "black", size = 20, vjust = 1.7, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 16, face = "bold"),
      legend.key.size = unit(2, "lines"),
      legend.key.height = unit(1.6, "lines"),
      strip.text = element_text(size = 14, face = "bold")
    )
  ggsave(paste0("graphs/multi_join/all_treatments/", tolower(var), "_all.pdf"),
         plot = p2, width = 12, height = 16, units = "in")
  
  message("Done: ", var)
}

###############################################################################
## Make Joint PDF
###############################################################################

# --- all treatments booklet ---
files_all <- list.files("graphs/multi_join/all_treatments/",
                        pattern = "\\.pdf$", full.names = TRUE)
pdf_combine(files_all, output = "graphs/multi_join/booklets/all_treatments_booklet.pdf")

# --- p-val booklet ---
files_pval <- list.files("graphs/multi_join/p-val/",
                         pattern = "\\.pdf$", full.names = TRUE)
pdf_combine(files_pval, output = "graphs/multi_join/booklets/p_val_booklet.pdf")