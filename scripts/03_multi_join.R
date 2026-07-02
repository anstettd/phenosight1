##################################################################################
## Plot out all mutispectral data for phenosight1 
## 
## Author Daniel Anstett
## 
## 
## Written July 2026
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(lmtest)
library(ggh4x)

############################################################################################################
############################################ Data Setup ####################################################
############################################################################################################
#Import Data

#Multispectral Data
PS2_early <- read_csv("data/Experiment01_EarlyPS2_Mimulus.csv") %>% filter(TOE>100)
PS2_late <- read_csv("data/Experiment01_LatePS2_Mimulus.csv")

#Table that gives the pop/year info
decode <- read_csv("data/decode.csv")

#Merge early and late
PS2_early2 <- merge(PS2_early, decode, by = "Population")
PS2_early2$Treatment <- factor(PS2_early2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

PS2_late2 <- merge(PS2_late, decode, by = "Population")
PS2_late2$Treatment <- factor(PS2_late2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

#Make a data frame that has line(pop) means
PS2_early3 <- PS2_early2 %>%
  group_by(Population, Treatment, TOE) %>%
  summarise(across(Fv.Fm:NDVI, ~mean(.x, na.rm = TRUE)), Site = first(Site), Year = first(Year),.groups = 'drop')

PS2_late3 <- PS2_late2 %>%
  group_by(Population, Treatment, TOE) %>%
  summarise(across(Fv.Fm:NDVI, ~mean(.x, na.rm = TRUE)), Site = first(Site), Year = first(Year),.groups = 'drop')

#Merge early and late into a single dataset
PS2_early3 <- PS2_early3 %>% mutate(Timepoint = "Early")
PS2_late3 <- PS2_late3 %>% mutate(Timepoint = "Late")

PS2_all <- bind_rows(PS2_early3, PS2_late3) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Early", "Late")))

write_csv(PS2_all, "data/PS2_all.csv")



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
           Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
    ggplot(aes(x = TOE, y = .data[[var]], color = factor(Year))) +
    geom_point() +
    scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"), name = "Year") +
    geom_smooth() +
    facet_grid(Treatment ~ Site,
               labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California"))) +
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
  
  # --- 20% treatment plot ---
  p3 <- PS2_all %>%
    dplyr::filter(Treatment == "20%") %>%
    mutate(Year = factor(Year, levels = c(2010, 2014)),
           Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
    ggplot(aes(x = TOE, y = .data[[var]], color = factor(Year))) +
    geom_point() +
    scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"), name = "Year") +
    geom_smooth() +
    facet_wrap(~ Site, labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California"))) +
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
  ggsave(paste0("graphs/multi_join/20_treatment/", tolower(var), "_20.pdf"),
         plot = p3, width = 8, height = 6, units = "in")
  
  message("Done: ", var)
}

###############################################################################
## Make Joint PDF
###############################################################################
library(pdftools)

# --- 20% treatment booklet ---
files_20 <- list.files("graphs/multi_join/20_treatment/",
                       pattern = "\\.pdf$", full.names = TRUE)
pdf_combine(files_20, output = "graphs/multi_join/booklets/20_treatment_booklet.pdf")

# --- all treatments booklet ---
files_all <- list.files("graphs/multi_join/all_treatments/",
                        pattern = "\\.pdf$", full.names = TRUE)
pdf_combine(files_all, output = "graphs/multi_join/booklets/all_treatments_booklet.pdf")

# --- p-val booklet ---
files_pval <- list.files("graphs/multi_join/p-val/",
                         pattern = "\\.pdf$", full.names = TRUE)
pdf_combine(files_pval, output = "graphs/multi_join/booklets/p_val_booklet.pdf")


