##################################################################################
## Predict which variables explain 2010 vs 2014
## For both sites
## Author Daniel Anstett
## 
## 
## Written July 2026
###################################################################################
#Import libraries
library(tidyverse)

#Import Data
PS2_all <- read_csv("data/PS2_all.csv")
PS2_all <- PS2_all %>% mutate(Treatment = factor(Treatment, levels = c("20%", "40%", "60%", "80%", "100%")),
    Timepoint = factor(Timepoint, levels = c("Early", "Late")))

###############################################################################
## Correlation matrix for PS2 variables
###############################################################################

corr_vars <- c("Fv.Fm", "Fq.Fm", "NPQ", "qP", "qN", "qL", "qI",
               "phiNO", "phiNPQ", "npq.t", "ChlIdx", "AriIdx", "NDVI")

r_corr <- PS2_all %>%
  dplyr::select(all_of(corr_vars)) %>%
  cor(use = "pairwise.complete.obs")

print(r_corr)


###############################################################################
## Correlation heatmap
###############################################################################
library(reshape2)

#Reshape r_corr (matrix) into long format for ggplot
r_corr_long <- melt(r_corr, varnames = c("Var1", "Var2"), value.name = "Correlation")

#Mask diagonal (self-correlations = 1) so they don't get colored
r_corr_long <- r_corr_long %>%
  mutate(Correlation_plot = ifelse(Var1 == Var2, NA, Correlation))

#Set factor levels: X axis left-to-right in corr_vars order, Y axis top-to-bottom in corr_vars order
r_corr_long <- r_corr_long %>%
  mutate(Var1 = factor(Var1, levels = corr_vars),
         Var2 = factor(Var2, levels = rev(corr_vars)))

p_heatmap <- ggplot(r_corr_long, aes(x = Var1, y = Var2, fill = Correlation_plot)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 3.5) +
  scale_fill_gradient2(low = "#378ADD", mid = "white", high = "#D85A30",
                       midpoint = 0, limit = c(-1, 1), name = "r",
                       na.value = "grey90") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Correlation matrix: PS2 variables") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  ) +
  coord_fixed()

p_heatmap

ggsave("graphs/multi_join/booklets/correlation_heatmap.pdf", plot = p_heatmap, width = 9, height = 8, units = "in")


###################### Predict year based on trait data #########################

#Make sure Year is a proper binary factor (2010 vs 2014)
PS2_all <- PS2_all %>%
  mutate(Year = factor(Year, levels = c(2010, 2014)))

#Break up Sites
PS2_cali <- PS2_all %>% filter(Site == "SCalifornia")
PS2_oregon <- PS2_all %>% filter(Site == "Oregon")

glm_cali <- glm(Year ~ Fv.Fm + Fq.Fm + NPQ + npq.t + ChlIdx,family=binomial,data=PS2_cali)
summary(glm_cali)

glm_oregon <- glm(Year ~ Fv.Fm + Fq.Fm + NPQ + npq.t + ChlIdx,family=binomial,data=PS2_oregon)
summary(glm_oregon)


###############################################################################
## Random forest variable importance — MeanDecreaseAccuracy, separate by Site
###############################################################################

#Extract importance values into tidy data frames
imp_cali <- as.data.frame(importance(rf_cali)) %>%
  mutate(Variable = rownames(.), Site = "California")

imp_oregon <- as.data.frame(importance(rf_oregon)) %>%
  mutate(Variable = rownames(.), Site = "Oregon")

#--- California plot ---
p_importance_cali <- imp_cali %>%
  ggplot(aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_col(fill = "#378ADD") +
  coord_flip() +
  labs(x = NULL, y = "Mean decrease in accuracy",
       title = "Random forest variable importance: California") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

p_importance_cali

ggsave("graphs/multi_join/booklets/rf_importance_california.pdf", plot = p_importance_cali, width = 7, height = 5, units = "in")

#--- Oregon plot ---
p_importance_oregon <- imp_oregon %>%
  ggplot(aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_col(fill = "#378ADD") +
  coord_flip() +
  labs(x = NULL, y = "Mean decrease in accuracy",
       title = "Random forest variable importance: Oregon") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

p_importance_oregon

ggsave("graphs/multi_join/booklets/rf_importance_oregon.pdf", plot = p_importance_oregon, width = 7, height = 5, units = "in")



###############################################################################
## XY Graphs
###############################################################################
#Specific xy plots, cholrophy, anthosyenin, on FV.Fm
lm1 <-lm(Fv.Fm~ChlIdx,data=PS2_all)
summary(lm1) #P<0.001

lm1 <-lm(Fv.Fm~AriIdx,data=PS2_all)
summary(lm1) # 0.00692

### Graphs 

#--- Fv.Fm ~ ChlIdx ---
p_chl <- PS2_all %>%
  ggplot(aes(x = ChlIdx, y = Fv.Fm)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Chlorophyll Index", y = "Fv/Fm") +
  theme_classic() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 24, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 24, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_blank(),
    strip.background = element_blank()
  )

p_chl

ggsave("graphs/multi_join/xy/fvfm_chlidx.pdf", plot = p_chl, width = 8, height = 6, units = "in")

#--- Fv.Fm ~ AriIdx ---
p_ari <- PS2_all %>%
  ggplot(aes(x = AriIdx, y = Fv.Fm)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Anthocyanin Index", y = "Fv/Fm") +
  theme_classic() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 24, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 24, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_blank(),
    strip.background = element_blank()
  )

p_ari

ggsave("graphs/multi_join/xy/fvfm_ariidx.pdf", plot = p_ari, width = 8, height = 6, units = "in")

###############################################################################
## XY Graphs faceted by Site x Year
###############################################################################
#Specific xy plots, chlorophyll, anthocyanin, on Fv.Fm
lm1 <- lm(Fv.Fm ~ ChlIdx+Site*Year, data = PS2_all)
Anova(lm1,type=3) #P<0.001
lm2 <- lm(Fv.Fm ~ AriIdx+Site*Year, data = PS2_all)
Anova(lm2,type=3) # 0.00692

#--- Fv.Fm ~ ChlIdx, faceted by Site x Year ---
p_chl <- PS2_all %>%
  ggplot(aes(x = ChlIdx, y = Fv.Fm, fill = factor(Year))) +
  geom_point(shape = 21, color = "black", size = 2) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  facet_grid(Year ~ Site) +
  labs(x = "Chlorophyll Index", y = "Fv/Fm") +
  theme_classic() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 24, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 24, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank()
  )
p_chl
ggsave("graphs/multi_join/xy/fvfm_chlidx_siteyear.pdf", plot = p_chl, width = 10, height = 8, units = "in")

#--- Fv.Fm ~ AriIdx, faceted by Site x Year ---
p_ari <- PS2_all %>%
  ggplot(aes(x = AriIdx, y = Fv.Fm, fill = factor(Year))) +
  geom_point(shape = 21, color = "black", size = 2) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  facet_grid(Year ~ Site) +
  labs(x = "Anthocyanin Index", y = "Fv/Fm") +
  theme_classic() + theme(
    axis.text.x = element_text(size = 20, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 24, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 24, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank()
  )
p_ari
ggsave("graphs/multi_join/xy/fvfm_ariidx_siteyear.pdf", plot = p_ari, width = 10, height = 8, units = "in")