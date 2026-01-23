##################################################################################
## Organize data for phenosight 1 experiment
## 
## Author Daniel Anstett
## 
## 
## Written Jan 22 2026
###################################################################################
#Import libraries
library(tidyverse)
library(car)

GR_early <- read_csv("data/Experiment01_EarlyGrowth_Mimulus.csv")
GR_late <- read_csv("data/Experiment01_LateGrowth_Mimulus.csv")

PS2_early <- read_csv("data/Experiment01_EarlyPS2_Mimulus.csv") %>% filter(TOE>100)
PS2_late <- read_csv("data/Experiment01_LatePS2_Mimulus.csv")
decode <- read_csv("data/decode.csv")

#Merge early and late
GR_early2 <- merge(GR_early, decode, by = "Population")
GR_early2$Treatment <- factor(GR_early2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

GR_late2 <- merge(GR_late, decode, by = "Population")
GR_late2$Treatment <- factor(GR_late2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

PS2_early2 <- merge(PS2_early, decode, by = "Population")
PS2_early2$Treatment <- factor(PS2_early2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

PS2_late2 <- merge(PS2_late, decode, by = "Population")
PS2_late2$Treatment <- factor(PS2_late2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

#Make a data frame that has line(pop) means
GR_early3 <- GR_early2 %>%
  group_by(Population, Treatment) %>%
  summarise(GR = mean(GR, na.rm = TRUE),Site = first(Site), Year = first(Year), .groups = 'drop')

GR_late3 <- GR_late2 %>%
  group_by(Population, Treatment) %>%
  summarise(GR = mean(GR, na.rm = TRUE),Site = first(Site), Year = first(Year), .groups = 'drop')

PS2_early3 <- PS2_early2 %>%
  group_by(Population, Treatment, TOE) %>%
  summarise(across(Fv.Fm:NDVI, ~mean(.x, na.rm = TRUE)), Site = first(Site), Year = first(Year),.groups = 'drop')

PS2_late3 <- PS2_late2 %>%
  group_by(Population, Treatment, TOE) %>%
  summarise(across(Fv.Fm:NDVI, ~mean(.x, na.rm = TRUE)), Site = first(Site), Year = first(Year),.groups = 'drop')



###################################################################################
###################################################################################
#Run Models
gr_early <- lm(GR ~ Site*Year+Treatment,data=GR_early3)
summary(gr_early)
Anova(gr_early,type=3)

gr_late <- lm(GR ~ Site*Year+Treatment,data=GR_late3)
summary(gr_late)
Anova(gr_late,type=3)

Fv.Fm_early <- lm(Fv.Fm ~ Site*Year+Treatment+TOE,data=PS2_early3) 
summary(Fv.Fm_early)
Anova(Fv.Fm_early)

Fv.Fm_late <- lm(Fv.Fm ~ Site*Year+Treatment+TOE,data=PS2_late3) 
summary(Fv.Fm_late)
Anova(Fv.Fm_late)




########## Make Plot ############

#Early Growth Rate All Data Points
ggplot(GR_early2, aes(x = Treatment, y = GR, fill = factor(Year))) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  facet_wrap(~ Site) +
  labs(x = "Water Holding Capacity", y = "Early Growth Rate") +
  theme_classic() + theme(
  axis.text.x = element_text(size=20, face="bold", hjust = 0.4),
  axis.text.y = element_text(size=20,face="bold"),
  axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 16,face="bold"),  # Increase the size of the legend text
  legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
  legend.key.height = unit(1.6, "lines"), #Reduce height
  strip.text = element_text(size = 14, face = "bold") #Make strip text larger
)+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))

ggsave("graphs/Site_year_gr_early_all.pdf",width=12, height = 6, units = "in")


#Early Growth Rate Line Means
ggplot(GR_early3, aes(x = Treatment, y = GR, fill = factor(Year))) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  scale_y_continuous(breaks = c(0,5,10,15,20))+
  facet_wrap(~ Site) +
  labs(x = "Water Holding Capacity", y = "Early Growth Rate") +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold", hjust = 0.4),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16,face="bold"),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"), #Reduce height
    strip.text = element_text(size = 14, face = "bold") #Make strip text larger
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))

ggsave("graphs/Site_year_gr_early_means.pdf",width=12, height = 6, units = "in")

######## Late ##########


#late Growth Rate All Data Points
ggplot(GR_late2, aes(x = Treatment, y = GR, fill = factor(Year))) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  facet_wrap(~ Site) +
  labs(x = "Water Holding Capacity", y = "Late Growth Rate") +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold", hjust = 0.4),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16,face="bold"),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"), #Reduce height
    strip.text = element_text(size = 14, face = "bold") #Make strip text larger
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))

ggsave("graphs/Site_year_gr_late_all.pdf",width=12, height = 6, units = "in")


#late Growth Rate Line Means
ggplot(GR_late3, aes(x = Treatment, y = GR, fill = factor(Year))) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  scale_y_continuous(breaks = c(0,5,10,15,20))+
  facet_wrap(~ Site) +
  labs(x = "Water Holding Capacity", y = "Late Growth Rate") +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold", hjust = 0.4),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16,face="bold"),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"), #Reduce height
    strip.text = element_text(size = 14, face = "bold") #Make strip text larger
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))

ggsave("graphs/Site_year_gr_late_means.pdf",width=12, height = 6, units = "in")
###################################################################################
#PS2 plots

#Fv/Fm
ggplot(PS2_early3, aes(x = TOE, y = Fv.Fm, color = factor(Year))) +
  geom_line(size = 1) +
  scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                     name = "Year") +
  facet_grid(Treatment ~ Site) +
  labs(x = "Time (h)", y = "Fv/Fm") +
  theme_bw() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size=12,hjust = 0.4),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
        axis.title.y = element_text(color="black", size=20,vjust = 1.7, face="bold",hjust=0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 14,face="bold"))  # Increase the size of the legend text)

ggsave("graphs/FvFm_early.pdf",width=6, height = 8, units = "in")

ggplot(PS2_late3, aes(x = TOE, y = Fv.Fm, color = factor(Year))) +
  geom_line(size = 1) +
  scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                     name = "Year") +
  facet_grid(Treatment ~ Site) +
  labs(x = "Time (h)", y = "Fv/Fm") +
  theme_bw() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size=12,hjust = 0.4),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
        axis.title.y = element_text(color="black", size=20,vjust = 1.7, face="bold",hjust=0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 14,face="bold"))  # Increase the size of the legend text)

ggsave("graphs/FvFm_late.pdf",width=8, height = 11, units = "in")



###################################################################################

