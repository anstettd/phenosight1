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
library(lme4)
library(lmtest)
#library(ggh4x)
###################################################################################
#Import Data
#Growth Rate Data
GR_early <- read_csv("data/Experiment01_EarlyGrowth_Mimulus.csv")
GR_late <- read_csv("data/Experiment01_LateGrowth_Mimulus.csv")
#Table that gives the pop/year info
decode <- read_csv("data/decode.csv")
#Merge early and late
GR_early2 <- merge(GR_early, decode, by = "Population")
GR_early2$Treatment <- factor(GR_early2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))
GR_late2 <- merge(GR_late, decode, by = "Population")
GR_late2$Treatment <- factor(GR_late2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))
#Make a data frame that has line(pop) means
GR_early3 <- GR_early2 %>%
  group_by(Population, Treatment) %>%
  summarise(GR = mean(GR, na.rm = TRUE),Site = first(Site), Year = first(Year), .groups = 'drop')
GR_late3 <- GR_late2 %>%
  group_by(Population, Treatment) %>%
  summarise(GR = mean(GR, na.rm = TRUE),Site = first(Site), Year = first(Year), .groups = 'drop')
###################################################################################
###################################################################################
#Run Models for Growth Rate for line means
gr_early <- lm(GR ~ Site*Year+Treatment,data=GR_early3)
#summary(gr_early)
Anova(gr_early,type=3)
gr_late <- lm(GR ~ Site*Year+Treatment,data=GR_late3)
#summary(gr_late)
Anova(gr_late,type=3)



#################################### Make Plot Growth Rates #####################################
########################################### For NSF Grant ##########################################
### Early Growth Rate Line Means ###
ggplot(GR_early3, aes(x = Treatment, y = GR, fill = factor(Year))) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"),
                    name = "Year") +
  scale_y_continuous(breaks = c(0,5,10,15,20))+
  facet_wrap(~ Site) +
  labs(x = "Water Holding Capacity", y = "Early Growth Rate") +
  theme_classic() + theme(
    axis.text.x = element_text(size=22, face="bold", hjust = 0.4),
    axis.text.y = element_text(size=22,face="bold"),
    axis.title.x = element_text(color="black", size=30, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=30,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20,face="bold"),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"), #Reduce height
    strip.text = element_text(size = 22, face = "bold"), #Make strip text larger
    # Move the legend into the whitespace in the top-left corner of the left
    # facet panel (Oregon), instead of it sitting off to the side.
    # >>> ADJUST THESE TWO NUMBERS (x, y) TO REPOSITION THE LEGEND <<<
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.99),   # <-- tweak x/y (0-1 range) here
    legend.justification = c(0, 1),           # anchors legend's top-left corner
    legend.background = element_rect(fill = "white", color = NA)
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
ggsave("graphs/growth/Site_year_gr_early_means.pdf",width=12, height = 6, units = "in")


#################################### Make Plot Growth Rates #####################################
########################################### For paper ##########################################
### Early Growth Rate Line Means ###
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
    strip.text = element_text(size = 14, face = "bold"), #Make strip text larger
    # Move the legend into the whitespace in the top-left corner of the left
    # facet panel (Oregon), instead of it sitting off to the side.
    # >>> ADJUST THESE TWO NUMBERS (x, y) TO REPOSITION THE LEGEND <<<
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.99),   # <-- tweak x/y (0-1 range) here
    legend.justification = c(0, 1),           # anchors legend's top-left corner
    legend.background = element_rect(fill = "white", color = NA)
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
ggsave("graphs/growth/Site_year_gr_early_means.pdf",width=12, height = 6, units = "in")





#####################################################################################################
#####################################################################################################
### late Growth Rate Line Means ###
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

ggsave("graphs/growth/Site_year_gr_late_means.pdf",width=12, height = 6, units = "in")


##### Merge early and late into one 2x2 graphic.


################ For Grants & Presentations ##########################

#Only S California
GR_early3_scal <- GR_early3 %>%
  filter(Site %in% c("SCalifornia"))

ggplot(GR_early3_scal, aes(x = Treatment, y = GR, fill = factor(Year))) +
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
    legend.text = element_text(size = 16,face="bold"),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"),
    strip.text = element_blank(),        # removes header text
    strip.background = element_blank()   # removes header box
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))

ggsave("graphs/growth/SCal_site_year_gr_early_means_cali.pdf", width=6, height=4, units="in")














###############################################################################################
###############################################################################################


######################################## All Data Points #######################################

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

#ggsave("graphs/growth/Site_year_gr_early_all.pdf",width=12, height = 6, units = "in")

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

#ggsave("graphs/growth/Site_year_gr_late_all.pdf",width=12, height = 6, units = "in")

