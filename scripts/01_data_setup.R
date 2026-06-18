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
library(lmerTest)
library(lmtest)
library(ggh4x)

###################################################################################
#Functions
PS2 <- function(name, data) {
  
  # Get unique TOE values from the data
  time_PS2 <- unique(data$TOE)
  
  # Initialize results data frame
  PS2_results <- data.frame()
  
  # Loop through each TOE value
  for(i in 1:length(time_PS2)){
    toe_select <- time_PS2[i] 
    PS2_select <- data %>% filter(TOE == toe_select)
    PS2_early <- lm(name ~ Site*Year+Treatment, data = PS2_select)
    
    PS2_results[i,1] <- toe_select
    PS2_results[i,2] <- Anova(PS2_early, type=3)[["Pr(>F)"]][2]
    PS2_results[i,3] <- Anova(PS2_early, type=3)[["Pr(>F)"]][3]
    PS2_results[i,4] <- Anova(PS2_early, type=3)[["Pr(>F)"]][4]
    PS2_results[i,5] <- Anova(PS2_early, type=3)[["Pr(>F)"]][5]
    PS2_results[i,6] <- summary(PS2_early)$r.squared 
  }
  
  # Add column names
  colnames(PS2_results) <- c("TOE", "Site", "Year", "Treatment", "Site*Year", "R2")
  
  # Reshape to long format and add Response column
  PS2_long <- PS2_results %>%
    pivot_longer(cols = -TOE, 
                 names_to = "Variable", 
                 values_to = "Value") %>%
    mutate(Response = name)
  
  return(PS2_long)
}

# Example usage:
# PS2_long <- PS2(name = "Fv.Fm", data = PS2_early3)


###################################################################################
#Import Data
#Growth Rate Data
GR_early <- read_csv("data/Experiment01_EarlyGrowth_Mimulus.csv")
GR_late <- read_csv("data/Experiment01_LateGrowth_Mimulus.csv")
#Multispectral Data
PS2_early <- read_csv("data/Experiment01_EarlyPS2_Mimulus.csv") %>% filter(TOE>100)
PS2_late <- read_csv("data/Experiment01_LatePS2_Mimulus.csv")

#Table that gives the pop/year info
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
#Run Models for Growth Rate for line means
gr_early <- lm(GR ~ Site*Year+Treatment,data=GR_early3)
#summary(gr_early)
Anova(gr_early,type=3)

gr_late <- lm(GR ~ Site*Year+Treatment,data=GR_late3)
#summary(gr_late)
Anova(gr_late,type=3)

#################################### Make Plot Growth Rates #####################################

########################################### Line Means ##########################################


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
    strip.text = element_text(size = 14, face = "bold") #Make strip text larger
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))

ggsave("graphs/Site_year_gr_early_means.pdf",width=12, height = 6, units = "in")


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

ggsave("graphs/Site_year_gr_late_means.pdf",width=12, height = 6, units = "in")


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

ggsave("graphs/SCal_site_year_gr_early_means_cali.pdf", width=6, height=4, units="in")

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

ggsave("graphs/Site_year_gr_early_all.pdf",width=12, height = 6, units = "in")

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


###################################################################################

#Single time point
toe_select<-unique(PS2_early3$TOE)[5] 
PS2_select <- PS2_early3 %>% filter(TOE==toe_select) #Select on TOE point
PS2_early <- lm(Fv.Fm ~ Site*Year+Treatment,data=PS2_select)
summary(PS2_early)
Anova(PS2_early,type=3)

summary(PS2_early)$r.squared 
Anova(PS2_early,type=3)[["Pr(>F)"]][2]
Anova(PS2_early,type=3)[["Pr(>F)"]][3]
Anova(PS2_early,type=3)[["Pr(>F)"]][4]
Anova(PS2_early,type=3)[["Pr(>F)"]][5]


#Loop
time_PS2<-unique(PS2_early3$TOE)
PS2_results <- data.frame()

for(i in 1:length(time_PS2)){
  toe_select<-time_PS2[i] 
  PS2_select <- PS2_early3 %>% filter(TOE==toe_select) #Select on TOE point
  PS2_early <- lm(Fv.Fm ~ Site*Year+Treatment,data=PS2_select)
  PS2_results[i,1]<-unique(PS2_early3$TOE)[i] 
  PS2_results[i,2]<-Anova(PS2_early,type=3)[["Pr(>F)"]][2]
  PS2_results[i,3]<-Anova(PS2_early,type=3)[["Pr(>F)"]][3]
  PS2_results[i,4]<-Anova(PS2_early,type=3)[["Pr(>F)"]][4]
  PS2_results[i,5]<-Anova(PS2_early,type=3)[["Pr(>F)"]][5]
  PS2_results[i,6]<-summary(PS2_early)$r.squared 
}
PS2_results
colnames(PS2_results)<-c("TOE","Site","Year","Treatment","Site*Year","R2")

PS2_long <- PS2_results %>%
  pivot_longer(cols = -TOE, 
               names_to = "Variable", 
               values_to = "Value") %>%
  mutate(Response = "Fv.Fm")
#################################################################################################
#################################################################################################



###################################### Multi-Spectral Data  ##########################################
### Single Variable sample analysis ###
###################################### Fv.Fm Early  ##########################################

  # Get unique TOE values from the data
  time_PS2 <- unique(PS2_early3$TOE)
  
  # Initialize results data frame
  PS2_Fv.Fm_early <- data.frame()
  
  # Loop through each TOE value
  for(i in 1:length(time_PS2)){
    toe_select <- time_PS2[i] 
    PS2_select <- PS2_early3 %>% dplyr::filter(TOE == toe_select)
    PS2_early <- lm(Fv.Fm ~ Site*Year+Treatment, data = PS2_select)
    
    PS2_Fv.Fm_early[i,1] <- toe_select
    PS2_Fv.Fm_early[i,2] <- Anova(PS2_early, type=3)[["Pr(>F)"]][2]
    PS2_Fv.Fm_early[i,3] <- Anova(PS2_early, type=3)[["Pr(>F)"]][3]
    PS2_Fv.Fm_early[i,4] <- Anova(PS2_early, type=3)[["Pr(>F)"]][4]
    PS2_Fv.Fm_early[i,5] <- Anova(PS2_early, type=3)[["Pr(>F)"]][5]
    PS2_Fv.Fm_early[i,6] <- summary(PS2_early)$r.squared 
  }
  
  # Add column names
  colnames(PS2_Fv.Fm_early) <- c("TOE", "Site", "Year", "Treatment", "Site*Year", "R2")
  
  
######################### p-value Graphs #########################   


  PS2_Fv.Fm_early %>%
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
    labs(x = "Hours of Imaging", y = "") +
    theme_bw() + theme(
      axis.text.x = element_text(size = 16, face = "bold", hjust = 0.4),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(color = "black", size = 20, vjust = 0.5, face = "bold"),
      axis.title.y = element_text(color = "black", size = 20, vjust = 1.7, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 14, face = "bold")
    )
  
  ggsave("graphs/multi-spectral/fv.fm_early_p_val.pdf",width=6, height = 8, units = "in")
  

  
######################### TOE Graphs ######################### 
  
#All Treatments
  #head(PS2_early3)
  #unique(PS2_early3$Site)
  #unique(PS2_early3$Year)
  #unique(PS2_early3$Treatment)  

PS2_early3 %>%
  mutate(Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
  ggplot(aes(x = TOE, y = Fv.Fm, color = factor(Year))) +
  geom_point() +
  scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"), name = "Year") +
  geom_smooth() +
  facet_grid(Treatment ~ Site, 
             labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California"))) +
  labs(x = "Time of Experiment", y = "Fv/Fm") +
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
ggsave("graphs/multi-spectral/fv.fm_early_all.pdf",width=10, height = 16, units = "in")

### Just 20% Treatment
PS2_early3_ <- PS2_early3 %>% filter(Treatment=="20%") %>% select(Population,TOE,Treatment,Fv.Fm,Site,Year)
#head(PS2_early3_20)
#unique(PS2_early3_20$Site)
#unique(PS2_early3_20$Year)

PS2_early3_20 %>%
  mutate(Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
  ggplot(aes(x = TOE, y = Fv.Fm, color = factor(Year))) +
  geom_point() +
  scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"), name = "Year") +
  geom_smooth() +
facet_wrap(~ Site, labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California")))+
  labs(x = "Time of Experiment", y = "Fv/Fm") +
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
ggsave("graphs/multi-spectral/fv.fm_early_20.pdf",width=8, height = 6, units = "in")



###################################### Fv.Fm late  ##########################################

# Get unique TOE values from the data
time_PS2 <- unique(PS2_late3$TOE)

# Initialize results data frame
PS2_Fv.Fm_late <- data.frame()

# Loop through each TOE value
for(i in 1:length(time_PS2)){
  toe_select <- time_PS2[i] 
  PS2_select <- PS2_late3 %>% dplyr::filter(TOE == toe_select)
  PS2_late <- lm(Fv.Fm ~ Site*Year+Treatment, data = PS2_select)
  
  PS2_Fv.Fm_late[i,1] <- toe_select
  PS2_Fv.Fm_late[i,2] <- Anova(PS2_late, type=3)[["Pr(>F)"]][2]
  PS2_Fv.Fm_late[i,3] <- Anova(PS2_late, type=3)[["Pr(>F)"]][3]
  PS2_Fv.Fm_late[i,4] <- Anova(PS2_late, type=3)[["Pr(>F)"]][4]
  PS2_Fv.Fm_late[i,5] <- Anova(PS2_late, type=3)[["Pr(>F)"]][5]
  PS2_Fv.Fm_late[i,6] <- summary(PS2_late)$r.squared 
}

# Add column names
colnames(PS2_Fv.Fm_late) <- c("TOE", "Site", "Year", "Treatment", "Site*Year", "R2")


######################### p-value Graphs #########################   


PS2_Fv.Fm_late %>%
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
  labs(x = "Hours of Imaging", y = "") +
  theme_bw() + theme(
    axis.text.x = element_text(size = 16, face = "bold", hjust = 0.4),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(color = "black", size = 20, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, vjust = 1.7, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave("graphs/multi-spectral/fv.fm_late_p_val.pdf",width=6, height = 8, units = "in")



######################### TOE Graphs ######################### 

#All Treatments
#head(PS2_late3)
#unique(PS2_late3$Site)
#unique(PS2_late3$Year)
#unique(PS2_late3$Treatment)  

PS2_late3 %>%
  mutate(Year = factor(Year, levels = c(2010, 2014)),
         Site = factor(Site, levels = c("Oregon", "SCalifornia"))) %>%
  ggplot(aes(x = TOE, y = Fv.Fm, color = factor(Year))) +
  geom_point() +
  scale_color_manual(values = c("2010" = "skyblue3", "2014" = "#FF7700"), name = "Year") +
  geom_smooth() +
  facet_grid(Treatment ~ Site, 
             labeller = labeller(Site = c("Oregon" = "Oregon", "SCalifornia" = "Southern California"))) +
  labs(x = "Time of Experiment", y = "Fv/Fm") +
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
ggsave("graphs/multi-spectral/fv.fm_late_all.pdf",width=10, height = 16, units = "in")



























#################################################################################################
# Run Function PS2
#Tabulate P-values for Site, Year, Treatment and Site*Year and R2

#Run Models for PS2 Early
early01<-PS2("Fv.Fm",PS2_early3)
early02<-PS2("Fq.Fm",PS2_early3)
early03<-PS2("NPQ",PS2_early3)
early04<-PS2("qP",PS2_early3)
early05<-PS2("qN",PS2_early3)
early06<-PS2("qL",PS2_early3)
early07<-PS2("qI",PS2_early3)
early08<-PS2("phiNO",PS2_early3)
early09<-PS2("phiNPQ",PS2_early3)
early10<-PS2("npq.t",PS2_early3)
early11<-PS2("ChlIdx",PS2_early3)
early12<-PS2("AriIdx",PS2_early3)
early13<-PS2("NDVI",PS2_early3)

#Run Models for PS2 Late

late01<-PS2("Fv.Fm",PS2_late3)
late02<-PS2("Fq.Fm",PS2_late3)
late03<-PS2("NPQ",PS2_late3)
late04<-PS2("qP",PS2_late3)
late05<-PS2("qN",PS2_late3)
late06<-PS2("qL",PS2_late3)
late07<-PS2("qI",PS2_late3)
late08<-PS2("phiNO",PS2_late3)
late09<-PS2("phiNPQ",PS2_late3)
late10<-PS2("npq.t",PS2_late3)
late11<-PS2("ChlIdx",PS2_late3)
late12<-PS2("AriIdx",PS2_late3)
late13<-PS2("NDVI",PS2_late3)

###################################################################################
###################################################################################
#Make Graphs
#Early


ggplot(PS2_long, aes(x = TOE, y = Value)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
  labs(x = "TOE", 
       y = "P-Value") +
  theme_bw() +
  theme(strip.text = element_text(size = 11, face = "bold"))




































###################################################################################
#Defunct
# Run PS2 mixed effect models
#Fv.Fm

Fv.Fm_3way <- lmer(Fv.Fm ~ Site*Year*Treatment+TOE+(1|Population),data=PS2_early3) 
Fv.Fm_early <- lmer(Fv.Fm ~ Site*Year+Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fv.Fm_3way,Fv.Fm_early) #drop interaciton

Fv.Fm_early <- lmer(Fv.Fm ~ Site*Year+Treatment+TOE+(1|Population),data=PS2_early3) 
Fv.Fm_early_no <- lmer(Fv.Fm ~ Site+Year+Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fv.Fm_early, Fv.Fm_early_no) #drop interaciton

Fv.Fm_early_no_year <- lmer(Fv.Fm ~ Site+Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fv.Fm_early_no,Fv.Fm_early_no_year) #drop year

Fv.Fm_early_Treatment <- lmer(Fv.Fm ~ Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fv.Fm_early_no_year,Fv.Fm_early_Treatment) #drop site

Fv.Fm_early_Time <- lmer(Fv.Fm ~ TOE+(1|Population),data=PS2_early3) 
lrtest(Fv.Fm_early_Treatment,Fv.Fm_early_Time) #Keep Treatment only model.

# Run PS2 mixed effect models
Fv.Fm_late <- lmer(Fv.Fm ~ Site*Year+Treatment+TOE+(1|Population),data=PS2_late3) 
Fv.Fm_late_no <- lmer(Fv.Fm ~ Site+Year+Treatment+TOE+(1|Population),data=PS2_late3) 
lrtest(Fv.Fm_late, Fv.Fm_late_no) #drop interaciton

Fv.Fm_late_no_year <- lmer(Fv.Fm ~ Site+Treatment+TOE+(1|Population),data=PS2_late3) 
lrtest(Fv.Fm_late_no,Fv.Fm_late_no_year) #drop year

Fv.Fm_late_Treatment <- lmer(Fv.Fm ~ Treatment+TOE+(1|Population),data=PS2_late3) 
lrtest(Fv.Fm_late_no_year,Fv.Fm_late_Treatment) #drop site

Fv.Fm_late_Time <- lmer(Fv.Fm ~ TOE+(1|Population),data=PS2_late3) 
lrtest(Fv.Fm_late_Treatment,Fv.Fm_late_Time) #Keep Treatment only model.





############
#Fq.Fm
Fq.Fm_early <- lmer(Fq.Fm ~ Site*Year+Treatment+TOE+(1|Population),data=PS2_early3) 
Fq.Fm_early_no <- lmer(Fq.Fm ~ Site+Year+Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fq.Fm_early, Fq.Fm_early_no) #drop interaciton

Fq.Fm_early_no_year <- lmer(Fq.Fm ~ Site+Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fq.Fm_early_no,Fq.Fm_early_no_year) #drop year

Fq.Fm_early_Treatment <- lmer(Fq.Fm ~ Treatment+TOE+(1|Population),data=PS2_early3) 
lrtest(Fq.Fm_early_no_year,Fq.Fm_early_Treatment) #drop site

Fq.Fm_early_Time <- lmer(Fq.Fm ~ TOE+(1|Population),data=PS2_early3) 
lrtest(Fq.Fm_early_Treatment,Fq.Fm_early_Time) #Keep Treatment only model.






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

