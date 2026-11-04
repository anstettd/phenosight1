##################################################################################
## Setup data frames for figure production & analysis
## Bring in year and site information into data frame
## 
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(ggh4x)


############################################ Data Setup ####################################################
############################################################################################################
#Import Data

#Multispectral Data
PS2_early <- read_csv("data/Experiment01_EarlyPS2_Mimulus.csv") %>% filter(TOE>140) %>% rename(PSA = Roi.Mask)
PS2_late <- read_csv("data/Experiment01_LatePS2_Mimulus.csv") %>% rename(PSA = Roi.Mask)

#Table that gives the pop/year info (SCalifornia relabeled to California)
decode <- read_csv("data/decode.csv") %>% mutate(Site = dplyr::recode(Site, SCalifornia = "California"))

#Add region and year to dataframe
PS2_early2 <- merge(PS2_early, decode, by = "Population")
PS2_early2$Treatment <- factor(PS2_early2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

PS2_late2 <- merge(PS2_late, decode, by = "Population")
PS2_late2$Treatment <- factor(PS2_late2$Treatment, levels = c("20%", "40%", "60%", "80%", "100%"))

#Calculate line means for Population/Treatment/TOE
PS2_early3 <- PS2_early2 %>%
  group_by(Population, Treatment, TOE) %>%
  summarise(across(PSA:NDVI, ~mean(.x, na.rm = TRUE)), Site = first(Site), Year = first(Year),.groups = 'drop')

PS2_late3 <- PS2_late2 %>%
  group_by(Population, Treatment, TOE) %>%
  summarise(across(PSA:NDVI, ~mean(.x, na.rm = TRUE)), Site = first(Site), Year = first(Year),.groups = 'drop')

#Merge early and late into a single dataset
PS2_early3 <- PS2_early3 %>% mutate(Timepoint = "Early")
PS2_late3 <- PS2_late3 %>% mutate(Timepoint = "Late")

PS2_all <- bind_rows(PS2_early3, PS2_late3) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Early", "Late"))) 

write_csv(PS2_all, "data/PS2_all.csv")

############################################################################################################
#Growth Rate Setup

#Magda this is what I currently have for growth setup. However it might be better to put in the calculation
#for growth rate directly from Roi.Mask in here, so its clear and either one of us can easily modify what
#the TOE cut of is. 

#So we should just calculate GR from PS2_early3$Roi.Mark and PS2_late3$Roi.Mark, then export the 
# GR_late3 files as 
#write_csv(GR_early3, "data/GR_early3.csv")
#write_csv(GR_late3, "data/GR_late3.csv")


#Growth Rate Data
#GR_early <- read_csv("data/Experiment01_EarlyGrowth_Mimulus.csv")
#GR_late <- read_csv("data/Experiment01_LateGrowth_Mimulus.csv")
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


