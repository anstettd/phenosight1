##################################################################################
## Setup data frames for figure production & analysis
## Bring in year and site information into data frame
## 
###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(ggh4x)

#Single place to change the time cut-off used everywhere below (MMJ, Nov 2026).
#144 h is the first grid point kept and is the minimum of projected shoot area in
#all four site x year groups; before it the splines are noisy on very small plants.
TOE_MIN <- 140


############################################ Data Setup ####################################################
############################################################################################################
#Import Data

#Multispectral Data
PS2_early <- read_csv("data/Experiment01_EarlyPS2_Mimulus.csv") %>% filter(TOE>TOE_MIN) %>% rename(PSA = Roi.Mask)
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
#
# Dan: done, and the cut-off is now the TOE_MIN variable at the top of this script (MMJ, Nov 2026).
#
# Two notes on what changed:
#
# 1. GR is fitted PER PLANT and only then averaged to line means, rather than
#    fitting one slope to the line-mean series. Same quantity, but it keeps the
#    within-line variation, and it lets the R2 filter below throw out individual
#    plants whose growth was not linear. Fitting the line means first would hide
#    both. The fit therefore runs off PS2_early2 / PS2_late2 (plant level), not
#    PS2_early3 / PS2_late3 (line means).
#
# 2. The old Experiment01_*Growth_Mimulus.csv files used TOE > 120, and only for
#    Early - the Late frame starts at 411 h, so the cut-off never applied there.
#    Both windows now use TOE_MIN, so the two phases are treated the same way.
#    That drops the 132 h point from Early, so the Early slopes shift slightly
#    from the old files.
#
# The filters (R2 > 0.7, GR > 0) are the same ones used for the original files.

gr_per_plant <- function(df, toe_min = TOE_MIN) {
  df %>%
    filter(TOE > toe_min, is.finite(PSA)) %>%
    group_by(Plant.ID, Population, Treatment, Site, Year) %>%
    filter(n() > 2) %>%
    summarise(GR = coef(lm(PSA ~ TOE))[2],
              R2 = summary(lm(PSA ~ TOE))$r.squared,
              .groups = "drop") %>%
    filter(R2 > 0.7, GR > 0)
}

GR_early <- gr_per_plant(PS2_early2)
GR_late  <- gr_per_plant(PS2_late2)

#Make a data frame that has line(pop) means
GR_early3 <- GR_early %>%
group_by(Population, Treatment) %>%
summarise(GR = mean(GR, na.rm = TRUE), n_plants = n(),
          Site = first(Site), Year = first(Year), .groups = 'drop')

GR_late3 <- GR_late %>%
group_by(Population, Treatment) %>%
summarise(GR = mean(GR, na.rm = TRUE), n_plants = n(),
          Site = first(Site), Year = first(Year), .groups = 'drop')

write_csv(GR_early3, "data/GR_early3.csv")
write_csv(GR_late3, "data/GR_late3.csv")

#Combined, in the same shape as PS2_all
GR_all <- bind_rows(GR_early3 %>% mutate(Timepoint = "Early"),
                    GR_late3  %>% mutate(Timepoint = "Late")) %>%
mutate(Timepoint = factor(Timepoint, levels = c("Early", "Late")))

write_csv(GR_all, "data/GR_all.csv")
