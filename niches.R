# this script takes averages for atlantis functional groups of the thermal niches that Szymon pulled from AquaMaps.
# If large differences are seen within a functional group, we should average with the weights from RACE-GAP, but probably not crucial as this is already high-level


library(tidyverse)

niches <- read.csv('niches.csv')
atlantis_groups <- read.csv('GOA_Groups.csv')

atlantis_groups <- atlantis_groups %>% select(Code,Name)

niches1 <- niches %>%
  select(Atlantis.group,Preferred.minimum.SST..ºC...10th.percentile.,Preferred.maximum.SST..ºC...90th.percentile.) %>%
  set_names(c('Name','Minimum','Maximum')) %>%
  left_join(atlantis_groups, by = 'Name') %>%
  drop_na()

# take averages by group

niches2 <- niches1 %>%
  group_by(Name, Code) %>%
  summarise(Minimum=mean(Minimum),
            Maximum=mean(Maximum))
