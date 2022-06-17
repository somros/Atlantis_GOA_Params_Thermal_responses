# This script explores surface and bottom temperature of each haul from the RACE-GAP data, to compare to the niches we derive from Aquamaps

library(tidyverse)
library(data.table)

# Read data ---------------------------------------------------------------

race <- read.csv('data/RACE/race_catch_by_haul.csv', skip = 5)
haul <- read.csv('data/RACE/Haul Descriptions.csv', fileEncoding = 'UTF-8-BOM')
fg <- read.csv('data/RACE/GOA_Groups.csv')
key <- read.csv('data/RACE/RACE_species_goa_Atlantis_Nov162021.csv')
flagdem <- read.csv('flagdem.csv')

t <- race %>%
  select(Haul.Join.ID, Species.Code, Weight..kg.) %>%
  filter(Weight..kg. > 0) %>%
  left_join((haul %>% select(Haul.Join.ID, Gear.Temperature..C., Surface.Temperature..C.)), by = 'Haul.Join.ID') %>%
  left_join((key %>% select(Species.Code, Atlantis.group)), by = 'Species.Code') %>%
  filter(Atlantis.group != '?') %>%
  left_join((fg %>% select(Code, Name)), by = c('Atlantis.group'='Code')) %>%
  left_join(flagdem, by = c('Atlantis.group'='atlantis_fg')) %>%
  rowwise() %>%
  mutate(Temp = ifelse(flagdemXXX == 1, Gear.Temperature..C., Surface.Temperature..C.)) %>%
  ungroup() %>%
  select(Name, Temp) %>%
  drop_na()

all_fg <- unique(t$Name)

temps_frame <- data.frame(fg = all_fg)

temps <- list()
for(i in 1:length(all_fg)){
  this_t <- t %>% filter(Name == all_fg[i])
  temps[[i]] <- as.data.frame(matrix(quantile(this_t$Temp, probs = c(0.05, 0.5, 0.95)), ncol = 3))
}

temps_frame <- cbind(temps_frame, rbindlist(temps))
temps_frame <- temps_frame %>%
  set_names(c('fg','min','med','max')) %>%
  arrange(fg)

# compare to aquamaps
aquamaps <- read.csv('thermal_niches_aquamaps_5_95_percentiles.csv')

shared <- intersect(temps_frame$fg, aquamaps$Name)

tt <- aquamaps %>%
  filter(Name %in% shared) %>%
  select(Name, min, max) %>%
  set_names(c('Name','min_aquamaps','max_aquamaps')) %>%
  left_join((temps_frame %>% filter(fg %in% shared) %>% select(fg, min, max) %>% set_names(c('Name','min_race','max_race'))), by = 'Name')
