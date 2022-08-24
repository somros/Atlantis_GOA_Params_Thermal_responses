# this code reads files of spatial distributions downloaded from AquaMaps.org 
# each file contains temperature data for the cells where the species has been found to occur
# These are intended to be used as high-level thermal windows for species in Atlantis
# Each data set reports both surface and bottom temperatures. Szymon seems to have used surface temperatures, which result in thermal ranges that are possibly too strict for some bottom-dwelling species

# here we pull surface and bottom and we get to decide for each species which to use depending on where they live
# Since in Atlantis these will determine whether or not a species will be present in a cell at all or not, I am inclined to relax the lower bounds (cold water)
# using 10th and 90th percentiles may be a good idea here to avod outliers

library(data.table)
library(tidyverse)

niche_files <- list.files('data/aquamaps', full.names = T)
key <- read.csv('species_list.csv')
atlantis_groups <- read.csv('GOA_Groups.csv')
flagdem <- read.csv('flagdem.csv')



make_niche <- function(this_file){
  
  this_species <- word(readLines(this_file,1), 4, 5)
  
  dat <- read.csv(this_file, skip=7) %>%
    filter(between(Surface.Water.Temp....C.,-10,40), Center.Lat > 0, Center.Long > 120 | Center.Long < -100) 
  
  min_bottom_000 <- dat %>% pull(Bottom.Water.Temp....C.) %>% min()
  min_bottom_005 <- dat %>% pull(Bottom.Water.Temp....C.) %>% quantile(probs = seq(0,1,.05)) %>% .[2]
  max_bottom_100 <- dat %>% pull(Bottom.Water.Temp....C.) %>% max()
  max_bottom_095 <- dat %>% pull(Bottom.Water.Temp....C.) %>% quantile(probs = seq(0,1,.05)) %>% .[20]
  
  min_surface_000 <- dat %>% pull(Surface.Water.Temp....C.) %>% min()
  min_surface_005 <- dat %>% pull(Surface.Water.Temp....C.) %>% quantile(probs = seq(0,1,.05)) %>% .[2]
  max_surface_100 <- dat %>% pull(Surface.Water.Temp....C.) %>% max()
  max_surface_095 <- dat %>% pull(Surface.Water.Temp....C.) %>% quantile(probs = seq(0,1,.05)) %>% .[20]
  
  this_niche <- data.frame(this_species, min_bottom_000, max_bottom_100, min_bottom_005, max_bottom_095, 
                           min_surface_000, max_surface_100, min_surface_005, max_surface_095, row.names = NULL)
  
}

all_niches <- rbindlist(lapply(niche_files,make_niche))

# bind species with Atlantis groups

# rename capelin
all_niches$this_species <- gsub('villosus','catervarius',all_niches$this_species)
all_niches$this_species <- gsub('Beringraja rhina','Raja rhina',all_niches$this_species)

all_niches <- all_niches %>%
  left_join(key, by = c('this_species'='Species'))

# take averages by group

all_niches_fg <- all_niches %>%
  select(Code, min_bottom_000:max_surface_095) %>%
  left_join((atlantis_groups %>% select(Code,Name)), by='Code') %>%
  group_by(Code,Name) %>%
  summarise(across(everything(), ~mean(.x))) %>%
  ungroup()

# make some plots#####################################################
# reformat so that it is convenient for ggplot
all_niches1 <- all_niches_fg %>%
  pivot_longer(!c(Code,Name), names_to = 'vartype', values_to = 't') %>%
  mutate(bound=substr(vartype,1,3),
         location=substr(vartype,5,5),
         percentile=substr(vartype,(nchar(vartype)-2),nchar(vartype))) %>%
  select(-vartype) %>%
  rowwise() %>%
  mutate(percentile = ifelse(percentile %in% c('000','100'), 'all','05-95')) %>%
  ungroup() %>%
  # mutate(vartype=paste(bound,percentile,sep='_')) %>% # make a new vartype
  # select(-bound,-percentile) %>%
  pivot_wider(names_from = bound, values_from = t) %>%
  rowwise() %>%
  mutate(opt = (min+max)/2) %>%
  ungroup()

# some renaming
all_niches1$location <- gsub('b','Bottom',all_niches1$location)
all_niches1$location <- gsub('s','Surface',all_niches1$location)
all_niches1$percentile <- gsub('all','Full range',all_niches1$percentile)
all_niches1$percentile <- gsub('05-95','5 and 95 percentiles',all_niches1$percentile)


# view
ggplot(all_niches1, aes(x=location, y=opt, color=percentile))+
  geom_linerange(aes(min=min,max=max), position = position_dodge(width = 0.3), size =1.5)+
  theme_bw()+
  labs(y='Temperature',x='',color='')+
  facet_wrap(~Name, ncol=10)

# here we decide for each fg whether we use bottom or surface temperature as limiting
pull_temp <- function(fg,flagdem,is_minimum,fullrange){
  if(isTRUE(fullrange)){
    if(flagdem>0){
      if(isTRUE(is_minimum)){
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Bottom',percentile=='Full range') %>% pull(min)
      } else {
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Bottom',percentile=='Full range') %>% pull(max)
      }
    } else {
      if(isTRUE(is_minimum)){
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Surface',percentile=='Full range') %>% pull(min)
      } else {
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Surface',percentile=='Full range') %>% pull(max)
      }
    }
  } else {
    if(flagdem>0){
      if(isTRUE(is_minimum)){
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Bottom',percentile=='5 and 95 percentiles') %>% pull(min)
      } else {
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Bottom',percentile=='5 and 95 percentiles') %>% pull(max)
      }
    } else {
      if(isTRUE(is_minimum)){
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Surface',percentile=='5 and 95 percentiles') %>% pull(min)
      } else {
        this_bound <- all_niches1 %>% filter(Code==fg,location=='Surface',percentile=='5 and 95 percentiles') %>% pull(max)
      }
    }
  }
  return(this_bound)
}

niches_final <- all_niches1 %>% select(Code, Name) %>% distinct() %>%
  left_join(flagdem, by = c('Code'='atlantis_fg')) %>%
  mutate(min = purrr::pmap(list(Code,flagdemXXX,T,T),pull_temp),
         max = purrr::pmap(list(Code,flagdemXXX,F,T),pull_temp)) %>%
  unnest(cols = c(min,max))

# view 
p <- niches_final %>%
  mutate(opt=(min+max)/2) %>% # this may not be true
  rowwise() %>%
  mutate(Type = ifelse(flagdemXXX == 1, 'Demersal', 'Pelagic')) %>%
  ungroup() %>%
  ggplot(aes(x=Name, y=opt))+
  geom_point(aes(color = factor(Type)))+
  geom_linerange(aes(min=min, max=max, color = factor(Type)))+
  scale_color_manual(values = c('orange','blue'))+
  scale_y_continuous(breaks = seq(0,30,2), labels = seq(0,30,2))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), text = element_text(size = 14))+
  labs(x='', y='Temperature (C)', color = '')
p

p <- niches_final %>%
  mutate(opt=(min+max)/2) %>%
  rowwise() %>%
  mutate(Type = ifelse(flagdemXXX == 1, 'Demersal', 'Pelagic')) %>%
  ungroup() %>%
  ggplot(aes(x=Name, y=opt))+
  #geom_point(aes(color = factor(Type)))+
  geom_linerange(aes(min=min, max=max, color = factor(Type)), size = 1.5)+
  scale_color_manual(values = c('orange','blue'))+
  scale_y_continuous(breaks = seq(0,30,2), labels = seq(0,30,2))+
  coord_flip()+
  theme_bw()+
  theme(text = element_text(size = 14))+
  labs(x='', y='Temperature (C)', color = '')
p

ggsave('niches_0_100.png', p, width = 8, height = 8)

# write out

niches_final %>% 
  left_join((atlantis_groups %>% select(Index, Code)), by = 'Code') %>%
  arrange (Index) %>%
  write.csv('thermal_niches_aquamaps_0_100_percentiles.csv', row.names = F)

  
