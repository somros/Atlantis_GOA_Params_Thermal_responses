# Alberto Rovellini
# 6/1/2023
# Construct default Q10 responses for all species that we have Aquamaps data for
# The goal will be to use some of these in some runs for the ICES paper
# Let's stay away from invertebrates and endotherms
# Let's apply them to all fish groups
# Expect major changes to model dynamics, something more realistic will be using them for a small subset of species

# here's a limitation to this: these are really high. We could use percentiles 

library(tidyverse)
library(viridis)

# read atlantis groups
atlantis_groups <- read.csv('GOA_Groups.csv')


dat100 <- read.csv('thermal_niches_aquamaps_0_100_percentiles.csv')
dat95 <- read.csv('thermal_niches_aquamaps_5_95_percentiles.csv')

# calculate ratio of optimum/maximum from CEATTLE
ceattle <- data.frame('Par' = c('Cq', 'Tc0', 'Tcm'),
                  'Pollock' = c(2.6, 10, 15),
                  'Cod' = c(2.41, 13.7, 21),
                  'ATF' = c(2.497, 20.512, 26),
                  'Halibut' = c(3.084, 12.97, 18)) 

rownames(ceattle) <- ceattle$Par
ceattle <- ceattle %>% dplyr::select(-Par) %>% t() %>% data.frame() %>% mutate(ratio = Tc0/Tcm)

mean_ratio <- mean(ceattle$ratio) # on average Tc0 = 0.7*Tcm for CEATTLE species

# start from dat100
bioen_pars <- dat100 %>%
  mutate(Cq = 2,
         Tcm = max,
         Tc0 = max * 0.7) %>%
  select(Name, Cq, Tc0, Tcm)

# throw away all that is MAMMAL and BIRD
bioen_pars <- bioen_pars %>%
  filter(Name %in% setdiff(atlantis_groups$Name, (atlantis_groups %>% 
                                                    filter(GroupType %in% c('MAMMAL', 'BIRD')) %>% 
                                                    pull(Name))))

# function to make the Wisconsin curve
make_curve <- function(Cq, Tc0, Tcm){
  
  Tamb <- seq(0.1, 30, 0.1)
  
  Y <- log(Cq) * (Tcm - Tc0 + 2)
  Z <- log(Cq) * (Tcm - Tc0)
  X <- (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
  V <- (Tcm - Tamb) / (Tcm - Tc0)
  
  Tcorr <- V^X * exp((X * (1 - V)))
  
  Tcorr_df <- data.frame(Tamb, Tcorr)
  
  return(Tcorr_df)
  
}


tcorr_frame <- bioen_pars %>%
  mutate(data = purrr::pmap(list(Cq, Tc0, Tcm), make_curve)) %>%
  unnest(data)

tcorr_frame$Tcorr[is.nan(tcorr_frame$Tcorr)] <- NA

# drop NA's but then add 0 as last data point
tcorr_frame <- tcorr_frame %>%
  drop_na()

tcorr_limit <- tcorr_frame %>%
  dplyr::select(-Tcorr) %>%
  group_by(Name,Cq,Tc0,Tcm) %>%
  arrange(Tamb) %>%
  slice_tail() %>%
  ungroup() %>%
  mutate(Tamb = Tamb + 0.1, Tcorr = 0)

tcorr_frame <- rbind(tcorr_frame, tcorr_limit) %>%
  arrange(Name, Tamb)
  

# plot for fish only
fish_grps <- atlantis_groups %>% filter(GroupType == 'FISH') %>% pull(Name)

# also remove ATF, POL, COD, and HAL

# add guilds for plot
guilds <- read.csv('../../Paper1/data/fg_to_guild.csv')
guilds$fg <- gsub('_N','',guilds$fg)

tcorr_frame <- tcorr_frame %>% left_join(guilds, by = c('Name'='fg'))

# add code
tcorr_frame <- tcorr_frame %>% left_join(atlantis_groups %>% select(Name, Code), by = 'Name')

p1 <- tcorr_frame %>%
  filter(Name %in% fish_grps) %>%
  filter(Name != 'Pollock', Name != 'Cod', Name != 'Arrowtooth_flounder', Name != 'Halibut') %>%
  ggplot(aes(x = Tamb, y = Tcorr, color = Name))+
  geom_line(linewidth = 2)+
  #scale_color_viridis_d(begin = 0.1, end = 0.9)+
  theme_bw()+
  labs(x = 'Temperature (C)', 'Tcorr')

p1

# save a table with the new values for biol_prm
# keep the defaults for everything tat is not a fish - these will still use method 0 anyway

# add functional group Code
bioen_pars <- bioen_pars %>% 
  filter(Name %in% fish_grps) %>%
  left_join(atlantis_groups %>% select(Code, Name), by = 'Name')

# codes to skip
skip <- data.frame(Code = c('POL','COD','HER','ATF','HAL','CAP','SAN','FOS','EUL'),
                   coefftA = c(2.6,2.41,2.6,2.497,3.084,2.6,2.6,2.6,2.6),
                   q10_optimal_temp = c(10,13.7,10,20.512,12.97,10,10,10,10),
                   q10_correction = c(15,21,15,26,18,15,15,15,15)) # do forage fish separate

# coefftA
file_params <- 'dome_shaped_parameters_Aquamaps.txt'
file.create(file_params)

for(i in atlantis_groups$Code){
  
  if(i %in% skip$Code){
    this_coefftA <- skip %>% filter(Code == i) %>% pull(coefftA)
  } else if(i %in% bioen_pars$Code){
    this_coefftA <- bioen_pars %>% filter(Code == i) %>% pull(Cq)
  } else {
    this_coefftA <- 0.851 # default
  }
  
  cat(paste0('temp_coefftA_',i, ' ', this_coefftA, '\n'), file = file_params, append = T)
}

cat('\n\n', file = file_params, append = T)

# q10_optimal_temp

for(i in atlantis_groups$Code){
  
  if(i %in% skip$Code){
    this_q10_optimal_temp <- skip %>% filter(Code == i) %>% pull(q10_optimal_temp)
  } else if(i %in% bioen_pars$Code){
    this_q10_optimal_temp <- bioen_pars %>% filter(Code == i) %>% pull(Tc0)
  } else {
    this_q10_optimal_temp <- 15 # default
  }
  
  cat(paste0('q10_optimal_temp_',i, ' ', this_q10_optimal_temp, '\n'), file = file_params, append = T)
}

cat('\n\n', file = file_params, append = T)

# q10_correction

for(i in atlantis_groups$Code){
  
  if(i %in% skip$Code){
    this_q10_correction <- skip %>% filter(Code == i) %>% pull(q10_correction)
  } else if(i %in% bioen_pars$Code){
    this_q10_correction <- bioen_pars %>% filter(Code == i) %>% pull(Tcm)
  } else {
    this_q10_correction <- 1000 # default
  }
  
  cat(paste0('q10_correction_',i, ' ', this_q10_correction, '\n'), file = file_params, append = T)
}

cat('\n\n', file = file_params, append = T)

# q10_method
for(i in atlantis_groups$Code){
  
  if(i %in% skip$Code){
    this_q10_method <- 2
  } else if(i %in% bioen_pars$Code){
    this_q10_method <- 2
  } else {
    this_q10_method <- 0 # default
  }
  
  cat(paste0('q10_method_',i, ' ', this_q10_method, '\n'), file = file_params, append = T)
}

cat('\n\n', file = file_params, append = T)

# flagq10eff

for(i in atlantis_groups$Code){
  
  if(i %in% skip$Code){
    this_flagq10eff <- 2
  } else if(i %in% bioen_pars$Code){
    this_flagq10eff <- 2
  } else {
    this_flagq10eff <- 0 # default
  }
  
  cat(paste0('flagq10eff',i, ' ', this_flagq10eff, '\n'), file = file_params, append = T)
}

cat('\n\n', file = file_params, append = T)