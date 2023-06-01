# Alberto Rovellini
# 12/02/2022
# Thermal niche plotter
# Using flagtempsensitive 2

library(tidyverse)
library(viridis)

dat <- read.csv('thermal_niches_aquamaps_0_100_percentiles.csv')

# added for AMSS 1/17/2023
# Change COD, POL, ATF, HAL so that max is same as the maximum in the bioenergetics
dat$max[dat$Code=='POL'] <- 15
dat$max[dat$Code=='COD'] <- 21
dat$max[dat$Code=='ATF'] <- 26
dat$max[dat$Code=='HAL'] <- 18

temp <- seq(-5,35,0.1)
k <- 2 # this is same for all for now

make_niche <- function(min_sp, max_sp, current_enviro, K_const_sp = k){
  
  step1 <- K_const_sp * exp(current_enviro - min_sp) / (K_const_sp + (exp(current_enviro - min_sp) - 1.0))
  step2 <- K_const_sp * exp(max_sp - current_enviro) / (K_const_sp + (exp(max_sp - current_enviro) - 1.0))
  
  # case sensitive_biologistic_window: // Gaussian shape
  numScalar <- 1.0
  step3 <- step1 / K_const_sp
  step4 <- step2 / K_const_sp
  if (step3 > step4) {
    numScalar <- step4
  } else {
    numScalar <- step3
  }
  
  return(numScalar)
}

dat1 <- dat %>%
  rename(mint = min, maxt = max)

# atach temperature as column
dat2 <- do.call("rbind", replicate(length(temp), dat1, simplify = FALSE)) %>%
  arrange(Index) %>%
  mutate(temp = rep(temp, nrow(dat)))

# apply function
to_show <- c('Arrowtooth_flounder','Cod','Halibut','Pollock')

dat3 <- dat2 %>%
  mutate(niche = purrr::pmap(list(min_sp = mint, max_sp = maxt, current_enviro = temp), make_niche)) %>%
  unnest(cols = c('niche')) %>%
  select(Code, Name, mint, maxt, temp, niche)

dat_vline <- dat3 %>%
  select(Code, Name, mint, maxt) %>%
  filter(Name %in% to_show)%>%
  pivot_longer(-c(Code, Name), names_to = 'edge', values_to = 'temp')
  
# plot
p1 <- dat3 %>%
  filter(Name %in% to_show)%>%
  ggplot()+
  geom_line(aes(x = temp, y = niche), linewidth = 1.5)+
  geom_vline(data = dat_vline, aes(xintercept = temp, color = edge), linewidth = 1.5, linetype = 'dashed')+
  # geom_vline(aes(xintercept = mint), color = 'blue', linewidth = 1.5, linetype = 'dashed')+
  # geom_vline(aes(xintercept = maxt), color = 'orange', linewidth = 1.5, linetype = 'dashed')+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  theme_bw()+
  labs(x = 'Temperature (C)', y = 'Scalar on abundance', color = '')+
  facet_wrap(~Name, ncol = 4)
p1
ggsave('../../../ECCWO/code/output/thermal_niche_eccwo.png', p1, width = 8.5, height = 2.5)
  