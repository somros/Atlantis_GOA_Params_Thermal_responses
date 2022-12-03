# Alberto Rovellini
# 12/02/2022
# Thermal niche plotter
# Using flagtempsensitive 2

library(tidyverse)

dat <- read.csv('thermal_niches_aquamaps_0_100_percentiles.csv')

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
dat3 <- dat2 %>%
  mutate(niche = purrr::pmap(list(min_sp = mint, max_sp = maxt, current_enviro = temp), make_niche)) %>%
  unnest(cols = c('niche')) %>%
  select(Code, Name, mint, maxt, temp, niche)

# plot
ggplot(dat3)+
  geom_line(aes(x = temp, y = niche), linewidth = 2)+
  geom_vline(aes(xintercept = mint), color = 'blue', linewidth = 2)+
  geom_vline(aes(xintercept = maxt), color = 'red', linewidth = 2)+
  theme_bw()+
  labs(x = 'Temperature (C)', y = 'Scalar')+
  facet_wrap(~Name)
  