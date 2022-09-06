# Bioenergetics from CEATTLE
library(tidyverse)

# Goal: fit Gary Griffith's 6-parameter Q10 function to the bioenergetic curves to calculate f(T) in Holsman et al. (2015, 2019)

# Pollock, Cod, and ATF from: Holsman, K. K., & Aydin, K. (2015). Comparative methods for evaluating climate change impacts on the foraging ecology of Alaskan groundfish. Marine Ecology Progress Series, 521, 217–235. https://doi.org/10.3354/meps11102
# Halibut from: Holsman, KK, Aydin, K, Sullivan, J, Hurst, T, Kruse, GH. Climate effects and bottom-up controls on growth and size-at-age of Pacific halibut (Hippoglossus stenolepis) in Alaska (USA). Fish Oceanogr. 2019; 28: 345– 358. https://doi.org/10.1111/fog.12416

dat <- data.frame('Par' = c('Cq', 'Tc0', 'Tcm'),
                  'Pollock' = c(2.6, 10, 15),
                  'Cod' = c(2.41, 13.7, 21),
                  'ATF' = c(2.497, 20.512, 26),
                  'Halibut' = c(3.084, 12.97, 18))

# TC0 is the temperature where laboratory consumption rates are highest, TCM is the maximum water
# temperature above which consumption ceases and CQ approximates the Q10 or the rate at which 
# the function increases over relatively low water temperatures.

make_curve <- function(species, dat, Tamb){
  
  Cq <- dat %>% filter(Par == 'Cq') %>% pull(species)
  Tc0 <- dat %>% filter(Par == 'Tc0') %>% pull(species)
  Tcm <- dat %>% filter(Par == 'Tcm') %>% pull(species)
  
  Y <- log(Cq) * (Tcm - Tc0 + 2)
  Z <- log(Cq) * (Tcm - Tc0)
  X <- (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
  V <- (Tcm - Tamb) / (Tcm - Tc0)
  
  Tcorr <- V^X * exp((X * (1 - V)))
  
  return(Tcorr)
  
}

tcorr_frame <- data.frame('Tamb' = seq(0, 30, 0.1)) %>%
  mutate(Pollock = make_curve('Pollock', dat, Tamb),
         Cod = make_curve('Cod', dat, Tamb),
         ATF = make_curve('ATF', dat, Tamb),
         Halibut = make_curve('Halibut', dat, Tamb)) 

tcorr_frame_long <- tcorr_frame %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr')


p1 <- tcorr_frame_long %>%
  ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('orange','blue','grey','black'))+
  theme_bw()+
  labs(x = 'Temperature (C)', 'Tcorr')

p1
ggsave('CEATTLE.png', p1, width = 6, height = 4)


# Atlantis Q10method = 1 --------------------------------------------------
# 6-parameter function

# thi is the original parametrization in the manual/code
coeffA <- 0.851
coeffB <- 1.066 # global
Tc0 <- 19 # optimum for each species
coeffC <- 1 # global
coeffD <- 3 # global
corr <- 1000

Tamb <- seq(0,30,0.1) # temperature

tcorr_at <- log(2) * coeffA * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0)^coeffD) / corr)

plot(Tamb, tcorr_at, 'l')

# Fit model for single species --------------------------------------------

set.seed(999)

all_species <- unique(tcorr_frame_long$Species)
all_Tc0 <- unlist(dat[2,2:5])

par_vector <- matrix(0, nrow = length(all_species), ncol = 5)

# as initial values for the parameters, use the defaults from GG
start_tcorr <- c(log(0.851), log(1.066), log(1), log(3), log(1000)) # these are the initial values in the Atlantis manual

for(i in 1:length(all_species)){
  
  this_data <- tcorr_frame %>% 
    select(Tamb, all_species[i]) %>% 
    set_names(c('Tamb','Tcorr')) %>% 
    filter(!is.nan(Tcorr)) 
  
  this_Tc0 <- all_Tc0[i]
  
  fit_at_tcorr <- function(par, dat, Tc0){
    
    coeffA <- exp(par[1])
    coeffB <- exp(par[2]) # global
    coeffC <- exp(par[3]) # global
    coeffD <- exp(par[4]) # global
    corr <- exp(par[5])
    
    #data
    Tamb <- dat$Tamb
    Tcorr <- dat$Tcorr
    
    #gg
    Tcorr_at <- log(2) * coeffA * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0)^coeffD) / corr)
    
    obj_fun <- sum((log1p(Tcorr)-log1p(Tcorr_at))^2)
    obj_fun
  }
  
  tcorr_fit <- optim(start_tcorr, fit_at_tcorr, dat = this_data, Tc0 = this_Tc0)
  
  par_vector[i,] <- exp(tcorr_fit$par)
}

# Global fixed ------------------------------------------------------------

# set values of the global parameters to some defaults: defaults chosen based on values that work with some of the single-species fitting above
# taking an average of them will not work, q10 1 is EXTREMELY sensitive to input values
# fitting all ot once will not work with optim(), probably way overparameterized

# fit the rest species by species

Tamb <- seq(0,30,0.1) # temperature

fit_at_tcorr <- function(par, dat, Tc0){
  
  # par <- start_tcorr_POL
  # dat <- this_data_POL
  # Tc0 <- 10
  
  coeffA <- exp(par[1])
  coeffB <- 0.95 # global
  coeffC <- 7 # global
  coeffD <- 1.2 # global
  corr <- exp(par[2])
  
  #data
  Tamb <- dat$Tamb
  Tcorr <- dat$Tcorr
  
  #gg
  Tcorr_at <- log(2) * coeffA * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0)^coeffD) / corr)
  
  obj_fun <- sum((log1p(Tcorr)-log1p(Tcorr_at))^2)
  obj_fun
}

# use as starting values the parameter estimates from the first fit
start_tcorr_POL <- c(log(par_vector[1,1]), log(par_vector[1,5]))
start_tcorr_COD <- c(log(par_vector[2,1]), log(par_vector[2,5]))
start_tcorr_ATF <- c(log(par_vector[3,1]), log(par_vector[3,5]))
start_tcorr_HAL <- c(log(par_vector[4,1]), log(par_vector[4,5]))

# make data
this_data_POL <- tcorr_frame %>% select(Tamb, Pollock) %>% filter(!is.nan(Pollock)) %>% set_names(c('Tamb','Tcorr'))
this_data_COD <- tcorr_frame %>% select(Tamb, Cod) %>% filter(!is.nan(Cod)) %>% set_names(c('Tamb','Tcorr'))
this_data_ATF <- tcorr_frame %>% select(Tamb, ATF) %>% filter(!is.nan(ATF)) %>% set_names(c('Tamb','Tcorr'))
this_data_HAL <- tcorr_frame %>% select(Tamb, Halibut) %>% filter(!is.nan(Halibut)) %>% set_names(c('Tamb','Tcorr'))

# fit
tcorr_fit_POL <- optim(start_tcorr_POL, fit_at_tcorr, dat = this_data_POL, Tc0 = 10)
tcorr_fit_COD <- optim(start_tcorr_COD, fit_at_tcorr, dat = this_data_COD, Tc0 = 13.7)
tcorr_fit_ATF <- optim(start_tcorr_ATF, fit_at_tcorr, dat = this_data_ATF, Tc0 = 20.512)
tcorr_fit_HAL <- optim(start_tcorr_HAL, fit_at_tcorr, dat = this_data_HAL, Tc0 = 12.97)

#pol
fit_A_POL <- exp(tcorr_fit_POL$par[1])
fit_corr_POL <- exp(tcorr_fit_POL$par[2])
#cod
fit_A_COD <- exp(tcorr_fit_COD$par[1])
fit_corr_COD <- exp(tcorr_fit_COD$par[2])
#atf
fit_A_ATF <- exp(tcorr_fit_ATF$par[1])
fit_corr_ATF <- exp(tcorr_fit_ATF$par[2])
#halibut
fit_A_HAL <- exp(tcorr_fit_HAL$par[1])
fit_corr_HAL <- exp(tcorr_fit_HAL$par[2])

# make predictions and see how they align
T_test <- seq(0,30,0.1)

fit_tcorr_POL <- log(2) * fit_A_POL * 0.95^T_test * exp(-7 * (abs(T_test - 10)^1.2) / fit_corr_POL)
fit_tcorr_COD <- log(2) * fit_A_COD * 0.95^T_test * exp(-7 * (abs(T_test - 13.7)^1.2) / fit_corr_COD)
fit_tcorr_ATF <- log(2) * fit_A_ATF * 0.95^T_test * exp(-7 * (abs(T_test - 20.512)^1.2) / fit_corr_ATF)
fit_tcorr_HAL <- log(2) * fit_A_HAL * 0.95^T_test * exp(-7 * (abs(T_test - 12.97)^1.2) / fit_corr_ATF)

fit_tcorr <- data.frame(T_test, fit_tcorr_POL, fit_tcorr_COD, fit_tcorr_ATF, fit_tcorr_HAL) %>%
  set_names(c('Tamb', 'Pollock', 'Cod', 'ATF', 'Halibut')) %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr')

fit_tcorr %>% ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('orange','blue','grey','black'))+
  theme_bw()+
  labs(x = 'Temperature (C)', 'Tcorr')

write.csv(data.frame('Parameter' = c('temp_coefftA_POL','q10_correction_POL',
  'temp_coefftA_COD','q10_correction_COD',
  'temp_coefftA_ATF','q10_correction_ATF',
  'temp_coefftA_HAL','q10_correction_HAL',
  'temp_coeffB',
  'temp_coeffC',
  'temp_exp'),
  'Value' = c(fit_A_POL, fit_corr_POL, fit_A_COD, fit_corr_COD, fit_A_ATF, fit_corr_ATF, fit_A_HAL, fit_corr_HAL, 0.95, 7, 1.2)),
  file = 'params_q10.csv', row.names = F)

# Merge with thermal maxima ----------------------------------------------

# It is interesting to note that the thermal optima and maximum temperatures we get from Aquamaps are 
# generally lower than the values used in Holsman and Aydin (2015)
# TODO: check those as it will be important to get those right

# Tcm is the temperature after which consumption ceased in the lab setting from Holsman and Aydin (2015)

fit_tcorr1 <- fit_tcorr %>%
  rowwise() %>%
  mutate(Tcm = ifelse(Species == 'Pollock', 15, ifelse(Species == 'Cod', 21, ifelse(Species == 'ATF', 26, 18)))) %>%
  mutate(Tcorr = ifelse(Tamb <= Tcm, Tcorr, 0))

fit_tcorr1 %>% ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue','black'))+
  theme_bw()

# Plot for methods section ------------------------------------------------
# make a plot for the methods
t <- tcorr_frame_long %>% mutate(Source = 'Wisconsin')
t1 <- fit_tcorr %>% mutate(Source = 'Atlantis')
t2 <- rbind(t,t1)

t2$Source <- factor(t2$Source, levels = c('Wisconsin','Atlantis'))

p2 <- t2 %>% ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('orange','blue','grey','black'))+
  theme_bw()+
  labs(x = 'Temperature (C)', y = expression(T[corr]))+ # chage this to a damn subscript
  facet_wrap(~Source)
p2

ggsave('fit.png', p2, width = 8, height = 4)


# Merge q10 curves with thermal windows -----------------------------------
niches <- read.csv('thermal_niches_aquamaps_0_100_percentiles.csv')
niches <- niches %>%
  filter(Name %in% c('Pollock', 'Cod', 'Arrowtooth_flounder', 'Halibut')) %>%
  select(Name, min, max)

niches$Name <- gsub('Arrowtooth_flounder', 'ATF', niches$Name)

fit_tcorr2 <- fit_tcorr %>%
  left_join(niches, by = c('Species'='Name')) %>%
  rowwise() %>%
  mutate(Tcorr1 = ifelse(Tamb >= min & Tamb <= max, Tcorr, 0)) %>%
  ungroup()

p3 <- fit_tcorr2 %>% ggplot(aes(x = Tamb, y = Tcorr1, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('orange','blue','grey','black'))+
  theme_bw()+
  labs(x = 'Temperature (C)', 'Tcorr')
p3

ggsave('window_and_bioenergetics.png', p3, width = 6, height = 4)
