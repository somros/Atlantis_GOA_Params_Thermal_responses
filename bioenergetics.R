# Bioenergetics from CEATTLE
library(tidyverse)

# from: Holsman, K. K., & Aydin, K. (2015). Comparative methods for evaluating climate change impacts on the foraging ecology of Alaskan groundfish. Marine Ecology Progress Series, 521, 217–235. https://doi.org/10.3354/meps11102

dat <- data.frame('Par' = c('Cq', 'Tc0', 'Tcm'),
                  'Pollock' = c(2.6, 10, 15),
                  'Cod' = c(2.41, 13.7, 21),
                  'ATF' = c(2.497, 20.512, 26))

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
         ATF = make_curve('ATF', dat, Tamb)) 

tcorr_frame_long <- tcorr_frame %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr')


tcorr_frame_long %>%
  ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()


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


# Fit model for single species --------------------------------------------

this_data <- tcorr_frame %>% select(Tamb, ATF) %>% filter(!is.nan(ATF)) %>% set_names(c('Tamb','Tcorr'))

fit_at_tcorr <- function(par, dat){
  coeffA <- exp(par[1])
  coeffB <- exp(par[2]) # global
  Tc0 <- 20.512 # optimum for each species
  coeffC <- exp(par[3]) # global
  coeffD <- exp(par[4]) # global
  corr <- exp(par[5])
  
  #data
  Tamb <- dat$Tamb
  Tcorr <- dat$Tcorr
  
  #gg
  Tcorr_at <- log(2) * coeffA * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0)^coeffD) / corr)
  
  NLL <- -sum(dnorm(Tcorr, Tcorr_at, log = T))
  
  NLL
}

start_tcorr <- c(log(0.851), log(1.066), log(1), log(3), log(1000))
tcorr_fit <- optim(tcorr_fit$par, fit_at_tcorr, dat = this_data)
exp(tcorr_fit$par)

fit_A <- exp(tcorr_fit$par[1])
fit_B <- exp(tcorr_fit$par[2])
fit_C <- exp(tcorr_fit$par[3])
fit_D <- exp(tcorr_fit$par[4])
fit_corr <- exp(tcorr_fit$par[5])

# make predictions and see how they align
T_test <- seq(0,30,0.1)
fit_tcorr <- log(2) * fit_A * fit_B^T_test * exp(-fit_C * (abs(T_test - 20.512)^fit_D) / fit_corr)

plot(T_test, fit_tcorr, 'l') # not awesome. Sensitive to some of the initial values (functional response is very sensititve to parameter changes)



# Global fixed ------------------------------------------------------------

# set values of the global parameters to some defaults
# defaults chosen based on single-species fits above
# fit the rest species by species

set.seed(999)

Tamb <- seq(0,30,0.1) # temperature

fit_at_tcorr <- function(par, dat, Tc0){
  
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
start_tcorr_POL <- c(log(2.58), log(83.50))
start_tcorr_COD <- c(log(2.88), log(107.06))
start_tcorr_ATF <- c(log(4.29), log(102.08))

# make data
this_data_POL <- tcorr_frame %>% select(Tamb, Pollock) %>% filter(!is.nan(Pollock)) %>% set_names(c('Tamb','Tcorr'))
this_data_COD <- tcorr_frame %>% select(Tamb, Cod) %>% filter(!is.nan(Cod)) %>% set_names(c('Tamb','Tcorr'))
this_data_ATF <- tcorr_frame %>% select(Tamb, ATF) %>% filter(!is.nan(ATF)) %>% set_names(c('Tamb','Tcorr'))

# fit
tcorr_fit_POL <- optim(start_tcorr_POL, fit_at_tcorr, dat = this_data_POL, Tc0 = 10)
tcorr_fit_COD <- optim(start_tcorr_COD, fit_at_tcorr, dat = this_data_COD, Tc0 = 13.7)
tcorr_fit_ATF <- optim(start_tcorr_ATF, fit_at_tcorr, dat = this_data_ATF, Tc0 = 20.512)

#pol
fit_A_POL <- exp(tcorr_fit_POL$par[1])
fit_corr_POL <- exp(tcorr_fit_POL$par[2])
#cod
fit_A_COD <- exp(tcorr_fit_COD$par[1])
fit_corr_COD <- exp(tcorr_fit_COD$par[2])
#atf
fit_A_ATF <- exp(tcorr_fit_ATF$par[1])
fit_corr_ATF <- exp(tcorr_fit_ATF$par[2])

# make predictions and see how they align
T_test <- seq(0,30,0.1)

fit_tcorr_POL <- log(2) * fit_A_POL * 0.95^T_test * exp(-7 * (abs(T_test - 10)^1.2) / fit_corr_POL)
fit_tcorr_COD <- log(2) * fit_A_COD * 0.95^T_test * exp(-7 * (abs(T_test - 13.7)^1.2) / fit_corr_COD)
fit_tcorr_ATF <- log(2) * fit_A_ATF * 0.95^T_test * exp(-7 * (abs(T_test - 20.512)^1.2) / fit_corr_ATF)

fit_tcorr <- data.frame(T_test, fit_tcorr_POL, fit_tcorr_COD, fit_tcorr_ATF) %>%
  set_names(c('Tamb', 'Pollock', 'Cod', 'ATF')) %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr')

fit_tcorr %>% ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()


# Merge with thermal maxima ----------------------------------------------

# It is interesting to note that the thermal optima and maximum temperatures we get from Aquamaps are 
# generally lower than the values used in Holsman and Aydin (2015)
# TODO: check those as it will be important to get those right

# Tcm is the temperature after which consumption ceased in the lab setting from Holsman and Aydin (2015)

fit_tcorr <- fit_tcorr %>%
  rowwise() %>%
  mutate(Tcm = ifelse(Species == 'Pollock', 15, ifelse(Species == 'Cod', 21, 26))) %>%
  mutate(Tcorr = ifelse(Tamb <= Tcm, Tcorr, 0))

fit_tcorr %>% ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()
