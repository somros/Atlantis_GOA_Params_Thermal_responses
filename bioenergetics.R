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

# can we fit a model here?

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
tcorr_fit <- optim(start_tcorr, fit_at_tcorr, dat = this_data)
exp(tcorr_fit$par)

fit_A <- exp(tcorr_fit$par[1])
fit_B <- exp(tcorr_fit$par[2])
fit_C <- exp(tcorr_fit$par[3])
fit_D <- exp(tcorr_fit$par[4])
fit_corr <- exp(tcorr_fit$par[5])

# make predictions and see how they align
T_test <- seq(0,30,0.1)
fit_tcorr <- log(2) * fit_A * fit_B^T_test * exp(-fit_C * (abs(T_test - 13.7)^fit_D) / fit_corr)

plot(T_test, fit_tcorr, 'l') # not awesome. Sensitive to some of the initial values (functional response is very sensititve to parameter changes)


# try all at once

tcorr_frame <- tcorr_frame %>% mutate(across(-Tamb, ~replace(., is.nan(.), 0)))

fit_at_tcorr <- function(par, dat){
  
  # par <- start_tcorr
  dat <- tcorr_frame

  # parameters
  coeffA_POL <- exp(par[1])
  coeffA_COD <- exp(par[2])
  coeffA_ATF <- exp(par[3])
  
  coeffB <- exp(par[4]) # global
  
  # optimum T - known (from Holsman and Aydin (2015))
  Tc0_POL <- 10 
  Tc0_COD <- 13.7 
  Tc0_ATF <- 20.512 
  
  coeffC <- exp(par[5]) # global
  coeffD <- exp(par[6]) 
  
  corr_POL <- exp(par[7]) # correction factor
  corr_COD <- exp(par[8])
  corr_ATF <- exp(par[9])
  
  #data - Tcorr from Holsman and Aydin (2015)
  Tamb <- dat$Tamb
  Tcorr_POL <- dat$Pollock
  Tcorr_COD <- dat$Cod
  Tcorr_ATF <- dat$ATF
  
  #calculate Tcorr vectors from q10method 1
  Tcorr_at_POL <- log(2) * coeffA_POL * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0_POL)^coeffD) / corr_POL)
  Tcorr_at_COD <- log(2) * coeffA_COD * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0_COD)^coeffD) / corr_COD)
  Tcorr_at_ATF <- log(2) * coeffA_ATF * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0_ATF)^coeffD) / corr_ATF)
  
  # likelihoods
  NLL_POL <- -sum(dnorm(Tcorr_POL, Tcorr_at_POL, log = T))
  NLL_COD <- -sum(dnorm(Tcorr_COD, Tcorr_at_COD, log = T))
  NLL_ATF <- -sum(dnorm(Tcorr_ATF, Tcorr_at_ATF, log = T))
  
  NLL <- NLL_POL + NLL_COD + NLL_ATF
  
  return(NLL)
}

# start from initial values from GG's function, equal for all 3
start_tcorr <- c(rep(log(0.851), 3), log(1.066), log(1), log(3), rep(log(1000), 3))
tcorr_fit <- optim(start_tcorr, fit_at_tcorr, dat = tcorr_frame)
exp(tcorr_fit$par)

fit_A_POL <- exp(tcorr_fit$par[1])
fit_A_COD <- exp(tcorr_fit$par[2])
fit_A_ATF <- exp(tcorr_fit$par[3])

fit_B <- exp(tcorr_fit$par[4])
fit_C <- exp(tcorr_fit$par[5])
fit_D <- exp(tcorr_fit$par[6])

fit_corr_POL <- exp(tcorr_fit$par[7])
fit_corr_COD <- exp(tcorr_fit$par[8])
fit_corr_ATF <- exp(tcorr_fit$par[9])


# make predictions and see how they align
Temp <- seq(0,30,0.1)

fit_tcorr_POL <- log(2) * fit_A_POL * fit_B^Temp * exp(-fit_C * (abs(Temp - 10)^fit_D) / fit_corr_POL)
fit_tcorr_COD <- log(2) * fit_A_COD * fit_B^Temp * exp(-fit_C * (abs(Temp - 13.7)^fit_D) / fit_corr_COD)
fit_tcorr_ATF <- log(2) * fit_A_ATF * fit_B^Temp * exp(-fit_C * (abs(Temp - 20.512)^fit_D) / fit_corr_ATF)

tcorr_frame_fit <- data.frame('Tamb' = Temp, 
                              'Pollock_fit' = fit_tcorr_POL,
                              'Cod_fit' = fit_tcorr_COD,
                              'ATF_fit' = fit_A_ATF)

tcorr_frame_fit_long <- tcorr_frame_fit %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr from fit')

tcorr_frame_fit_long %>% 
  ggplot(aes(x = Tamb, y = `Tcorr from fit`, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()

# Loop --------------------------------------------------------------------

fit_at_tcorr <- function(par, dat){
  
  # dat <- tcorr_frame
  # par <- start_tcorr
  
  # parameters
  coeffA <- c(exp(par[1]), exp(par[2]), exp(par[3]))
  coeffB <- exp(par[4]) # global
  Tc0 <- c(10, 13.7, 20.512) # optimum T - known (from Holsman and Aydin (2015))
  coeffC <- exp(par[5]) # global
  coeffD <- exp(par[6]) # global
  corr <- c(exp(par[7]), exp(par[8]), exp(par[9])) # correction factor
  
  species <- colnames(dat)[-1]
  
  NLL_vec <- rep(0, length(species))
  
  for(i in 1:length(species)){

    #data - Tcorr from Holsman and Aydin (2015)
    this_dat <- dat %>% select(Tamb, species[i]) %>% drop_na() %>% set_names(c('Tamb','Tcorr'))
    Tamb <- this_dat$Tamb
    Tcorr <- this_dat$Tcorr
    
    #calculate Tcorr vectors from q10method 1
    Tcorr_at <- log(2) * coeffA[i] * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0[i])^coeffD) / corr[i])
    
    #likelihood - most likely wrong
    NLL <- -sum(dnorm(Tcorr, Tcorr_at, log = T))
    
    NLL_vec[i] <- NLL
  }
  
  NLL_all <- sum(NLL_vec)
  
  return(NLL_all)
}

# start from initial values from GG's function, equal for all 3
start_tcorr <- c(rep(log(0.851), 3), log(1.066), log(1), log(3), rep(log(1000), 3))
# start_tcorr <- tcorr_fit$par
tcorr_fit <- optim(start_tcorr, fit_at_tcorr, dat = tcorr_frame)
exp(tcorr_fit$par)


fit_A_POL <- exp(tcorr_fit$par[1])
fit_A_COD <- exp(tcorr_fit$par[2])
fit_A_ATF <- exp(tcorr_fit$par[3])

fit_B <- exp(tcorr_fit$par[4])
fit_C <- exp(tcorr_fit$par[5])
fit_D <- exp(tcorr_fit$par[6])

fit_corr_POL <- exp(tcorr_fit$par[7])
fit_corr_COD <- exp(tcorr_fit$par[8])
fit_corr_ATF <- exp(tcorr_fit$par[9])

# make predictions and see how they align
Temp <- seq(0,30,0.1)

fit_tcorr_POL <- log(2) * fit_A_POL * fit_B^Temp * exp(-fit_C * (abs(Temp - 10)^fit_D) / fit_corr_POL)
fit_tcorr_COD <- log(2) * fit_A_COD * fit_B^Temp * exp(-fit_C * (abs(Temp - 13.7)^fit_D) / fit_corr_COD)
fit_tcorr_ATF <- log(2) * fit_A_ATF * fit_B^Temp * exp(-fit_C * (abs(Temp - 20.512)^fit_D) / fit_corr_ATF)

tcorr_frame_fit <- data.frame('Tamb' = Temp, 
                              'Pollock_fit' = fit_tcorr_POL,
                              'Cod_fit' = fit_tcorr_COD,
                              'ATF_fit' = fit_A_ATF)

tcorr_frame_fit_long <- tcorr_frame_fit %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr from fit')

tcorr_frame_fit_long %>% 
  ggplot(aes(x = Tamb, y = `Tcorr from fit`, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()

# Global fixed ------------------------------------------------------------

# set values of the global parameters to the defaults
# fit the rest species by species

set.seed(999)

Tamb <- seq(0,30,0.1) # temperature

fit_at_tcorr <- function(par, dat, Tc0){
  
  coeffA <- exp(par[1])
  coeffB <- 1.066 # global
  coeffC <- 1 # global
  coeffD <- 3 # global
  corr <- exp(par[2])
  
  #data
  Tamb <- dat$Tamb
  Tcorr <- dat$Tcorr
  
  #gg
  Tcorr_at <- log(2) * coeffA * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0)^coeffD) / corr)
  
  obj_fun <- sum((log1p(Tcorr)-log1p(Tcorr_at))^2)
  obj_fun
}

start_tcorr <- c(log(0.851), log(100))

# make data
this_data_POL <- tcorr_frame %>% select(Tamb, Pollock) %>% filter(!is.nan(Pollock)) %>% set_names(c('Tamb','Tcorr'))
this_data_COD <- tcorr_frame %>% select(Tamb, Cod) %>% filter(!is.nan(Cod)) %>% set_names(c('Tamb','Tcorr'))
this_data_ATF <- tcorr_frame %>% select(Tamb, ATF) %>% filter(!is.nan(ATF)) %>% set_names(c('Tamb','Tcorr'))

# fit
tcorr_fit_POL <- optim(start_tcorr, fit_at_tcorr, dat = this_data_POL, Tc0 = 10)
tcorr_fit_COD <- optim(start_tcorr, fit_at_tcorr, dat = this_data_COD, Tc0 = 13.7)
tcorr_fit_ATF <- optim(start_tcorr, fit_at_tcorr, dat = this_data_ATF, Tc0 = 20.512)

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

fit_tcorr_POL <- log(2) * fit_A_POL * 1.066^T_test * exp(-1 * (abs(T_test - 10)^3) / fit_corr_POL)
fit_tcorr_COD <- log(2) * fit_A_COD * 1.066^T_test * exp(-1 * (abs(T_test - 13.7)^3) / fit_corr_COD)
fit_tcorr_ATF <- log(2) * fit_A_ATF * 1.066^T_test * exp(-1 * (abs(T_test - 20.512)^3) / fit_corr_ATF)

fit_tcorr <- data.frame(T_test, fit_tcorr_POL, fit_tcorr_COD, fit_tcorr_ATF) %>%
  set_names(c('Tamb', 'Pollock', 'Cod', 'ATF')) %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr')

fit_tcorr %>% ggplot(aes(x = Tamb, y = Tcorr, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()

# nope. Probably too constrained now.


# Fix C and D -------------------------------------------------------------------

# estimate species-specific together and B

fit_at_tcorr <- function(par, dat){
  
  # dat <- tcorr_frame
  # par <- start_tcorr
  
  # parameters
  coeffA <- c(exp(par[1]), exp(par[2]), exp(par[3]))
  coeffB <- exp(par[4]) # global
  Tc0 <- c(10, 13.7, 20.512) # optimum T - known (from Holsman and Aydin (2015))
  coeffC <- 1 # global
  coeffD <- 3 # global
  corr <- c(exp(par[5]), exp(par[6]), exp(par[7])) # correction factor
  
  species <- colnames(dat)[-1]
  
  obj_vec <- rep(0, length(species))
  
  for(i in 1:length(species)){
    
    #data - Tcorr from Holsman and Aydin (2015)
    this_dat <- dat %>% select(Tamb, species[i]) %>% drop_na() %>% set_names(c('Tamb','Tcorr'))
    Tamb <- this_dat$Tamb
    Tcorr <- this_dat$Tcorr
    
    #calculate Tcorr vectors from q10method 1
    Tcorr_at <- log(2) * coeffA[i] * coeffB^Tamb * exp(-coeffC * (abs(Tamb - Tc0[i])^coeffD) / corr[i])
    
    #likelihood - most likely wrong
    obj_fun <- sum((log1p(Tcorr)-log1p(Tcorr_at))^2)

    obj_vec[i] <- obj_fun
  }
  
  obj_all <- sum(obj_vec)
  
  return(obj_all)
}

# start from initial values from GG's function, equal for all 3
start_tcorr <- c(rep(log(0.851), 3), log(1.066), rep(log(1000), 3))
start_tcorr <- tcorr_fit$par
tcorr_fit <- optim(start_tcorr, fit_at_tcorr, dat = tcorr_frame)
exp(tcorr_fit$par)


fit_A_POL <- exp(tcorr_fit$par[1])
fit_A_COD <- exp(tcorr_fit$par[2])
fit_A_ATF <- exp(tcorr_fit$par[3])

fit_B <- exp(tcorr_fit$par[4])

fit_corr_POL <- exp(tcorr_fit$par[5])
fit_corr_COD <- exp(tcorr_fit$par[6])
fit_corr_ATF <- exp(tcorr_fit$par[7])

# make predictions and see how they align
Temp <- seq(0,30,0.1)

fit_tcorr_POL <- log(2) * fit_A_POL * fit_B^Temp * exp(-1 * (abs(Temp - 10)^3) / fit_corr_POL)
fit_tcorr_COD <- log(2) * fit_A_COD * fit_B^Temp * exp(-1 * (abs(Temp - 13.7)^3) / fit_corr_COD)
fit_tcorr_ATF <- log(2) * fit_A_ATF * fit_B^Temp * exp(-1 * (abs(Temp - 20.512)^3) / fit_corr_ATF)

tcorr_frame_fit <- data.frame('Tamb' = Temp, 
                              'Pollock_fit' = fit_tcorr_POL,
                              'Cod_fit' = fit_tcorr_COD,
                              'ATF_fit' = fit_A_ATF)

tcorr_frame_fit_long <- tcorr_frame_fit %>%
  pivot_longer(-Tamb, names_to = 'Species', values_to = 'Tcorr from fit')

tcorr_frame_fit_long %>% 
  ggplot(aes(x = Tamb, y = `Tcorr from fit`, color = Species))+
  geom_line(size = 2)+
  scale_color_manual(values = c('red','green','blue'))+
  theme_bw()

# also no

