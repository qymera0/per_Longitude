
# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(rethinking)
library(fs)

stanPath <- 'C:/Users/Administrator/Documents/.cmdstan/cmdstan-2.30.1'

set_cmdstan_path(stanPath)

# SIMULATE MODEL ----------------------------------------------------------

t <- seq(0.1, 30, length = 100) #simulating 30 days, don't start at 0 to avoid 0/inf in plot

alpha <- 20; 

beta <- 2; #just some values to show shape

mu <- log(t^alpha*exp(-beta*t)) #log virus load

plot(t, mu, type = "l", ylim = c(0,30)) #looks somewhat like virus load in acute infections

# SIMULATE DATA -----------------------------------------------------------

set.seed(123)

# days at which we assume outcome is measured

timevec <- c(0.1,1,3,5,7,10,14,21,28,35,42)

#different number of individuals per dose to make it clearer which is which
#also, that's the structure of the data which motivated the tutorial

Nlow <- 7

Nmed <- 8

Nhigh <- 9

filename <- "simdat.Rds"

#if you want to explore how model fitting changes if you increase sample size
#turn on this line of code
#this is used in part 4 of the tutorial
#Nlow = 70; Nmed = 80; Nhigh = 90; filename = "simdat_big.Rds"

Ntot <- Nlow + Nmed + Nhigh #total number of individuals

# Set values for dose
# since we only consider dose on a log scale
# we'll log transform right here and then always use it in those log units

high_dose <- log(1000)

med_dose <- log(100)

low_dose <- log(10)

dosevec <- c(rep(low_dose,Nlow),rep(med_dose,Nmed),rep(high_dose,Nhigh))

# we are also creating a version of the dose variable
# that consists of ordered categories instead of numeric values
# we'll use that mostly for plotting

dosevec_cat <- 
        ordered(
                c(
                        rep("low", Nlow),
                        rep("medium",Nmed),
                        rep("high",Nhigh)),
                levels = c("low","medium","high")
)

## Setting parameters level ------------------------------------------------

sigma <- 1

a1 <- 0.1

b1 <- -0.1


### Model 01 --------------------------------------------------------------

m1_mua <- 3

m1_mub <- 1

m1_sigmaa <- 1

m1_sigmab <- 1

m1_a0 <- rnorm(n = Ntot, m1_mua, m1_sigmaa)

m1_b0 <- rnorm(n = Ntot, m1_mub, m1_sigmaa)

m1pars <-
        c(sigma = sigma, 
          a1 = a1, 
          b1 = b1,
          a0_mu = m1_mua, 
          b0_mu = m1_mub
)

### Model 02 --------------------------------------------------------------

m2_mua <- 3

m2_mub <- 1

m2_sigmaa <- 0.0001

m2_sigmab <- 0.0001

m2_a0 <- rnorm(n = Ntot, m2_mua, m2_sigmaa)

m2_b0 <- rnorm(n = Ntot, m2_mub, m2_sigmab)

m2pars <- 
        c(sigma = sigma, 
          a1 = a1, 
          b1 = b1,
          a0_mu = m2_mua, 
          b0_mu = m2_mub
)


### Model 3 ---------------------------------------------------------------

m3_mua <- 3

m3_mub <- 1

m3_sigmaa <- 0.1

m3_sigmab <- 0.1

m3_a0 <- rnorm(n = Ntot, m3_mua, m3_sigmaa)

m3_b0 <- rnorm(n = Ntot, m3_mub, m3_sigmaa)

m3pars <-
        c(sigma = sigma, 
          a1 = a1, 
          b1 = b1,
          a0_mu = m3_mua, 
          b0_mu = m3_mub)


## Simulate data ----------------------------------------------------------


### Model 1 ---------------------------------------------------------------

m1_alpha <- m1_a0 + a1*(dosevec - med_dose)

m1_beta <- m1_b0 + b1*(dosevec - med_dose)

#doing matrix multiplication to get time-series for each individual
#for that to work, the timevec vector needs to be transposed

m1_mu <- exp(m1_alpha) %*% t(log(timevec)) - exp(m1_beta) %*% t(timevec)

# apply variation following a normal distribution to each data point

m1_y <- rnorm(length(m1_mu),m1_mu, sigma)

# in a final step, we reorganize the data into a long data frame with
# columns id, time, dose, model,
# the deterministic mean mu, and the normally distributed outcome.
# We store dose in 3 versions, the original (log transformed one),
# the one that has the middle value subtracted, and a categorical one.
# Note that trick using sort to get time in the right order.
# Not a robust way of doing things, but works here

m1_dat <- 
        data.frame(
                id = rep(1:Ntot,length(timevec)),
                dose = rep(dosevec,length(timevec)),
                dose_adj = rep(dosevec,length(timevec)) - med_dose,
                dose_cat =  rep(dosevec_cat,length(timevec)),
                time = sort(rep(timevec,Ntot)),
                mu = as.vector(m1_mu),
                outcome = as.vector(m1_y),
                model = "m1"
)

### Model 2 ---------------------------------------------------------------

m2_alpha <- m2_a0 + a1*(dosevec - med_dose)

m2_beta <- m2_b0 + b1*(dosevec - med_dose)

m2_mu <-  exp(m2_alpha) %*% t(log(timevec)) - exp(m2_beta) %*% t(timevec)

m2_y <- rnorm(length(m2_mu),m2_mu, sigma)

m2_dat <- 
        data.frame(
                id = rep(1:Ntot,length(timevec)),
                dose = rep(dosevec,length(timevec)),
                dose_adj = rep(dosevec,length(timevec)) - med_dose,
                dose_cat =  rep(dosevec_cat,length(timevec)),
                time = sort(rep(timevec,Ntot)),
                mu = as.vector(m2_mu),
                outcome = as.vector(m2_y),
                model = "m2"
)

#model 3

m3_alpha <- m3_a0 + a1*(dosevec - med_dose)

m3_beta <- m3_b0 + b1*(dosevec - med_dose)

m3_mu <-  exp(m3_alpha) %*% t(log(timevec)) - exp(m3_beta) %*% t(timevec)

m3_y <- rnorm(length(m3_mu),m3_mu, sigma)

m3_dat <- 
        data.frame(
                id = rep(1:Ntot,length(timevec)),
                dose = rep(dosevec,length(timevec)),
                dose_adj = rep(dosevec,length(timevec)) - med_dose,
                dose_cat =  rep(dosevec_cat,length(timevec)),
                time = sort(rep(timevec,Ntot)),
                mu = as.vector(m3_mu),
                outcome = as.vector(m3_y),
                model = "m3"
)

## PLOT SIMULATED DATA ----------------------------------------------------

p1 <- ggplot(m1_dat) +
        geom_line(aes(x = time, y = mu, col = dose_cat, group = id)) +
        geom_point(aes(x = time, y = outcome, col = dose_cat)) +
        scale_y_continuous(limits = c(-30,200)) +
        labs(y = "Outcome (log virus load)",  x = "Time (days)") +
        theme_minimal()


p2 <- ggplot(m2_dat) +
        geom_line(aes(x = time, y = mu, col = dose_cat, group = id)) +
        geom_point(aes(x = time, y = outcome, col = dose_cat)) +
        scale_y_continuous(limits = c(-30,50)) +
        labs(y = "Outcome (log virus load)",  x = "Time (days)") +
        theme_minimal()

p3 <- ggplot(m3_dat) +
        geom_line(aes(x = time, y = mu, col = dose_cat, group = id)) +
        geom_point(aes(x = time, y = outcome, col = dose_cat)) +
        scale_y_continuous(limits = c(-30,50)) +
        labs(y = "Outcome (log virus load)",  x = "Time (days)") +
        theme_minimal()


## Combine simulated data -------------------------------------------------

#save a plot so we can use it in the blog post

simdat <- 
        list(m1 = m1_dat, 
             m2 = m2_dat, 
             m3 = m3_dat, 
             m1pars = m1pars, 
             m2pars = m2pars, 
             m3pars = m3pars, 
             Nlow = Nlow, 
             Nmed = Nmed, 
             Nhigh = Nhigh
)

saveRDS(simdat, file = filename)

ggsave(file = paste0("featured.png"), p3, dpi = 300, units = "in", width = 6, height = 6)

# MODELING WITH RETHINC ----------------------------------------------------

## Load data --------------------------------------------------------------

simdat <- readRDS("simdat.Rds")

#using dataset 3 for fitting
#also removing anything in the dataframe that's not used for fitting
#makes the ulam/Stan code more robust

fitdat <-
        list(id = simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            dose_adj = simdat[[3]]$dose_adj,
            time = simdat[[3]]$time)

#pulling out number of observations

Ntot = length(unique(simdat$m3$id))

## Model 1 ----------------------------------------------------------------

#wide-prior, no-pooling model

#separate intercept for each individual/id

#2x(N+1)+1 parameters

m1 <- alist(
        # distribution of outcome
        outcome ~ dnorm(mu, sigma),
        
        # main equation for time-series trajectory
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        
        #equations for alpha and beta
        alpha <-  a0[id] + a1*dose_adj,
        beta <-  b0[id] + b1*dose_adj,
        
        #priors
        a0[id] ~ dnorm(2,  10),
        b0[id] ~ dnorm(0.5, 10),
        
        a1 ~ dnorm(0.3, 1),
        b1 ~ dnorm(-0.3, 1),
        sigma ~ cauchy(0,1)
        
)


## Model 2 ----------------------------------------------------------------

#narrow-prior, full-pooling model
#2x(N+2)+1 parameters

m2 <- alist(
        outcome ~ dnorm(mu, sigma),
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        alpha <-  a0[id] + a1*dose_adj,
        beta <-  b0[id] + b1*dose_adj,
        a0[id] ~ dnorm(mu_a,  0.0001),
        b0[id] ~ dnorm(mu_b, 0.0001),
        mu_a ~ dnorm(2, 1),
        mu_b ~ dnorm(0.5, 1),
        a1 ~ dnorm(0.3, 1),
        b1 ~ dnorm(-0.3, 1),
        sigma ~ cauchy(0,1)
)

## Model 3 ----------------------------------------------------------------

#regularizing prior, partial-pooling model

m3 <- alist(
        outcome ~ dnorm(mu, sigma),
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        alpha <-  a0[id] + a1*dose_adj,
        beta <-  b0[id] + b1*dose_adj,
        a0[id] ~ dnorm(2,  1),
        b0[id] ~ dnorm(0.5, 1),
        a1 ~ dnorm(0.3, 1),
        b1 ~ dnorm(-0.3, 1),
        sigma ~ cauchy(0,1)
)

## Model 4 ---------------------------------------------------------------

#adaptive priors, partial-pooling model
#2x(N+2)+1 parameters

m4 <- alist(
        outcome ~ dnorm(mu, sigma),
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        alpha <-  a0[id] + a1*dose_adj,
        beta <-  b0[id] + b1*dose_adj,
        a0[id] ~ dnorm(mu_a,  sigma_a),
        b0[id] ~ dnorm(mu_b, sigma_b),
        mu_a ~ dnorm(2, 1),
        mu_b ~ dnorm(0.5, 1),
        sigma_a ~ cauchy(0, 1),
        sigma_b ~ cauchy(0, 1),
        a1 ~ dnorm(0.3, 1),
        b1 ~ dnorm(-0.3, 1),
        sigma ~ cauchy(0, 1)
)

## Model 2a ---------------------------------------------------------------

#full-pooling model, population-level parameters only
#2+2+1 parameters

m2a <- alist(
        outcome ~ dnorm(mu, sigma),
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        alpha <-  a0 + a1*dose_adj,
        beta <-  b0 + b1*dose_adj,
        a0 ~ dnorm(2,  0.1),
        b0 ~ dnorm(0.5, 0.1),
        a1 ~ dnorm(0.3, 1),
        b1 ~ dnorm(-0.3, 1),
        sigma ~ cauchy(0,1)
)

## Model 4a ---------------------------------------------------------------

#adaptive priors, partial-pooling model
#non-centered

m4a <- alist(
        outcome ~ dnorm(mu, sigma),
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        #rewritten to non-centered
        alpha <-  mu_a + az[id]*sigma_a + a1*dose_adj,
        beta  <-  mu_b + bz[id]*sigma_b + b1*dose_adj,
        #rewritten to non-centered
        az[id] ~ dnorm(0, 1),
        bz[id] ~ dnorm(0, 1),
        mu_a ~ dnorm(2, 1),
        mu_b ~ dnorm(0.5, 1),
        sigma_a ~ cauchy(0, 1),
        sigma_b ~ cauchy(0, 1),
        a1 ~ dnorm(0.3, 1),
        b1 ~ dnorm(-0.3, 1),
        sigma ~ cauchy(0, 1)
        
)

## Model 5 ----------------------------------------------------------------

#no dose effect
#separate intercept for each individual/id
#2xN+1 parameters

m5 <- alist(
        # distribution of outcome
        outcome ~ dnorm(mu, sigma),
        
        # main equation for time-series trajectory
        mu <- exp(alpha)*log(time) - exp(beta)*time,
        
        #equations for alpha and beta
        alpha <-  a0[id],
        beta <-  b0[id],
        
        #priors
        a0[id] ~ dnorm(2,  10),
        b0[id] ~ dnorm(0.5, 10),
        
        sigma ~ cauchy(0,1)
)

## Starting value --------------------------------------------------------

#starting values for model 1
startm1 <- list(a0 = rep(2,Ntot), b0 = rep(0.5,Ntot), a1 = 0.3 , b1 = -0.3, sigma = 1)

#starting values for model 2
startm2 <- list(a0 = rep(2,Ntot), b0 = rep(0.5,Ntot), mu_a = 2, mu_b = 1, a1 = 0.3 , b1 = -0.3, sigma = 1)

#starting values for model 3
startm3 <- startm1

#starting values for models 4 and 4a
startm4 <- list(mu_a = 2, sigma_a = 1, mu_b = 1, sigma_b = 1, a1 = 0.3 , b1 = -0.3, sigma = 1)
startm4a <- startm4

#starting values for model 2a
startm2a <- list(a0 = 2, b0 = 0.5, a1 = 0.3, b1 = -0.3, sigma = 1)

#starting values for model 5
startm5 <- list(a0 = rep(2,Ntot), b0 = rep(0.5,Ntot), sigma = 1)

#put different starting values in list
#need to be in same order as models below
startlist <- list(startm1,startm2,startm3,startm4,startm2a,startm4,startm5)

## Model fitting ----------------------------------------------------------

#general settings for fitting
#you might want to adjust based on your computer

warmup = 6000

iter = warmup + floor(warmup/2)

max_td = 18 #tree depth

adapt_delta = 0.9999

chains = 4

cores  = 4

seed = 4321

#stick all models into a list
modellist = list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, m2a = m2a, m4a = m4a, m5 = m5)

# set up a list in which we'll store our results
fl = vector(mode = "list", length = length(modellist))

#setting for parameter constraints
constraints = list(sigma = "lower=0", sigma_a = "lower=0", sigma_b = "lower=0")

# fitting models
#loop over all models and fit them using ulam
for (n in 1:length(modellist))
{
        
        cat('************** \n')
        cat('starting model', names(modellist[n]), '\n')
        
        tstart = proc.time(); #capture current time
        
        #run model fit
        fl[[n]]$fit <- ulam(flist = modellist[[n]],
                            data = fitdat,
                            start = startlist[[n]],
                            constraints = constraints,
                            log_lik = TRUE, 
                            cmdstan = TRUE,
                            control = list(adapt_delta=adapt_delta,
                                         max_treedepth = max_td),
                            chains = chains, cores = cores,
                            warmup = warmup, iter = iter,
                            seed = seed
        )# end ulam
        
        #capture time taken for fit
        tdiff = proc.time()-tstart;
        
        runtime_minutes = tdiff[[3]]/60;
        
        cat('model fit took this many minutes:', runtime_minutes, '\n')
        cat('************** \n')
        
        #add some more things to the fit object
        fl[[n]]$runtime = runtime_minutes
        fl[[n]]$model = names(modellist)[n]
        
} #end fitting of all models

# saving the list of results so we can use them later
# the file is too large for GitHub
# thus I am saving here to a local folder
# adjust accordingly for your setup
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits", ext="Rds")
saveRDS(fl,filepath)