# example script to compile and sample stan models for devil ray life history analyses
# for length-weight regression, growth models, and length/age-maturity logistic regression
# all stan models used available in 'stan_models' directory
# all stan model fits available in 'model_fits' directory

# load required packages
library(tidyverse)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# read in data
# data for Spinetail devil ray (Mobula mobular) and Bentfin devil ray (M. thurstoni) from three locations (Indonesia, Pakistan, and Kenya)
d <- read.csv("data/devil-ray-data.csv", header = TRUE, stringsAsFactors = FALSE)

# Length-weight regression------------------------------
#subset data to remove rows missing weight (kg) and length (cm)
d1 <- subset(d, !is.na(Wgt)&!is.na(DW))
d1 <- subset(d1, Sp == "M. mobular") # subset for each species
d1 <- d1[, c("DW", "Wgt")] # select rows
rownames(d1) <- NULL # reset row names 

# integer for DW and Wgt
d1$DW <- ceiling(d1$DW) 
d1$Wgt <- ceiling(d1$Wgt) 
  
# data passed to Stan needs to be a list of named objects
DW_wgt_dat <- with(d1, list(DW = DW, Wgt = Wgt, N = length(DW)))

# Compiling Stan model
m <- stan_model("stan_models/DW-wgt-str.stan") # stronger priors

# Sampling Stan model i.e. drawing from the posterior distribution
# same model with same priors but using data for each species
fit.LW.str.st <- sampling(object = m, data=DW_wgt_dat, iter = 5000, chains = 3) 

# Model output
print(fit.LW.str.st)

# Growth models---------------------
#subset data to remove rows missing age and length
d2 <- subset(d, !is.na(Age)&!is.na(DW))
d2 <- subset(d2, Sp == "M. mobular") #subset for each species
d2 <- d2[, c("DW", "Age")] # select row
d2$DW <- d2$DW * 10 # convert DW in cm to mm
d2$Age <- ceiling(d2$Age) # model specifies an integer for DW and Age

# turning data frame into list for use by Stan
growth_dat <- with(d2, list(DW = DW, age = Age, N = length(DW)))

# compile stan model
m2 <- stan_model("stan_models/VGBM-ST-str.stan")

# Sampling Stan model
fit.v.str.st <- sampling(object = m2, data=growth_dat, iter = 5000, chains = 3)

# Model output
print(fit.v.str.st)

# Length-maturity logistic regression--------------------
d4 <- subset(d, !is.na(Mat)&!is.na(DW)&Sp == "M. mobular") # for length at maturity
d5 <- subset(d, !is.na(Mat)&!is.na(Age)&Sp == "M. mobular") # for age-at-maturity
d4$DW <- ceiling(d4$DW) # whole integer for stan model
d5$Age <- ceiling(d5$Age) # whole integer for stan model

#subset to select rows
d4 <- d4[, c("DW", "Mat")] #subset to select rows 
d5 <- d5[, c("Age", "Mat")] #subset to select rows 

# Maturity status needs to be binary
d4$Mat [d4$Mat  == "N"] <- 0
d4$Mat [d4$Mat  == "Y"] <- 1
d4$Mat <- as.numeric(d4$Mat)
rownames(d4) <- NULL #reset row names 

d5$Mat [d5$Mat  == "N"] <- 0
d5$Mat [d5$Mat  == "Y"] <- 1
d5$Mat <- as.numeric(d5$Mat)
rownames(d5) <- NULL #reset row names

# set number of predictions and for which DWs
# add in known Lmat
#M. mobular
DW_mat_dat <- with(d4, list(DW = DW, mat = Mat, N = length(DW), Lmat = 225,
                            K = 100, DW_pred = seq(from = 0, to = 300, length.out = 100)))

# compile stan model
m3 <- stan_model("stan_models/DW-mat-str.stan")

# Sampling Stan model
# same models with same priors but fit to different datasets
fit.LM.str.st <- sampling(object = m3, data=DW_mat_dat, iter = 5000, chains = 3)

# Model output
print(fit.LM.str.st)

# set number of predictions and for which  Ages
# add in known Amat (only known for M. mobular) so same for both species (with different datasets)
Age_mat_dat <- with(d5, list(Age = Age, mat = Mat, N = length(Age), Amat = 5.5,
                             K = 30, Age_pred = seq(from = 0, to = 30, length.out = 30)))

# compile stan model
m4 <- stan_model("stan_models/Age-mat-str.stan")

# Sampling Stan model
# same models with same priors but fit to different datasets
fit.AM.str.st <- sampling(object = m4, data=Age_mat_dat, iter = 5000, chains = 3)

# Model output
print(fit.AM.str.st)
