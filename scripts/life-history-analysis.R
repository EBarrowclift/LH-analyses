# script to analyse devil ray length-at-age dataset

# load required packages
library(tidyverse)
library(cowplot)
library(ggmcmc) # for ggs function for plotting diagnostics to transform stan output into longformat tibble for plotting
library(loo) # to compare models using looic (equivalent to aic)
library(FSA) # has chapmanRobson function
library(rstan)
library(triangle) # for triangular distribution

# load data--------------------
# data for Spinetail devil ray (Mobula mobular) and Bentfin devil ray (M. thurstoni) from three locations (Indonesia, Pakistan, and Kenya)
d <- read.csv("data/devil-ray-data.csv", header = TRUE, stringsAsFactors = FALSE)
# add in common name
d$Sp2 <- d$Sp
d[d$Sp2 == "M. mobular", "Sp2"] <- "Spinetail devil ray" 
d[d$Sp2 == "M. thurstoni", "Sp2"] <- "Bentfin devil ray"

# set themes for plots
set_theme <- theme_classic() +
  theme(legend.title=element_blank(), 
        legend.position=c(0.01,0.99),
        legend.justification=c(0,1), 
        text=element_text(size=15), #change text font size
        legend.text=element_text(face="italic"))

set_theme2 <- theme_classic() +
  theme(legend.title=element_blank(), 
        legend.position=c(0.75,0.99),
        legend.justification=c(0,1), 
        text=element_text(size=15), #change text font size
        legend.text=element_text(face="italic"))

set_theme3 <- theme_classic() +
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.42),
        legend.justification=c(0,1), 
        text=element_text(size=15), #change text font size
        legend.text=element_text(face="italic"))

set_theme4 <- theme_classic() +
  theme(legend.title=element_blank(), 
        legend.position=c(0.85,0.2),
        legend.justification=c(0,1), 
        text=element_text(size=15), #change text font size
        legend.text=element_text(face="italic"))

# Disc Width-Weight Relationship----------------
# subset data to remove rows missing weight and length
# DW in cm and weight in kg
d1 <- subset(d, !is.na(Wgt)&!is.na(DW)) # remove rows without weight data
rownames(d1) <- NULL # reset row names

d1$lnWgt <- log(d1$Wgt) #find natural log of weight and length
d1$lnDW <- log(d1$DW)

# loads 4 LW regression stanfit objects for both species
load("model_fits/stanfits-LW-reg.hyp") 

# Figure 3a M. mobular LW plot------------
# extract posterior estimates (loga and b) for model fit
fit1_trans <- ggs(fit.LW.str.st) %>% # the ggs function transforms stan output into a longformat tibble, that we can use to make different types of plots
  # inc_warmup = FALSE by default in ggs function so does not include warmup samples
  filter(Parameter %in% c("loga", "b")) %>% 
  as.data.frame
ests <- spread(fit1_trans, Parameter, value)
ests <- ests[,3:4]

DW2plot <- 50:400 #set log DW
LCI = UCI = Fit <- numeric(length(DW2plot))
for (i in 1:length(DW2plot)) {
  pv = ests$loga + (ests$b * log(DW2plot[i])) # log(Wi) = log(a) + b*log(L) + Ei
  #pv = exp(ests$loga) * DW2plot^ests$b
  LCI[i] <- quantile(pv,0.025)
  UCI[i] <- quantile(pv,0.975)
  Fit[i] <- quantile(pv,0.5) # this is our log(Wt)
}
CI <- data.frame(Fit, LCI, UCI) # credible intervals to plot
#CI$DW <- DW2plot
CI$lnDW <- log(DW2plot)
head(CI)

# plot
p2 <- ggplot(CI, aes(x=lnDW,y=Fit)) +
  geom_smooth(method=loess, se=FALSE, colour="black") +
  geom_ribbon(aes(ymin=LCI,ymax=UCI), colour="black", fill="NA", linewidth=0.5) +
  geom_point(data=subset(d1, Sp == "M. mobular"), 
             aes(x=lnDW, y=lnWgt, colour=Sex),
             size =2) +
  scale_colour_manual(values=c("#FF9900", "#009933")) +
  # if wanted to label in cm and kg but because natural log scale, looks a bit weird
  # scale_x_continuous(breaks=log(seq(5000, 30000, by = 5000)),  # Breaks in actual cm values
  #   labels=seq(50, 300, by = 50),  # Labels in cm
  #  expand=c(0.01, 0.01)) +
  scale_x_continuous(breaks=seq(4.5, 6, by=0.5), expand=c(0.01, 0.01)) +
  scale_y_continuous(breaks=seq(1.5, 6, by=0.5), expand=c(0.01, 0.01)) +
  coord_cartesian(ylim=c(1.5, 6), xlim=c(4.5, 6)) +
  xlab("ln(disc width, cm)") +
  ylab("ln(weight, kg)") +
  facet_grid(. ~ Sp2) + # so have species label on each plot
  set_theme
p2

p2 <- p2 + annotate("text", x = 5.15, y = 5.5, size = 4, 
                    label = "ln(weight) = 2.51 x ln(disc width) - 9.26",
                    fontface = 'italic') # to make entire text italic
p2

# Figure 3b M. thurstoni LW plot------------
# extract posterior estimates (loga and b) for model fit
fit1_trans <- ggs(fit.LW.str.bf) %>% # the ggs function transforms stan output into a longformat tibble, that we can use to make different types of plots
  # inc_warmup = FALSE by default in ggs function so does not include warmup samples
  filter(Parameter %in% c("loga", "b")) %>% 
  as.data.frame
ests <- spread(fit1_trans, Parameter, value)
ests <- ests[,3:4]

DW2plot <- 50:400 #set log DW
LCI = UCI = Fit <- numeric(length(DW2plot))
for (i in 1:length(DW2plot)) {
  pv = ests$loga + (ests$b * log(DW2plot[i])) # log(Wi) = log(a) + b*log(L) + Ei
  #pv = exp(ests$loga) * DW2plot^ests$b
  LCI[i] <- quantile(pv,0.025)
  UCI[i] <- quantile(pv,0.975)
  Fit[i] <- quantile(pv,0.5) # this is our log(Wt)
}
CI <- data.frame(Fit, LCI, UCI) # credible intervals to plot
#CI$DW <- DW2plot
CI$lnDW <- log(DW2plot)
head(CI)

p3 <- ggplot(CI, aes(x=lnDW,y=Fit)) +
  geom_smooth(method=loess, se=FALSE, colour="black") +
  geom_ribbon(aes(ymin=LCI,ymax=UCI), colour="black", fill="NA", linewidth=0.5) +
  geom_point(data=subset(d1, Sp == "M. thurstoni"), 
             aes(x=lnDW, y=lnWgt, colour=Sex),
             size =2) +
  scale_colour_manual(values=c("#FF9900", "#009933", "grey33")) + # NA doesn't show anyway
  scale_x_continuous(breaks=seq(4, 5.5,by=0.5), expand=c(0.01,0.01)) + 
  scale_y_continuous(breaks=seq(0,5, by=0.5), expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(0,5), xlim=c(4,5.5)) +
  xlab("ln(disc width, cm)") +
  ylab("ln(weight, kg)") +
  facet_grid(. ~ Sp2) +
  set_theme
p3

p3 <- p3 + annotate("text", x = 4.65, y = 4.5, size = 4,
                    label = "ln(weight) = 2.84 x ln(disc width) - 10.74",
                    fontface = 'italic')
p3

# combine plots
plot_grid(p2, p3, ncol=1,labels=c('a','b'), label_size = 15)

# Growth Models---------------------
# loads VGBM, Lester, Gompertz, Logistic growth model stanfit objects for each species
load("model_fits/stanfits-growth-mod-ST.hyp") # M. mobular models
load("model_fits/stanfits-growth-mod-BF.hyp") # M. thurstoni models

# Stan does not automatically compute and store the log-likelihood. 
# has to be incorporated into the Stan program if it is to be extracted after fitting the model.
# calculate loo for M. mobular models
loo1 <- loo(fit.v.str.st, pars = "log_lik") # computes loo for the model fit
loo2 <- loo(fit.v.wk.st, pars = "log_lik")
loo3 <- loo(fit.g.str.st, pars = "log_lik")
loo4 <- loo(fit.g.wk.st, pars = "log_lik")
loo5 <- loo(fit.l.str.st, pars = "log_lik")
loo6 <- loo(fit.l.wk.st, pars = "log_lik")

model_list1 <- list(loo1, loo2, loo3, loo4, loo5, loo6)
loo_compare(model_list1) 

# calculate loo for M. thurstoni models
loo7 <- loo(fit.v.str.bf, pars = "log_lik") # computes loo for the model fit
loo8 <- loo(fit.v.wk.bf, pars = "log_lik")
loo9 <- loo(fit.g.str.bf, pars = "log_lik")
loo10 <- loo(fit.g.wk.bf, pars = "log_lik")
loo11 <- loo(fit.l.str.bf, pars = "log_lik")
loo12 <- loo(fit.l.wk.bf, pars = "log_lik")

# Warning message:
# Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
# pareto_k_table(loo9) # says pareto values are ok

model_list2 <- list(loo7, loo8, loo9, loo10, loo11, loo12)
loo_compare(model_list2) 

# extract looic
#loo1 # can see output
#loo1$estimates[3] # gives looic estimate
#loo1$estimates[6] # gives looic sd estimate

# stacking weights
wts1 <- loo::loo_model_weights(model_list1, method = "stacking")
print(wts1)
wts1 <- as.vector(wts1)
wts2 <- loo::loo_model_weights(model_list2, method = "stacking")
print(wts2)
wts2 <- as.vector(wts2)

# create table of looic and weights-----------
looics1 <- sapply(model_list1, function (x) x$estimates[3]) # extract looics from models
looics_sd1 <- sapply(model_list1, function (x) x$estimates[6]) # extract looics sds from models
looics2 <- sapply(model_list2, function (x) x$estimates[3]) # extract looics from models
looics_sd2 <- sapply(model_list2, function (x) x$estimates[6]) # extract looics sds from models

LOOIC_table1 <- data.frame(Species = rep(c("M. mobular"), 6),
                           Model = rep(c("VGBM", "Gompertz", "Logistic"), each = 2),
                           Prior = rep(c("Strong", "Weak"), 3),
                           LOOIC = round(looics1, 2), # round to 2dp
                           LOOIC_SD = round(looics_sd1, 2),
                           Weights = signif(wts1, 3)) # round to 3dp 

LOOIC_table2 <- data.frame(Species = rep(c("M. thurstoni"), 6),
                           Model = rep(c("VGBM", "Gompertz", "Logistic"), each = 2),
                           Prior = rep(c("Strong", "Weak"), 3),
                           LOOIC = round(looics2, 2), # round to 2dp
                           LOOIC_SD = round(looics_sd2, 2),
                           Weights = signif(wts2, 3)) # round to 3 sig figs 

# create table of mean parameter estimates (95% credible intervals)-------------
# M. mobular
# function to extract parameters for each model
# Define a function to perform the common transformation
extract_fit <- function(data, model, prior) {
  ggs(data) %>%
    filter(Parameter %in% c("Linf", "k", "L0", "sigma")) %>%
    as.data.frame() %>%
    group_by(Parameter) %>%
    summarise(
      Model = model,
      Prior = prior,
      mean = round(mean(value), 2),
      LWR = round(quantile(value, probs = 0.025), 2),
      UPR = round(quantile(value, probs = 0.975), 2),
      est = paste0(mean, " (", LWR, " , ", UPR, ")")
    )
}

# Apply the function to each model for M. mobular
fit1_trans <- extract_fit(fit.v.str.st, "VBGM", "Strong")
fit2_trans <- extract_fit(fit.v.wk.st, "VBGM", "Weaker")
fit3_trans <- extract_fit(fit.g.str.st, "Gompertz", "Strong")
fit4_trans <- extract_fit(fit.g.wk.st, "Gompertz", "Weaker")
fit5_trans <- extract_fit(fit.l.str.st, "Logistic", "Strong")
fit6_trans <- extract_fit(fit.l.wk.st, "Logistic", "Weaker")

# combine tables
Param_ests_table <- rbind(fit1_trans, fit2_trans, fit3_trans, fit4_trans, fit5_trans, fit6_trans)

# convert to long format
Param_ests_table <-  Param_ests_table[, c(1:3, 7)] # select columns
Param_ests_table <- spread(Param_ests_table, Prior, est)

# Apply the function to each model for M. thurstoni
fit1_trans <- extract_fit(fit.v.str.bf, "VBGM", "Strong")
fit2_trans <- extract_fit(fit.v.wk.bf, "VBGM", "Weaker")
fit3_trans <- extract_fit(fit.g.str.bf, "Gompertz", "Strong")
fit4_trans <- extract_fit(fit.g.wk.bf, "Gompertz", "Weaker")
fit5_trans <- extract_fit(fit.l.str.bf, "Logistic", "Strong")
fit6_trans <- extract_fit(fit.l.wk.bf, "Logistic", "Weaker")

# combine tables
Param_ests_table2 <- rbind(fit1_trans, fit2_trans, fit3_trans, fit4_trans, fit5_trans, fit6_trans)

# convert to long format
Param_ests_table2 <-  Param_ests_table2[, c(1:3, 7)] # select columns
Param_ests_table2 <- spread(Param_ests_table2, Prior, est)

# predict growth model parameters M. mobular--------------------
# Function to extract and reshape parameter estimates for each model
extract_params <- function(data) {
  ggs(data) %>%
    filter(Parameter %in% c("Linf", "k", "L0")) %>%
    as.data.frame() %>%
    spread(Parameter, value) %>%
    select(Linf, k, L0) # Select only the necessary columns
}

# Extract parameters for each M. mobular model
ests_1 <- extract_params(fit.v.str.st)
ests_2 <- extract_params(fit.v.wk.st)
ests_3 <- extract_params(fit.g.str.st)
ests_4 <- extract_params(fit.g.wk.st)
ests_5 <- extract_params(fit.l.str.st)
ests_6 <- extract_params(fit.l.wk.st)

# Function to generate model fits
compute_fit <- function(ages, params, formula_func) {
  sapply(ages, function(age) {
    pv <- formula_func(params, age)
    pv <- exp(pv) # Transform back from log scale
    quantile(pv, 0.5) # Return median (0.5 quantile)
  })
}

# Define each growth model formula
vbgm_formula <- function(params, age) {
  log(params$Linf * (1 - ((1 - (params$L0 / params$Linf)) * exp(-params$k * age))))
}

gompertz_formula <- function(params, age) {
  log(params$L0 * exp(log(params$Linf / params$L0) * (1 - exp(-params$k * age))))
}

logistic_formula <- function(params, age) {
  log((params$Linf * params$L0 * exp(params$k * age)) / 
        (params$Linf + params$L0 * (exp(params$k * age) - 1)))
}

# set ages to plot
ages2plot <- 0:30 #set appropriate ages 
LCI = UCI = Fit1 = Fit2 = Fit3 = Fit4 = Fit5 = Fit6  <- numeric(length(ages2plot)) # create empty vectors for data

# Compute fits for each model
Fit1 <- compute_fit(ages2plot, ests_1, vbgm_formula)
LCI <- sapply(ages2plot, function(age) quantile(exp(vbgm_formula(ests_1, age)), 0.025))
UCI <- sapply(ages2plot, function(age) quantile(exp(vbgm_formula(ests_1, age)), 0.975))
Fit2 <- compute_fit(ages2plot, ests_2, vbgm_formula)
Fit3 <- compute_fit(ages2plot, ests_3, gompertz_formula)
Fit4 <- compute_fit(ages2plot, ests_4, gompertz_formula)
Fit5 <- compute_fit(ages2plot, ests_5, logistic_formula)
Fit6 <- compute_fit(ages2plot, ests_6, logistic_formula)

# Combine the results into a data frame
Spinetail_fit <- data.frame(Fit1, LCI, UCI, Fit2, Fit3, Fit4, Fit5, Fit6)
Spinetail_fit$age <- 0:(nrow(Spinetail_fit) - 1)

# Display the first few rows of the result
head(Spinetail_fit)

# predict growth model parameters M. thurstoni--------------------
# Extract parameters for each M. thurstoni model
ests_1 <- extract_params(fit.v.str.bf)
ests_2 <- extract_params(fit.v.wk.bf)
ests_3 <- extract_params(fit.g.str.bf)
ests_4 <- extract_params(fit.g.wk.bf)
ests_5 <- extract_params(fit.l.str.bf)
ests_6 <- extract_params(fit.l.wk.bf)

# set ages to plot
ages2plot <- 0:30 #set appropriate ages 
LCI = UCI = Fit1 = Fit2 = Fit3 = Fit4 = Fit5 = Fit6  <- numeric(length(ages2plot)) # create empty vectors for data

# Compute fits for each model
Fit1 <- compute_fit(ages2plot, ests_1, vbgm_formula)
Fit2 <- compute_fit(ages2plot, ests_2, vbgm_formula)
Fit3 <- compute_fit(ages2plot, ests_3, gompertz_formula)
Fit4 <- compute_fit(ages2plot, ests_4, gompertz_formula)
Fit5 <- compute_fit(ages2plot, ests_5, logistic_formula)
LCI <- sapply(ages2plot, function(age) quantile(exp(logistic_formula(ests_5, age)), 0.025))
UCI <- sapply(ages2plot, function(age) quantile(exp(logistic_formula(ests_5, age)), 0.975))
Fit6 <- compute_fit(ages2plot, ests_6, logistic_formula)

Bentfin_fit <- data.frame(Fit1, Fit2, Fit3, Fit4, Fit5, Fit6, LCI, UCI)
Bentfin_fit$age <- 1:nrow(Bentfin_fit)-2
head(Bentfin_fit)

# Plot top growth model for each species--------------
d2 <- subset(d, !is.na(Age)&!is.na(DW)) # subset for age data
d2$DW <- d2$DW * 10 # convert DW in cm to mm # makes whole number
d2$Age <- ceiling(d2$Age) # plot in whole numbers

# spinetail asymtoptic DW for top model = 3509
# bentfin asymtoptic DW for top model = 2022

# create df with max size info for each species
MaxSize <- data.frame(Sp = rep(c("M. mobular", "M. thurstoni"), each =3),
                      Prior = rep(c("TopModel","MaxKnown", "Observed"), 2),
                      MaxDW = c(3509, 3500, 2600, 2180, 1970, 1900)) # IUCN RLA - note, M. thurstoni 197 (possibly 220cm DW)

# Figure 4a: multiple growth models M. mobular--------------------
# edit colours so have dark and light of three colours for strong and weaker priors for each model

p4 <- ggplot(Spinetail_fit, aes(x=age, y=Fit1)) +
  geom_ribbon(data=Spinetail_fit, aes(ymin=predict(loess(LCI~age)), ymax=predict(loess(UCI~age))), colour="black", fill="NA", linewidth=0.5) + # plot CIs for top model
  # note, order affects plotting order so plot top model last
  geom_smooth(data=Spinetail_fit, aes(x=age, y=Fit3),  # fit model 2 (gompertz)
              colour="#330099", method=loess, se=FALSE, linetype = "dashed") + # Dark blue
  geom_smooth(data=Spinetail_fit, aes(x=age, y=Fit4),# Light blue
              colour="#6699FF", method=loess, se=FALSE, linetype = "dashed") + 
  geom_smooth(data=Spinetail_fit, aes(x=age, y=Fit5), # fit model 3 (logistic)
              colour="#990000", method=loess, se=FALSE, linetype = "dashed") + # Dark red
  geom_smooth(data=Spinetail_fit, aes(x=age, y=Fit6), # Light red
              colour="#FF3333", method=loess, se=FALSE, linetype = "dashed") + 
  # fit top model (VGBM) with solid line
  geom_smooth(data=Spinetail_fit, aes(x=age, y=Fit1),# fit model 1 (VGBM)
              colour="#660066", method=loess, se=FALSE) + # Dark purple (strong priors)
  geom_smooth(data=Spinetail_fit, aes(x=age, y=Fit2), # Light purple (weak priors)
              colour="#CC66CC", method=loess, se=FALSE, linetype = "dashed") + 
  geom_point(data=subset(d2, Sp == "M. mobular"), 
             aes(x=Age,y=DW, colour = Sex, shape = Loc), alpha = 0.6) + 
  scale_x_continuous(breaks=seq(0,25,by=2), expand=c(0.01,0.01)) + 
  scale_y_continuous(breaks=seq(800,3700,by=200), expand=c(0.01,0.01),
                     labels = seq(80, 360, by = 20)) + # update y axis labels
  # scale_alpha_manual(values=c(0.3,0.2,0.7,1, 0.7,1)) + # can't seem to change line transparency
  coord_cartesian(ylim=c(800,3700), xlim=c(0,25)) +
  scale_colour_manual(values=c("#FF9900", "#009933")) +
  xlab("Age (years)") +
  ylab("Disc width (cm)") +
  facet_grid(.~ Sp2) +
  set_theme3
p4

# add in dotted lines for asymptotic size (top model), max. known and observed size
p4 <- p4 + geom_hline(data = subset(MaxSize, Sp == "M. mobular"), 
                      aes(yintercept = MaxDW),
                      colour = c("#660066", "#999999", "chocolate4"),
                      linetype = "dotted",
                      linewidth= 1) 
p4

# label lines with annotate function
p4 <- p4 + annotate("text", x = c(3,3,3), y = c(3600, 3400, 2700), label = c("Asymptotic size", "Species maximum size", "Study maximum size"),
                    colour = c("#660066", "#999999", "chocolate4"))
p4

p4 <- p4 + annotate("text", x = c(22,22,10,10,18,18), y = c(2200,2100, 2900,2800, 3100,3000), label = c("VBGM: strong priors", "weaker priors", "Gompertz: strong priors", "weaker priors", "Logistic: strong priors", "weaker priors"),
                    colour = c("#660066", "#CC66CC", "#330099", "#6699FF", "#990000", "#FF3333"))
p4

# # Figure 4b: multiple growth models M. thurstoni--------------------
p5 <- ggplot(Bentfin_fit, aes(x=age, y=Fit1)) +
  geom_ribbon(data=Bentfin_fit, aes(ymin=predict(loess(LCI~age)), ymax=predict(loess(UCI~age))), colour="black", fill="NA", linewidth=0.5) + # plot CIs for top model
  # note, order affects plotting order so plot top model last
  geom_smooth(data=Bentfin_fit, aes(x=age, y=Fit1), # fit model 1 (VGBM) 
              colour="#660066", method=loess, se=FALSE, linetype = "dashed") + # Dark purple (strong priors)
  geom_smooth(data=Bentfin_fit, aes(x=age, y=Fit2), # Light purple (weak priors)
              colour="#CC66CC", method=loess, se=FALSE, linetype = "dashed") + 
  geom_smooth(data=Bentfin_fit, aes(x=age, y=Fit3), # fit model 2 (gompertz)
              colour="#330099", method=loess, se=FALSE, linetype = "dashed") + # Dark blue
  geom_smooth(data=Bentfin_fit, aes(x=age, y=Fit4),# Light blue
              colour="#6699FF", method=loess, se=FALSE, linetype = "dashed") + 
  # fit top model (logistic) with solid line
  geom_smooth(data=Bentfin_fit, aes(x=age, y=Fit5), # fit model 3 (logistic)
              colour="#990000", method=loess, se=FALSE) + # Dark red
  geom_smooth(data=Bentfin_fit, aes(x=age, y=Fit6), # Light red
              colour="#FF3333", method=loess, se=FALSE, linetype = "dashed") + # fit model 3 +
  geom_point(data=subset(d2, Sp == "M. thurstoni"), 
             aes(x=Age,y=DW, colour = Sex), alpha = 0.6) +
  scale_x_continuous(breaks=seq(0,20,by=2), expand=c(0.01,0.01)) + 
  scale_y_continuous(breaks=seq(600,2600,by=200), expand=c(0.01,0.01),
                     labels = seq(60, 260, by = 20)) + # update y axis labels
  coord_cartesian(ylim=c(600,2600), xlim=c(0,20)) +
  scale_colour_manual(values=c("#FF9900", "#009933")) +
  xlab("Age (years)") +
  ylab("Disc width (cm)") +
  facet_grid(. ~ Sp2) +
  set_theme4
p5

# add in dotted lines for asymptotic size (top model), max. known and observed size
p5 <- p5 + geom_hline(data = subset(MaxSize, Sp == "M. thurstoni"), 
                      aes(yintercept = MaxDW),
                      colour = c("#660066", "#999999", "chocolate4"),
                      linetype = "dotted",
                      linewidth= 1) 
p5

# label lines with annotate function
p5 <- p5 + annotate("text", x = c(3,3,3), y = c(2250, 2025, 1850), label = c("Asymptotic size", "Species maximum size", "Study maximum size"),
                    colour = c("#660066", "#999999", "chocolate4"))
p5

p5 <- p5 + annotate("text", x = c(18,18,12,12,12,12), y = c(1400,1300, 1400,1300, 2500,2400), label = c("VBGM: strong priors", "weaker priors", "Gompertz: strong priors", "weaker priors", "Logistic: strong priors", "weaker priors"),
                    colour = c("#660066", "#CC66CC", "#330099", "#6699FF", "#990000", "#FF3333"))
p5

# combine plots
plot_grid(p4, p5, ncol=1,labels=c('a','b'), label_size = 15)

# Predicting length and age at maturity-------------------
# loads 4 length-maturity and 4 age-maturity logistic regression model stanfit objects for both species
load("model_fits/stanfits-LMat-log-reg.hyp")

#subset data to remove rows missing size and maturity
d4 <- subset(d, !is.na(Mat)&!is.na(DW))
d5 <- subset(d, !is.na(Mat)&!is.na(Age)) # if do age and mat
d4$DW <- ceiling(d4$DW) # needs to be whole integer for stan model
d5$Age <- ceiling(d5$Age) # needs to be whole integer for stan model

# Maturity status needs to be binary
d4$Mat[d4$Mat  == "N"] <- 0
d4$Mat[d4$Mat  == "Y"] <- 1
d4$Mat <- as.numeric(d4$Mat)
rownames(d4) <- NULL #reset row names 

d5$Mat[d5$Mat  == "N"] <- 0
d5$Mat[d5$Mat  == "Y"] <- 1
d5$Mat <- as.numeric(d5$Mat)
rownames(d5) <- NULL #reset row names

# generate required data to predict from M. mobular-----------
# plot for strong priors for DW-mat and Age-mat
# Extract posterior samples
posterior_samples1 <- as.data.frame(extract(fit.LM.str.st))
posterior_samples2 <- as.data.frame(extract(fit.AM.str.st))

# Generate a range of DW values for prediction
d4 <- subset(d4, Sp == "M. mobular")
d5 <- subset(d5, Sp == "M. mobular")
DW2Plot <- seq(from = min(d4$DW), to = max(d4$DW), length.out = 100)
Age2Plot <- seq(from = min(d5$Age), to = max(d5$Age), length.out = 100)

# Function to calculate predicted probabilities for a given DW value
predict_probs1 <- function(DW_value) {
  # Calculate log-odds (logit scale)
  log_odds <- posterior_samples1$alpha + posterior_samples1$beta * DW_value
  
  # Transform log-odds to probabilities using the inverse logit function
  probs <- 1 / (1 + exp(-log_odds))
  
  # Return the mean and 95% credible interval of the predicted probabilities
  c(mean_prob = mean(probs), ci_low = quantile(probs, 0.025), ci_high = quantile(probs, 0.975))
}

# Function to calculate predicted probabilities for a given Age value
predict_probs2 <- function(Age_value) {
  # Calculate log-odds (logit scale)
  log_odds <- posterior_samples2$alpha + posterior_samples2$beta * Age_value
  
  # Transform log-odds to probabilities using the inverse logit function
  probs <- 1 / (1 + exp(-log_odds))
  
  # Return the mean and 95% credible interval of the predicted probabilities
  c(mean_prob = mean(probs), ci_low = quantile(probs, 0.025), ci_high = quantile(probs, 0.975))
}

# Apply the function to each DW and agevalue
pred_mod_fit1 <- t(sapply(DW2Plot, predict_probs1))
pred_mod_fit2 <- t(sapply(Age2Plot, predict_probs2))

# Create a data frame for plotting
Spinetail_fit1 <- data.frame(
  DW = DW2Plot,
  mean_prob = pred_mod_fit1[, "mean_prob"],
  ci_low = pred_mod_fit1[, grep("ci_low", colnames(pred_mod_fit1))], # because for some reason above code for quantiles generates slightly different column headers
  ci_high = pred_mod_fit1[, grep("ci_high", colnames(pred_mod_fit1))]
)

Spinetail_fit2 <- data.frame(
  Age = Age2Plot,
  mean_prob = pred_mod_fit2[, "mean_prob"],
  ci_low = pred_mod_fit2[, grep("ci_low", colnames(pred_mod_fit2))], # because for some reason above code for quantiles generates slightly different column headers
  ci_high = pred_mod_fit2[, grep("ci_high", colnames(pred_mod_fit2))]
)

# determine DW-at-50% mat and Age-at-50% mat so can add dashed lines to plot
#add column to dataframe that finds value closest to 0.5
Spinetail_fit1$Mature.5 <- sqrt((0.5-Spinetail_fit1$mean_prob)^2)
Spinetail_fit1$upr.5 <- sqrt((0.5-Spinetail_fit1$ci_high)^2)
Spinetail_fit1$lwr.5 <- sqrt((0.5-Spinetail_fit1$ci_low)^2)

Spinetail_fit2$Mature.5 <- sqrt((0.5-Spinetail_fit2$mean_prob)^2)
Spinetail_fit2$upr.5 <- sqrt((0.5-Spinetail_fit2$ci_high)^2)
Spinetail_fit2$lwr.5 <- sqrt((0.5-Spinetail_fit2$ci_low)^2)

#sort dataframe by Mature value closest to 0.5 and extract DW value for 50% maturity
Spinetail_fit1 <- Spinetail_fit1[order(Spinetail_fit1$Mature.5),]
head(Spinetail_fit1)
DW50 <- Spinetail_fit1[1,1]

Spinetail_fit <- Spinetail_fit2[order(Spinetail_fit2$Mature.5),]
head(Spinetail_fit2)
A50 <- Spinetail_fit2[1,1]

# Calculate length and age at 50% maturity (alternative)
# uniroot function finds the DW value where predicted probability is closest to 0.5.
DW50 <- uniroot(function(x) {
  predict_probs1(x)["mean_prob"] - 0.5
}, interval = range(d4$DW))$root

A50 <- uniroot(function(x) {
  predict_probs2(x)["mean_prob"] - 0.5
}, interval = range(d5$Age))$root

#sort dataframe by upr value closest to 0.5 and extract DW value for first row
Spinetail_fit1 <- Spinetail_fit1[order(Spinetail_fit1$upr.5),]
DW50u <- Spinetail_fit1[1,1]

Spinetail_fit2 <- Spinetail_fit2[order(Spinetail_fit2$upr.5),]
A50u <- Spinetail_fit2[1,1]

DW50u <- uniroot(function(x) {
  predict_probs1(x)[grep("ci_high", colnames(pred_mod_fit1))] - 0.5
}, interval = range(d4$DW))$root

A50u <- uniroot(function(x) {
  predict_probs2(x)[grep("ci_high", colnames(pred_mod_fit2))] - 0.5
}, interval = range(d5$Age))$root

#sort dataframe by lwr value closest to 0.5 and extract DW value for first row
Spinetail_fit1 <- Spinetail_fit1[order(Spinetail_fit1$lwr.5),]
DW50l <- Spinetail_fit1[1,1]

Spinetail_fit2 <- Spinetail_fit2[order(Spinetail_fit2$lwr.5),]
A50l <- Spinetail_fit2[1,1]

DW50l <- uniroot(function(x) {
  predict_probs1(x)[grep("ci_low", colnames(pred_mod_fit1))] - 0.5
}, interval = range(d4$DW))$root

A50l <- uniroot(function(x) {
  predict_probs2(x)[grep("ci_low", colnames(pred_mod_fit2))] - 0.5
}, interval = range(d5$Age))$root

# Fig. 7b DW-Mat curve M. mobular--------------------
p6 <- ggplot(Spinetail_fit1, aes(x = DW, y = mean_prob)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), colour = "black", fill = "NA", linewidth=0.5) + # plots line for CIs
  geom_point(data=subset(d4, Sp == "M. mobular"), 
             aes(x=DW, y=Mat, colour = Sex), alpha = 0.6) +
  #DW at 50% mature lines
  geom_segment(aes(x=0,xend=DW50l,y=0.5,yend=0.5), linetype=2, colour="black", linewidth=0.5) +
  geom_segment(aes(x=DW50,xend=DW50,y=0.5,yend=0), linetype=2, colour="black", linewidth=0.5) +
  geom_segment(aes(x=DW50u,xend=DW50u,y=0.5,yend=0), linetype=2, colour="black", linewidth=0.5) +
  geom_segment(aes(x=DW50l,xend=DW50l,y=0.5,yend=0), linetype=2, colour="black", linewidth=0.5) +
  scale_colour_manual(values=c("#FF9900", "#009933")) + # colour by sex
  scale_y_continuous(breaks=seq(0,1,by=0.5), expand=c(0.01,0.01)) + 
  scale_x_continuous(breaks=seq(60,260,by=20), expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(0,1), xlim=c(60,260)) +
  ylab("Maturity (0=Immature, 1=Mature)")+
  xlab("Disc Width (cm)") +
  set_theme
p6

p6 <- p6 + annotate("text", x = c(120, 120), y = c(0.65, 0.6), 
                    label = c("DW at 50% maturity = 203.7cm", 
                              "(95% CI 191.3cm, 218.8cm)"),
                    fontface = 'italic', # to make entire text italic
                    colour = "black")
p6


# Fig. 7a Age-mat curve for M. mobular-------------------------
p7 <- ggplot(Spinetail_fit2, aes(x = Age, y = mean_prob)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), colour = "black", fill = "NA", linewidth=0.5) + # plots line for CIs
  geom_point(data=subset(d5, Sp == "M. mobular"), 
             aes(x=Age, y=Mat, colour = Sex), alpha = 0.6) +
  #DW at 50% mature lines
  geom_segment(aes(x=0,xend=A50l,y=0.5,yend=0.5), linetype=2, colour="black", linewidth=0.5) +
  geom_segment(aes(x=A50,xend=A50,y=0.5,yend=0), linetype=2, colour="black", linewidth=0.5) +
  geom_segment(aes(x=A50u,xend=A50u,y=0.5,yend=0), linetype=2, colour="black", linewidth=0.5) +
  geom_segment(aes(x=A50l,xend=A50l,y=0.5,yend=0), linetype=2, colour="black", linewidth=0.5) +
  scale_colour_manual(values=c("#FF9900", "#009933")) + # colour by sex
  scale_y_continuous(breaks=seq(0,1,by=0.5), expand=c(0.01,0.01)) + 
  scale_x_continuous(breaks=seq(0,20,by=2), expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,20)) +
  ylab("Maturity (0=Immature, 1=Mature)")+
  xlab("Age (years)") +
  set_theme
p7

p7 <- p7 + annotate("text", x = c(14, 14), y = c(0.5, 0.45), 
                    label = c("Age at 50% maturity = 8.2 yrs", 
                              "(95% CI 6.8 yrs, 10.3 yrs)"),
                    fontface = 'italic', # to make entire text italic
                    colour = "black")
p7

# combine plots
plot_grid(p7, p6, ncol=1,labels=c('a','b'), label_size = 15)

# Predict age given width-------------------
d3 <- subset(d, Sp == "M. mobular") # subset dataset for species
d4 <- subset(d, Sp == "M. thurstoni") # subset dataset for species

# function to predict Age given DW based on the inverse VBGM for M. mobular
# used mean values from posterior distributions
IvbC <- function(DW) {
  Linf <- 350.295 #define values supplied by the previous analyses top model # done in same unit as DW
  K <- 0.05
  L0 <- 116.653 
  
  Age <- (1/K)*log((Linf-L0)/(Linf-DW))
  
  return(Age)
}

#write a function predicting Age given DW based on the inverse logistic for M. thurstoni
IvLC <- function(DW) {
  Linf <- 202.211 #define values supplied by the previous analysestop model
  g <- 0.19
  L0 <- 87.021 
  
  Age <- log((DW * Linf - DW * L0) / (Linf * L0 - DW * L0 )) / g
  
  return(Age)
}

#apply function to all data
d3$Est_Age <- IvbC(d3$DW)
d4$Est_Age <- IvLC(d4$DW)

#all NaN change to NA 
#although none in this case
d3$Est_Age[is.nan(d3$Est_Age)] <- NA
d4$Est_Age[is.nan(d4$Est_Age)] <- NA

#all estimates <0 change to 0
d3$Est_Age[d3$Est_Age <= 0] <- 0
d4$Est_Age[d4$Est_Age <= 0] <- 0

#round all estimates to the closest whole number
d3$Est_Age <- round(d3$Est_Age,digits=0)
d4$Est_Age <- round(d4$Est_Age,digits=0)

#extract animals with estimated ages
# all in this case
d3 <- subset(d3, !is.na(Est_Age))
d4 <- subset(d4, !is.na(Est_Age))

#table of age frequencies - M. mobular
tblAge <- table(d3$Est_Age)
tblAge <- as.data.frame(tblAge)
names(tblAge)[1] <- 'Age'
tblAge$Age <- as.numeric(levels(tblAge$Age)) [tblAge$Age]
tblAge

#table of age frequencies - M. thurstoni
tblAge2 <- table(d4$Est_Age)
tblAge2 <- as.data.frame(tblAge2)
names(tblAge2)[1] <- 'Age'
tblAge2$Age <- as.numeric(levels(tblAge2$Age)) [tblAge2$Age]
tblAge2

# Catch Curve-----------------
# generate catch curve
#plot natural log of frequency
tblAge$logFreq <- log(tblAge$Freq)
plot(logFreq~Age, data=tblAge, ylab="log(Catch)", pch=20)

tblAge2$logFreq <- log(tblAge2$Freq)
plot(logFreq~Age, data=tblAge2, ylab="log(Catch)", pch=20)

#3 appears to be age at full selection to the fishery for M. mobular, 19 is the max age estimated
#Therefore run curve for ages 4-19 (chapman robson approach)
#run catch curve regression
CCurve <- chapmanRobson(Freq~Age, data=tblAge, ages2use=4:19)
#examine outputs for instantenous mortality (Z) estimates, extract mean, CI and SD
summary(CCurve) #mean and SE
confint(CCurve) # UCI and LCI

#2 appears to be age at full selection to the fishery for M. thurstoni, 13 is the max age estimated
#Therefore run curve for ages 3-13 (chapman robson approach)
#run catch curve regression
CCurve2 <- chapmanRobson(Freq~Age, data=tblAge2, ages2use=3:13)
#examine outputs for instantenous mortality (Z) estimates, extract mean, CI and SD
summary(CCurve2) #mean and SE
confint(CCurve2) # UCI and LCI

#then extract Z data as follows:
Z_st <- summary(CCurve)[2,1]
Zse_st <- summary(CCurve)[2,2]
Zlwr_st <- confint(CCurve)[2,1]
Zupr_st <- confint(CCurve)[2,2]

Z_bf <- summary(CCurve2)[2,1]
Zse_bf <- summary(CCurve2)[2,2]
Zlwr_bf <- confint(CCurve2)[2,1]
Zupr_bf <- confint(CCurve2)[2,2]

#estimate Annual mortality (A) from Z (including CIs)
A_st <- 1-exp(-Z_st) ##mean
Alwr_st <- 1-exp(-Zlwr_st) ##95%LCI
Aupr_st <- 1-exp(-Zupr_st) ##95%UCI

A_bf <- 1-exp(-Z_bf) ##mean
Alwr_bf <- 1-exp(-Zlwr_bf) ##95%LCI
Aupr_bf <- 1-exp(-Zupr_bf) ##95%UCI

# Fig. 9 plot the catch curves-----------------
tblAge$Sp <- "Spinetail devil ray"
#plot curve in ggplot for M. mobular
cc1 <- ggplot(data=tblAge, aes(x=Age,y=logFreq)) +
  ## subset smooth to only assess 4:19
  geom_smooth(data=subset(tblAge,Age>=4&Age<=19), method="lm", se=FALSE, colour="black") +
  ##subset data for used
  geom_point(data=subset(tblAge,Age>=4&Age<=19), size=1, shape=21, fill="black") +
  ##subset data for unused
  geom_point(data=subset(tblAge, Age<=3), size=1, shape=21, fill="white") +
  scale_y_continuous(breaks=seq(0,3,by=0.5), expand=c(0.01,0.01)) + 
  scale_x_continuous(breaks=seq(0,19,by=1), expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(0,3), xlim=c(0,20)) +
  ylab("ln(Catch Frequency)")+
  xlab("Age (Years)") +
  facet_grid(. ~ Sp) +
  set_theme
cc1

# add text for parameters
Zlab <- expression(paste(italic("Z = 0.215 "), italic(year^-1), italic(" (95% CI 0.157, 0.272)")))
cc1 <- cc1 + annotate("text",  x=Inf, y=Inf, label=Zlab, hjust=1, vjust=1, fontface = "italic")
cc1  

#plot curve in ggplot for M. mobular
tblAge2$Sp <- "Bentfin devil ray"
cc2 <- ggplot(data=tblAge2, aes(x=Age,y=logFreq)) +
  ## subset smooth to only assess 4:13
  geom_smooth(data=subset(tblAge2,Age>=3&Age<=13), method="lm", se=FALSE, colour="black") +
  ##subset data for used
  geom_point(data=subset(tblAge2,Age>=3&Age<=13), size=1, shape=21, fill="black") +
  ##subset data for unused
  geom_point(data=subset(tblAge2, Age<=2), size=1, shape=21, fill="white") +
  scale_y_continuous(breaks=seq(0,2,by=0.5), expand=c(0.01,0.01)) + 
  scale_x_continuous(breaks=seq(0,13,by=1), expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(0,3), xlim=c(0,13)) +
  ylab("ln(Catch Frequency)")+
  xlab("Age (Years)") +
  facet_grid(. ~ Sp) +
  set_theme
cc2

# add text for parameters
Zlab2 <- expression(paste(italic("Z = 0.232 "), italic(year^-1), italic(" (95% CI 0.136, 0.328)")))
cc2 <- cc2 + annotate("text",  x=Inf, y=Inf, label=Zlab2, hjust=1, vjust=1, fontface = "italic")
cc2  

# combine plots
plot_grid(cc1, cc2, ncol=1,labels=c('a','b'), label_size = 15)

# Estimating Rmax and Exploitation Ratio-------------
# 1. Annual reproductive output (b): we have assumed litters are 1 or 2 pups every 1-3 years, assuming a 1:1 sex ratio and annual reproduction 
# this means an annual rate of 0.17-1 female pup per year for both species
iter <- 10000 # 10,000 iterations sufficient
#reproductive output
b <- runif(iter, 0.17, 1)
mean(b)
hist(b)

# 2. age at 50% maturity-------------------
# is estimated based on the size at 50% maturity of combined sex, assuming males and females share a similar size at maturity. 

# for M. mobular, we directly estimated age-at-mat
# 8.2 years (6.8 , 10.3 years)
amatupr <- 10.3 # A50u
amatlwr <- 6.8 # A50l

#treat as uniform (conservative)
amat_st <- runif(iter, amatlwr, amatupr)
mean(amat_st)
hist(amat_st)

#treat as uniform (conservative)
amat <- runif(iter, amatlwr, amatupr)
mean(amat)
hist(amat)

# For M. thurstoni, using range of 7-10.4 years for amat

#treat as uniform (conservative)
amatupr <- 10.4
amatlwr <- 7
amat_bf <- runif(iter, amatlwr, amatupr)
mean(amat_bf)
hist(amat_bf)

# 3. Maximum Age---------------
# Maximum age can be estimated by age at 95%DW (5LN(2)/K) and 99%DW (7LN(2)/K) for an exponential function (VGBM)
#need to generate a distribution for k in order to predict the distribution of max age

# can calculate max age at 95%DW and 99%DW for M. mobular as VGBM (exponential function) is top growth curve
# gives high ages for spinetail
K_st <- 0.05
amax95_st <- 5*log(2)/K_st
amax95_st
#amax99_st <- 7*log(2)/K_st

# for sigmoid logistic function
# amax can be estimated as 99% of Linf
# (((Linf x 0.99) x (Linf x L0)) / L0 (Linf - (Linf * 0.99))) / g
g_bf <- 0.19
Linf_bf <- 2022.11 # in mm as used in function
L0_bf <- 870.21 # in mm
amax99_bf <- log(((Linf_bf * 0.99) * (Linf_bf * L0_bf)) / (L0_bf * (Linf_bf - (Linf_bf * 0.99)))) / g_bf
amax99_bf

#however the max age for females observed directly was 17.5 and 6 for M. mobular and M. thurstoni so this should be used as the lower bound instead
# given 6 likely too low, used 17.5 for both

# also 95% and 99%DW for M.mobular and M. thurstoni of 69 & 64 years seems unrealistically high
# tested with upper bound of 26 for both
amaxD <- 17.5
#amaxC <- 26
amaxC_st <- 69 # M. mobular
amaxC_bf <- 64 # M. thurstoni

#take these to be the limits of max age estimates and generate a max age distribution assuming a uniform distribution between the two ages
#amax_st <- runif(iter, amaxD, amaxC)# ran with 26 years as upper bound
amax_st <- runif(iter, amaxD, amaxC_st) # runif assumes a uniform distribution
mean(amax_st)
hist(amax_st)

#amax_bf <- runif(iter, amaxD, amaxC)
amax_bf <- runif(iter, amaxD, amaxC_bf) # runif assumes a uniform distribution
mean(amax_bf)
hist(amax_bf)

# 4. get distribution of M----------
M_st <- ((amax_st+amat_st)/2)^-1 # based on distributions of amax and amat
mean(M_st)
median(M_st)
hist(M_st)
Min_st <- quantile(M_st,0.025) # then use 95% quanitles for min and max
Max_st <- quantile(M_st,0.975)
Mode_st <- (Min_st + Max_st) / 2 # and mode for triangular distribution

#get distribution of M
M_bf <- ((amax_bf+amat_bf)/2)^-1
mean(M_bf)
median(M_bf)
hist(M_bf)
Min_bf <- quantile(M_bf,0.025)
Max_bf <- quantile(M_bf,0.975)
Mode_bf <- (Min_bf + Max_bf) / 2

# uniform distribution for M
# M_st <- runif(iter, Min_st, Max_st)
# M_bf <- runif(iter, Min_bf, Max_bf)

# triangular distribution for M (less confidence in M being at the minimum/maximum of the range than the mid-point)
M_st <- rtriangle(iter, Min_st, Max_st, Mode_st)
mean(M_st)
median(M_st)
quantile(M_st,0.025) 
quantile(M_st,0.975)
hist(M_st)
M_bf <- rtriangle(iter, Min_bf, Max_bf, Mode_bf)
mean(M_bf)
median(M_bf)
quantile(M_bf,0.025)
quantile(M_bf,0.975)
hist(M_bf)

# 4. Calculate *Rmax*----------------
# We now need to draw data from these ranges to create a data set for estimating rmax from
# We assume a uniform distribution. Define the equation

# function to calculate rmax based on lifespan estimator using input of M
rmaxnlm <- function(x, mortality="lifespan") {
  amat <- as.numeric(x["amat"])
  b <- as.numeric(x["b"])
  k <- as.numeric(x["k"])
  amax <- as.numeric(x["amax"])
  M <- as.numeric(x["M"])
  
  l_alpha <- exp(-M*amat)
  alpha_til <- b * l_alpha # new method accounting for survival to maturity
  
  #This is the optimization function, which uses the nlminb() function in R
  func <- function(rmax) (exp(rmax)^amat - exp(-M)*(exp(rmax)^(amat - 1)) - (alpha_til))^2
  nlminb(1/(0.4*amat), func, lower = 0, upper = 5)$par
}

#create dataframe to run rmax from including M
newdat_st <- as.data.frame(cbind(amat_st, b, amax_st, M_st)) # b the same for both spp
names(newdat_st) <- c("amat", "b", "amax", "M") # has to match function

newdat_bf <- as.data.frame(cbind(amat_bf, b, amax_bf, M_bf)) 
names(newdat_bf) <- c("amat", "b", "amax", "M") # has to match function

#run function and bind rmax values to the data frame
rmax_st <- apply(newdat_st, 1, rmaxnlm)
rmax_bf <- apply(newdat_bf, 1, rmaxnlm)
##mean and quantiles
mean(rmax_st)
median(rmax_st)
quantile(rmax_st,0.025)
quantile(rmax_st,0.975)
hist(rmax_st)

mean(rmax_bf)
median(rmax_bf)
quantile(rmax_bf,0.025)
quantile(rmax_bf,0.975)
hist(rmax_bf)

# calculate exploitation ratio E---------------------
#distribution of Z, assumed normal distribution, required n=10000, mean and sd calculated from the catch curve
#SD <- SQRT(n)*SE
#SD <- SQRT(n)*(UCI-LCI)/((tval at df)*2)

Zn_st <- nrow(subset(tblAge, Age>=4))
Zn_bf <- nrow(subset(tblAge2, Age>=3))
#Zsd <- Zse*(sqrt(Zn))
Zsd_st <- sqrt(Zn_st)*(Zupr_st-Zlwr_st)/(qt(0.025, df=Zn_st-1, lower.tail=FALSE)*2)
Zsd_bf <- sqrt(Zn_bf)*(Zupr_bf-Zlwr_bf)/(qt(0.025, df=Zn_bf-1, lower.tail=FALSE)*2)

#create Distribution of Z assuming a normal distribution
ZDist_st <- rnorm(iter, mean=Z_st, sd=Zsd_st)
mean(ZDist_st)
hist(ZDist_st) #distribution of Z

ZDist_bf <- rnorm(iter, mean=Z_bf, sd=Zsd_bf)
mean(ZDist_bf)
hist(ZDist_bf) #distribution of Z

#Fisheries mortality = Z - M
F_st <- ZDist_st-M_st #calc F for M.mobular
mean(F_st)
median(F_st)
quantile(F_st,0.025)
quantile(F_st,0.975)
hist(F_st) #distribution of F

F_bf <- ZDist_bf-M_bf #calc F for M.thurstoni
mean(F_bf)
median(F_bf)
quantile(F_bf,0.025)
quantile(F_bf,0.975)
hist(F_bf) #distribution of F

E_st <- F_st/(ZDist_st) #calc E for M. mobular
mean(E_st)
median(E_st)
quantile(E_st,0.025)
quantile(E_st,0.975)

E_bf <- F_bf/(ZDist_bf) #calc E for M. thurstoni
mean(E_bf)
median(E_bf)
quantile(E_bf,0.025)
quantile(E_bf,0.975)

#proportion >0.5
length(which(E_st>0.5))/length(E_st) # M. mobular
length(which(E_bf>0.5))/length(E_bf) # M. thurstoni

# Fig.8 Plot Rmax and E distributions-----------
Edf_st <- as.data.frame(cbind(rmax_st,ZDist_st, F_st, M_st, E_st))
Edf_bf <- as.data.frame(cbind(rmax_bf,ZDist_bf, F_bf, M_bf, E_bf))

#Rmax
Edf_st$Sp <- "Spinetail devil ray"
St_rm1 <- ggplot(Edf_st, aes(x=rmax_st)) +
  geom_density(alpha=1, fill="grey50") +
  geom_vline(aes(xintercept=median(rmax_st, na.rm=T)), colour="black", linetype="longdash") +
  scale_x_continuous(limits=c(0,0.2), breaks=seq(0,0.2,by=0.05), expand=c(0.01,0.01)) +
  scale_y_continuous(breaks=seq(0,0,by=0.1), expand=c(0.01,0.01)) +
  annotate("text", x=Inf, y=Inf, label="Median = 0.109", hjust=1, vjust=2, fontface = "italic") +
  ylab("Relative density") +
  labs(x=expression(paste("Population growth rate ","(","r"[max],")"))) +
  facet_grid(. ~Sp) +
  set_theme

Edf_bf$Sp <- "Bentfin devil ray"
Bf_rm1 <- ggplot(Edf_bf, aes(x=rmax_bf)) +
  geom_density(alpha=1, fill="grey50") +
  geom_vline(aes(xintercept=median(rmax_bf, na.rm=T)), colour="black", linetype="longdash") +
  scale_x_continuous(limits=c(0,0.2), breaks=seq(0,0.2,by=0.05), expand=c(0.01,0.01)) +
  scale_y_continuous(breaks=seq(0,0,by=0.1), expand=c(0.01,0.01)) +
  annotate("text", x=Inf, y=Inf, label="Median = 0.107", hjust=1, vjust=2, fontface = "italic") +
  ylab("Relative density") +
  labs(x=expression(paste("Population growth rate ","(","r"[max],")"))) +
  facet_grid(. ~Sp) +
  set_theme

#E
St_e1 <- ggplot(Edf_st, aes(x=E_st)) +
  geom_density(alpha=1, fill="grey50") +
  geom_vline(aes(xintercept=median(E_st, na.rm=T)), colour="black", linetype="longdash") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,by=0.2), expand=c(0.01,0.01)) +
  scale_y_continuous(breaks=seq(0,0,by=0.1), expand=c(0.01,0.01)) +
  annotate("text", x=Inf, y=Inf, label="Median = 0.77", hjust=3, vjust=1, fontface= "italic") +
  ylab("Relative density") +
  xlab("Exploitation ratio (E)") +
  set_theme2

Bf_e1 <- ggplot(Edf_bf, aes(x=E_bf)) +
  geom_density(alpha=1, fill="grey50") +
  geom_vline(aes(xintercept=median(E_bf, na.rm=T)), colour="black", linetype="longdash") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,by=0.2), expand=c(0.01,0.01)) +
  scale_y_continuous(breaks=seq(0,0,by=0.1), expand=c(0.01,0.01)) +
  annotate("text", x=Inf, y=Inf, label="Median = 0.80", hjust=3, vjust=1, fontface= "italic") +
  ylab("Relative density") +
  xlab("Exploitation ratio (E)") +
  set_theme

#plot r and e distributions
St_redist1 <- plot_grid(St_rm1, St_e1, ncol=1, labels=c('a','b'))
St_redist1

Bf_redist1 <- plot_grid(Bf_rm1, Bf_e1, ncol=1, labels=c('c','d'))
Bf_redist1

combine_all <- plot_grid(St_redist1, Bf_redist1, nrow = 1, hjust=-0.6, label_fontface = "italic")
combine_all
