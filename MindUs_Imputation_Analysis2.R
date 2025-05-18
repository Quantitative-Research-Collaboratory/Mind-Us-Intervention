### Imputation Analysis ###

#remotes::install_github('blimp-stats/rblimp')
#remotes::install_github("bkeller2/fdir")

library(tidyr)
library(dplyr)
library(readr)
library(rblimp)
library(fdir)
library(nlme)
library(mitml)
library(naniar)

# Make wide version of data
longdat <- read.csv("Mind-Us_Long_12.26.23.csv", sep = ",", header = TRUE)
names(longdat)[1] <- "id"
longdat$Time.reorder <- factor(longdat$Time, levels = c("1", "2", "0"))
widedat <- longdat %>% pivot_wider(id_cols = c(id, Group, Age, EnoughMoney, Sex),
                                   names_from = Time,
                                   values_from = c(PHQ9,GAD7,PSS, MAAS, ERQ_ES, ERQ_CR, RRS, SCSSF, BEAQ),
)

# Missing Data Patterns
arrange(widedat[!complete.cases(widedat),], PHQ9_1, PHQ9_2)

## 13 cases with missing data
## 8 completed baseline of all measures, no follow-up on any measure
## 4 completed baseline and first follow-up on all measures, no second follow-up
## 1 person missed first follow-up (all measures), completed second follow-up (all measures)

##########################################################
# Descriptive Analyses of Missingness
##########################################################

# Treatment Condition
table(widedat[!complete.cases(widedat),]$Group)
prop.table(table(widedat[!complete.cases(widedat),]$Group))
table(widedat[complete.cases(widedat),]$Group)
prop.table(table(widedat[complete.cases(widedat),]$Group))

# Age
mean(widedat[!complete.cases(widedat),]$Age)
mean(widedat[complete.cases(widedat),]$Age)

# Enough Money
table(widedat[!complete.cases(widedat),]$EnoughMoney)
prop.table(table(widedat[!complete.cases(widedat),]$EnoughMoney))
table(widedat[complete.cases(widedat),]$EnoughMoney)
prop.table(table(widedat[complete.cases(widedat),]$EnoughMoney))

# Sex
table(widedat[!complete.cases(widedat),]$Sex)
prop.table(table(widedat[!complete.cases(widedat),]$Sex))
table(widedat[complete.cases(widedat),]$Sex)
prop.table(table(widedat[complete.cases(widedat),]$Sex))

##########################################################
# Missing Completely at Random
##########################################################

mcar_test(widedat)

## Depression

summary(lme(PHQ9 ~ as.factor(Time) * Group + Sex + Age + EnoughMoney, 
    random = ~ 1 | id, 
    correlation = corAR1(form = ~ 1 | id), data = longdat[complete.cases(longdat),]))

summary(lme(PHQ9 ~ Time.reorder * Group + Sex + Age + EnoughMoney, 
            random = ~ 1 | id, 
            correlation = corAR1(form = ~ 1 | id), data = longdat[complete.cases(longdat),]))

## Stress

summary(lme(PSS ~ as.factor(Time) * Group + Sex + Age + EnoughMoney, 
            random = ~ 1 | id, 
            correlation = corAR1(form = ~ 1 | id), data = longdat[complete.cases(longdat),]))

summary(lme(PSS ~ Time.reorder * Group + Sex + Age + EnoughMoney, 
            random = ~ 1 | id, 
            correlation = corAR1(form = ~ 1 | id), data = longdat[complete.cases(longdat),]))



## Anxiety

summary(lme(GAD7 ~ as.factor(Time) * Group + Sex + Age + EnoughMoney, 
            random = ~ 1 | id, 
            correlation = corAR1(form = ~ 1 | id), data = longdat[complete.cases(longdat),]))

summary(lme(GAD7 ~ Time.reorder * Group + Sex + Age + EnoughMoney, 
            random = ~ 1 | id, 
            correlation = corAR1(form = ~ 1 | id), data = longdat[complete.cases(longdat),]))

##########################################################
# Single Depression - MAR imputation model
##########################################################

mar.dep <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept.dep',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept.dep ~ 1; 
    measurement.model:
    PHQ9_0 ~ 1@0 icept.dep@1 Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 icept.dep@1 (PHQ9_0)@lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 icept.dep@1 (PHQ9_1)@lagprior Group Sex Age EnoughMoney;',
  parameters = 
    'lagprior ~ truncate(-1,1);',
  seed = 1016287,
  burn = 10000,
  iter = 10000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all')

# reshape data to long format
mar.dep_implist <- lapply(mar.dep@imputations, (function(x) 
  reshape(x, 
          varying = c("PHQ9_0", "PHQ9_1", "PHQ9_2"), 
          v.names = "PHQ9",
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mar.dep_implist <- as.mitml.list(mar.dep_implist)

# fit model and pool estimates
fit_mar.dep <- with(mar.dep_implist, 
                lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                    random = ~ 1 | id, 
                    correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mar.dep, extra.pars = T)
confint(testEstimates(fit_mar.dep, extra.pars = T))

fit_mar.dep2 <- with(mar.dep_implist, 
                    lme(PHQ9 ~ factor(wave) * Group + Sex + Age + EnoughMoney, 
                        random = ~ 1 | id, 
                        correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mar.dep2, extra.pars = T)
confint(testEstimates(fit_mar.dep2, extra.pars = T))

##########################################################
# Single Stress - MAR imputation model
##########################################################

mar.str <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept.str',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept.str ~ 1; 
    measurement.model:
    PSS_0 ~ 1@0 icept.str@1 Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 icept.str@1 (PSS_0)@lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 icept.str@1 (PSS_1)@lagprior Group Sex Age EnoughMoney;',
  parameters = 
    'lagprior ~ truncate(-1,1);',
  seed = 1016287,
  burn = 10000,
  iter = 10000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all')

# reshape data to long format
mar.str_implist <- lapply(mar.str@imputations, (function(x) 
  reshape(x, 
          varying = c("PSS_0", "PSS_1", "PSS_2"), 
          v.names = "PSS",
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mar.str_implist <- as.mitml.list(mar.str_implist)

# fit model and pool estimates
fit_mar.str <- with(mar.str_implist, 
                    lme(PSS ~ factor(wave) * Group + Sex + Age + EnoughMoney, 
                        random = ~ 1 | id, 
                        correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mar.str, extra.pars = T)
confint(testEstimates(fit_mar.str, extra.pars = T))

fit_mar.str2 <- with(mar.str_implist, 
                    lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                        random = ~ 1 | id, 
                        correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mar.str2, extra.pars = T)
confint(testEstimates(fit_mar.str2, extra.pars = T))

##########################################################
# Single Anxiety - MAR imputation model
##########################################################

mar.anx <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept.anx',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept.anx ~ 1; 
    measurement.model:
    GAD7_0 ~ 1@0 icept.anx@1 Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 icept.anx@1 (GAD7_0)@lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 icept.anx@1 (GAD7_1)@lagprior Group Sex Age EnoughMoney;',
  parameters = 
    'lagprior ~ truncate(-1,1);',
  seed = 1016287,
  burn = 10000,
  iter = 10000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all')

# reshape data to long format
mar.anx_implist <- lapply(mar.anx@imputations, (function(x) 
  reshape(x, 
          varying = c("GAD7_0", "GAD7_1", "GAD7_2"), 
          v.names = "GAD7",
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mar.anx_implist <- as.mitml.list(mar.anx_implist)

# fit model and pool estimates
fit_mar.anx <- with(mar.anx_implist, 
                    lme(GAD7 ~ factor(wave) * Group + Sex + Age + EnoughMoney, 
                        random = ~ 1 | id, 
                        correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mar.anx, extra.pars = T)
confint(testEstimates(fit_mar.anx, extra.pars = T))

fit_mar.anx2 <- with(mar.anx_implist, 
                     lme(GAD7 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                         random = ~ 1 | id, 
                         correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mar.anx2, extra.pars = T)
confint(testEstimates(fit_mar.anx2, extra.pars = T))

#############################################################################
# Multiple (depression, stress, anxiety) correlated shared parameter, MAR
#############################################################################

mar.3shared <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'dep.icept str.icept anx.icept',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    dep.icept ~ 1;
    str.icept ~ 1;
    anx.icept ~ 1;
    dep.icept ~~ str.icept;
    dep.icept ~~ anx.icept;
    str.icept ~~ anx.icept;
    measurement.model:
    PHQ9_0 ~ 1@0 dep.icept@1 Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 dep.icept@1 (PHQ9_0)@dep.lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 dep.icept@1 (PHQ9_1)@dep.lagprior Group Sex Age EnoughMoney;
    PSS_0 ~ 1@0 str.icept@1 Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 str.icept@1 (PSS_0)@str.lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 str.icept@1 (PSS_1)@str.lagprior Group Sex Age EnoughMoney;
    GAD7_0 ~ 1@0 anx.icept@1 Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 anx.icept@1 (GAD7_0)@anx.lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 anx.icept@1 (GAD7_1)@anx.lagprior Group Sex Age EnoughMoney;',
  parameters = 
    'dep.lagprior ~ truncate(-1,1); 
    str.lagprior ~ truncate(-1,1); 
    anx.lagprior ~ truncate(-1,1)',
  seed = 6907734,
  burn = 100000,
  iter = 100000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all')

# reshape data to long format
mar.3shared_implist <- lapply(mar.3shared@imputations, (function(x) 
  reshape(x, 
          varying = list(c("PSS_0", "PSS_1", "PSS_2"),
                         c("PHQ9_0", "PHQ9_1", "PHQ9_2"),
                         c("GAD7_0", "GAD7_1", "GAD7_2")),
          v.names = c("PSS", "PHQ9", "GAD7"),
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mar.3shared_implist <- as.mitml.list(mar.3shared_implist)

# Depression Model: fit model and pool estimates
dep.mar.3shared <- with(mar.3shared_implist, 
                            lme(PHQ9 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mar.3shared, extra.pars = T)
confint(testEstimates(dep.mar.3shared, extra.pars = T))

dep.mar.3shared2 <- with(mar.3shared_implist, 
                            lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mar.3shared2, extra.pars = T)
confint(testEstimates(dep.mar.3shared2, extra.pars = T))

# Stress Model: fit model and pool estimates
str.mar.3shared <- with(mar.3shared_implist, 
                            lme(PSS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mar.3shared, extra.pars = T)
confint(testEstimates(str.mar.3shared, extra.pars = T))

str.mar.3shared2 <- with(mar.3shared_implist, 
                        lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mar.3shared2, extra.pars = T)
confint(testEstimates(str.mar.3shared2, extra.pars = T))

# Anxiety Model: fit model and pool estimates
anx.mar.3shared <- with(mar.3shared_implist, 
                            lme(GAD7 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mar.3shared, extra.pars = T)
confint(testEstimates(anx.mar.3shared, extra.pars = T))

anx.mar.3shared2 <- with(mar.3shared_implist, 
                        lme(GAD7 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mar.3shared2, extra.pars = T)
confint(testEstimates(anx.mar.3shared2, extra.pars = T))

#############################################################################
# Multiple (Primary & Secondary) correlated shared parameter, MAR
#############################################################################

mar.9shared <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'dep.icept str.icept anx.icept maas.icept erqes.icept erqcr.icept rrs.icept scssf.icept beaq.icept',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    dep.icept ~ 1;
    str.icept ~ 1;
    anx.icept ~ 1;
    maas.icept ~ 1;
    erqes.icept ~ 1;
    erqcr.icept ~ 1;
    rrs.icept ~ 1;
    scssf.icept ~ 1;
    beaq.icept ~ 1;
    dep.icept str.icept anx.icept maas.icept erqes.icept erqcr.icept rrs.icept scssf.icept beaq.icept ~~ dep.icept str.icept anx.icept maas.icept erqes.icept erqcr.icept rrs.icept scssf.icept beaq.icept;

    measurement.model:
    PHQ9_0 ~ 1@0 dep.icept@1 Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 dep.icept@1 (PHQ9_0)@dep.lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 dep.icept@1 (PHQ9_1)@dep.lagprior Group Sex Age EnoughMoney;
    
    PSS_0 ~ 1@0 str.icept@1 Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 str.icept@1 (PSS_0)@str.lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 str.icept@1 (PSS_1)@str.lagprior Group Sex Age EnoughMoney;
    
    GAD7_0 ~ 1@0 anx.icept@1 Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 anx.icept@1 (GAD7_0)@anx.lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 anx.icept@1 (GAD7_1)@anx.lagprior Group Sex Age EnoughMoney;
  
    MAAS_0 ~ 1@0 maas.icept@1 Group Sex Age EnoughMoney; 
    MAAS_1 ~ 1 maas.icept@1 (MAAS_0)@maas.lagprior Group Sex Age EnoughMoney; 
    MAAS_2 ~ 1 maas.icept@1 (MAAS_1)@maas.lagprior Group Sex Age EnoughMoney;
    
    ERQ_ES_0 ~ 1@0 erqes.icept@1 Group Sex Age EnoughMoney; 
    ERQ_ES_1 ~ 1 erqes.icept@1 (ERQ_ES_0)@erqes.lagprior Group Sex Age EnoughMoney; 
    ERQ_ES_2 ~ 1 erqes.icept@1 (ERQ_ES_1)@erqes.lagprior Group Sex Age EnoughMoney;
      
    ERQ_CR_0 ~ 1@0 erqcr.icept@1 Group Sex Age EnoughMoney; 
    ERQ_CR_1 ~ 1 erqcr.icept@1 (ERQ_CR_0)@erqcr.lagprior Group Sex Age EnoughMoney; 
    ERQ_CR_2 ~ 1 erqcr.icept@1 (ERQ_CR_1)@erqcr.lagprior Group Sex Age EnoughMoney;
        
    RRS_0 ~ 1@0 rrs.icept@1 Group Sex Age EnoughMoney; 
    RRS_1 ~ 1 rrs.icept@1 (RRS_0)@rrs.lagprior Group Sex Age EnoughMoney; 
    RRS_2 ~ 1 rrs.icept@1 (RRS_1)@rrs.lagprior Group Sex Age EnoughMoney;
          
    SCSSF_0 ~ 1@0 scssf.icept@1 Group Sex Age EnoughMoney; 
    SCSSF_1 ~ 1 scssf.icept@1 (SCSSF_0)@scssf.lagprior Group Sex Age EnoughMoney; 
    SCSSF_2 ~ 1 scssf.icept@1 (SCSSF_1)@scssf.lagprior Group Sex Age EnoughMoney;
            
    BEAQ_0 ~ 1@0 beaq.icept@1 Group Sex Age EnoughMoney; 
    BEAQ_1 ~ 1 beaq.icept@1 (BEAQ_0)@beaq.lagprior Group Sex Age EnoughMoney; 
    BEAQ_2 ~ 1 beaq.icept@1 (BEAQ_1)@beaq.lagprior Group Sex Age EnoughMoney;'
  ,
  parameters = 
    'dep.lagprior ~ truncate(-1,1); 
    str.lagprior ~ truncate(-1,1); 
    anx.lagprior ~ truncate(-1,1);
    maas.lagprior ~ truncate(-1,1);
    erqes.lagprior ~ truncate(-1,1);
    erqcr.lagprior ~ truncate(-1,1);
    rrs.lagprior ~ truncate(-1,1);
    scssf.lagprior ~ truncate(-1,1);
    beaq.lagprior ~ truncate(-1,1)',
  seed = 6907734,
  burn = 100000,
  iter = 100000, 
  chains = 20,
  nimps = 20,
  print_output = 'all')

#save(mar.9shared, file = "mar9shared.R")
load(file = "mar9shared.R")

# reshape data to long format
mar.9shared_implist <- lapply(mar.9shared@imputations, (function(x) 
  reshape(x, 
          varying = list(c("PSS_0", "PSS_1", "PSS_2"),
                         c("PHQ9_0", "PHQ9_1", "PHQ9_2"),
                         c("GAD7_0", "GAD7_1", "GAD7_2"),
                         c("MAAS_0", "MAAS_1", "MAAS_2"),
                         c("ERQ_ES_0", "ERQ_ES_1", "ERQ_ES_2"),
                         c("ERQ_CR_0", "ERQ_CR_1", "ERQ_CR_2"),
                         c("RRS_0", "RRS_1", "RRS_2"),
                         c("SCSSF_0", "SCSSF_1", "SCSSF_2"),
                         c("BEAQ_0", "BEAQ_1", "BEAQ_2")
                         ),
          v.names = c("PSS", "PHQ9", "GAD7", "MAAS", "ERQ_ES", "ERQ_CR", "RRS", "SCSSF", "BEAQ"),
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mar.9shared_implist <- as.mitml.list(mar.9shared_implist)

# MAAS Model: fit model and pool estimates
maas.mar.9shared <- with(mar.9shared_implist, 
                         lme(MAAS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(maas.mar.9shared, extra.pars = T)
confint(testEstimates(maas.mar.9shared, extra.pars = T))

maas.mar.9shared2 <- with(mar.9shared_implist, 
                          lme(MAAS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(maas.mar.9shared2, extra.pars = T)
confint(testEstimates(maas.mar.9shared2, extra.pars = T))

maas.mar.9shared3 <- with(mar.9shared_implist, 
                         lme(MAAS ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(maas.mar.9shared3, extra.pars = T)
confint(testEstimates(maas.mar.9shared3, extra.pars = T))

maas.mar.9shared4 <- with(mar.9shared_implist, 
                          lme(MAAS ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(maas.mar.9shared4, extra.pars = T)
confint(testEstimates(maas.mar.9shared4, extra.pars = T))

# ERQ_ES Model: fit model and pool estimates
erqes.mar.9shared <- with(mar.9shared_implist, 
                          lme(ERQ_ES ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(erqes.mar.9shared, extra.pars = T)
confint(testEstimates(erqes.mar.9shared, extra.pars = T))

erqes.mar.9shared2 <- with(mar.9shared_implist, 
                           lme(ERQ_ES ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(erqes.mar.9shared2, extra.pars = T)
confint(testEstimates(erqes.mar.9shared2, extra.pars = T))

erqes.mar.9shared3 <- with(mar.9shared_implist, 
                          lme(ERQ_ES ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(erqes.mar.9shared3, extra.pars = T)
confint(testEstimates(erqes.mar.9shared3, extra.pars = T))

erqes.mar.9shared4 <- with(mar.9shared_implist, 
                          lme(ERQ_ES ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(erqes.mar.9shared4, extra.pars = T)
confint(testEstimates(erqes.mar.9shared4, extra.pars = T))

# ERQ_CR Model: fit model and pool estimates
erqcr.mar.9shared <- with(mar.9shared_implist, 
                          lme(ERQ_CR ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(erqcr.mar.9shared, extra.pars = T)
confint(testEstimates(erqcr.mar.9shared, extra.pars = T))

erqcr.mar.9shared2 <- with(mar.9shared_implist, 
                           lme(ERQ_CR ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(erqcr.mar.9shared2, extra.pars = T)
confint(testEstimates(erqcr.mar.9shared2, extra.pars = T))

# SCSSF Model: fit model and pool estimates
scssf.mar.9shared <- with(mar.9shared_implist, 
                          lme(SCSSF ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(scssf.mar.9shared, extra.pars = T)
confint(testEstimates(scssf.mar.9shared, extra.pars = T))

scssf.mar.9shared2 <- with(mar.9shared_implist, 
                           lme(SCSSF ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(scssf.mar.9shared2, extra.pars = T)
confint(testEstimates(scssf.mar.9shared2, extra.pars = T))

scssf.mar.9shared3 <- with(mar.9shared_implist, 
                           lme(SCSSF ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(scssf.mar.9shared3, extra.pars = T)
confint(testEstimates(scssf.mar.9shared3, extra.pars = T))

scssf.mar.9shared4 <- with(mar.9shared_implist, 
                           lme(SCSSF ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(scssf.mar.9shared4, extra.pars = T)
confint(testEstimates(scssf.mar.9shared4, extra.pars = T))

# BEAQ Model: fit model and pool estimates
beaq.mar.9shared <- with(mar.9shared_implist, 
                          lme(BEAQ ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(beaq.mar.9shared, extra.pars = T)
confint(testEstimates(beaq.mar.9shared, extra.pars = T))

beaq.mar.9shared2 <- with(mar.9shared_implist, 
                           lme(BEAQ ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(beaq.mar.9shared2, extra.pars = T)
confint(testEstimates(beaq.mar.9shared2, extra.pars = T))

beaq.mar.9shared3 <- with(mar.9shared_implist, 
                           lme(BEAQ ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(beaq.mar.9shared3, extra.pars = T)
confint(testEstimates(beaq.mar.9shared3, extra.pars = T))

beaq.mar.9shared4 <- with(mar.9shared_implist, 
                           lme(BEAQ ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                               random = ~ 1 | id, 
                               correlation = corAR1(form = ~ 1 | id)))
testEstimates(beaq.mar.9shared4, extra.pars = T)
confint(testEstimates(beaq.mar.9shared4, extra.pars = T))

# RRS Model: fit model and pool estimates
rrs.mar.9shared <- with(mar.9shared_implist, 
                         lme(RRS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(rrs.mar.9shared, extra.pars = T)
confint(testEstimates(rrs.mar.9shared, extra.pars = T))

rrs.mar.9shared2 <- with(mar.9shared_implist, 
                          lme(RRS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(rrs.mar.9shared2, extra.pars = T)
confint(testEstimates(rrs.mar.9shared2, extra.pars = T))

rrs.mar.9shared3 <- with(mar.9shared_implist, 
                          lme(RRS ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(rrs.mar.9shared3, extra.pars = T)
confint(testEstimates(rrs.mar.9shared3, extra.pars = T))

rrs.mar.9shared4 <- with(mar.9shared_implist, 
                          lme(RRS ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                              random = ~ 1 | id, 
                              correlation = corAR1(form = ~ 1 | id)))
testEstimates(rrs.mar.9shared4, extra.pars = T)
confint(testEstimates(rrs.mar.9shared4, extra.pars = T))

#########################################################################
# Multiple (depression, stress, anxiety) single shared parameter, MAR
#########################################################################

mar.1shared <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept ~ 1@0;
    icept ~~ icept@1;
    measurement.model:
    PHQ9_0 ~ 1 icept@depscale Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 icept@depscale (PHQ9_0)@dep.lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 icept@depscale (PHQ9_1)@dep.lagprior Group Sex Age EnoughMoney;
    PSS_0 ~ 1@strmean icept@strscale Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 icept@strscale (PSS_0)@str.lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 icept@strscale (PSS_1)@str.lagprior Group Sex Age EnoughMoney;
    GAD7_0 ~ 1@anxmean icept@anxscale Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 icept@anxscale (GAD7_0)@anx.lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 icept@anxscale (GAD7_1)@anx.lagprior Group Sex Age EnoughMoney;',
  parameters = 
    'dep.lagprior ~ truncate(-1,1); 
    str.lagprior ~ truncate(-1,1); 
    anx.lagprior ~ truncate(-1,1)',
  seed = 1016287,
  burn = 10000,
  iter = 10000, 
  chains = 20,
  nimps = 20,
  print_output = 'all')

#save(mar.1shared, file = "mar1shared.R")
load(file = "mar1shared.R")

# reshape data to long format
mar.1shared_implist <- lapply(mar.1shared@imputations, (function(x) 
  reshape(x, 
          varying = list(c("PSS_0", "PSS_1", "PSS_2"),
                         c("PHQ9_0", "PHQ9_1", "PHQ9_2"),
                         c("GAD7_0", "GAD7_1", "GAD7_2")),
          v.names = c("PSS", "PHQ9", "GAD7"),
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mar.1shared_implist <- as.mitml.list(mar.1shared_implist)

# Depression Model: fit model and pool estimates
dep.mar.1shared <- with(mar.1shared_implist, 
                        lme(PHQ9 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mar.1shared, extra.pars = T)
confint(testEstimates(dep.mar.1shared, extra.pars = T))

dep.mar.1shared2 <- with(mar.1shared_implist, 
                         lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mar.1shared2, extra.pars = T)
confint(testEstimates(dep.mar.1shared2, extra.pars = T))

dep.mar.1shared3 <- with(mar.1shared_implist, 
                        lme(PHQ9 ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mar.1shared3, extra.pars = T)
confint(testEstimates(dep.mar.1shared, extra.pars = T))

dep.mar.1shared4 <- with(mar.1shared_implist, 
                         lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mar.1shared4, extra.pars = T)
confint(testEstimates(dep.mar.1shared4, extra.pars = T))

# Stress Model: fit model and pool estimates
str.mar.1shared <- with(mar.1shared_implist, 
                        lme(PSS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mar.1shared, extra.pars = T)
confint(testEstimates(str.mar.1shared, extra.pars = T))

str.mar.1shared2 <- with(mar.1shared_implist, 
                         lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mar.1shared2, extra.pars = T)
confint(testEstimates(str.mar.1shared2, extra.pars = T))

str.mar.3shared <- with(mar.1shared_implist, 
                        lme(PSS ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mar.3shared, extra.pars = T)
confint(testEstimates(str.mar.3shared, extra.pars = T))

str.mar.1shared4 <- with(mar.1shared_implist, 
                         lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mar.1shared4, extra.pars = T)
confint(testEstimates(str.mar.1shared4, extra.pars = T))

# Anxiety Model: fit model and pool estimates
anx.mar.1shared <- with(mar.1shared_implist, 
                        lme(GAD7 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mar.1shared, extra.pars = T)
confint(testEstimates(anx.mar.1shared, extra.pars = T))

anx.mar.1shared2 <- with(mar.1shared_implist, 
                         lme(GAD7 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mar.1shared2, extra.pars = T)
confint(testEstimates(anx.mar.1shared2, extra.pars = T))

anx.mar.3shared <- with(mar.1shared_implist, 
                        lme(GAD7 ~ as.factor(wave) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                            random = ~ 1 | id, 
                            correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mar.3shared, extra.pars = T)
confint(testEstimates(anx.mar.3shared, extra.pars = T))

anx.mar.1shared4 <- with(mar.1shared_implist, 
                         lme(GAD7 ~ factor(wave, levels = c("1", "2", "0")) * factor(Group, levels = c("1", "0")) + Sex + Age + EnoughMoney, 
                             random = ~ 1 | id, 
                             correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mar.1shared4, extra.pars = T)
confint(testEstimates(anx.mar.1shared4, extra.pars = T))

##########################################################
# Single Depression - MNAR imputation model
##########################################################

sharedpar.dep <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept.dep',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept.dep ~ 1; 
    measurement.model:
    PHQ9_0 ~ 1@0 icept.dep@1 Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 icept.dep@1 (PHQ9_0)@lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 icept.dep@1 (PHQ9_1)@lagprior Group Sex Age EnoughMoney;
    #PHQ9_0 PHQ9_1 PHQ9_2~~ PHQ9_0 PHQ9_1 PHQ9_2;
    missingness.model:
    PHQ9_2.missing ~ icept.dep Group EnoughMoney;',
  parameters = 
    'lagprior ~ truncate(-1,1); 
    #dep1mean = t0mean;
    #dep2mean = t0mean + t1diff;
    #dep3mean = t0mean + t2diff;
   ',
  seed = 6907733,
  burn = 20000,
  iter = 20000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all',
  options = 'labels')

# reshape data to long format
sharedpar.dep_implist <- lapply(sharedpar.dep@imputations, (function(x) 
  reshape(x, 
          varying = c("PHQ9_0", "PHQ9_1", "PHQ9_2"), 
          v.names = "PHQ9",
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
sharedpar.dep_implist <- as.mitml.list(sharedpar.dep_implist)

# fit model and pool estimates
fit_shared.dep <- with(sharedpar.dep_implist, 
                   lme(PHQ9 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                       random = ~ 1 | id, 
                       correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_shared.dep, extra.pars = T)
confint(testEstimates(fit_shared.dep, extra.pars = T))

fit_mnar.dep2 <- with(sharedpar.dep_implist, 
                     lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                         random = ~ 1 | id, 
                         correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mnar.dep2, extra.pars = T)
confint(testEstimates(fit_mnar.dep2, extra.pars = T))


##########################################################
# Single Stress - MNAR imputation model
##########################################################

sharedpar.str <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept.str',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept.str ~ 1; 
    measurement.model:
    PSS_0 ~ 1@0 icept.str@1 Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 icept.str@1 (PSS_0)@lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 icept.str@1 (PSS_1)@lagprior Group Sex Age EnoughMoney;
    #PHQ9_0 PHQ9_1 PHQ9_2~~ PHQ9_0 PHQ9_1 PHQ9_2;
    missingness.model:
    PSS_2.missing ~ icept.str Group EnoughMoney;',
  parameters = 
    'lagprior ~ truncate(-1,1); 
   ',
  seed = 6907733,
  burn = 20000,
  iter = 20000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all',
  options = 'labels')

# reshape data to long format
sharedpar.str_implist <- lapply(sharedpar.str@imputations, (function(x) 
  reshape(x, 
          varying = c("PSS_0", "PSS_1", "PSS_2"), 
          v.names = "PSS",
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
sharedpar.str_implist <- as.mitml.list(sharedpar.str_implist)

# fit model and pool estimates
fit_mnar.str <- with(sharedpar.str_implist, 
                       lme(PSS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                           random = ~ 1 | id, 
                           correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mnar.str, extra.pars = T)
confint(testEstimates(fit_mnar.str, extra.pars = T))

fit_mnar.str2 <- with(sharedpar.str_implist, 
                      lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                          random = ~ 1 | id, 
                          correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mnar.str2, extra.pars = T)
confint(testEstimates(fit_mnar.str2, extra.pars = T))

##########################################################
# Single Anxiety - MNAR imputation model
##########################################################

sharedpar.anx <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept.anx',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept.anx ~ 1; 
    measurement.model:
    GAD7_0 ~ 1@0 icept.anx@1 Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 icept.anx@1 (GAD7_0)@lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 icept.anx@1 (GAD7_1)@lagprior Group Sex Age EnoughMoney;
    #PHQ9_0 PHQ9_1 PHQ9_2~~ PHQ9_0 PHQ9_1 PHQ9_2;
    missingness.model:
    GAD7_2.missing ~ icept.anx Group EnoughMoney;',
  parameters = 
    'lagprior ~ truncate(-1,1); 
   ',
  seed = 6907733,
  burn = 20000,
  iter = 20000, 
  #chains = 20,
  #nimps = 20,
  print_output = 'all',
  options = 'labels')

# reshape data to long format
sharedpar.anx_implist <- lapply(sharedpar.anx@imputations, (function(x) 
  reshape(x, 
          varying = c("GAD7_0", "GAD7_1", "GAD7_2"), 
          v.names = "GAD7",
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
sharedpar.anx_implist <- as.mitml.list(sharedpar.anx_implist)

# fit model and pool estimates
fit_mnar.anx <- with(sharedpar.anx_implist, 
                       lme(GAD7 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                           random = ~ 1 | id, 
                           correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mnar.anx, extra.pars = T)
confint(testEstimates(fit_mnar.anx, extra.pars = T))

fit_mnar.anx2 <- with(sharedpar.anx_implist, 
                      lme(GAD7 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                          random = ~ 1 | id, 
                          correlation = corAR1(form = ~ 1 | id)))
testEstimates(fit_mnar.anx2, extra.pars = T)
confint(testEstimates(fit_mnar.anx2, extra.pars = T))


##########################################################################
# Multiple (depression, stress, anxiety) correlated shared 
# MNAR
##########################################################################

mult.3sharedpar <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'dep.icept str.icept anx.icept',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    dep.icept ~ 1;
    str.icept ~ 1;
    anx.icept ~ 1;
    dep.icept ~~ str.icept;
    dep.icept ~~ anx.icept;
    str.icept ~~ anx.icept;
    measurement.model:
    PHQ9_0 ~ 1@0 dep.icept@1 Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 dep.icept@1 (PHQ9_0)@dep.lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 dep.icept@1 (PHQ9_1)@dep.lagprior Group Sex Age EnoughMoney;
    PSS_0 ~ 1@0 str.icept@1 Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 str.icept@1 (PSS_0)@str.lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 str.icept@1 (PSS_1)@str.lagprior Group Sex Age EnoughMoney;
    GAD7_0 ~ 1@0 anx.icept@1 Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 anx.icept@1 (GAD7_0)@anx.lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 anx.icept@1 (GAD7_1)@anx.lagprior Group Sex Age EnoughMoney;
    missingness.model:
    PHQ9_2.missing ~ dep.icept str.icept anx.icept Group EnoughMoney;',
  parameters = 
    'dep.lagprior ~ truncate(-1,1); 
    str.lagprior ~ truncate(-1,1); 
    anx.lagprior ~ truncate(-1,1)',
  seed = 6907734,
  burn = 100000,
  iter = 100000, 
  chains = 20,
  nimps = 20,
  print_output = 'all')

# reshape data to long format
mult.3sharedpar_implist <- lapply(mult.3sharedpar@imputations, (function(x) 
  reshape(x, 
          varying = list(c("PSS_0", "PSS_1", "PSS_2"),
                         c("PHQ9_0", "PHQ9_1", "PHQ9_2"),
                         c("GAD7_0", "GAD7_1", "GAD7_2")),
          v.names = c("PSS", "PHQ9", "GAD7"),
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mult.3sharedpar_implist <- as.mitml.list(mult.3sharedpar_implist)

# Depression Model: fit model and pool estimates
dep.mult.fit_shared <- with(mult.3sharedpar_implist, 
                            lme(PHQ9 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mult.fit_shared, extra.pars = T)
confint(testEstimates(dep.mult.fit_shared, extra.pars = T))

dep.mult.fit_shared2 <- with(mult.3sharedpar_implist, 
                            lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mult.fit_shared2, extra.pars = T)
confint(testEstimates(dep.mult.fit_shared2, extra.pars = T))

# Stress Model: fit model and pool estimates
str.mult.fit_shared <- with(mult.3sharedpar_implist, 
                            lme(PSS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mult.fit_shared, extra.pars = T)
confint(testEstimates(str.mult.fit_shared, extra.pars = T))

str.mult.fit_shared2 <- with(mult.3sharedpar_implist, 
                            lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mult.fit_shared2, extra.pars = T)
confint(testEstimates(str.mult.fit_shared2, extra.pars = T))

# Anxiety Model: fit model and pool estimates
anx.mult.fit_shared <- with(mult.3sharedpar_implist, 
                            lme(GAD7 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mult.fit_shared, extra.pars = T)
confint(testEstimates(anx.mult.fit_shared, extra.pars = T))

anx.mult.fit_shared2 <- with(mult.3sharedpar_implist, 
                            lme(GAD7 ~ factor(wave, levels = c("1", "2", "0"))* Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mult.fit_shared2, extra.pars = T)
confint(testEstimates(anx.mult.fit_shared2, extra.pars = T))

##########################################################################
# Multiple (depression, stress, anxiety) single parameter imputation model
# MNAR
##########################################################################

mult.sharedpar <- rblimp(
  data = widedat,
  nominal = 'Group Sex',
  latent = 'icept',
  fixed = 'Group Sex Age EnoughMoney',
  center = 'Group Sex Age EnoughMoney',
  model = '
    random.intercept.model:
    icept ~ 1@0;
    icept ~~ icept@1;
    
    measurement.model:
    PHQ9_0 ~ 1 icept@depscale Group Sex Age EnoughMoney; 
    PHQ9_1 ~ 1 icept@depscale (PHQ9_0)@dep.lagprior Group Sex Age EnoughMoney; 
    PHQ9_2 ~ 1 icept@depscale (PHQ9_1)@dep.lagprior Group Sex Age EnoughMoney;
    PSS_0 ~ 1@strmean icept@strscale Group Sex Age EnoughMoney; 
    PSS_1 ~ 1 icept@strscale (PSS_0)@str.lagprior Group Sex Age EnoughMoney; 
    PSS_2 ~ 1 icept@strscale (PSS_1)@str.lagprior Group Sex Age EnoughMoney;
    GAD7_0 ~ 1@anxmean icept@anxscale Group Sex Age EnoughMoney; 
    GAD7_1 ~ 1 icept@anxscale (GAD7_0)@anx.lagprior Group Sex Age EnoughMoney; 
    GAD7_2 ~ 1 icept@anxscale (GAD7_1)@anx.lagprior Group Sex Age EnoughMoney;
    missingness.model:
    PHQ9_2.missing ~ icept Group EnoughMoney;',
  parameters = 
    'dep.lagprior ~ truncate(-1,1); 
    str.lagprior ~ truncate(-1,1); 
    anx.lagprior ~ truncate(-1,1)',
  seed = 6907733,
  burn = 20000,
  iter = 20000, 
  chains = 20,
  nimps = 20,
  print_output = 'all')

# reshape data to long format
mult.sharedpar_implist <- lapply(mult.sharedpar@imputations, (function(x) 
  reshape(x, 
          varying = list(c("PSS_0", "PSS_1", "PSS_2"),
                         c("PHQ9_0", "PHQ9_1", "PHQ9_2"),
                         c("GAD7_0", "GAD7_1", "GAD7_2")),
          v.names = c("PSS", "PHQ9", "GAD7"),
          timevar = "wave", 
          times = c(0,1,2), 
          new.row.names = 1:(nrow(widedat)*3),
          direction = "long")
))

# convert list of files to mitml list
mult.sharedpar_implist <- as.mitml.list(mult.sharedpar_implist)

# Depression Model: fit model and pool estimates
dep.mult.fit_shared <- with(mult.sharedpar_implist, 
                            lme(PHQ9 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mult.fit_shared, extra.pars = T)
confint(testEstimates(dep.mult.fit_shared, extra.pars = T))

dep.mult.fit_shared2 <- with(mult.sharedpar_implist, 
                             lme(PHQ9 ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                                 random = ~ 1 | id, 
                                 correlation = corAR1(form = ~ 1 | id)))
testEstimates(dep.mult.fit_shared2, extra.pars = T)
confint(testEstimates(dep.mult.fit_shared2, extra.pars = T))

# Stress Model: fit model and pool estimates
str.mult.fit_shared <- with(mult.sharedpar_implist, 
                            lme(PSS ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mult.fit_shared, extra.pars = T)
confint(testEstimates(str.mult.fit_shared, extra.pars = T))

str.mult.fit_shared2 <- with(mult.sharedpar_implist, 
                             lme(PSS ~ factor(wave, levels = c("1", "2", "0")) * Group + Sex + Age + EnoughMoney, 
                                 random = ~ 1 | id, 
                                 correlation = corAR1(form = ~ 1 | id)))
testEstimates(str.mult.fit_shared2, extra.pars = T)
confint(testEstimates(str.mult.fit_shared2, extra.pars = T))

# Anxiety Model: fit model and pool estimates
anx.mult.fit_shared <- with(mult.sharedpar_implist, 
                            lme(GAD7 ~ as.factor(wave) * Group + Sex + Age + EnoughMoney, 
                                random = ~ 1 | id, 
                                correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mult.fit_shared, extra.pars = T)
confint(testEstimates(anx.mult.fit_shared, extra.pars = T))

anx.mult.fit_shared2 <- with(mult.sharedpar_implist, 
                             lme(GAD7 ~ factor(wave, levels = c("1", "2", "0"))* Group + Sex + Age + EnoughMoney, 
                                 random = ~ 1 | id, 
                                 correlation = corAR1(form = ~ 1 | id)))
testEstimates(anx.mult.fit_shared2, extra.pars = T)
confint(testEstimates(anx.mult.fit_shared2, extra.pars = T))
