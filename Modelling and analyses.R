#*******************************************************************************
# R codes for modelling and statistical analyses
# Author: Shengman Lyu
# E-mail:shengman.lyu@gmail.com
# More information can be found in Lyu, S. and Alexander, J. (2022) Compensatory responses of vital rates attenuate impacts of competition on population growth and promote coexistence
# Date updated: 19.12.2022
#*******************************************************************************
rm(list=ls())

library(readxl)
library(tidyverse)
library(car)
library(GGally)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(lme4)
library(lmerTest)

#*****************************************************************************
# 1. Responses and contributions of individual vital rates ----
#*****************************************************************************

#*****************************************************************************
# ** - 1.1 Perturb intrinsic and invasion IPMs ----
#*****************************************************************************
# vital rates
# Note: This should be modified accordingly.
vr <- read_excel("Estimated_vital_rates.xlsx", col_names = TRUE, na="NA")
vr <- vr[-1,]
vr

# reshuffle for each competition-dependent vital rates
vr.reshuffle <- vr
vr.reshuffle[,10:44] <- as.numeric(NA)

vr.survival <- 
  vr.growth <- 
  vr.flowering <- 
  vr.fecundity <- 
  vr.establish <- 
  vr.size.lower <- 
  vr.size.upper <- vr.reshuffle

for(a in 1:nrow(vr.reshuffle)) {
  # a = 3154
  da <- vr.reshuffle[a, ]
  fc <- as.character(da$focal.species)
  bg <- as.character(da$background.species)
  st <- as.character(da$site)
  
  # skip when no pair or no competitor
  if(is.na(da$pair) | bg == "site" | bg == "site_none" | bg == "species" | bg == "species_none") next
  
  # print
  print(a)
  print(paste(fc,st,bg, sep="_"))
  
  # intrinsic vital rates
  vr.intrinsic <- filter(vr, focal.species == fc & site == st  & background.species == "none")
  
  # vital rates in the presence of competition
  vr.comp <- filter(vr, focal.species == fc & site == st  & background.species == bg)
  
  # copy intrinsic vital rates to competition pair
  vr.survival[a, 10:44] <- 
    vr.growth[a, 10:44] <- 
    vr.flowering[a, 10:44] <- 
    vr.fecundity[a, 10:44] <- 
    vr.establish[a, 10:44] <- 
    vr.size.lower[a, 10:44] <- 
    vr.size.upper[a, 10:44] <- vr.intrinsic[10:44] 
  
  #*******************************************************************
  # replace intrinsic vital rates with invading vital rates
  #*******************************************************************
  
  ## **** 1.1.1 survival--------------------------------------------------------
  vr.survival[a, "surv.int"] = vr.comp[, "surv.int"]
  vr.survival[a, "surv.slope"] =  vr.comp[, "surv.slope"]
  
  ## **** 1.1.2 growth----------------------------------------------------------
  vr.growth[a, "growth.int"] = vr.comp[, "growth.int"]
  vr.growth[a, "growth.slope"] = vr.comp[, "growth.slope"]
  vr.growth[a,"growth.sd"] = vr.comp[,"growth.sd"]
  
  ## **** 1.1.3 flowering-------------------------------------------------------
  vr.flowering[a, "flowering.int"] = vr.comp[, "flowering.int"]
  vr.flowering[a, "flowering.slope"] = vr.comp[, "flowering.slope"]
  
  # **** 1.1.4 fecundity_linear ------------------------------------------------
  vr.fecundity[a, "fecundity.int"] = vr.comp[, "fecundity.int"]
  vr.fecundity[a, "fecundity.slope"] = vr.comp[, "fecundity.slope"]
  
  # **** 1.1.5 germination (competition-independent) ---------------------------
  
  # **** 1.1.6.1 seedling establishment of FuNiche (competition independent) ---
  
  # **** 1.1.6.2 seedling establishment with competition------------------------
  vr.establish[a,"establishment.comp.int"] = vr.comp[,"establishment.comp.int"]
  
  # **** 1.1.7 recruit size (competition independent) --------------------------
  
  # **** 1.1.8.1 size range: L -------------------------------------------------
  vr.size.lower[a, "L"] = vr.comp[, "L"]
  
  # **** 1.1.8.2 size range: U -------------------------------------------------
  vr.size.upper[a, "U"] = vr.comp[, "U"]
  
  rm(vr.comp)
}

# reshape to the long format
vr$vital.rate <- "Original"
vr.survival$vital.rate <- "Survival"
vr.growth$vital.rate <- "Growth"
vr.flowering$vital.rate <- "Flowering"
vr.fecundity$vital.rate <- "Fecundity"
vr.establish$vital.rate <- "Establishment"
vr.size.lower$vital.rate <- "Size.lower"
vr.size.upper$vital.rate <- "Size.upper"

vr.reorganize <- bind_rows(vr, vr.survival, vr.growth, vr.flowering, vr.fecundity, vr.establish, vr.size.lower, vr.size.upper) %>%
  filter(!background.species %in% c("site", "site_none", "species", "species_none"))

vr.reorganize$pgr == as.numeric(NA)
vr.reorganize$vital.rate

# output
vr_perturbed_20220125 <- vr.reorganize

#*****************************************************************************
# ** - 1.2 Calculate lambda of perturbed IPMs using cluster ----
#*****************************************************************************
library(parallel)

# load functions
# Note: This should be modified accordingly.
source("Functions.R")

vr.pert <- vr_perturbed_20220125
n.cores <- 2
n.seq <- 1:10
# Note: If you want to compute the whole data set. This can take very long to compute on laptops.
# n.seq <- 1:nrow(vr.pert)

vital.rate_competition <- function(a) {
  # a = 3153
  da <- vr.pert[a,]
  
  if(is.na(da$pair) | 
     da$pair == "no" | 
     da$background.species %in% c("site","site_none","species", "species_none")) {
    lambda.a <- list(lambda = "skip")
    print(paste(a, "_skip"))
  }
  else {
    ipm <- mk.kernel(n=3000, L=as.numeric(da$L), U=as.numeric(da$U), par=da)
    if(sum(is.na(ipm$K))>0) { 
      lambda.a <- list(lambda = "null")
      print(paste(a, "_null")) 
    } 
    else {
      cal.a <- lambda.eigen(ipm$K)
      lam.a <- cal.a$lambda
      lambda.a <- list(lambda = lam.a)
      print(a)
    }
  }
  return(lambda.a)
}

# n.seq = 1
lambda.pert <- mclapply(n.seq, vital.rate_competition, mc.cores = n.cores)

# output
# save(lambda.pert, file="lambda_perturbed.RData")

#*****************************************************************************
# 2. Vital rates contribution ----
#*****************************************************************************

#*****************************************************************************
# ** - 2.1 Calculate vital rate contributions ----
#*****************************************************************************
# perturbed vital rates
vr.pert <- vr_perturbed_20220125

# perturbed lambdas
# Note: You can have this data set by running Section 1.2.
load("lambda_perturbed.RData")
length(lambda.pert)

# combine vital rates with perturbed lambdas
vr.pert$pgr.pert <- as.numeric(NA)

for(i in 1:nrow(vr.pert)) {
  # i = 1
  if(lambda.pert[[i]]$lambda != "skip" & lambda.pert[[i]]$lambda != "null") {
    vr.pert[i,"pgr.pert"] = lambda.pert[[i]]$lambda
  }
}

# divide the perturbated and real lambdas
vr.orig <- filter(vr.pert, vital.rate == "Original")
vr.orig$pgr.pert
vr.pert <- filter(vr.pert, vital.rate != "Original")
vr.pert$pgr.pert
unique(vr.pert$vital.rate)

#*****************************************************************
# Calculate changes in PGR
vr.pert$pgr.intrinsic <- as.numeric(NA) 
vr.pert$pgr.delta.vr <- as.numeric(NA)
vr.pert$pgr.delta.full <- as.numeric(NA)
vr.pert$pgr.delta.rii <- as.numeric(NA)
vr.pert$pgr.delta.logRR <- as.numeric(NA)
vr.pert$pgr.delta.percent <- as.numeric(NA)

for(i in 1:nrow(vr.pert)) {
  # i = 2523
  vr.i <- vr.pert[i,]
  pgr.pert.i <- vr.i$pgr.pert
  if(is.na(pgr.pert.i) | vr.i$background.species == "none") next
  
  pgr.intrinsic.i <- filter(vr.orig, focal.species == vr.i$focal.species & site == vr.i$site & background.species == "none")$pgr.pert
  pgr.original.i <- filter(vr.orig, focal.species == vr.i$focal.species & site == vr.i$site & background.species == vr.i$background.species)$pgr.pert
  vr.pert[i, "pgr.intrinsic"] = pgr.intrinsic.i
  
  # using raw PGR
  vr.pert[i, "pgr.delta.vr"] = round(pgr.pert.i,8) - round(pgr.intrinsic.i,8)
  vr.pert[i, "pgr.delta.full"] = round(pgr.original.i,8) - round(pgr.intrinsic.i,8)
  vr.pert[i, "pgr.delta.rii"] = (round(pgr.original.i,8) - round(pgr.intrinsic.i,8))/(round(pgr.original.i,8) + round(pgr.intrinsic.i,8))
  vr.pert[i, "pgr.delta.percent"] = (round(pgr.original.i,8) - round(pgr.intrinsic.i,8))/(round(pgr.intrinsic.i,8))
  vr.pert[i, "pgr.delta.logRR"] = log(round(pgr.original.i,8)) - log(round(pgr.intrinsic.i,8))
}

# combine delta pgr
pgr.delta <- tibble(focal.species = vr.orig$focal.species,
                    background.species = vr.orig$background.species,
                    site = vr.orig$site,
                    elevation = vr.orig$elevation,
                    origin.focal = vr.orig$origin.focal,
                    origin.background = vr.orig$origin.background,
                    pgr.original = vr.orig$pgr.pert,
                    pgr.delta.full = filter(vr.pert, vital.rate == "Survival")$pgr.delta.full,
                    pgr.delta.rii = filter(vr.pert, vital.rate == "Survival")$pgr.delta.rii,
                    pgr.delta.logRR = filter(vr.pert, vital.rate == "Survival")$pgr.delta.logRR,
                    pgr.delta.percent = filter(vr.pert, vital.rate == "Survival")$pgr.delta.percent,
                    pgr.delta.survival = filter(vr.pert, vital.rate == "Survival")$pgr.delta.vr,
                    pgr.delta.growth = filter(vr.pert, vital.rate == "Growth")$pgr.delta.vr,
                    pgr.delta.flowering = filter(vr.pert, vital.rate == "Flowering")$pgr.delta.vr,
                    pgr.delta.fecundity = filter(vr.pert, vital.rate == "Fecundity")$pgr.delta.vr,
                    pgr.delta.establishment = filter(vr.pert, vital.rate == "Establishment")$pgr.delta.vr,
                    pgr.delta.size.lower = filter(vr.pert, vital.rate == "Size.lower")$pgr.delta.vr,
                    pgr.delta.size.upper = filter(vr.pert, vital.rate == "Size.upper")$pgr.delta.vr
)

#*****************************************************************
# standardize vital rate contributions
pgr.delta$Survival <-
  pgr.delta$Growth <-
  pgr.delta$Flowering <-
  pgr.delta$Fecundity <-
  pgr.delta$Establishment <- 
  pgr.delta$Size.lower <- 
  pgr.delta$Size.upper <- as.numeric(NA)

for(i in 1:nrow(pgr.delta)) {
  # i = 1
  di <- pgr.delta[i,]
  if(is.na(di[8])[1,1]) { next; print("skip") }
  
  di.standarsized <- ltre.standardise(as.numeric(di[,c("pgr.delta.survival", "pgr.delta.growth", "pgr.delta.flowering", "pgr.delta.fecundity", "pgr.delta.establishment")]))
  pgr.delta[i,"Survival"] <- di.standarsized[1]
  pgr.delta[i,"Growth"] <- di.standarsized[2]
  pgr.delta[i,"Flowering"] <- di.standarsized[3]
  pgr.delta[i,"Fecundity"] <- di.standarsized[4]
  pgr.delta[i,"Establishment"] <- di.standarsized[5]
  #vr.orig[i,"Size.lower"] <- di.standarsized[6]
  #vr.orig[i,"Size.upper"] <- di.standarsized[7]
}

# output
# write.csv(pgr.delta, "vital.rate_contribution_logRR_20221020.csv")
vr.cont <- pgr.delta

#*****************************************************************************
# ** - 2.2 Competition effects on population growth ----
#*****************************************************************************
# vital rate contributions
vr.cont

# In total 312 pairs, 284 pair after exluding Armo and Trba (too few data)
# 20/284 = 7% facilitation exlucding Armo and Trba
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(!is.na(pgr.delta.rii)) %>%
  filter(pgr.delta.rii > 0)

# Competition reduce population growth rate by 55.1%
# excluding Armo and Trba: -0.5511067 -0.9094631 -0.1009579
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(!is.na(pgr.delta.rii)) %>%
  filter(pgr.delta.rii < 0) %>%
  .$pgr.delta.percent %>%
  mean_ci_quantile()

# Competition intensity
fig.logRR <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  ggplot(aes(x=pgr.delta.logRR)) +
  geom_histogram(col="grey") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_x_continuous(name = "ln(Response ratio)") +
  scale_y_continuous(name="Number of species pairs")
fig.logRR

# Competition across site and species
fig.logRR_site_focal <- vr.cont %>%
  filter(!(focal.species == "Armo" & site == "Les Posses")) %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low","Middle","High"))) %>% 
  ggplot(aes(x=site, y = pgr.delta.logRR, group=focal.species)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_point(col="grey") +
  stat_summary(col="black") +
  stat_summary(fun="mean", geom="line") +
  geom_vline(xintercept = 0, linetype="dotted") +
  facet_wrap(~focal.species, ncol=7) +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name="ln(Response ratio)")
fig.logRR_site_focal

# lmer test
vr.cont$ID.plot <- as.factor(paste(vr.cont$background.species, vr.cont$site, sep="_"))
vr.cont %>% 
  filter(!(focal.species == "Armo" & site == "Les Posses")) %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii <= 0) %>% # include only competitive interactions
  lmer(pgr.delta.logRR ~ site*origin.focal + (site|focal.species) + (site|background.species) + (1|ID.plot), data=.) %>%
  Anova()
  #summary()
  #rand()
# site: F2,297 = 1.8505 ,  P = 0.396
# focal species x site: F5, 297=15.100, P = 0.009 **

#*****************************************************************************
# ** - 2.3 PGR delta vs. sum of each vital rate ----
#*****************************************************************************
# Vital rate replacement works very well! It could be useful for other studies.
fig.pgr.delta_compare <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  mutate(pgr.delta.sum = 
           pgr.delta.survival + 
           pgr.delta.growth +
           pgr.delta.flowering +
           pgr.delta.fecundity +
           pgr.delta.establishment) %>%
  ggplot(aes(x=pgr.delta.full, y = pgr.delta.sum, label=paste(focal.species, background.species, site))) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_abline(intercept = 0, slope=1) +
  #geom_text() +
  geom_point()
fig.pgr.delta_compare

#*****************************************************************************
# ** - 2.4 Contributions across vital rates ----
#*****************************************************************************
# plot
fig.contribution_vital.rate <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii <= 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  ggplot(aes(x=vital.rate, y=relative.contribution)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_jitter(col="grey", alpha=0.5, height=0, width=0.1) +
  #stat_summary(geom="errorbar",col="grey30", fun.data = ggplot2::mean_se, position = position_dodge(width=0.3), width=0.2)+
  geom_boxplot(alpha=0.5, outlier.alpha = 0, width=0.5) +
  scale_x_discrete(name="Vital rate") +
  scale_y_continuous(name="Normalized contribution")
fig.contribution_vital.rate

# mean and SE of each vital rates
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  group_by(vital.rate) %>%
  summarise(mean_ci_se(relative.contribution))

#*****************************************************************************
# ** - 2.5 Contributions across focal species and site ----
#*****************************************************************************

#***********************************
# focal species
#***********************************
fig.contribution_focal_vital.rate <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  mutate(origin.focal=dplyr::recode(origin.focal, "Highland"= "Highland species", "Lowland" = "Lowland species")) %>%
  filter(relative.contribution != 0) %>%
  ggplot(aes(x=focal.species, y=relative.contribution, col=origin.focal)) +
  #ggplot(aes(x=focal.species, y=relative.contribution)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.1,linetype="dotted") +
  geom_hline(yintercept = -0.1,linetype="dotted") +
  geom_jitter(col="grey", alpha=0.5, height=0, width=0) +
  stat_summary(geom="pointrange",fun.data = mean_sd.ggplot)  +
  scale_color_manual(values=c("Lowland species" = "#FB9A06FF", "Highland species" = "#3E4A89FF")) +
  facet_wrap(~vital.rate, nrow = 1) +
  scale_x_discrete(name="Focal species") +
  scale_y_continuous(name="Normalized contribution")
fig.contribution_focal_vital.rate

# mean and SE
mean.se_focal.species <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  filter(relative.contribution != 0) %>%
  group_by(vital.rate, focal.species) %>%
  summarise(mean_ci_quantile(relative.contribution))

# growth
mean.se_focal.species %>%
  filter(vital.rate == "Growth")

# vital rate contributions  > 10%
mean.se_focal.species %>%
  filter(y < -0.1)

#***********************************
# sites
#***********************************
fig.contribution_site_vital.rate <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low","Middle","High"))) %>% 
  filter(relative.contribution != 0) %>%
  ggplot(aes(x=site, y=relative.contribution)) +
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_jitter(col="grey", alpha=0.5, height=0, width=0) +
  stat_summary(geom="pointrange",fun.data = mean_sd.ggplot) +
  facet_wrap(~vital.rate, nrow=1, scale="free_y", ) +
  scale_x_discrete(name="Site") + 
  scale_y_continuous(name="Normalized contribution")
fig.contribution_site_vital.rate

# mean +_ SE
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low","Middle","High"))) %>% 
  filter(relative.contribution != 0) %>%
  group_by(vital.rate, site) %>%
  summarise(mean_ci_quantile(relative.contribution))

#***********************************
# lmer test
# overall model
vr.cont$ID.plot <- as.factor(paste(vr.cont$background.species, vr.cont$site, sep="_"))

vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii <= 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  filter(relative.contribution != 0) %>%  # include only non-zero vital rate contributions
  lmer(relative.contribution ~ site*vital.rate + origin.focal + site:origin.focal + vital.rate:origin.focal  +(vital.rate|focal.species) + (vital.rate|background.species) + (1|ID.plot), data=.) %>%
  #lmer(relative.contribution ~ site*vital.rate*origin.focal + (vital.rate|focal.species) + (vital.rate|background.species) + (1|ID.plot), data=.) %>%
  Anova()
  #rand() # test for focal and competitor species 
  #summary()

# individual vital rates
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii <= 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%# include only non-zero vital rate contributions
  split(.$vital.rate) %>%
  purrr::map(~lmer(relative.contribution ~ site*origin.focal + (1|focal.species) + (1|background.species) + (1|ID.plot), data=.)) %>%
  purrr::map(Anova)
  #purrr::map(rand) # test for focal and competitor species 
  #purrr::map(summary)

#*****************************************************************************
# 3. Effects of demographic compensation on competition ----
#*****************************************************************************

#*****************************************************************************
# ** - 3.1 The presence of compensatory responses ----
#*****************************************************************************

#***********************************
# how many positive responses
# 373/(264*5) = 28.3% of vital rates
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low","Middle","High"))) %>% 
  filter(relative.contribution != 0) %>%
  filter(relative.contribution >0)

# 243/264 = 92.1% of IPMs
vr.cont$n.positive <- rowSums(vr.cont[,c("Survival", "Growth", "Flowering", "Fecundity", "Establishment")] > 0)
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  filter(!is.na(n.positive)) %>%
  .$n.positive %>% table()

#***********************************
# which vital rates?
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  mutate(if.positive = relative.contribution > 0) %>%
  group_by(vital.rate, if.positive) %>%
  summarise(table(if.positive))

# survial: 100/(100+113) = 46.9%
# growth: 12/(12+252) = 4.5%
# flowering: 178/(178+46) = 79.4%
# fecundity: 80/(80+94) = 45.97%
# establishment: 3/(47+3) = 6%

#***********************************
# the presence of positve effects
fig.positive_focal_vital.rate <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  mutate(if.positive = ifelse(relative.contribution >0, 1, 0)) %>%
  ggplot(aes(x=focal.species, y=if.positive)) +
  geom_jitter(col="grey", alpha=0.5, height=0.1, width=0) +
  #stat_summary(geom="pointrange",fun.data = ggplot2::mean_se)  +
  scale_color_manual(values=c("Lowland" = "#FB9A06FF", "Highland" = "#3E4A89FF")) +
  facet_wrap(~vital.rate, nrow = 1) +
  scale_x_discrete(name="Focal species") +
  scale_y_continuous(name="Positive response or not")
fig.positive_focal_vital.rate

# glmer test
vr.cont$ID.plot <- as.factor(paste(vr.cont$background.species, vr.cont$site, sep="_"))

model1 <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  mutate(if.positive = ifelse(relative.contribution >0, 1, 0)) %>%
  glmer(if.positive ~ site*vital.rate + (1|focal.species) + (1|background.species) + (1|ID.plot), data=., family=binomial)

model2 <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  mutate(if.positive = ifelse(relative.contribution >0, 1, 0)) %>%
  glmer(if.positive ~ site*vital.rate + (vital.rate|focal.species) + (1|background.species) + (1|ID.plot), data=.,  family=binomial)

model3 <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  mutate(if.positive = ifelse(relative.contribution >0, 1, 0)) %>%
  glmer(if.positive ~ site*vital.rate + (1|focal.species) + (vital.rate|background.species) + (1|ID.plot), data=., family=binomial)

model4 <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  pivot_longer(cols=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"), values_to="relative.contribution", names_to="vital.rate") %>%
  mutate(vital.rate = factor(vital.rate, levels=c("Survival", "Growth", "Flowering", "Fecundity", "Establishment"))) %>%
  filter(relative.contribution != 0) %>%
  mutate(if.positive = ifelse(relative.contribution >0, 1, 0)) %>%
  glmer(if.positive ~ site*vital.rate + origin.focal + site:origin.focal + vital.rate:origin.focal  +(vital.rate|focal.species) + (vital.rate|background.species) + (1|ID.plot), data=., family=binomial)

# test for fixed effects
Anova(model4)

# test for focal species
anova(model1, model2)

# test for competitor species
anova(model1, model3)

#*****************************************************************************
# ** - 3.3 Effects of positive responses on competition intensity ----
#*****************************************************************************

#***********************************
# The strength of positive responses at the IPM level
vr.cont$strength.postive <- as.numeric(NA)
for(i in 1:nrow(vr.cont)) {
  # i = 1
  di <- vr.cont[i,]
  ci <- c(di$Survival,di$Growth, di$Flowering,di$Fecundity,di$Establishment)
  ci.sum <- sum((ci>0) * ci)
  vr.cont[i,"strength.postive"] = ci.sum
  ci.sum = NA
}
hist(vr.cont$strength.postive)

fig.positive_strength_logRR <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  filter(strength.postive != 0) %>%
  ggplot(aes(x=pgr.delta.logRR, y = log(strength.postive))) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x, col="black", se=FALSE) +
  #facet_wrap(~site) +
  scale_x_continuous(name=("ln(Response ratio)")) +
  scale_y_continuous(name=("ln(Proportion of compensatory response)"))
fig.positive_strength_logRR

# lmer test
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(strength.postive != 0) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  lmer(log(strength.postive) ~ pgr.delta.logRR + (1|focal.species) + (1|background.species), data=.) %>%
  Anova()
  #summary()

#*****************************************************************************
# 4. Effects of demographic compensation on population growth ----
#*****************************************************************************

#*****************************************************************************
# ** - 4.1 Replace vital rate showing compensatory responses ----
#*****************************************************************************
# vital rate contribution
vr.cont

# original and perturbated vital rates
vr.pert <- vr_perturbed_20220125
vr.orig <- filter(vr.pert, vital.rate == "Original")

# check if they are the same
sum(paste(vr.cont$focal.species,  vr.cont$background.species, vr.cont$site) != 
      paste(vr.orig$focal.species, vr.orig$background.species, vr.orig$site))

# data frame to hold replaced IPMs
vr.replace <- vr.orig

# Replace positive vital rates
for(a in 1:nrow(vr.cont)) {
  # a = 1
  da = vr.cont[a,]
  fc <- as.character(da$focal.species)
  bg <- as.character(da$background.species)
  st <- as.character(da$site)
  
  # skip when no pair or no competitor
  if(is.na(da$pgr.original) | bg == "none") next
  
  # print
  print(a)
  print(paste(fc,st,bg, sep="_"))
  
  # intrinsic vital rates
  vr.intrinsic.i <- vr.orig[vr.orig$focal.species == fc & vr.orig$site == st & vr.orig$background.species == "none", ]
  
  # **** - 4.1.1 survival---------------------------------------------------
  # if invasion vital rate is "improved", replace it with intrinsic vital rate
  if(da$Survival >0) {
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"surv.int"] = vr.intrinsic.i[, "surv.int"]
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"surv.slope"] = vr.intrinsic.i[, "surv.slope"]
    print("Survival")
  }
  
  # **** - 4.1.2 growth------------------------------------------------------
  if(da$Growth >0) {
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"growth.int"] = vr.intrinsic.i[, "growth.int"]
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"growth.slope"] = vr.intrinsic.i[, "growth.slope"]
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"growth.sd"] = vr.intrinsic.i[, "growth.sd"]
    print("Growth")
  }
  
  # **** - 4.1.3 flowering---------------------------------------------------
  if(da$Flowering >0) {
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"flowering.int"] = vr.intrinsic.i[, "flowering.int"]
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"flowering.slope"] = vr.intrinsic.i[, "flowering.slope"]
    print("Flowering")
  }
  
  # **** - 4.1.4 fecundity---------------------------------------------------
  if(da$Fecundity >0) {
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"fecundity.int"] = vr.intrinsic.i[, "fecundity.int"]
    vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"fecundity.slope"] = vr.intrinsic.i[, "fecundity.slope"]
    print("Fecundity")
  }
  
  # **** - 4.1.5 establish---------------------------------------------------
    if(da$Establishment >0) {
      vr.replace[vr.replace$focal.species == fc & vr.replace$background.species == bg & vr.replace$site == st,"establishment.comp.int"] = vr.intrinsic.i[, "establishment.comp.int"]
      print("Establishment")
    }
  rm(vr.intrinsic.i)
}

# combine original and replaced IPMs
vr.orig$Positive.response = "Presence"
vr.replace$Positive.response = "Absence"

# output
# write.csv(bind_rows(vr.orig, vr.replace), "vr_replaced_20210129.csv")
vr_replaced_20210129 <- bind_rows(vr.orig, vr.replace)

#*****************************************************************************
# ** - 4.2 Calculate PGR of original and replaced IPM ----
#*****************************************************************************
vr.pert <- vr_replaced_20210129
n.cores <- 2
n.seq <- 1:4

# Note: If you want to compute the whole data set. This can take very long to compute on laptops.
# n.seq <- 1:nrow(vr.pert)

vital.rate_competition <- function(a) {
  # a = 3153
  da <- vr.pert[a,]
  
  if(is.na(da$pair) | 
     da$pair == "no" | 
     da$background.species %in% c("site","site_none","species", "species_none")) {
    lambda.a <- list(lambda = "skip")
    print(paste(a, "_skip"))
  }
  else {
    ipm <- mk.kernel(n=3000, L=as.numeric(da$L), U=as.numeric(da$U), par=da)
    if(sum(is.na(ipm$K))>0) { 
      lambda.a <- list(lambda = "null")
      print(paste(a, "_null")) 
    } 
    else {
      cal.a <- lambda.eigen(ipm$K)
      lam.a <- cal.a$lambda
      lambda.a <- list(lambda = lam.a)
      print(a)
    }
  }
  return(lambda.a)
}

# n.seq = 1
lambda.replace <- mclapply(n.seq, vital.rate_competition, mc.cores = n.cores)
lambda.replace

# output
# save(lambda.replace, file="lambda_replaced.RData")

#*****************************************************************************
# ** - 4.3 Population growth rates when positive response present vs absent ----
#*****************************************************************************
# vital rate contribution
vr.cont
vr.cont$n.positive <- rowSums(vr.cont[,c("Survival", "Growth", "Flowering", "Fecundity", "Establishment")] > 0)
vr.cont 

# vital rate replaced
vr.repl <- vr_replaced_20210129
vr.repl

# replaced lambdas
# Note: This should be modified accordingly. You can also have this data set by running Section 4.2
load("lambda_replaced.RData")
length(lambda.replace)

# combine lambda with replaced vital rates
vr.repl$pgr.replace <- as.numeric(NA)
for(i in 1:nrow(vr.repl)) {
  # i = 1
  if(lambda.replace[[i]]$lambda != "skip" & lambda.replace[[i]]$lambda != "null") {
    vr.repl[i,"pgr.replace"] = lambda.replace[[i]]$lambda
  }
}

#*****************************************************************
# Combine original and replaced PGR
vr.cont$pgr_positive_presence <- filter(vr.repl, Positive.response  == "Presence")$pgr.replace
vr.cont$pgr_positive_absence <- filter(vr.repl, Positive.response  == "Absence")$pgr.replace

# calculate if.intra and percent.change
vr.cont$pgr_percent_change <- as.numeric(NA)
vr.cont$col_pgr_positive <- as.character(NA)
for(i in 1:nrow(vr.cont)) {
  # i = 1
  di <- vr.cont[i,]
  fc <- di$focal.species
  bg <- di$background.species
  
  # intra vs. inter
  if(bg %in% c("none", "site", "site.none", "species", "species.none")) {
    vr.cont[i,"intra"] = "NA"
  }
  else if(fc == bg) {
    vr.cont[i,"intra"] = "intraspecific"
  }
  else if(fc!= bg) {
    vr.cont[i,"intra"] = "interspecific"
  }
  
  # change in lambda with/without compensation
  pgr.orig.i <- di$pgr_positive_presence
  pgr.repl.i <- di$pgr_positive_absence
  if(!is.na(pgr.orig.i) & !is.na(pgr.repl.i)) {
    # log RR 
    vr.cont[i,"pgr_percent_change"] = log(pgr.orig.i/pgr.repl.i)
    # proportional change
    #vr.cont[i,"pgr_percent_change"] = (pgr.orig.i - pgr.repl.i)/pgr.repl.i
  }
  # color
  if(!is.na(pgr.orig.i) & !is.na(pgr.repl.i)) {
    if(pgr.orig.i > 1 & pgr.repl.i >1) {vr.cont[i,"col_pgr_positive"] = "green"}
    else if (pgr.orig.i < 1 & pgr.repl.i < 1) {vr.cont[i,"col_pgr_positive"] = "blue"}
    else {vr.cont[i,"col_pgr_positive"] = "orange"}
  }
}

#*****************************************************************
# PGR in the presence and absence of positive responses
fig.positive_pgr <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  ggplot(aes(x = log(pgr_positive_presence), y = log(pgr_positive_absence), fill=col_pgr_positive)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(intercept = 0, slope=1) +
  geom_point(alpha=1, shape=21) +
  #scale_color_manual(values=c("orange" = "white", "green" = "black", "blue" = "grey")) +
  scale_fill_manual(values=c("orange" = "black", "green" = "grey65", "blue" = "white")) +
  scale_x_continuous(name=("ln(λ, compensatory responses present)")) +
  scale_y_continuous(name=("ln(λ, compensatory responses removed)"))
fig.positive_pgr

# lmer test
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  filter(pgr_percent_change > 0) %>%
  pivot_longer(cols=c("pgr_positive_absence", "pgr_positive_presence"), values_to = "Population.growth.rate", names_to="invasion.vs.replace") %>%
  lmer(log(Population.growth.rate) ~ invasion.vs.replace*site + (invasion.vs.replace|focal.species) + (invasion.vs.replace|background.species), data=.) %>%
  Anova()
  #rand()
  #summary()
  
# how many pair
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  .$col_pgr_positive %>%
  table()

#*****************************************************************
# distribution of percent changes
fig.positive_percent <- vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  ggplot(aes(y = pgr_percent_change, x=NA)) +
  geom_jitter(width=0.05, col="grey", size=0.5) +
  stat_summary(alpha=0.6, fun.data = mean_sd.ggplot, size=0.8)  +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(name=("")) +
  scale_y_continuous(name=("Proportional change in λ"))
fig.positive_percent

#*****************************************************************
# Positive response promote PGR by 22.5% 
# mean_se:  0.2256997 0.1801036 0.2712959
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  filter(pgr_percent_change > 0) %>% # include only competitive interactions
  .$pgr_percent_change %>% mean_ci_se()

#*****************************************************************
# 14/239 = 5% of PGR changed persistence
vr.cont %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(pgr.delta.rii < 0) %>% # include only competitive interactions
  filter(pgr_percent_change > 0) %>% # include only competitive interactions
  filter(!is.na(pgr_positive_presence) & !is.na(pgr_positive_absence)) %>%
  filter(pgr_positive_presence > 1 & pgr_positive_absence <1)

#*****************************************************************************
# 5. Effects of positive responses on coexistence ----
#*****************************************************************************

#*****************************************************************
# ** - 5.1 Sensitivity ---- 
#*****************************************************************
vr.cont
vr.cont$pgr.intrinsic.original <- as.numeric(NA)
vr.cont$pgr.intrinsic.replaced <- as.numeric(NA)
vr.cont$sensitivity.original <- as.numeric(NA)
vr.cont$sensitivity.replaced <- as.numeric(NA)

for(a in 1:nrow(vr.cont)) {
  # a = 1
  di <- vr.cont[a, ]
  st <- as.character(di$site)
  fc <- as.character(di$focal.species)
  bg <- as.character(di$background.species)
  
  # original PGRs
  pgr.invasion.orig.i <- di$pgr_positive_presence
  pgr.intrinsic.orig.i <- filter(vr.cont, focal.species == fc & site == st & background.species == "none")$pgr_positive_presence
  
  # PGRs after replacement
  pgr.invasion.repl.i <- di$pgr_positive_absence
  pgr.intrinsic.repl.i <- filter(vr.cont, focal.species == fc & site == st & background.species == "none")$pgr_positive_absence
  
  # skip when pgr is NA or lambda is NA
  if(is.na(pgr.intrinsic.orig.i) | is.na(pgr.invasion.orig.i) | 
     bg == "none" | 
     bg == "site" | bg == "site_none" |
     bg == "species" | bg == "species_none") next
  
  # Original PGRs
  vr.cont[a, "pgr.intrinsic.original"] <- pgr.intrinsic.orig.i
  vr.cont[a, "sensitivity.original"] <- sensitivity.igr(pgr.intrinsic.orig.i, pgr.invasion.orig.i)
  
  # Replaced PGRs
  vr.cont[a, "pgr.intrinsic.replaced"] <- pgr.intrinsic.repl.i
  vr.cont[a, "sensitivity.replaced"] <- sensitivity.igr(pgr.intrinsic.repl.i, pgr.invasion.repl.i)
  
  print(a)
}

#*****************************************************************
# ** - 5.2 Calculate competitive outcomes ---- 
#*****************************************************************
# data frame to hold outcomes
# Note: This should be modified accordingly.
outcome.original <- read_excel("Coexistence_null.xlsx")
outcome.original <- outcome.original[-1,]

# add two columns
outcome.original$igr12.percent.change <- as.numeric(NA)
outcome.original$igr21.percent.change <- as.numeric(NA)
outcome.original

# outcome using replaced IPMs
outcome.replaced <- outcome.original
outcome.replaced

# predict outcomes
for(a in 1:nrow(outcome.original)) {
  # a = 3
  di <- outcome.original[a, ]
  sps1 <- di$sps1
  sps2 <- di$sps2
  st <- as.character(di$site)
  type <- di$type
  origin1 <- di$origin.sps1
  origin2 <- di$origin.sps2
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # original IPMs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # percernt change of PGR
  pgr.percent12 <- filter(vr.cont, focal.species == sps1 & background.species == sps2 & site == st)$pgr_percent_change
  pgr.percent21 <- filter(vr.cont, focal.species == sps2 & background.species == sps1 & site == st)$pgr_percent_change
  
  # intrinsic PGRs
  igr1.orig <- filter(vr.cont, focal.species == sps1 & background.species == "none" & site == st)$pgr_positive_presence
  igr2.orig <- filter(vr.cont, focal.species == sps2 & background.species == "none" & site == st)$pgr_positive_presence
  # invaion PGRs
  igr12.orig <- filter(vr.cont, focal.species == sps1 & background.species == sps2 & site == st)$pgr_positive_presence
  igr21.orig <- filter(vr.cont, focal.species == sps2 & background.species == sps1 & site == st)$pgr_positive_presence
  # sensitivities
  s12.orig <- filter(vr.cont, focal.species == sps1 & background.species == sps2 & site == st)$sensitivity.original
  s21.orig <- filter(vr.cont, focal.species == sps2 & background.species == sps1 & site == st)$sensitivity.original
  
  outcome.original[a,"igr12.percent.change"] <- pgr.percent12
  outcome.original[a,"igr21.percent.change"] <- pgr.percent21
  outcome.original[a,"lambda1"] <- igr1.orig
  outcome.original[a,"lambda2"] <- igr2.orig
  outcome.original[a,"igr12"] <- igr12.orig
  outcome.original[a,"igr21"] <- igr21.orig
  outcome.original[a,"sensitivity.12"] <- s12.orig
  outcome.original[a,"sensitivity.21"] <- s21.orig
  
  # skip when sensitivity is NA
  if(is.na(s12.orig) | is.na(s21.orig)) { next; print("skip missing sensitivity")}
  
  # skip facilitation
  if((s12.orig < 0) | (s21.orig < 0)) { next; print("skip facilitation") }
  
  # skip when IGR is NA
  if(is.na(igr1.orig) | is.na(igr2.orig) | is.na(igr12.orig) | is.na(igr21.orig)) next
  
  # score negative sensitivty
  #if(s12.orig < 0) { s12.orig = 0.1; print("original facilitation scored") }
  #if(s21.orig < 0) { s21.orig = 0.1; print("original facilitation scored") }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # coexistence by IGR
  si <- coex.igr(sps1, sps2, igr1.orig, igr2.orig, igr12.orig, igr21.orig)
  
  outcome.original[a, "outcome.igr"] <- si$outcome[1]
  outcome.original[a, "winner.igr"] <- si$winner[1]
  outcome.original[a, "superior.igr"] <- si$superior[1]
  si <- NULL
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # coexistence by NDFD
  ni <- coex.ndfd(s12.orig, s21.orig)
  outcome.original[a,"nd"] <- ni$ndfd[1]
  outcome.original[a,"fd"] <- ni$ndfd[2]
  outcome.original[a, "coexistence.metric"] <- ni$ndfd[3] # coexistence metric
  
  outcome.original[a, "outcome.ndfd"] <- ni$outcome[1]
  outcome.original[a, "winner.ndfd"] <- ni$winner[1]
  outcome.original[a, "superior.ndfd"] <- ni$superior[1]
  
  # FD highland/lowland
  if(origin1 == "Lowland" & origin2 == "Highland") outcome.original[a, "fd.highlow"] <- ni$ndfd[2]
  else outcome.original[a, "fd.highlow"] <- 1/ni$ndfd[2]
  
  # superior origin
  if(ni$superior[1] == "sps1") outcome.original[a, "origin.superior"] <- origin1
  else outcome.original[a, "origin.superior"] <- origin2
  
  ni <- NULL
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # replaced IPMs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # intrinsic PGRs
  igr1.repl <- filter(vr.cont, focal.species == sps1 & background.species == "none" & site == st)$pgr_positive_absence
  igr2.repl <- filter(vr.cont, focal.species == sps2 & background.species == "none" & site == st)$pgr_positive_absence
  # invaion PGRs
  igr12.repl <- filter(vr.cont, focal.species == sps1 & background.species == sps2 & site == st)$pgr_positive_absence
  igr21.repl <- filter(vr.cont, focal.species == sps2 & background.species == sps1 & site == st)$pgr_positive_absence
  # sensitivities
  s12.repl<- filter(vr.cont, focal.species == sps1 & background.species == sps2 & site == st)$sensitivity.replaced
  s21.repl <- filter(vr.cont, focal.species == sps2 & background.species == sps1 & site == st)$sensitivity.replaced
  
  outcome.replaced[a,"igr12.percent.change"] <- pgr.percent12
  outcome.replaced[a,"igr21.percent.change"] <- pgr.percent21
  outcome.replaced[a,"lambda1"] <- igr1.repl
  outcome.replaced[a,"lambda2"] <- igr2.repl 
  outcome.replaced[a,"igr12"] <- igr12.repl
  outcome.replaced[a,"igr21"] <- igr21.repl
  outcome.replaced[a,"sensitivity.12"] <- s12.repl
  outcome.replaced[a,"sensitivity.21"] <- s21.repl

  # score negative sensitivity
  #if(s12.repl < 0.1) { s12.repl = 0.1; print("replaced facilitation scored") }
  #if(s21.repl < 0.1) { s21.repl = 0.1; print("replaced facilitation scored") }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # coexistence by IGR
  si <- coex.igr(sps1, sps2, igr1.repl, igr2.repl, igr12.repl, igr21.repl)
  
  outcome.replaced[a, "outcome.igr"] <- si$outcome[1]
  outcome.replaced[a, "winner.igr"] <- si$winner[1]
  outcome.replaced[a, "superior.igr"] <- si$superior[1]
  si <- NULL
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # coexistence by NDFD
  ni <- coex.ndfd(s12.repl, s21.repl)
  outcome.replaced[a,"nd"] <- ni$ndfd[1]
  outcome.replaced[a,"fd"] <- ni$ndfd[2]
  outcome.replaced[a, "coexistence.metric"] <- ni$ndfd[3] # coexistence metric
  
  outcome.replaced[a, "outcome.ndfd"] <- ni$outcome[1]
  outcome.replaced[a, "winner.ndfd"] <- ni$winner[1]
  outcome.replaced[a, "superior.ndfd"] <- ni$superior[1]
  
  # FD highland/lowland
  if(origin1 == "Lowland" & origin2 == "Highland") outcome.replaced[a, "fd.highlow"] <- ni$ndfd[2]
  else outcome.replaced[a, "fd.highlow"] <- 1/ni$ndfd[2]
  
  # superior origin
  if(ni$superior[1] == "sps1") outcome.replaced[a, "origin.superior"] <- origin1
  else outcome.replaced[a, "origin.superior"] <- origin2
  
  ni <- NULL
  
  print(a)
}

# combine data
outcome.original$Compensatory.response = "Present"
outcome.replaced$Compensatory.response = "Absent"
outcome.long <- bind_rows(outcome.original, outcome.replaced)

#*****************************************************************
# ** - 5.3 Compare IGR outcomes ---- 
#*****************************************************************
# plot
fig.outcomes_igr <- outcome.long %>%
  filter(!(igr12.percent.change == 0 & igr21.percent.change == 0)) %>%
  mutate(`Compensatory response` = Compensatory.response)  %>%
  ggplot(aes(x=log(igr12), y=log(igr21), col=`Compensatory response`)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data=data.frame(igr12.original = log(outcome.original$igr12),
                               igr21.original = log(outcome.original$igr21),
                               igr12.replaced = log(outcome.replaced$igr12),
                               igr21.replaced = log(outcome.replaced$igr21)),
               aes(x=igr12.replaced, y=igr21.replaced, xend=igr12.original, yend=igr21.original), size=0.4, col="grey40", arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_manual(values=c("Absent" = "#FB9A06FF", "Present" = "#3E4A89FF")) +
  scale_x_continuous(name="ln(Invasion growth rate of sp. 1)") +
  scale_y_continuous(name="ln(Invasion growth rate of sp. 2)") 
fig.outcomes_igr

# how many pairs shift outcomes?
# 8/107 pairs = 7.4% pairs change the outcomes
data.frame(outcome.igr.original = filter(outcome.original, !is.na(outcome.igr) & !(igr12.percent.change == 0 & igr21.percent.change == 0))$outcome.igr,
           outcome.igr.replaced = filter(outcome.replaced, !is.na(outcome.igr) & !(igr12.percent.change == 0 & igr21.percent.change == 0))$outcome.igr) %>%
  filter(!is.na(outcome.igr.original) & !is.na(outcome.igr.replaced)) %>%
  filter(outcome.igr.original != outcome.igr.replaced)

#*****************************************************************
# ** - 5.4 Compare NDFD outcomes ---- 
#*****************************************************************

#*****************************************************************
# ** - 5.4.1 NDFD outcomes ---- 
#*****************************************************************
# plot
fig.outcomes_ndfd <- outcome.long %>%
  mutate(`Compensatory response` = Compensatory.response)  %>%
  filter(!(igr12.percent.change == 0 & igr21.percent.change == 0)) %>% 
  filter(sensitivity.12 > 0 & sensitivity.21 > 0) %>%
  ggplot(aes(x=log(1-nd), y=log(fd.highlow), col=`Compensatory response`)) +
  geom_polygon(data=polygan.coex.no(x1=-2.7), aes(x=x,y=y), col="black",fill="grey90", inherit.aes = FALSE) + 
  geom_polygon(data=polygan.prio.no(x2=2.7), aes(x=x,y=y), col="black", fill="grey80", inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", col="grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", col="grey30") +
  geom_point(alpha=0.8) +
  geom_segment(data=data.frame(nd.original = log(1-outcome.original$nd),
                               fd.original = log(outcome.original$fd.highlow),
                             nd.replaced = log(1-outcome.replaced$nd),
                             fd.replaced = log(outcome.replaced$fd.highlow)),
           aes(x=nd.replaced, y=fd.replaced, xend=nd.original, yend=fd.original), size=0.4, col="black", arrow = arrow(length = unit(0.2, "cm")), alpha=0.8) +
  scale_color_manual(values=c("Absent" = "#FB9A06FF", "Present" = "#3E4A89FF")) +
  coord_cartesian(xlim=c(-2.7,2.5), ylim=c(-1,2)) +
  scale_x_continuous(name="ln(Niche overlap)", expand = c(0,0)) +
  scale_y_continuous(name="ln(Relative fitness difference)", expand = c(0, 0)) 
fig.outcomes_ndfd

# how many pairs shift outcomes?
# 8/107 pairs = 7.4% pairs change the outcomes
data.frame(outcome.igr.original = outcome.original$outcome.igr,
           outcome.igr.replaced = outcome.replaced$outcome.igr) %>%
  filter(!is.na(outcome.igr.original) & !is.na(outcome.igr.replaced)) %>%
  filter(outcome.igr.original != outcome.igr.replaced)

#*****************************************************************
# ** - 5.4.2 Change in NO ---- 
#*****************************************************************
# changes in NO
no.original <- 1 - outcome.original$nd
no.replaced <- 1 -outcome.replaced$nd
dif.nd.log <- log(no.original) - log(no.replaced)
dif.nd.percent <- (no.replaced - no.original)/no.original

mean_sd.ggplot(dif.nd.percent)
mean(dif.nd.percent, na.rm=TRUE)
mean_se(dif.nd.percent)
sqrt(var(dif.nd.percent, na.rm=TRUE))

# 0.2575333 0.1860385 0.3290282
mean_ci_se(dif.nd.percent[!is.na(dif.nd.percent)])

# plot
fig.no <- data.frame(no.original = log(no.original),
                     no.replaced = log(no.replaced)) %>%
  filter(!is.na(no.original)) %>%
  mutate(no.change = no.original - no.replaced) %>%
  #filter(no.change > 0)
  ggplot(aes(x=no.change, y=NA)) +
  geom_vline(xintercept = 0, linetype = "dashed", col="grey30") +
  geom_jitter(alpha=0.8, height=0.1, width=0, col="grey50") +
  geom_boxplot(outlier.alpha = 0, width=0.6, alpha=0.6) +
  scale_x_continuous(name=("Change in niche overlap")) +
  scale_y_discrete(name=(" "))
fig.no

# t.test
# NO
data.frame(dif = dif.nd.log, coexistence.index = "Niche overlap") %>%
  filter(!is.na(dif)) %>%
  .$ dif %>%
  t.test()

# percent change in NO
-0.2036854/mean(log(no.original), na.rm=TRUE)
-0.1609136/mean(log(no.original), na.rm=TRUE)
-0.2464571/mean(log(no.original), na.rm=TRUE)

#*****************************************************************
# ** - 5.4.3 Change in RFD ----
#*****************************************************************
fd.original <- abs(log(outcome.original$fd.highlow))
fd.replaced <- abs(log(outcome.replaced$fd.highlow))
dif.fd.log <- fd.replaced - fd.original
dif.fd.percent <- (fd.replaced - fd.original)/fd.original
mean_ci_se(dif.fd.percent[!is.na(dif.fd.percent)])

# 0.2575333 0.1860385 0.3290282

fig.fd <- data.frame(fd.original = abs(log(outcome.original$fd.highlow)),
                     fd.replaced = abs(log(outcome.replaced$fd.highlow))) %>%
  mutate(fd.change = fd.original - fd.replaced) %>%
  ggplot(aes(y=fd.change, x=NA)) +
  geom_hline(yintercept = 0, linetype = "dashed", col="grey30") +
  geom_jitter(alpha=0.8, height=0, width=0.1, col="grey50") +
  geom_boxplot(outlier.alpha = 0, width=0.6, alpha=0.6) +
  scale_x_discrete(name=(" ")) +
  scale_y_continuous(name=("Change in relative fitness difference"))
fig.fd

# t.test
data.frame(fd.original = abs(log(outcome.original$fd.highlow)),
           fd.replaced = abs(log(outcome.replaced$fd.highlow))) %>%
  mutate(fd.change = fd.original - fd.replaced) %>%
  filter(fd.change != 0) %>% 
  .$fd.change %>%
  t.test()

# percent change in FD
0.04727371/mean(fd.original, na.rm=TRUE)
0.01566570/mean(fd.original, na.rm=TRUE)
0.07888172/mean(fd.original, na.rm=TRUE)

# how many pairs change domiance
sum(log(outcome.original$fd.highlow) > 0 & log(outcome.replaced$fd.highlow) < 0, na.rm=TRUE)
sum(log(outcome.original$fd.highlow) < 0 & log(outcome.replaced$fd.highlow) > 0, na.rm=TRUE)

#*****************************************************************************
# Figures ----
#*****************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ** - Figures 1 ----
# Vital rate controbution across focal species
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# overall
p <- fig.contribution_vital.rate + 
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none") 
p

ggsave(p, filename = "fig.1.pdf", device = cairo_pdf, 
       width = 4.63, height = 4.35, units = "in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ** - Figures 2 ----
# PGR and outcomes of competition after replacing “enhanced” vital rates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# by focal species
p <- fig.contribution_focal_vital.rate + 
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text.x=element_text(size=9, angle = 90, vjust=0.3, hjust=1),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=11),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=10, hjust=0)) 
p
ggsave(p, filename = "fig.2a.pdf", device = cairo_pdf, 
       width = 8.41, height = 2.51, units = "in")

# by site
p <- fig.contribution_site_vital.rate + 
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=11),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=10, hjust=0)) 
p
ggsave(p, filename = "fig.2b.pdf", device = cairo_pdf, 
       width = 8.41, height = 2.11, units = "in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ** - Figures 3 ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# outcomes
p <- fig.outcomes_ndfd +
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    axis.line = element_line(colour = "black"),

    legend.position = c(0.35,0.93),
    legend.key.size = unit(1,"lines"),
    legend.text = element_text(size=9),
    legend.title = element_text(size=10),
    
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=12, hjust=0)) 

p 
ggsave(p, filename = "fig.3.pdf", device = cairo_pdf, 
       width = 5, height = 5, units = "in")
# NO
p <- fig.no + 
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=12, hjust=0)) 
p

ggsave(p, filename = "fig.3_no.pdf", device = cairo_pdf, 
       width = 5, height = 1.2, units = "in")

# FD
p <- fig.fd + 
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=12, hjust=0)) 
p
ggsave(p, filename = "fig.3_fd.pdf", device = cairo_pdf, 
       width = 1.4, height = 5, units = "in")