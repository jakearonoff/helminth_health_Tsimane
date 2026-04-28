library("tidyverse")
library("plm")
library("jtools")
library("stargazer")
library("mediation")

# ldl = lipids data
# glu = glucose data 
# bp = blood pressure data
# bmi = BMI data
# weight = weight data
# hb = hemoglobin data

# pid = person id 


#####################################################################

## Panel models of infection predicting cardiometabolic health

#####################################################################

m1 <- plm(ldl ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = ldl)
m2 <- plm(ldl ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = ldl)
m3 <- plm(hdl ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(hdl ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(trig ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = ldl)
m6 <- plm(trig ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = ldl)
m7 <- plm(glucose ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = glu)
m8 <- plm(glucose ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = glu)


stargazer(m1, m2, m3, m4, m5, m6, m7, m8, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_lipglu_within.html")



m1 <- plm(SBP ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = bp)
m2 <- plm(SBP ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = bp)
m3 <- plm(DBP ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = bp)
m4 <- plm(DBP ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = bp)
m5 <- plm(BMI ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = bmi)
m6 <- plm(BMI ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = bmi)
m7 <- plm(weight ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = weight)
m8 <- plm(weight ~ years + helminth,
          index = c("pid"),
          model = "within",
          data = weight)
m9 <- plm(hb ~ years + hookworm,
          index = c("pid"),
          model = "within",
          data = hb)
m10 <- plm(hb ~ years + helminth,
           index = c("pid"),
           model = "within",
           data = hb)

stargazer(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10,
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_bpbmihb_within.html")



#####################################################################

## Figure 1

#####################################################################

ldl <- ungroup(ldl) %>% mutate(ldlw = ldl/mean(ldl)*100)
ldl <- ungroup(ldl) %>% mutate(hdlw = hdl/mean(hdl)*100)
ldl <- ungroup(ldl) %>% mutate(trigw = trig/mean(trig)*100)
glu <- ungroup(glu) %>% mutate(glucosew = glucose/mean(glucose)*100)
bp <- ungroup(bp) %>% mutate(SBPw = SBP/mean(SBP)*100)
bp <- ungroup(bp) %>% mutate(DBPw = DBP/mean(DBP)*100)
bmi <- ungroup(bmi) %>% mutate(BMIw = BMI/mean(BMI)*100)
hb <- ungroup(hb) %>% mutate(hbw = hb/mean(hb)*100)
weight <- ungroup(weight) %>% mutate(weightw = weight/mean(weight)*100)


ldl$infection <- ldl$hookworm
glu$infection <- glu$hookworm

m1 <- plm(ldlw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = ldl)
m3 <- plm(hdlw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(trigw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = ldl)
m7 <- plm(glucosew ~ years + infection,
          index = c("pid"),
          model = "within",
          data = glu)

bp$infection <- bp$hookworm
hb$infection <- hb$hookworm

m9 <- plm(SBPw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = bp)
m11 <- plm(DBPw ~ years + infection,
           index = c("pid"),
           model = "within",
           data = bp)
m13 <- plm(hbw ~ years + infection,
           index = c("pid"),
           model = "within",
           data = hb)


ldl$infection <- ldl$helminth
glu$infection <- glu$helminth

m2 <- plm(ldlw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(hdlw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = ldl)
m6 <- plm(trigw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = ldl)
m8 <- plm(glucosew ~ years + infection,
          index = c("pid"),
          model = "within",
          data = glu)


bp$infection <- bp$helminth
hb$infection <- hb$helminth

m10 <- plm(SBPw ~ years + infection,
           index = c("pid"),
           model = "within",
           data = bp)
m12 <- plm(DBPw ~ years + infection,
           index = c("pid"),
           model = "within",
           data = bp)
m14 <- plm(hbw ~ years + infection,
           index = c("pid"),
           model = "within",
           data = hb)


# plot settings
apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Helvetica'),
        legend.title=element_blank(), 
        axis.text=element_text(size=15),
        axis.title=element_text(size = 15),
        legend.text = element_text(size = 15))


coef_names <- c(
  "Infection" = "infection")


plot <- plot_summs(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14,
                   coefs = coef_names,
                   model.names = c("LDL (hookworm)", "LDL (any helminth)", "HDL (hookworm)", "HDL (any helminth)", 
                                   "Trig (hookworm)", "Trig (any helminth)", "Glu (hookworm)", "Glu (any helminth)",
                                   "SBP (hookworm)", "SBP (any helminth)", "DBP (hookworm)", "DBP (any helminth)", 
                                   "Hb (hookworm)", "Hb (any helminth)"),
                   position = position_dodge(width = 1.3),
                   colors = c("brown1", "brown4", "orchid", "orchid4", "orange", "orange4",
                              "dodgerblue", "dodgerblue4", "green", "green4", "gray", "gray1", "brown2", "red4"),
                   point.shape = c(21,22,23,24,25,22,16,17,18,3,4,8,18,20))
plot <- plot  + apatheme + 
  labs(x = "% change with infection", y = "Infection") + 
  theme_classic() + theme(axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 16),
                          axis.text.y = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 14),
                          axis.text.x = element_text(size = 14),
                          legend.key.spacing.y = unit(0.175, "cm")) + 
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.15))) + 
  guides(linetype = "none")
plot 

ggsave("/filepath/hookhelminthresults.png", plot = plot, dpi = "retina",
       height = 5, width = 8)



#####################################################################

## Between-person models, testing associations between 
##   average infection and average BMI/weight

#####################################################################

bmim <- group_by(bmi, pid) %>% mutate(helminth = mean(helminth),
                                      hookworm = mean(hookworm),
                                      age = mean(age),
                                      years = mean(years),
                                      BMI = mean(BMI)) %>% filter(!duplicated(pid))

weightm <- group_by(weight, pid) %>% mutate(helminth = mean(helminth),
                                            hookworm = mean(hookworm),
                                            age = mean(age),
                                            years = mean(years),
                                            weight = mean(weight)) %>% filter(!duplicated(pid))


m1 <- lm(BMI ~ age + female + years + route_distance_townp + hookworm, data = bmim)
m2 <- lm(BMI ~ age + female + years + route_distance_townp + helminth, data = bmim)

m3 <- lm(weight ~ age + female + years + route_distance_townp + hookworm, data = weightm)
m4 <- lm(weight ~ age + female + years + route_distance_townp + helminth, data = weightm)

stargazer(m1, m2, m3, m4,
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_bmiweight_between.html")



#####################################################################

## Figure 2

#####################################################################

bmim <- ungroup(bmim) %>% mutate(BMIw = BMI/mean(BMI)*100)
weightm <- ungroup(weightm) %>% mutate(weightw = weight/mean(weight)*100)


bmi$infection <- bmi$hookworm
weight$infection <- weight$hookworm

bmim$infection <- bmim$hookworm
weightm$infection <- weightm$hookworm


m1 <- plm(BMIw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = bmi)
m5 <- plm(weightw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = weight)

m3 <- lm(BMIw ~ age + male + time + route_distance_townp + infection, data = bmim)
m7 <- lm(weightw ~ age + male + time + route_distance_townp + infection, data = weightm)


bmi$infection <- bmi$helminth
weight$infection <- weight$helminth

bmim$infection <- bmim$helminth
weightm$infection <- weightm$helminth


m2 <- plm(BMIw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = bmi)
m6 <- plm(weightw ~ years + infection,
          index = c("pid"),
          model = "within",
          data = weight)

m4 <- lm(BMIw ~ age + male + time + route_distance_townp + infection, data = bmim)
m8 <- lm(weightw ~ age + male + time + route_distance_townp + infection, data = weightm)


plot <- plot_summs(m1, m2, m3, m4, m5, m6, m7, m8,
                   coefs = coef_names,
                   model.names = c("BMI (hookworm, \n    within-person)", "BMI (any helminth, \n    within-person)", 
                                   "BMI (hookworm, \n     between-person)", "BMI (any helminth, \n    between-person)", 
                                   "Weight (hookworm, \n     within-person)", "Weight (any helminth, \n     within-person)", 
                                   "Weight (hookworm, \n     between-person)", "Weight (any helminth, \n     between-person)"),
                   position = position_dodge(width = 1.3),
                   colors = c("dodgerblue", "dodgerblue4", "brown1","brown4",
                              "gray", "gray1", "green", "green4"),
                   point.shape = c(16,15,19,17,4,3,16,15))
plot <- plot  + apatheme + 
  labs(x = "% change with infection", y = "Infection") + 
  theme_classic() + theme(axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 16),
                          axis.text.y = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 14),
                          axis.text.x = element_text(size = 14),
                          legend.key.spacing.y = unit(0.38, "cm")) + 
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.15))) + 
  guides(linetype = "none")
plot 

ggsave("/filepath/hookhelminthresults_bmiweight.png", plot = plot, dpi = "retina",
       height = 5, width = 8)



#####################################################################

## Panel models of infection predicting cardiometabolic health,
##  including eosinophils 

#####################################################################

m1 <- plm(ldl ~ years + hookworm + eosinophils,
          index = c("pid"),
          model = "within",
          data = ldl)
m2 <- plm(ldl ~ years + helminth + eosinophils,
          index = c("pid"),
          model = "within",
          data = ldl)
m3 <- plm(hdl ~ years + hookworm + eosinophils,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(hdl ~ years + helminth + eosinophils,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(trig ~ years + hookworm + eosinophils,
          index = c("pid"),
          model = "within",
          data = ldl)
m6 <- plm(trig ~ years + helminth + eosinophils,
          index = c("pid"),
          model = "within",
          data = ldl)
m7 <- plm(glucose ~ years + hookworm + eosinophils,
          index = c("pid"),
          model = "within",
          data = glu)
m8 <- plm(glucose ~ years + helminth + eosinophils,
          index = c("pid"),
          model = "within",
          data = glu)

stargazer(m1, m2, m3, m4, m5, m6, m7, m8, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_lipglueos_within.html")



#####################################################################

## Panel models of infection predicting systolic blood pressure,
##  including lipids and hemoglobin 

#####################################################################

m1 <- plm(SBP ~ years + hookworm + ldl + trig,
          index = c("pid"),
          model = "within",
          data = bp)
m2 <- plm(SBP ~ years + helminth + ldl + trig,
          index = c("pid"),
          model = "within",
          data = bp)
m3 <- plm(SBP ~ years + hookworm + hb,
          index = c("pid"),
          model = "within",
          data = bp)
m4 <- plm(SBP ~ years + helminth + hb,
          index = c("pid"),
          model = "within",
          data = bp)

stargazer(m1, m2, m3, m4,
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_sbp_lipids_hb_within.html")



#####################################################################

## Panel models of infection predicting cardiometabolic health,
##  interaction by APOE-e4 allele carrier (car4 = e4 carrier)

#####################################################################

m1 <- plm(ldl ~ years + hookworm*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m2 <- plm(ldl ~ years + helminth*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m3 <- plm(hdl ~ years + hookworm*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(hdl ~ years + helminth*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(trig ~ years + hookworm*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m6 <- plm(trig ~ years + helminth*car4,
          index = c("pid"),
          model = "within",
          data = ldl)


stargazer(m1, m2, m3, m4, m5, m6, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_lipidscar4_within.html")



#####################################################################

## Figure 3

#####################################################################

helminthldl4 <- filter(ldl, car4 == 1)
helminthldl3 <- filter(ldl, car4 == 0)

helminthldl3$infection <- helminthldl3$hookworm
helminthldl4$infection <- helminthldl4$hookworm

hookldl3 <- plm(ldl ~ years + infection,
                index = c("pid"),
                model = "within",
                data = helminthldl3)
hookldl4 <- plm(ldl ~ years + infection,
                index = c("pid"),
                model = "within",
                data = helminthldl4)


helminthldl3$infection <- helminthldl3$helminth
helminthldl4$infection <- helminthldl4$helminth

helmldl3 <- plm(ldl ~ years + infection,
                index = c("pid"),
                model = "within",
                data = helminthldl3)
helmldl4 <- plm(ldl ~ years + infection,
                index = c("pid"),
                model = "within",
                data = helminthldl4)


plot <- plot_summs(hookldl3, helmldl3, hookldl4, helmldl4, 
                   coefs = coef_names,
                   model.names = c("ε3/ε3 (Hookworm)", "ε3/ε3 (Any Helminth)", "ε4 Carrier (Hookworm)", "ε4 Carrier (Any Helminth)"),
                   colors = c("brown1", "brown4", "dodgerblue", "dodgerblue4"))
plot <- plot  + apatheme + 
  labs(x = "              LDL (mg/dL) change with infection", y = "Infection") + 
  theme_classic() + theme(axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 16),
                          axis.text.y = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 14),
                          axis.text.x = element_text(size = 14),
                          legend.key.spacing.y = unit(.55, "cm"))
plot

ggsave("/filepath/hookhelminthldlapoe.png", plot = plot, dpi = "retina",
       height = 5, width = 8)


#####################################################################

## Conditional logit models assessing within-person changes in
##  infection odds over time

#####################################################################

m1 <- clogit(hookworm ~ years + strata(pid), data = ldl)
m2 <- clogit(helminth ~ years + strata(pid), data = ldl)
m3 <- clogit(hookworm ~ years + strata(pid), data = glu)
m4 <- clogit(helminth ~ years + strata(pid), data = glu)
m5 <- clogit(hookworm ~ years + strata(pid), data = bp)
m6 <- clogit(helminth ~ years + strata(pid), data = bp)

stargazer(m1, m2, m3, m4, m5, m6,
          coef = list(exp(coef(m1)), exp(coef(m2)), exp(coef(m3)), exp(coef(m4)), exp(coef(m5)), exp(coef(m6))),
          ci = TRUE, ci.level = 0.95, 
          ci.custom = list(exp(confint(m1)), exp(confint(m2)), exp(confint(m3)), exp(confint(m4)), 
                           exp(confint(m5)), exp(confint(m6))), 
          p = list(summary(m1)$coefficients[, 5], summary(m2)$coefficients[, 5], summary(m3)$coefficients[, 5], 
                   summary(m4)$coefficients[, 5], summary(m5)$coefficients[, 5], summary(m6)$coefficients[, 5]),
          align = T, single.row = T, digits = 2, star.cutoffs = c(0.10, 0.05, 0.01),
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelminth_time_within.html")


# restricting to first decade of the study period (2004-2015)
m1 <- clogit(hookworm ~ years + strata(pid), data = ldl[ldl$decade == 0,])
m2 <- clogit(helminth ~ years + strata(pid), data = ldl[ldl$decade == 0,])
m3 <- clogit(hookworm ~ years + strata(pid), data = glu[glu$decade == 0,])
m4 <- clogit(helminth ~ years + strata(pid), data = glu[glu$decade == 0,])
m5 <- clogit(hookworm ~ years + strata(pid), data = bp[bp$decade == 0,])
m6 <- clogit(helminth ~ years + strata(pid), data = bp[bp$decade == 0,])

stargazer(m1, m2, m3, m4, m5, m6, 
          coef = list(exp(coef(m1)), exp(coef(m2)), exp(coef(m3)), exp(coef(m4)), exp(coef(m5)), exp(coef(m6))),
          ci = TRUE, ci.level = 0.95, 
          ci.custom = list(exp(confint(m1)), exp(confint(m2)), exp(confint(m3)), exp(confint(m4)), 
                           exp(confint(m5)), exp(confint(m6))), 
          p = list(summary(m1)$coefficients[, 5], summary(m2)$coefficients[, 5], summary(m3)$coefficients[, 5], 
                   summary(m4)$coefficients[, 5], summary(m5)$coefficients[, 5], summary(m6)$coefficients[, 5]),
          align = T, single.row = T, digits = 2, star.cutoffs = c(0.10, 0.05, 0.01),
          star.char = c("t", "*", "**"), 
          out = "/filepath/helminth_time_firstdecade_within.html")


# restricting to last decade of the study period (2016-2025)
m1 <- clogit(hookworm ~ years + strata(pid), data = ldl[ldl$decade == 1,])
m2 <- clogit(helminth ~ years + strata(pid), data = ldl[ldl$decade == 1,])
m3 <- clogit(hookworm ~ years + strata(pid), data = glu[glu$decade == 1,])
m4 <- clogit(helminth ~ years + strata(pid), data = glu[glu$decade == 1,])
m5 <- clogit(hookworm ~ years + strata(pid), data = bp[bp$decade == 1,])
m6 <- clogit(helminth ~ years + strata(pid), data = bp[bp$decade == 1,])

stargazer(m1, m2, m3, m4, m5, m6, 
          coef = list(exp(coef(m1)), exp(coef(m2)), exp(coef(m3)), exp(coef(m4)), exp(coef(m5)), exp(coef(m6))),
          ci = TRUE, ci.level = 0.95, 
          ci.custom = list(exp(confint(m1)), exp(confint(m2)), exp(confint(m3)), exp(confint(m4)), 
                           exp(confint(m5)), exp(confint(m6))), 
          p = list(summary(m1)$coefficients[, 5], summary(m2)$coefficients[, 5], summary(m3)$coefficients[, 5], 
                   summary(m4)$coefficients[, 5], summary(m5)$coefficients[, 5], summary(m6)$coefficients[, 5]),
          align = T, single.row = T, digits = 2, star.cutoffs = c(0.10, 0.05, 0.01),
          star.char = c("t", "*", "**"), 
          out = "/filepath/helminth_time_lastdecade_within.html")



#####################################################################

## Within-person mediation analysis

#####################################################################

# creating variables for within-person changes 
#  (subtracting person-specific means)
ldl$ldlw <- ldl$ldl - ave(ldl$ldl, ldl$pid)
ldl$trigw <- ldl$trig - ave(ldl$trig, ldl$pid)
ldl$timew <- ldl$years - ave(ldl$years, ldl$pid)
ldl$helminthw <- ldl$helminth - ave(ldl$helminth, ldl$pid)
ldl$hookwormw <- ldl$hookworm - ave(ldl$hookworm, ldl$pid)

glu$glucosew <- glu$glucose - ave(glu$glucose, glu$pid)
glu$timew <- glu$years - ave(glu$years, glu$pid)
glu$helminthw <- glu$helminth - ave(glu$helminth, glu$pid)
glu$hookwormw <- glu$hookworm - ave(glu$hookworm, glu$pid)

bp$SBPw <- bp$SBP - ave(bp$SBP, bp$pid)
bp$timew <- bp$years - ave(bp$years, bp$pid)
bp$helminthw <- bp$helminth - ave(bp$helminth, bp$pid)
bp$hookwormw <- bp$hookworm - ave(bp$hookworm, bp$pid)


# mediation models 
# time --> reduced hookworm infection --> increased LDL
out <- lm(ldlw ~ timew + hookwormw,
          data = ldl[ldl$car4 == 0 & ldl$decade == 1,]) # data = no APOE-e4 carriers, and last decade (2016-2025)
med <- lm(hookwormw ~ timew, 
          data = ldl[ldl$car4 == 0 & ldl$decade == 1,]) # data = no APOE-e4 carriers, and last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "hookwormw")
summary(medresult)

# time --> reduced helminth infection --> increased LDL
out <- lm(ldlw ~ timew + helminthw,
          data = ldl[ldl$car4 == 0 & ldl$decade == 1,]) # data = no APOE-e4 carriers, and last decade (2016-2025)
med <- lm(helminthw ~ timew, 
          data = ldl[ldl$car4 == 0 & ldl$decade == 1,]) # data = no APOE-e4 carriers, and last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "helminthw")
summary(medresult)


# time --> reduced hookworm infection --> increased triglycerides 
out <- lm(trigw ~ timew + hookwormw,
          data = ldl[ldl$decade == 1,]) # data = last decade (2016-2025)
med <- lm(hookwormw ~ timew, 
          data = ldl[ldl$decade == 1,]) # data = last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "hookwormw")
summary(medresult)

# time --> reduced helminth infection --> increased triglycerides
out <- lm(trigw ~ timew + helminthw,
          data = ldl[ldl$decade == 1,]) # data = last decade (2016-2025)
med <- lm(helminthw ~ timew, 
          data = ldl[ldl$decade == 1,]) # data = last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "helminthw")
summary(medresult)


# time --> reduced hookworm infection --> increased glucose 
out <- lm(glucosew ~ timew + hookwormw,
          data = glu[glu$decade == 1,]) # data = last decade (2016-2025)
med <- lm(hookwormw ~ timew, 
          data = glu[glu$decade == 1,]) # data = last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "hookwormw")
summary(medresult)

# time --> reduced helminth infection --> increased glucose
out <- lm(glucosew ~ timew + helminthw,
          data = glu[glu$decade == 1,]) # data = last decade (2016-2025)
med <- lm(helminthw ~ timew, 
          data = glu[glu$decade == 1,]) # data = last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "helminthw")
summary(medresult)


# time --> reduced hookworm infection --> increased SBP 
out <- lm(SBPw ~ timew + hookwormw,
          data = bp[bp$decade == 1,]) # data = last decade (2016-2025)
med <- lm(hookwormw ~ timew, 
          data = bp[bp$decade == 1,]) # data = last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "hookwormw")
summary(medresult)


# time --> reduced helminth infection --> increased SBP
out <- lm(SBPw ~ timew + helminthw,
          data = bp[bp$decade == 1,]) # data = last decade (2016-2025)
med <- lm(helminthw ~ timew, 
          data = bp[bp$decade == 1,]) # data = last decade (2016-2025)
medresult <- mediation::mediate(med, out, treat = "timew", mediator = "helminthw")
summary(medresult)
