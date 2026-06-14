library("tidyverse")
library("plm")
library("jtools")
library("stargazer")
library("mediation")
library("ggeffects")

# ldl = lipids data
# glu = glucose data 
# bp = blood pressure data
# bmi = BMI data
# weight = weight data
# hb = hemoglobin data

# pid = person id 



# function for estimating LDL with the NIH2 equation 
calculate_LDL_NIH2 <- function(TC, HDL, TG) {
  # Calculate Non-HDL cholesterol
  NonHDL <- TC - HDL
  
  # NIH 2 equation formula
  # LDL-C = TC/0.948 - HDL-C/0.971 - (TG/8.56 + [TG * NonHDL-C]/2140 - TG^2/16100) - 9.44
  LDL_C <- (TC / 0.948) - (HDL / 0.971) - (TG / 8.56 + (TG * NonHDL) / 2140 - (TG^2) / 16100) - 9.44
  
  return(LDL_C)
}

ldl$nih2ldl <- calculate_LDL_NIH2(ldl$col, ldl$hdl, ldl$trig)

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
m8 <- plm(glucose ~ years + helminth,
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

## Panel models of infection predicting cardiometabolic health,
##  including eosinophils 

#####################################################################

# eosinophils were standardized: as.numeric(scale(eosinophils))

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
          out = "/filepath/hookhelm_liphbsbp_within.html")



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
m3 <- plm(trig ~ years + hookworm*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(trig ~ years + helminth*car4,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(glucose ~ years + hookworm*car4,
          index = c("pid"),
          model = "within",
          data = glu)
m6 <- plm(glucose ~ years + helminth*car4,
          index = c("pid"),
          model = "within",
          data = glu)
m7 <- plm(SBP ~ years + hookworm*car4,
          index = c("pid"),
          model = "within",
          data = bp)
m8 <- plm(SBP ~ years + helminth*car4,
          index = c("pid"),
          model = "within",
          data = bp)


stargazer(m1, m2, m3, m4, m5, m6, m7, m8, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_car4_within.html")



#####################################################################

## Figure 2

#####################################################################

ldl4 <- filter(ldl, car4 == 1)
ldl3 <- filter(ldl, car4 == 0)

hookldl3 <- plm(ldl ~ years + hookworm,
                index = c("pid"),
                model = "within",
                data = ldl3)
hookldl4 <- plm(ldl ~ years + hookworm,
                index = c("pid"),
                model = "within",
                data = ldl4)

pred1 <- ggpredict(
  hookldl3,
  terms = c("hookworm")
)
pred1$group <- as.factor("ε3/ε3")

pred2 <- ggpredict(
  hookldl4,
  terms = c("hookworm")
)
pred2$group <- as.factor("ε4 carrier")

pred <- rbind(pred1, pred2)


df <- tribble(
  ~infection, ~Genotype, ~LDL,
  "No infection", "ε3/ε3", 92.82277,  
  "Hookworm", "ε3/ε3", 85.25156,    
  "No infection", "ε4 carrier", 92.59464,  
  "Hookworm", "ε4 carrier", 93.51044  
)

df$infection <- factor(
  df$infection,
  levels = c("No infection", "Hookworm")
)

apoeplothook <- ggplot(df,
                       aes(x = infection,
                           y = LDL,
                           group = Genotype,
                           shape = Genotype,
                           color = Genotype)) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  
  labs(
    y = "Predicted LDL (mg/dL)",
    x = NULL,
    title = expression(
      paste(
        "APOE-", epsilon, "4 carriers maintain LDL during infection"
      )
    )
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )
apoeplothook

ggsave("filepath/apoeplothook.png", plot = apoeplothook, dpi = "retina",
       height = 4, width = 7.2)



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

## Figure 3

#####################################################################

ldl$ldlz <- as.numeric(scale(ldl$ldl))
ldl$trigz <- as.numeric(scale(ldl$trig))
glu$glucosez <- as.numeric(scale(glu$glucose))
glu$hookwormz <- as.numeric(scale(glu$hookworm))

# because R calculates time in reference to 1970, adding 1970 gives the actual years (e.g., 2010, 2011, etc.)
ldl$yearsmod <- ldl$years + 1970
glu$yearsmod <- glu$years + 1970

ldlz <- ldl[,c("yearsmod", "ldlz", "pid")]
trigz <- ldl[,c("yearsmod", "trigz", "pid")]
glucosez <- glu[,c("yearsmod", "glucosez", "pid")]
hookwormz <- glu[,c("yearsmod", "hookwormz", "pid")]


ldlmod <- plm(ldlz ~ yearsmod,
              index = c("pid"),
              model = "within",
              data = ldlz)

ldlz$values <- predict(ldlmod, ldlz)


trigmod <- plm(trigz ~ yearsmod,
               index = c("pid"),
               model = "within",
               data = trigz)

trigz$values <- predict(trigmod, trigz)


glumod <- plm(glucosez ~ yearsmod,
              index = c("pid"),
              model = "within",
              data = glucosez)

glucosez$values <- predict(glumod, glucosez)


hookmod <- plm(hookwormz ~ yearsmod,
               index = c("pid"),
               model = "within",
               data = hookwormz)

hookwormz$values <- predict(hookmod, hookwormz)


ldlz <- ldlz[,c(1,4)]
trigz <- trigz[,c(1,4)]
glucosez <- glucosez[,c(1,4)]
hookwormz <- hookwormz[,c(1,4)]


ldlz$data <- "LDL"
trigz$data <- "Triglycerides"
glucosez$data <- "Glucose"
hookwormz$data <- "Hookworms"


plotdat <- rbind(ldlz, trigz, glucosez, hookwormz)

plotdat$data <- factor(
  plotdat$data,
  levels = c("LDL", "Triglycerides", "Glucose", "Hookworms")
)

epitrans <- ggplot(plotdat[plotdat$yearsmod >= 2016,], aes(x = yearsmod, y = values, color = data)) + geom_smooth(method = "lm") + 
  scale_color_manual(values = c("brown1", "gold3", "dodgerblue", "green4")) + 
  labs(
    y = "Standardized values",
    x = "Time",
    title = "An ongoing epidemiological transition"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top",
        legend.title = element_blank()
  )


ggsave("filepath/epitrans.png", plot = epitrans, dpi = "retina",
       height = 5, width = 7.2)



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



#####################################################################

## Panel models of infection predicting cardiometabolic health,
##  interaction by town distance (proxy for market integration) 

#####################################################################

m1 <- plm(ldl ~ years + hookworm*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = ldl)
m2 <- plm(ldl ~ years + helminth*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = ldl)
m3 <- plm(hdl ~ years + hookworm*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(hdl ~ years + helminth*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(trig ~ years + hookworm*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = ldl)
m6 <- plm(trig ~ years + helminth*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = ldl)

stargazer(m1, m2, m3, m4, m5, m6, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_distanceinteraction1_within.html")



m1 <- plm(glucose ~ years + hookworm*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = glu)
m2 <- plm(glucose ~ years + helminth*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = glu)
m3 <- plm(SBP ~ years + hookworm*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = bp)
m4 <- plm(SBP ~ years + helminth*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = bp)
m5 <- plm(hb ~ years + hookworm*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = hb)
m6 <- plm(hb ~ years + helminth*route_distance_townp,
          index = c("pid"),
          model = "within",
          data = hb)

stargazer(m1, m2, m3, m4, m5, m6, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "/filepath/hookhelm_distanceinteraction2_within.html")



#####################################################################

## Panel models of infection predicting cardiometabolic health,
##  interaction by person-average cardiometabolic health 

#####################################################################

mldl <- filter(ldl, male == 1)
mldl <- group_by(mldl, age) %>% mutate(ldlz = scale(ldlcalc)) %>% ungroup() %>% group_by(pid) %>% mutate(ldlzave = mean(ldlz, na.rm = T)) %>% 
  ungroup() 
mldl <- group_by(mldl, age) %>% mutate(hdlz = scale(hdl)) %>% ungroup() %>% group_by(pid) %>% mutate(hdlzave = mean(hdlz, na.rm = T)) %>% 
  ungroup()
mldl <- group_by(mldl, age) %>% mutate(trigz = scale(trig)) %>% ungroup() %>% group_by(pid) %>% mutate(trigzave = mean(trigz, na.rm = T)) %>% 
  ungroup() 

fldl <- filter(ldl, male == 0)
fldl <- group_by(fldl, age) %>% mutate(ldlz = scale(ldlcalc)) %>% ungroup() %>% group_by(pid) %>% mutate(ldlzave = mean(ldlz, na.rm = T)) %>% 
  ungroup() 
fldl <- group_by(fldl, age) %>% mutate(hdlz = scale(hdl)) %>% ungroup() %>% group_by(pid) %>% mutate(hdlzave = mean(hdlz, na.rm = T)) %>% 
  ungroup() 
fldl <- group_by(fldl, age) %>% mutate(trigz = scale(trig)) %>% ungroup() %>% group_by(pid) %>% mutate(trigzave = mean(trigz, na.rm = T)) %>% 
  ungroup()

ldl <- rbind(mldl,fldl)



m1 <- plm(ldl ~ years + hookworm*ldlzave,
          index = c("pid"),
          model = "within",
          data = ldl)
m2 <- plm(ldl ~ years + helminth*ldlzave,
          index = c("pid"),
          model = "within",
          data = ldl)
m3 <- plm(hdl ~ years + hookworm*hdlzave,
          index = c("pid"),
          model = "within",
          data = ldl)
m4 <- plm(hdl ~ years + helminth*hdlzave,
          index = c("pid"),
          model = "within",
          data = ldl)
m5 <- plm(trig ~ years + hookworm*trigzave,
          index = c("pid"),
          model = "within",
          data = ldl)
m6 <- plm(trig ~ years + helminth*trigzave,
          index = c("pid"),
          model = "within",
          data = ldl)

stargazer(m1, m2, m3, m4, m5, m6, 
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "filepath/hookhelm_highlipids_within.html")



mglu <- filter(glu, male == 1)
mglu <- group_by(mglu, age) %>% mutate(glucosez = scale(glucose)) %>% ungroup() %>% group_by(pid) %>% 
  mutate(glucosezave = mean(glucosez, na.rm = T)) %>% ungroup()

fglu <- filter(glu, male == 0)
fglu <- group_by(fglu, age) %>% mutate(glucosez = scale(glucose)) %>% ungroup() %>% group_by(pid) %>% 
  mutate(glucosezave = mean(glucosez, na.rm = T)) %>% ungroup()

glu <- rbind(mglu,fglu)


mbp <- filter(bp, male == 1)
mbp <- group_by(mbp, age) %>% mutate(sbpz = scale(SBP)) %>% ungroup() %>% group_by(pid) %>% mutate(sbpzave = mean(sbpz, na.rm = T)) %>% 
  ungroup() 
fbp <- filter(bp, male == 0)
fbp <- group_by(fbp, age) %>% mutate(sbpz = scale(SBP)) %>% ungroup() %>% group_by(pid) %>% mutate(sbpzave = mean(sbpz, na.rm = T)) %>% 
  ungroup() 

bp <- rbind(mbp,fbp)


mhb <- filter(hb, male == 1)
mhb <- group_by(mhb, age) %>% mutate(hbz = scale(c_hb)) %>% ungroup() %>% group_by(pid) %>% mutate(hbzave = mean(hbz, na.rm = T)) %>% 
  ungroup() 
fhb <- filter(hb, male == 0)
fhb <- group_by(fhb, age) %>% mutate(hbz = scale(c_hb)) %>% ungroup() %>% group_by(pid) %>% mutate(hbzave = mean(hbz, na.rm = T)) %>% 
  ungroup() 

hb <- rbind(mhb,fhb)



m1 <- plm(glucose ~ years + hookworm*glucosezave,
          index = c("pid"),
          model = "within",
          data = glu)
m2 <- plm(glucose ~ years + helminth*glucosezave,
          index = c("pid"),
          model = "within",
          data = glu)
m3 <- plm(SBP ~ years + hookworm*sbpzave,
          index = c("pid"),
          model = "within",
          data = bp)
m4 <- plm(SBP ~ years + helminth*sbpzave,
          index = c("pid"),
          model = "within",
          data = bp)
m5 <- plm(hb ~ years + hookworm*hbzave,
          index = c("pid"),
          model = "within",
          data = hb)
m6 <- plm(hb ~ years + helminth*hbzave,
          index = c("pid"),
          model = "within",
          data = hb)

stargazer(m1, m2, m3, m4, m5, m6,
          align = T, single.row = T, digits = 2, ci = T, star.cutoffs = c(0.10, 0.05, 0.01), 
          star.char = c("t", "*", "**"), 
          out = "filepath/hookhelm_highglusbphb_within.html")












