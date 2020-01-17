rm(list = ls())

resp.dat <- read.csv("data/2.metabolic_rates.csv") %>%
  rename(SELECTION = X.SELECTION)

str(resp.dat) 

library(nlme)
library(tidyverse)
library(ggthemes)

# check balance of experimental design
xtabs(~ SEX + LINE, resp.dat)




# prelim models
# CO2
act.CO2 <- lme(VCO2 ~ ACTIVITY, random = ~1|LINE/SEX/SAMPLE, data = resp.dat, method = 'REML')
anova(act.CO2)
summary(act.CO2) # CO2 increases with activity

bod.CO2 <- lme(VCO2 ~ BODY_WEIGHT, random = ~1|LINE/SEX/SAMPLE, data = resp.dat, method = 'REML')
anova(bod.CO2)
summary(bod.CO2) # CO2 increases with body_Weight

# O2
act.O2 <- lme(VO2 ~ ACTIVITY, random = ~1|LINE/SEX/SAMPLE, data = resp.dat, method = 'REML')
anova(act.O2)
summary(act.O2)

bod.O2 <- lme(VO2 ~ BODY_WEIGHT, random = ~1|LINE/SEX/SAMPLE, data = resp.dat, method = 'REML')
anova(bod.O2)
summary(bod.O2) 

# RQ
act.RQ <- lme(RQ ~ ACTIVITY, random = ~1|LINE/SEX/SAMPLE, data = resp.dat, method = 'REML')
anova(act.RQ)
summary(act.RQ) 

bod.RQ <- lme(RQ ~ BODY_WEIGHT, random = ~1|LINE/SEX/SAMPLE, data = resp.dat, method = 'REML')
anova(bod.RQ)
summary(bod.RQ)



# MEAN CENTRE COVARIATES --------------------------------------------------
act_lme <- lme(ACTIVITY ~ SELECTION * SEX * BODY_WEIGHT,
               random = ~1|LINE/SEX/SAMPLE,
               data = resp.dat)

anova(act_lme) # selection effect on activity
summary(act_lme)



bod_lme <- lme(BODY_WEIGHT ~ SELECTION * SEX,
               random = ~1|LINE/SEX,
               data = resp.dat)

anova(bod_lme) # sex effect on body weight


resp.dat <- resp.dat %>%
  mutate(SELECTION = recode_factor(SELECTION, Mono = 'Monogamy', Poly = "Polyandry")) %>% 
  group_by(SEX) %>% # mean centre BODY_WEIGHT by SEX
  mutate(BODY_WEIGHT = scale(BODY_WEIGHT, scale = FALSE)) %>%
  group_by(SELECTION) %>% # mean centre ACTIVITY by SELECTION
  mutate(ACTIVITY = scale(ACTIVITY, scale = FALSE))



# TAKE MEAN ACROSS CYCLE -------------------------------------------------
resp.2 <- resp.dat %>% 
  group_by(SAMPLE, LINE, SELECTION, SEX) %>% 
  summarise(meanVCO2 = mean(VCO2),
            meanVO2 = mean(VO2),
            ACTIVITY = mean(ACTIVITY),
            BODY_WEIGHT = mean(BODY_WEIGHT),
            meanRQ = mean(RQ))

# correlation between VCO2 and VO2
plot(resp.2$meanVCO2, resp.2$meanVO2)
cor.test(resp.2$meanVCO2, resp.2$meanVO2)



# MEAN LME MODELS ---------------------------------------------------------
# metabolic rates models
O2_lme <- lme(meanVO2 ~ 
                SELECTION + 
                SEX + 
                ACTIVITY + 
                BODY_WEIGHT + 
                SELECTION:SEX + 
                SELECTION:ACTIVITY + 
                SELECTION:BODY_WEIGHT +
                SEX:ACTIVITY + 
                SEX:BODY_WEIGHT +
                SELECTION:SEX:ACTIVITY + 
                SELECTION:SEX:BODY_WEIGHT,
              random = ~1|LINE/SEX,
              data = resp.2, method = "REML")

anova(O2_lme)


CO2_lme <- lme(meanVCO2 ~ 
                 SELECTION + 
                 SEX + 
                 ACTIVITY + 
                 BODY_WEIGHT + 
                 SELECTION:SEX + 
                 SELECTION:ACTIVITY + 
                 SELECTION:BODY_WEIGHT +
                 SEX:ACTIVITY + 
                 SEX:BODY_WEIGHT +
                 SELECTION:SEX:ACTIVITY + 
                 SELECTION:SEX:BODY_WEIGHT,
               random = ~1|LINE/SEX,
               data = resp.2, method = "REML")

anova(CO2_lme)
summary(CO2_lme)
VarCorr(CO2_lme)

plot(fitted(CO2_lme), residuals(CO2_lme))
qqnorm(CO2_lme, abline = c(0, 1))



RQ_lme <- lme(meanRQ ~ 
                SELECTION + 
                SEX + 
                ACTIVITY + 
                BODY_WEIGHT + 
                SELECTION:SEX + 
                SELECTION:ACTIVITY + 
                SELECTION:BODY_WEIGHT +
                SEX:ACTIVITY + 
                SEX:BODY_WEIGHT +
                SELECTION:SEX:ACTIVITY + 
                SELECTION:SEX:BODY_WEIGHT,
              random = ~1|LINE/SEX,
              data = resp.2, method = "REML")

anova(RQ_lme)
summary(RQ_lme)
VarCorr(RQ_lme)


lsmeans::lsmeans(RQ_lme, pairwise ~ SELECTION | SEX)
lsmeans::lsmeans(RQ_lme, pairwise ~ SEX | SELECTION)





# ANOVA table output
# capture.output(anova(CO2_lme), file = "couch_potato_tables/ANOVA_CO2_lme.csv")
# capture.output(anova(RQ_lme), file = "couch_potato_tables/ANOVA_RQ_lme.csv")



# FIGURE 2 ----------------------------------------------------------------
prediction_data <- expand.grid(ACTIVITY = seq(from = min(resp.2$ACTIVITY), 
                                              to = max(resp.2$ACTIVITY), 
                                              length = 10),
                               BODY_WEIGHT = seq(from = min(resp.2$BODY_WEIGHT),
                                                 to = max(resp.2$BODY_WEIGHT),
                                                 length = 10), 
                               LINE = levels(resp.2$LINE),
                               SEX = levels(resp.2$SEX),
                               SELECTION = levels(resp.2$SELECTION))

newY_fixed <- predict(CO2_lme, newdata = prediction_data)

plot_data <- data.frame(prediction_data, newY_fixed)
head(plot_data)
colnames(plot_data)[6] <- 'meanVCO2'

head(plot_data)

facet_names <- c(F = "Female", M = "Male")

resp.2 %>% 
  ggplot(aes(x = BODY_WEIGHT, y = meanVCO2, colour = SELECTION)) +
  stat_smooth(method = "lm") + 
 # geom_line(aes(group = SELECTION, lty = SELECTION), lwd = 1) +
  geom_point(data = resp.2, aes(y = meanVCO2, fill = SELECTION), pch = 21, size = 5) +
  labs(x = "Body weight", y = expression(paste("Metabolic rate (µL CO"[2], ")"))) +
  facet_wrap(~SEX, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = c(0.09, .9),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(2, 'lines'),
        strip.background = element_rect(fill = "grey"),
        strip.text = element_text(face = 'bold', size = 20),
        plot.background = element_rect(colour = NA))

#ggsave(filename = "couch_potato_figs/meanCO2_predicted.pdf", width = 14, height = 7)
#ggsave(filename = "couch_potato_figs/meanCO2_predicted.tiff", width = 14, height = 7)


resp.2 %>% 
  group_by(SELECTION, SEX) %>% 
  summarise(mRQ = mean(meanRQ),
            se = sd(meanRQ)/sqrt(n()),
            CI.95 = qt(.95, df = length(meanRQ) - 1)*se,
            u95 = mRQ + CI.95,
            l95 = mRQ - CI.95) %>% 
  ggplot(aes(x = SEX, y = mRQ)) +
  geom_errorbar(aes(ymin = mRQ - se, ymax = mRQ + se), width = .2) +
  geom_line(aes(group = SELECTION)) +
  geom_point(size = 10, pch = 21, aes(fill = SELECTION)) +
  #scale_fill_manual(values = c("black", "white")) +
  scale_x_discrete(labels = c("Female", "Male")) +
  scale_y_continuous(limits = c(0.84, 1)) +
  labs(y = "Respiratory Quotient (mean ± SE)") +
  theme_bw() +
  theme(legend.position = c(.2, .85),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = NA))

#ggsave(filename = "couch_potato_figs/meanRQ_SE.pdf", width = 5.5, height = 5)
