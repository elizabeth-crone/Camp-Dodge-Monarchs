rm(list = ls())

# analysis of 2020 field experiment

setwd("C:/Users/ecrone01/Box/Documents/Monarch IA/final R code and data")

mon = read.csv("mon2020 cage summary.csv")

head(mon)
mon$fail = mon$Max.of.larvae - mon$Sum.of.adults

library(lme4)
library(emmeans)
library(car)
library(plotrix)

# effects model for significance testing
m0 = glmer(cbind(Sum.of.adults, fail) ~ as.factor(batch) + (1|cage), family = binomial, data = mon)
summary(m0)
Anova(m0)
(estimates = emmeans(m0, ~as.factor(batch)))
pairs(estimates)

# means model for graphing
m0a = glmer(cbind(Sum.of.adults, fail) ~ -1 + as.factor(batch) + (1|cage), family = binomial, data = mon)
plogis(confint(m0a))
(errs = plogis(confint(m0a, level = 0.834)[2:4,]))
(yvals = plogis(fixef(m0a) ))
#358 x 496
plotCI(1:3, 100*yvals, ui = 100*errs[,2], li = 100*errs[,1], xlim = c(0.5, 3.5), xaxt = "n", ylab = "", xlab = "", pch = 21, cex = 1.35, pt.bg = "indianred")  
mtext(c("earlier", "current", "later"), side = 1, line = 0, at = 1:3)
mtext("timing of release", side = 1, line = 1)
mtext("(compared to 1st natural cohort)", side = 1, line = 2)
mtext(side = 2, line = 2, "% survival to ecolosion")

mon2 = read.csv("summary with pupae and dates.csv")
head(mon2)
table(mon2$batch, mon2$date_pupa)

mon2$dev.time1 = mon2$date_pupa - mon2$date0
mon2$dev.time2 = mon2$date_adults - mon2$date0

# this analysis (time from release to pupa) not included in the manuscript
m1 = lmer(dev.time1 ~ as.factor(batch) + (1|cage), data = mon2)
Anova(m1)
hist(resid(m1))
summary(m1)

# development time - time from release to adult
# effects model for significance testing
m2 = lmer(dev.time2 ~ as.factor(batch) + (1|cage), data = mon2)
Anova(m2)
hist(resid(m2))
summary(m2)
(estimates2 = emmeans(m2, ~as.factor(batch)))
pairs(estimates2)


# means model for graphing 
m2a = lmer(dev.time2 ~ 0 + as.factor(batch) + (1|cage), data = mon2)
fixef(m2a)
confint(m2a)
(errs2 = (confint(m2a, level = 0.834)[3:5,]))
(yvals2 = (fixef(m2a) ))      
plotCI(1:3, yvals2, ui = errs2[,2], li = errs2[,1], xlim = c(0.5, 3.5), xaxt = "n", ylab = "", xlab = "", pch = 21, cex = 1.35, pt.bg = "indianred")  
mtext(c("earlier", "current", "later"), side = 1, line = 0, at = 1:3)
mtext("timing of release", side = 1, line = 1)
mtext("(compared to 1st natural cohort)", side = 1, line = 2)
mtext(side = 2, line = 2, "development time (days)")

