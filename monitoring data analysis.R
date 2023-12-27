rm(list = ls())

library(car)
library(quantreg)
library(bbmle)
library(mgcv)
library(lmtest)

setwd("C:/Users/ecrone01/Box/Documents/Monarch IA/final R code and data")
mon = read.csv("ratcliff_03_19.csv")

head(mon)
summary(mon)
mon$effort = 3.2
mon$effort[mon$year > 2015] = 2.4


# looking at some basic statistics
table(mon$count)
table(mon$year)
pal <- colorRampPalette(c("blue", "green"))
nyrs = length(unique(mon$year))
pal2 = pal(nyrs)
mon$syear = mon$year - min(mon$year)
plot(jitter(mon$DOY), jitter(mon$count), col = pal2[mon$syear+1], ylab = "monarchs per survey", xlab = "Day of year", xaxt = "n")
axis(side = 1, at = c(90, 120, 151, 181, 212, 243, 273, 304), labels = c("Apr-1", "May-1", "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1", "Nov-1"), cex = 0.8)
legend("topleft", lwd = 1, col = pal2, legend = 2003:2019, cex = 0.65)

# growing degree days, for comparison to day of year (supplemental analyses)
gdd = read.csv("GDD4R.csv")
head(gdd)

mon2 = merge(mon, gdd)
head(mon2)

##################################################
# analyses to illustrate the way GAMs work
# not included in the results section of the paper
##################################################
# warm-up example - using GAMs to characterize the phenology across all years
# because monarchs are multivoltine, we are fitting the expected number vs. day of year using GAMs to fit flexible functions
tapply(mon2$count, mon2$year, mean)
tapply(mon2$count, mon2$year, max)
m1 = gam(count ~ s(DOY), family = poisson, data = mon2[mon2$year == 2007,], niterPQL = 200)
# plot of the number of butterflies per transect walk, as a function of day of year
# the "rug" at the bottom of the plot shows whether there was at least one count on each day
# we may not want to keep this in final graphs, but it is a GAM
plot(m1, scale = 0, trans = exp, shift = coef(m1)[1], ylab = "butterflies / transect", xlab = "day of year", main = "2007", ylim = c(0,max(mon$count[mon2$year == 2007])), rug = F, lwd = 2)
# you can see that there is a pulse at the beginning of the season, then a couple generations in which the monarch population grows, then a pulse at the end (possibly due to migrants coming through, but a 4-fold increase would not be surprising just due to onsite breeding)



# tests for phenology vs. DOY and year
# these are the main results in the paper
m3a1 = gam(count ~ s(DOY), family = nb(), offset = log(effort), data = mon2) # no effect of year
m3a2 = gam(count ~ s(year), family = nb(), offset = log(effort), data = mon2) # no effect of year
m3b = gam(count ~ s(DOY) + s (year), family = nb(), offset = log(effort), data = mon2) # abundance, not phenology, changes among years
m3c = gam(count ~ s(DOY,year), family = nb(), offset = log(effort), data = mon2) # phenology (and abundance) differ among years

# marginal hypothesis tests
# test main effect of DOY
lrtest(m3a2, m3b)
#test main effect of year
lrtest(m3a1, m3b)
# test of interaction effect
lrtest(m3b, m3c)

# plot of the full model with the interaction (since all effects are statistically supported)
years = 2003:2019
DOYs = 100:300
DOYpreds = predict(m3c, se.fit = T, newdata = data.frame(year = years[1], DOY = DOYs))
plot(DOYs, 3.2*exp(DOYpreds$fit), type = "l", ylim = c(0,2.5), ylab = "", xlab = "", xaxt = "n", col = pal2[1])
mtext(side = 1, line = 2, "Day of year")
mtext(side = 2, line = 2, "monarchs per km")
for(i in 2:nyrs){
  DOYpreds = predict(m3c, se.fit = T, newdata = data.frame(year = years[i], DOY = DOYs))
  points(DOYs, 3.2*exp(DOYpreds$fit), type = "l", col = pal2[i])
}
legend("topleft", lwd = 1, col = pal2, legend = c(2003, "", 2005, "", 2007, "", 2009, "", 2011, "", 2013, "", 2015, "", 2017, "", 2019), cex = 0.5)
axis(side = 1, at = c(120, 151, 181, 212, 243, 273,304), labels = c("May-1", "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1", "Nov-1"))
mtext(side = 3, line = 0.5, "Monarch butterflies, Camp Dodge IA")
##mtext(side = 3, line = 0.5, "DATA: H. Ratcliff", cex = 0.75)


# comparing DOY to GDD as predictors of phenology (for supplemental analysis)
m3d = gam(count ~ s(cumGDD), family = nb(), offset = log(effort), data = mon2) # no effect of year
m3e = gam(count ~ s(cumGDD) + s (year), family = nb(), offset = log(effort), data = mon2) # abundance, not phenology, changes among years
m3f = gam(count ~ s(cumGDD,year), family = nb(), offset = log(effort), data = mon2) # phenology (and abundance) differ among years

ICtab(m3a1, m3a2, m3b, m3c, m3d, m3e, m3f) # interaction DOY model still wins



# analysis of derived metrics from the fitted GAM
years = 2003:2019
DOYs = 100:300
output = array(NA, dim = c(nyrs, 6))
ndays = length(DOYs)
for(i in 1:nyrs){
  DOYpreds = predict(m3c, se.fit = T, newdata = data.frame(year = years[i], DOY = DOYs))
  count.mean = sum(exp(DOYpreds$fit))
  count.slope = DOYpreds$fit[1:(ndays-1)] - DOYpreds$fit[2:ndays]
  change.slope = sign(count.slope)[2:(ndays-1)] - sign(count.slope)[1:(ndays-2)]
  peaks = DOYs[which(change.slope == 2)]
  peak.count = exp(DOYpreds$fit[which(change.slope == 2)])
  output[i,] = c(years[i], count.mean, min(peaks), max(peaks), peak.count[1], peak.count[length(peak.count)])
  
}

output = data.frame(output)
names(output) = c("year", "abundance", "DOY_onset", "DOY_end", "peak1", "peak3")
output
output$summer.grow = (output$peak3/output$peak1)
output


ndat = dim(mon2)[1]
nreps = 5000
outboot = array(NA, dim = c(nreps, nyrs, 6))
for(j in 1:nreps){
  monboot = mon2[sample(ndat, replace = T),]
  m.new = gam(count ~ s(DOY,year), family = nb(), offset = log(effort), data = monboot) # phenology (and abundance) differ among years
  for(i in 1:nyrs){
    DOYpreds = predict(m.new, se.fit = T, newdata = data.frame(year = years[i], DOY = DOYs))
    count.mean = sum(exp(DOYpreds$fit))
    count.slope = DOYpreds$fit[1:(ndays-1)] - DOYpreds$fit[2:ndays]
    change.slope = sign(count.slope)[2:(ndays-1)] - sign(count.slope)[1:(ndays-2)]
    peaks = DOYs[which(change.slope == 2)]
    peak.count = exp(DOYpreds$fit[which(change.slope == 2)])
    outboot[j,i,] = c(years[i], count.mean, min(peaks), max(peaks), peak.count[1], peak.count[length(peak.count)])
    
  }}

outboot[1,,]
for(i in 1:nyrs){
  print(quantile(outboot[,i,2], probs = c(0.025, 0.975)))
  
}

coef.N = array(NA, dim = c(nreps, 3))
year.N = array(NA, dim = c(nreps, nyrs))
for(i in 1:nreps){
  usedat = outboot[i,,]
  N.lm = lm(usedat[,2] ~ usedat[,1]) # abundance through time
  coef.N[i,1:2] = coef(N.lm)
  coef.N[i,3] = coefficients(summary(N.lm))[2,4]
  year.N[i,] = predict(N.lm)
}
quantile(coef.N[,3], probs = c(0.025, 0.5, 0.975))
sum(coef.N[,3] < 0.05)/nreps

count.out = array(NA, dim = c(nyrs, 3))
for(i in 1:nyrs){
  count.out[i,1:3] = quantile(outboot[,i,2], probs = c(0.16, 0.5, 0.84))
}
count.out



# outboot matrices summarize  
#1. years - 2003 through 2019
#2. count.mean - sum of the back-transformed area under the GAM
#3. min(peaks) - date of first GAM peak (metric of arrival phenology)
#4. max(peaks) - date of the last GAM peak (metric of departure phenology)
#5. peak.count[1] - height of first peak (abundance at start of summer)
#6. peak.count[length(peak.count) - height of last peak (abundance at end of summer)

usecol = 6 # choose the metric you want to summarize and plot

# code for making plots
coef.N = array(NA, dim = c(nreps, 3))
year.N = array(NA, dim = c(nreps, nyrs))
for(i in 1:nreps){
  usedat = outboot[i,,]
  N.lm = lm(usedat[,usecol] ~ usedat[,1]) # abundance through time
  coef.N[i,1:2] = coef(N.lm)
  coef.N[i,3] = coefficients(summary(N.lm))[2,4]
  year.N[i,] = predict(N.lm)
}
quantile(coef.N[,3], probs = c(0.025, 0.5, 0.975))
sum(coef.N[,3] < 0.05)/nreps

count.out = array(NA, dim = c(nyrs, 3))
for(i in 1:nyrs){
  count.out[i,1:3] = quantile(outboot[,i,usecol], probs = c(0.16, 0.5, 0.84), na.rm = T)
}
sig.pos = sum(coef.N[,2]>0 & coef.N[,3] < 0.05)/sum(!is.na(coef.N[,3]))
sig.neg = sum(coef.N[,2]<0 & coef.N[,3] < 0.05)/sum(!is.na(coef.N[,3]))

head(output)
#293 x 369
plot(output$year, output[, usecol], ylim = c(min(count.out), max(count.out)), type = "l", xlab = "", ylab = "")
mtext(side = 1, line = 2, "year")
mtext(side = 2, line = 2, "height of last peak") # need to change the name for plotting different metrics
polygon(c(output$year, rev(output$year)), c(count.out[,1], rev(count.out[,3])), col = "gray90", border = F)
points(output$year, output[, usecol], type = "l", lwd = 3)
text(x = 2015, y = (min(count.out) + 0.95*(max(count.out) - min(count.out))), paste(round(100*sig.pos,1), "% sig +"), col = "red", cex = 0.8)
text(x = 2015, y = (min(count.out) + 0.875*(max(count.out) - min(count.out))), paste(round(100*sig.neg,1), "% sig -"), col = "red", cex = 0.8)


#growth rate plot needs special code since it is a function of two metrics
coef.N = array(NA, dim = c(nreps, 3))
year.N = array(NA, dim = c(nreps, nyrs))
for(i in 1:nreps){
  usedat = outboot[i,,]
  N.lm = lm(usedat[,6]/usedat[,5] ~ usedat[,1]) # abundance through time
  coef.N[i,1:2] = coef(N.lm)
  coef.N[i,3] = coefficients(summary(N.lm))[2,4]
  year.N[i,] = predict(N.lm)
}
quantile(coef.N[,3], probs = c(0.025, 0.5, 0.975))
sum(coef.N[,3] < 0.05)/nreps

count.out = array(NA, dim = c(nyrs, 3))
for(i in 1:nyrs){
  count.out[i,1:3] = quantile(outboot[,i,6]/outboot[,i,5], probs = c(0.16, 0.5, 0.84), na.rm = T)
}
sig.pos = sum(coef.N[,2]>0 & coef.N[,3] < 0.05)/sum(!is.na(coef.N[,3]))
sig.neg = sum(coef.N[,2]<0 & coef.N[,3] < 0.05)/sum(!is.na(coef.N[,3]))

head(output)
#293 x 369
plot(output$year, output$summer.grow, ylim = c(min(count.out), max(count.out)), type = "l", xlab = "", ylab = "")
mtext(side = 1, line = 2, "year")
mtext(side = 2, line = 2, "summer growth (last/1st)")
polygon(c(output$year, rev(output$year)), c(count.out[,1], rev(count.out[,3])), col = "gray90", border = F)
points(output$year, output$summer.grow, type = "l", lwd = 3)
text(x = 2015, y = (min(count.out) + 0.95*(max(count.out) - min(count.out))), paste(round(100*sig.pos,1), "% sig +"), col = "red", cex = 0.8)
text(x = 2015, y = (min(count.out) + 0.875*(max(count.out) - min(count.out))), paste(round(100*sig.neg,1), "% sig -"), col = "red", cex = 0.8)


# supplemental analyses further exploring GDD effects

# changes in GDD through time
gmod = lm(cumGDD ~ scale(year)*I(DOY-100), data = mon2)
summary(gmod)
Anova(gmod)

gmod2 = lm(cumGDD ~ as.factor(year) + as.factor(year):I(DOY-100), data = mon2)
betas = coef(gmod2)[18:34]
lims = confint(gmod2)[18:34,]
library(plotrix)
#361x371
plotCI(2003:2019, betas, li = lims[,1], ui = lims[,2], pch = 19, col = pal2, cex = 1.2, xlab = "", ylab = "")
mtext(side = 2, line = 2, "degrees per day")
mtext(side = 1, line = 2, "year")


dat = mon2[mon2$year == min(mon2$year),]
with(dat, plot(DOY, cumGDD, type = "l", col = pal2[1], xlab = "", ylab = "", xaxt = "n"))
mtext(side = 1, line = 2, "day of year")
mtext(side = 2, line = 2, "cumulative degree days")
for(i in 2:nyrs){
  dat = mon2[mon2$year == years[i],]
  with(dat, points(DOY, cumGDD, type = "l", col = pal2[i]))
}
legend("bottomright", lwd = 1, col = pal2, legend = c(2003, "", 2005, "", 2007, "", 2009, "", 2011, "", 2013, "", 2015, "", 2017, "", 2019), cex = 0.45)
axis(side = 1, at = c(120, 151, 181, 212, 243, 273, 304), labels = c("May-1", "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1", "Nov-1"))


# yearly curves, colored by GDD
pal3 <- colorRampPalette(c("yellow", "red"))
pal4 = pal3(nyrs)

#291x371
yrGDDs = tapply(mon2$cumGDD, mon2$year, max)
plot(DOYs, 3.2*exp(DOYpreds$fit), type = "l", ylim = c(0,2.5), ylab = "", xlab = "", xaxt = "n", col = pal4[order(yrGDDs)][1])
mtext(side = 1, line = 2, "day of year")
mtext(side = 2, line = 2, "monarchs per survey mile")
for(i in 2:nyrs){
  DOYpreds = predict(m3c, se.fit = T, newdata = data.frame(year = years[i], DOY = DOYs))
  points(DOYs, 3.2*exp(DOYpreds$fit), type = "l", col = pal4[order(yrGDDs)][i])
}
legend("topleft", lwd = 1, col = pal4[order(yrGDDs)], legend = c(2003, "", 2005, "", 2007, "", 2009, "", 2011, "", 2013, "", 2015, "", 2017, "", 2019), cex = 0.45)
axis(side = 1, at = c(120, 151, 181, 212, 243, 273,304), labels = c("May-1", "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1", "Nov-1"))
mtext(side = 3, line = 0.5, "Ordered by total GDD")

yrGDDs = tapply(mon2$cumGDD, mon2$year, mean)
plot(DOYs, 3.2*exp(DOYpreds$fit), type = "l", ylim = c(0,2.5), ylab = "", xlab = "", xaxt = "n", col = pal4[order(yrGDDs)][1])
mtext(side = 1, line = 2, "day of year")
mtext(side = 2, line = 2, "monarchs per survey mile")
for(i in 2:nyrs){
  DOYpreds = predict(m3c, se.fit = T, newdata = data.frame(year = years[i], DOY = DOYs))
  points(DOYs, 3.2*exp(DOYpreds$fit), type = "l", col = pal4[order(yrGDDs)][i])
}
legend("topleft", lwd = 1, col = pal4[order(yrGDDs)], legend = c(2003, "", 2005, "", 2007, "", 2009, "", 2011, "", 2013, "", 2015, "", 2017, "", 2019), cex = 0.45)
axis(side = 1, at = c(120, 151, 181, 212, 243, 273,304), labels = c("May-1", "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1", "Nov-1"))
mtext(side = 3, line = 0.5, "Ordered by average GDD")

yrGDDs = betas
plot(DOYs, 3.2*exp(DOYpreds$fit), type = "l", ylim = c(0,2.5), ylab = "", xlab = "", xaxt = "n", col = pal4[order(yrGDDs)][1])
mtext(side = 1, line = 2, "day of year")
mtext(side = 2, line = 2, "monarchs per survey mile")
for(i in 2:nyrs){
  DOYpreds = predict(m3c, se.fit = T, newdata = data.frame(year = years[i], DOY = DOYs))
  points(DOYs, 3.2*exp(DOYpreds$fit), type = "l", col = pal4[order(yrGDDs)][i])
}
legend("topleft", lwd = 1, col = pal4[order(yrGDDs)], legend = c(2003, "", 2005, "", 2007, "", 2009, "", 2011, "", 2013, "", 2015, "", 2017, "", 2019), cex = 0.45)
axis(side = 1, at = c(120, 151, 181, 212, 243, 273,304), labels = c("May-1", "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1", "Nov-1"))
mtext(side = 3, line = 0.5, "Ordered by slope of GDD vs DOY")



