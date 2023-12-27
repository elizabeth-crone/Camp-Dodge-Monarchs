rm(list = ls())


library(quantreg)
setwd("C:/Users/ecrone01/Box/Documents/Monarch IA/final R code and data")
ASSY = read.csv("ASSY.csv") # Asclepias syriaca records from IA, WI and MN, 1970-2018

##############################
# additional analysis of milkweed specimen data, for comparison with changes in monarch phenology

table(ASSY$stateProvince)
sum(!is.na(ASSY$startDayOfYear)) # sample size

taus = c(0.1, 0.5, 0.9)
m1.assy = rq(startDayOfYear ~ I(year-2000), tau = taus, data = ASSY)
summary(m1.assy, se = "boot")

myyears = 1970:2019
mypreds = predict(m1.assy, newdata = data.frame(year = myyears))
with(ASSY, plot(year, startDayOfYear, xlab = "", ylab = "", pch = 19, cex = 0.8))
for(i in 1:length(taus)){points(myyears, mypreds[,i], type = "l", col = "indianred", lwd = 3)}
mtext(side = 1, line = 2, "Year")
mtext(side = 2, line = 2, "Collection Date (DOY)")
