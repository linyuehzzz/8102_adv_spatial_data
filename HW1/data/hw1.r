library(rgdal)
library(maps)
library(mapproj)

# set directory
setwd("C:\\Yue\\OneDrive - The Ohio State University\\osu\\SP20\\courses\\GEOG 8102\\homework\\Lab1\\data")
getwd()

# read shapefile
columbus.poly <- readOGR(dsn = ".", layer = "columbus")
columbus.df = as.data.frame(columbus.poly)
names(columbus.df)

# histogram
par(mfrow=c(1,3), pty="s")
hist(columbus.df$CRIME, nclass = 10, main = "Crime", xlab = "Crime (cases/1000 households)")
hist(columbus.df$HOVAL, nclass = 10, main = "Housing Value", xlab = "Housing Value ($1000)")
hist(columbus.df$INC, nclass = 10, main = "Household Income", xlab = "Household Income ($1000)")

# boxplot & qq-plot
par(mfrow=c(2,3), pty="s")
boxplot(columbus.df$CRIME, main="Crime")
boxplot(columbus.df$HOVAL, main="Housing Value")
boxplot(columbus.df$INC, main="Household Income")
qqnorm(columbus.df$CRIME)
qqline(columbus.df$CRIME)
qqnorm(columbus.df$HOVAL)
qqline(columbus.df$HOVAL)
qqnorm(columbus.df$INC)
qqline(columbus.df$INC)

# scatter plot matrix
plot(columbus.df[, c(9,7,8)], pch = 20) 

# regression
crime.lm = lm(columbus.df$CRIME ~ columbus.df$HOVAL + columbus.df$INC, data = columbus.df) 
summary(crime.lm)

# plot residuals
resids = resid(crime.lm)
fits = fitted(crime.lm)

graphics.off()
windows()
par(mfrow=c(2,2), pty = 's')
plot(columbus.df$HOVAL, resids, ylab = "Residuals (HOVAL)", xlab = "Housing Value ($1000)")
abline(h=0)
plot(columbus.df$INC, resids, ylab = "Residuals (INC)", xlab = "Household Income ($1000)")
abline(h=0)
hist(resids, xlab = "Residuals")
qqnorm(resids, ylab = "Residuals")
qqline(resids)