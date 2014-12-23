# Loading packages
library(lattice)       # Fancy graphics
library(nlme)          # Generalized linear mixed models 

# Reading data
setwd('~/Dropbox/quantumforest')  # Sets default working directory
un = read.csv('nzunemployment.csv', header = TRUE)

# Plotting first youth data and then adding adults
# as time series
with(un,
{
  plot(youth ~ q, type = 'l', ylim = c(0,30), col='red',
       xlab = 'Quarter', ylab = 'Percentage unemployment')
  lines(q, adult, lty=2, col='blue')
  legend('topleft', c('Youth', 'Adult'), lty=c(1, 2), col=c('red', 'blue'))
  abline(v = 90)
})


# Creating minimum wage policy factor
un$minwage = factor(ifelse(un$q < 90, 'Different', 'Equal'))

# And a scatterplot
xyplot(youth ~ adult, group=minwage, data = un, 
       type=c('p', 'r'), auto.key = TRUE)


# Linear regression accounting for change of policy
mod1 = lm(youth ~ adult*minwage, data = un)
summary(mod1)

# Centering continuous predictor
un$cadult = with(un, adult - mean(adult))
mod2 = lm(youth ~ cadult*minwage, data = un)
summary(mod2)

plot(mod2)    # Plots residuals for the model fit
acf(mod2$res) # Plots autocorrelation of the residuals

# Now we move to use nlme

# gls() is an nlme function when there are no random effects
mod3 = gls(youth ~ cadult*minwage, data = un)
summary(mod3)

# Adding autocorrelation
mod4 = gls(youth ~ cadult*minwage, correlation = corAR1(form=~1), data = un)
summary(mod4)

anova(mod3,mod4)
