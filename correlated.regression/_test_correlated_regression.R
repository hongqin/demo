#test of AIC

require(nlme)
set.seed(2014)

N = 500
x = rnorm(N)
y = rnorm(N)
z = x + y + 2* x*y + rnorm(100)*0.05
tb=data.frame(cbind(x,y,z))

m1 = gls( z ~ x + y, data=tb)
m2 = gls( z ~ x + y, correlation=corAR1(), data=tb)
anova(m1,m2)
summary(m1)
summary(m2)


set.seed(2014)
N=500
x = rnorm(N)
error = rnorm(N)
rho = sqrt(0.5)
y= rho*x + sqrt(1-rho^2)*error
summary(lm(y ~ x ))


#z= x + y + 5*x*y + rnorm(N)/100 #Wrong, this is not correlated residues
z= x + y + 5*x*y + rnorm(N)/20
tb=data.frame(cbind(x,y,z))
m0 = lm( z ~ x + y, data=tb)
m1 = lm( z ~ x * y, data=tb)
anova(m0,m1)

m1 = gls( z ~ x * y, data=tb)
plot(m1)
m2 = gls( z ~ x * y, correlation=corAR1(form=~1), data=tb)
plot(m2)

anova(m1,m2)
summary(m1)
summary(m2)




