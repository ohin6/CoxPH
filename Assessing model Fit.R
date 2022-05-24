
# This script assessing model fit by an example from Regression modelling
# strategies book page 236.
# 
# Where the a logistic regression model has cofactors X1 (age) and X2(sex). Where
# the model assumes there is no interaction between age and sex this can be tested
# by "stratfying samples into deciles". 
# 
# Requirement- Minimum sample size
# 
# Because there needs to be at least three Groups with at least 20 samples
# 
# population needs to be 2*3*20
# (2 cofactors --> this may increase accordingly)
# 


# install package
require(rms)

#downlaod data from source
getHdata(acath)

### model data using linear regression model

#convert sex from 0 and 1 to male and female and strored as a factor datatype
acath$sex = factor(acath$sex, 0:1, c('male','female'))

dd = datadist(acath)
options(datadist = 'dd')

# logistic regression  model with response (signz) and predictors sex and age(4 restricted cubic splines knots)
f = lrm(sigdz ~ rcs(age, 4) * sex, data = acath)
sumarise(f)

### model data for interactions
w = function (...)
  with(acath, {
  plsmo(age , sigdz , group=sex , fun=qlogis , lty = 'dotted', add=TRUE, grid = TRUE)
  af = cut2(age ,g=10, levels.mean = TRUE)
  prop = qlogis(tapply(sigdz, list(af, sex), mean, na.rm=TRUE))
  agem = as.numeric(row.names(prop))
  lpoints(agem, prop[, 'female'], pch = 4, col= 'purple')
  lpoints(agem, prop[, 'male'], pch=2, col= 'darkblue')
}) # Figure 10.6

plot(Predict(f, age, sex), ylim = c(-2, 4), addpanel=w,
          label.curve = list(offset = unit (0.5, 'cm')))

## Interpretation of graph
# There seems to be a slight interaction between age and sex for females as the lines is slightly curved
# whereas there seems to be a linear relationship between age and sex for males.
