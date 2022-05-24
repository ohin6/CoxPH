require(rms)
require(ISLR2)
require(tidyverse)


attach(Publication)
Publication

Publication = Publication %>%
  mutate(logbudget = log(budget)) %>%
  mutate(logimpact = log(impact)) %>%
  mutate(logimpact = replace(logimpact, logimpact == -Inf, 0))


#survival
S = Surv(time, status)

# kaplan meier
fit.posres = npsurv(S ~ posres, type = 'fleming')

survplot(fit.posres, n.risk = TRUE, conf = 'bands',
         label.curves = list(keys='lines'),
         levels.only = TRUE, col = c(4,2), col.fill = c(3,5))

# are these significant? Plot a log rank

logrank.posres = survdiff(S ~ posres)
print(logrank.posres)

#### P value of > 0.05 means we cannot reject the null hypothesis, therefore 
# there is evidence to suggest there is no difference in publication date  
# between postive and non postive results 

# However when we consider this same parameter when including all predictors
# using the the cox hazard proportional model. Here we see that posres is 
# significant with a p-value of 0.0025 and a copeficient of 0.55. Therefore 
# a publication with a postive result has a e^0.55 = 1.74 (times) more likely 
# to published than a negative result. 

f.full = coxph(Surv(time,status) ~ . -budget -impact, data = Publication)
f.nearfull = coxph(Surv(time,status) ~ . - mech -budget -impact, data = Publication)

############
# QUestion #
############

# 1) because there is 19 degrees of freedom would significance follow bonoforoni
# where 0.05/19 = 0.0026?


####################
# shrinkage methods #
####################
# this is a method of getting best fit of model by restricting predictor 
# coefficents to 0. How this works or what this means is a mystery. Not sure how to do this

###############
# Concordance #
###############
# similar to area under the curve AUC/ROC analysis Concordance indec determines 
# the performance of two classifiers but allows for censored data.
# The concordance index can be calculated by summarising the coxph output.

summary(f.full)
summary(f.nearfull)

# concordance = 0.805 whch indicates a 80% accuracy


########################
# Selecting variables #
########################
# does adding variables improve model or are there interactions between variables?

# Will dropping mechanisms variable signicantly reduce the predictive value of
# the model?
anova(f.nearfull, f.full, test = "LRT")

##Answer###
# Because the of high P-value we can reject null hypothesis, therefore we can 
# remove mech as a predictor without losing predictive power.


#####################################
# can we remove other predictors??? #
#####################################

#### Remove multi #####
f.nearfull2 = coxph(Surv(time, status) ~ . - multi -mech -impact -budget, data = Publication)

# compare predictive power of two models using anova
anova(f.nearfull2,f.nearfull, test = 'LRT')


##Answer###
# Because the of high P-value we can reject null hypothesis, therefore we can 
# remove multi as a predictor without losing predictive power.

##### Remove sampsize #####
f.nearfull3 = coxph(Surv(time, status) ~ . - multi -mech -impact -budget-sampsize , data = Publication)

# compare predictive power of two models using anova
anova(f.nearfull3,f.nearfull2, test = 'LRT')

##Answer###
# Because the of high P-value we can reject null hypothesis, therefore we can 
# remove sampsize as a predictor without losing predictive power. 

##### Remove clinend? #####
f.nearfull4 = coxph(Surv(time, status) ~ . - multi -mech -sampsize -budget -impact -clinend, data = Publication)

# compare predictive power of two models using anova
anova(f.nearfull4,f.nearfull3, test = 'LRT')


##Answer###
# Because the of low P-value we can not reject null hypothesis, therefore 
# there is evidence to suggest that removing clinend as a predictor would reduce
# the predictive power of this 

##################
# Interpretation #
##################
# By iteratively removing predictors and comparing predictive power of the new
# model using anova LRT test we have decerned that posres, budget and
# impact are significant predictors in this model.


################################################################################
#################### Checking assumptions of model #############################
################################################################################
# Need to check linearity and proportionality of model are met

###################
# check linearity #
###################
# to check that each variable is making a linear contribution. This can be 
# visualized by looking at residual plots such as 'martingale' or 'deviance'.

## Create function to look at goodness of fit 

smoothSEcurve <- function(yy, xx) {
  # use after a call to "plot"
  # fit a lowess curve and 95% confidence interval curve
  # make list of x values
  xx.list <- min(xx) + ((0:100)/100)*(max(xx) - min(xx))
  # Then fit loess function through the points (xx, yy)
  # at the listed values
  yy.xx <- predict(loess(yy ~ xx), se=T,
                   newdata=data.frame(xx=xx.list))
  lines(yy.xx$fit ~ xx.list, lwd=2)
  lines(yy.xx$fit -
          qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  lines(yy.xx$fit +
          qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
}


result.0.coxph <- coxph(Surv(time, status) ~ 1, data = Publication)
rr.0 <- residuals(result.0.coxph, type="martingale")

###### budget
plot(rr.0 ~ budget, data = Publication) + 
  smoothSEcurve(rr.0, budget) +
  title("Martingale residuals\nversus budget") +
  abline(h=0, col = 2)

### budget data needs transforming to log form ####
plot(rr.0 ~ logbudget, data = Publication) + 
smoothSEcurve(rr.0, Publication$logbudget) +
title("Martingale residuals\nversus budget") +
abline(h=0, col = 2)

#######################################
# budget needs further transformation #
#######################################


###### impact ######
plot(rr.0 ~ impact, data = Publication) + 
  smoothSEcurve(rr.0, impact) +
  title("Martingale residuals\nversus budget") +
  abline(h=0, col = 2)


### budget data needs transforming to log form ####
plot(rr.0 ~ logimpact, data = Publication) + 
  smoothSEcurve(rr.0, Publication$impact) +
  title("Martingale residuals\nversus budget") +
  abline(h=0, col = 2)


#################
# same method ? #
#################

plot(predict(f.nearfull3), residuals(f.nearfull3, type = 'martingale'),
     xlab = 'fitted values', ylab = 'Martingale residuals',
     main = 'Residual Plot', las = 1) +
  abline(h=0) +
  lines(smooth.spline(predict(f.nearfull3),
                      residuals(f.nearfull3, type = 'martingale')), col = 'red')


plot(predict(coxph(Surv(time, status) ~ log(budget))), residuals(coxph(Surv(time, status) ~ log(budget)), type = 'martingale'),
     xlab = 'fitted values', ylab = 'Martingale residuals',
     main = 'Residual Plot', las = 1) +
  abline(h=0) +
  lines(smooth.spline(predict(coxph(Surv(time, status) ~ log(budget))),
                      residuals(coxph(Surv(time, status) ~ log(budget)), type = 'martingale')), col = 'red')


#################
# very confused #
#################
# Do I check per variable or as a whole model?
# why is there so many knots and how do i interpret it?
# what do i do with categorical data?
# budget needs to be transformed to log






###########################################
# Checking proportional hazard assumption #
###########################################
# test prop hazard using schoenfld test.
# the null hypothesis is that the hazards are not proportional

cox.zph(f.nearfull3)

# The results show that for each variable P value is greater than 0.05. As such
# the null hypothesis can be rejected and therefore the is evidence to suggest
# that the variables are proportional.

#Posres - doesn't appear proportional
plot(cox.zph(f.nearfull3)[1])
abline(h=0, col = 2)

#Cline - appears proportional
plot(cox.zph(f.nearfull3)[2])
abline(h=0, col = 2)

#budget - appears proportional
plot(cox.zph(f.nearfull3)[3])
abline(h=0, col = 2)

#impact - doesnt appears proportional
plot(cox.zph(f.nearfull3)[4])
abline(h=0, col = 2)







plot(rr.0 ~ posres)
smoothSEcurve(rr.0, 1)
title("Martingale residuals\nversus budget")
