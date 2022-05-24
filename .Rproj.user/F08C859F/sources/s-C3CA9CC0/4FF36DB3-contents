# Survival analysis from Intro to statistical learning 

install.packages('ISLR2')
require(ISLR2)
require(tidyverse)
require(rms)



# look at Brain cancer data
attach(BrainCancer)
BrainCancer= tibble(BrainCancer)
print(BrainCancer)

table(sex)
table(diagnosis)
table(status)

# fit survival model

S = Surv(time, status)

# plot Kaplan-meier survival curve

plot(S , xlab = " Months ",
      ylab = " Estimated Probability of Survival ")

### stratify data by sex ####
fit.sex = npsurv(S ~ sex, type = 'fleming')

#fit kaplan-meier
plot(fit.sex , xlab = " Months ",
      ylab = "Estimated Probability of Survival", col = c(2,4)) +
      legend("bottomleft", levels(sex), col = c(2,4), lty = 1)

#alternative way to fit kaplen meier
survplot(fit.sex, n.risk = TRUE, conf = 'bands',
         label.curves = list(keys='lines'),
         levels.only = TRUE, col = c(4,2), col.fill = c(3,5))

### stratify data by diagnosis ###

fit.diagnosis = npsurv(S ~ diagnosis, type = 'fleming')

# fit kaplen meier
plot(fit.diagnosis, ylab = "Estimated Probability of Survival", col = c(2,4, 6,8)) +
legend("bottomleft", levels(diagnosis), col = c(2,4,6,8), lty = 1)


#alternative way to fit kaplen meier
survplot(fit.diagnosis, n.risk = TRUE, conf = 'bands',
         label.curves = list(keys='lines'),
         levels.only = TRUE, col = c(2,4,6,8))

#################
# log rank test #
#################
# How to ensure the observed variation in survival between male and female or
# between different diagnosis significant? Cant use t-test which looks are means
# survival times between different groups at time t, however this is invalid due
# to censored data. Instead a log rank test is used.

# compare survival of males and females using log rank test
log.rank= survdiff(S ~ sex)
print(log.rank)
##################
# Interpretation #
##################
# P value of 0.2. This means we cannot reject the null hypothesis that there is no
# difference between male and female survival and therefore there is evidence to
# suggest male and female have the same survival.
##################

# Compare survival of different diagnosis using log rank test
log.rank.diagnosis= survdiff(S ~ diagnosis)
print(log.rank.diagnosis)

##################
# Interpretation #
##################
# P value of <0.05. This means we can reject the null hypothesis therefore 
# indicating there is some evidence to suggest patients with different diagnosis
# have different survival rates.
##################

################################################################################
####################### fitting a regression model #############################
################################################################################
# Fitting a regression model to the data. We cant use nromal regression models
# because of censoring data. Instead we will use the COX PH model.

#########################
####### Cox model #######
#########################

##### What is hazard function f(h) #####
# This is the probability of an event happening (death) instantneously at time 
# t.
#  Ideally we would want to model survival as a function of a co-variate h(t|Xi)

# fitting cox model to different sex (h(t|x(sex)))
f.sex = coxph(S ~ sex)
summary(f.sex)

##################
# Interpretation #
##################


# no clear evidence for difference in survival

#### use additinal predictors ###
S_full = cph(S ~ sex + diagnosis + loc + ki + gtv + stereo)
S_full

#######
# results show HG glioma have a higher risk e^2.1546  = 8.62 (times more likely)
# while males are e^0.18 = 1.197 times higher hazard than females
# while one unit increase in Ki levels has a reduced risk of multiplier of 
# exp(-0.05) = 0.95. Therefore increased survival
########


S_estimate = cph(S ~diagnosis + ki)
S_estimate


modaldata = data.frame(
  diagnosis = levels(diagnosis),
  sex = rep("Female", 4),
  loc = rep("Supratentorial", 4),
  ki = rep(mean(ki), 4),
  gtv = rep(mean(gtv), 4),
  stereo = rep("SRT", 4)
)




