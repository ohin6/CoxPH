require(rms)
require(tidyverse)


##################
#import test data#
##################
getHdata(prostate)


###########
#Tidy data#
###########

#combine two levels in EKG column
levels(prostate$ekg)[levels(prostate$ekg) %in% c('old MI', 'recent MI')] = 'MI'
#new column that converts pf data into intger values (1 to 4)
prostate$pf.coded = as.integer(prostate$pf)
# Merge last two levels (levels 3 and 4 as 3)
levels(prostate$pf.coded) = c(levels(prostate$pf.coded)[1:3], levels(prostate$pf.coded)[3])

####################
# explore the data #
####################
# total number of patients 
nrow(prostate)

# number alive and dead
prostate %>% count(status)
sum(prostate$status != 'alive', na.rm=TRUE)

# is there missing values?????
tibble(prostate) %>% 
  filter(if_any(everything(), is.na)) # 27 rows with missing data


##########################
# predict missing values #
##########################
# Because there are missing data we can predict these values by imputing 
# them using transcan method
###########################

# predict missing values based on other scores
w = transcan(~ sz + sg + ap + sbp + dbp + age + wt + hg + ekg + pf + bm + hx,
             imputed = TRUE, data = prostate, pl = FALSE, pr = FALSE)

attach(prostate)
sz = impute(w, sz, data = prostate)
sg = impute(w, sg, data = prostate)
age = impute(w, age, data = prostate)
wt = impute(w, wt, data = prostate)
ekg = impute(w, ekg, data = prostate)


dd = datadist(prostate); options(datadist = 'dd')
units(dtime) = 'Month'

#############
# Questions #
#############
# 1) Does the transcan function predict missing values, by using the other 
# variables as a proxy. If so how do do we ensure these are sensible 
# predictions?
# 2) Is there a way to check if the transcan function worked, as I am unsure
# if this is the case.
# 3) what exactly is the datadist() function doing and when making models in the step below
# does this automatically refer to data stored in the datadist?
#############

#####################################################
# test whether full model (survival) is appropriate #
#####################################################

# create survival object based of time of death
S = Surv(dtime, status != 'alive')

# perform a cox proportional hazard model on all variables, use restricted cubic
# spline for certain variables
#
f = cph(S ~ rx + rcs(age, 4) + rcs(wt, 4) + pf + hx + rcs(sbp, 4) + rcs(dbp, 4) +
          ekg + rcs(hg, 4) + rcs(sg, 4) + rcs(sz, 4) + rcs(log(ap), 4) + bm)

print(f, latex = TRUE, coefs = FALSE)



###################
## Interpretation #
###################
# The likelihood statistic of 137.02 with a df of 37 is highly significant
# meaning that modelling is appropriate.
# 
# AIC for predicting model appropriateness is LR-2*df = 137.02-2*37 = 63.02
##################

#############
# Questions #
#############
# 1) How to interpret likelihood statistic? why is 137.02 good? what is bad?
# 2) book talks about shrinkage estimate is 0.74 by calculating  LS/100.2
#   where did they get the 100.2 value?
# 3) would the next step be to investigate which variables are significant in
#   this model and remove those which are not accordingly?




############################################################
# Next step  - data reduction to reduce degrees of freedom #
############################################################
# summarising data to reduce Df either by combining variables or reducing the
#  number of restricted cubic spline knots from 4 to 3
#
############################################################

############
# Question #
############
# 1) is the purpose of this to prevent over fitting???
############

heart = hx + ekg %nin% c('normal', 'benign') ###combining variables
label(heart) = 'Heart Disease Code'

map = (2*dbp+sbp)/3
label(map) = 'Mean Arterial Pressure/10' ##comnining vairables

dd = datadist(dd, heart, map) # add combined variables to datadist(dd)

################################################################################
## Compare reduced models with survival model and compare AIC with full model ##
################################################################################
# new model
f2 = cph(S ~ rx + rcs(age,3) + rcs(wt, 3) + pf.coded + heart + rcs(map,3) +
         rcs(hg, 4) + rcs(sg, 3) + rcs(sz, 3) + rcs(log(ap), 5) + bm,
         x = TRUE, y = TRUE, surv = TRUE, time.inc = 5*12)
print(f2, latex = TRUE, coefs = TRUE)
 
## Test that all variables are linear ###
latex(anova(f2), file='', label = 'tax:coccase-anova1')
#cant get above code to work

#########################
##### Interpretation ####
# reduction in number of df from 37 to 24 (-13) so with a new liklihood 
# ratio of 120.1
# AIC = 120.1-2*24 = 72.1
# This is an improvement  from the full model by 9.08, which had an AIC of 63.02
#
# However there are some variables associated with non linear effects 
# (note couldn't get the test to show this to work)
#########################


############################################
##### Checking Proportional Hazard #########
############################################
# Note there are multiple regression coefficients per predictor
# (for example different knots)
# As such compare model where all variables have only 1 df by using predict()


z = predict(f2, type = 'terms')
f.short = cph(S ~ z, x = TRUE, y= TRUE)

## Compute shoenfields residuals for each variable
# test the proportional hazard assumption using the cox.zph() function
phtest = cox.zph(f.short, transform = 'identity')
phtest
plot(phtest, var = 'rx')
# couldn't get to work - some reason phtest is not comparing different variables

######################
# Questions          #
######################
# 1) why are we not excluding non significant variables from model?
# 2) Why are they reducing degrees of freedom, does this make the test more powerful if so why?
# 3) How come they are looking at a single coefficient for each variable and
#    how come this approach is valid when a model has been built using numerous coefficients????
# 4) Why can't I plot residuals? If so what am I looking for?


#######################
# testing interactions
#######################
# this will look at whether dosage interacts with other variables
#######################

z.dose = z[, 1]  # z[,1] get column rx (dose)
z.other = z[, -1] #all other columns

## compare model
#
f.ia = cph(S ~ z.dose * z.other)
latex(anova(f2), file='', label = 'tab:coxcase-anova1')
# cant get above code to work

#############
# Questions #
#############
# 1) Why cant I get the above code to work?
# 2) Should I also be looking for interactions between other predictor 
#   variables???? Because surely there are numerous predictors that are not  
#   independent of each other.
#############

################################
# Describing predictor effects #
################################

# 
ggplot(predict(f2), sepdiscrete='vertical', nlevels=4, vnames='names')

