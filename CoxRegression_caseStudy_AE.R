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

attach(prostate)
####################
# explore the data #
####################
# total number of patients 
nrow(prostate)

# number alive and dead
prostate %>% count(status)
sum(prostate$status != 'alive', na.rm=TRUE)

table(status)

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

w = transcan(~ sz + sg + ap + sbp + dbp + age + wt + hg + ekg + pf + bm + hx,
             imputed = TRUE, trantab=TRUE, data = prostate, pl = FALSE, pr = FALSE)

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
##**Yes; the details of the algorithm can be found in Harrell's book 
##** if you add the trantab=TRUE argument to the transcan function:
w = transcan(~ sz + sg + ap + sbp + dbp + age + wt + hg + ekg + pf + bm + hx,
imputed = TRUE, trantab=TRUE,data = prostate, pl = FALSE, pr = FALSE)
##** you can then plot the transformations to see the strength of the associations of 
##*each variable with all the others. If the R^2 is small, it means the variable is mostly predicted
##*by their median or modal values from the marginal distribution
windows();ggplot(w,scale=TRUE)+theme(axis.text.x=element_text(size=6))
##**

# 2) Is there a way to check if the transcan function worked, as I am unsure
# if this is the case.
##* Look at both the above graph, and the output of summary(w); predictions with high
##* adjusted R^2 tend to be better (as they exploit information from other variables)

# 3) what exactly is the datadist() function doing and when making models in the step below
# does this automatically refer to data stored in the datadist?
##* datadist stores all the information on a given data set which are used by
##* other rms function (for example levels of a factor variable, units, minima, maxima etc.)
##* You may want to check Chapter 3 of Harrel's book for more details on missing data issues
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
##* The likelihood statistic tells us how good is the model at fitting the data.
##* It's distributed as chi2 with df degrees of freedom. The larger, the better the fit.
##* Note the output above provides Pr(> chi2) 0.0000, which is the p-value associated
##* to the LR; in this case, it's less than 0.00004, so strong rejection of the null hypothesis
##* of no effect of covariates on event outcome (death in this case)

# 2) book talks about shrinkage estimate is 0.74 by calculating  LS/100.2
#   where did they get the 100.2 value?
##* Note that the book's output doesn't match the output of R.
##* In general, shrinkage is calculated as (LR-df)/LR, so in the case of the book it's
##* (136.2-36)/136.2 = 100.2/136.2
##* In our case, it's (137.02-37)/137.02=0.73

###* this would mean that we would expect 27% of data is noise!


# 3) would the next step be to investigate which variables are significant in
#   this model and remove those which are not accordingly?
##* you could issue this command to see the strength of association of each predictor
##* and whether nonlinearities are warranted. Overall, there is evidence of nonlinearities
##* of some predictors, and also no evidence of association for some others (e.g. sbp, sg)
anova(f)



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
##** Yes. The aim should always build to be the most parsimonious model that attains
##*a reasonable predictive power. Combining predictors, simplifying transformations (e.g. 
##*reducing the number of spline knots) are possible strategies. 
##*However, Harrell himself underlines the fact that simplifying a model after fitting a more 
##*complex one is risky, as the statistics are not "honest" anymore (i.e. the statistics of the simpler model
##*do not account for the fact that OTHER statistics have already been performed, i.e. information has already been
##*extracted from the data!)
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
latex(anova(f2), file='', label = 'tab:coccase-anova1')
#cant get above code to work
##* This is a known bug. This instruction should be issued before the command:
##* However note that issuing this command will change ALL of the outputs to LaTeX
##* so beware!
options(prType='latex') # options(prType='plain') to restore the default
options(prType='plain')
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
##* To make it work, the phtest variable should be defined with:
phtest = cox.zph(f.short, transform = 'identity',terms=FALSE)
plot(phtest, var = 'rx')
abline(h=0, col = 2)
# couldn't get to work - some reason phtest is not comparing different variables

######################
# Questions          #
######################
# 1) why are we not excluding non significant variables from model?
##* Excluding variables would mean fitting another model. As previously commented, 
##* this will invalidate the test statistics

# 2) Why are they reducing degrees of freedom, does this make the test more powerful if so why?
##* if you reduce the df you make the model simpler, which may alleviate problems due
##* to overfitting, but rarely if ever improve predictive ability. But as discussed previously, it's dangerous
##* 

# 3) How come they are looking at a single coefficient for each variable and
#    how come this approach is valid when a model has been built using numerous coefficients????
##* The Schoenfeld residuals are designed for linear models, where each predictor has associated a single model parameter
##* so if a model's predictor requires more than one parameter (e.g. due to spline expansions) then the method
##* cannot be applied as-is

# 4) Why can't I plot residuals? If so what am I looking for?
##* The book forgot an argument (see above). What you are looking for are deviations
##* from an horizontal line, which means that the coefficient may be time-dependent, thus
##* violating the PH assumption 

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
options(prType='latex') # options(prType='plain') to restore the default
latex(anova(f2), file='', label = 'tab:coxcase-anova1')
options(prType='plain')
# cant get above code to work
##** See above for the fix

#############
# Questions #
#############
# 1) Why cant I get the above code to work?
##** See above
##*

# 2) Should I also be looking for interactions between other predictor 
#   variables???? Because surely there are numerous predictors that are not  
#   independent of each other.
##** In this case, since the goal of the analysis is to assess whether treatment is 
##* effective, it makes sense to test interaction of treatment with other predictors.
##* In general however, testing ALL possible interactions is not feasible (the number of parameters would be huge)
##* so it's advisable to use a priori knowledge on whether interactions make sense or not.
##* Note also that with regression analysis is more accurate to talk about "correlations" rather than "dependence"
##* (as statistical independence is very difficult if not impossible to assess.)
#############

################################
# Describing predictor effects #
################################

# 
ggplot(predict(f2), sepdiscrete='vertical', nlevels=4, vnames='names')

##################
# validate model #
##################

set.seed(1)
v = validate(f2, B=300)
latex(v, file = '')
print(v)

#############
# Questions #
#############

#* is validate showing the same thing as concordance index?
#* 

#############
# calibrate #
#############

cal = calibrate(f2, B=300, u=5*12, maxdim = 4)
plot(cal, subtitles = FALSE)
cal = calibrate(f2, cmethod = 'KM', u = 5*12, B = 120)
plot(cal, add = TRUE)

#############
# Questions #
#############

#* How do I interpret this?
#* 
