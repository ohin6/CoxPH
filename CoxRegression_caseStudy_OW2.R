###################
# import libaries #
###################

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

#################################
# make prostate data set global #
#################################
attach(prostate)

################
# explore data #
################

describe(prostate)

# deaths
table(status)

# deaths caused by prostate
table(status %in% 'dead - prostatic ca')
table(status == 'dead - prostatic ca')


##########################
# resolving missing data #
##########################

# Determine missing data
prostate %>%
  filter(if_any(everything(), is.na)) # 27 rows with missing data

# impute missing data
w = transcan(~ sz + sg + ap + sbp + dbp + age + wt + hg + ekg + pf + bm + hx,
             imputed = TRUE, trantab=TRUE, data = prostate, pl = FALSE, pr = FALSE)

# explore reliabilty of imputed values
summary(w)

# Create variable for imputed values
sz = impute(w, sz, data = prostate)
sg = impute(w, sg, data = prostate)
age = impute(w, age, data = prostate)
wt = impute(w, wt, data = prostate)
ekg = impute(w, ekg, data = prostate)

#  bind data to datadist object, this is a utility function to direct functions from to this data object
dd = datadist(prostate); options(datadist = 'dd')
units(dtime) = 'Month'
dd

############################################
# Create Kaplan-Meier survival plot for rx #
############################################

S = Surv(dtime, status == 'dead - prostatic ca' )
km = npsurv(S ~ rx, type='fleming')
     
survplot(km, conf = 'bands',
         label.curves = list(keys='lines'),
         levels.only = TRUE)

# Is there a significant difference? - logrank test
survdiff(S~rx)

##* There is evidence to suggest there is a significant difference between 
##* treatment types with 1mg and 5mg different to placebo.



#####################################
# Create model using more variables #
#####################################
## estimate how many parameters can be used in model before over fitting

events = sum(status %in% 'dead - prostatic ca')
parameters = events/15

#* As there are 130 events we can estimate the number of parameters we can use 
#* before model is overfitting data by dividing number of events by 15,where we
#* get around 9 parameters.

# Number of parameters 
ncol(select(prostate, everything(), -patno, -status, -pf))

#* Before fitting model we need to reduce the number of variables, as the number
#* of variables measured exceeds the estimated 9 parameter limit. This can be
#* achieved by data reduction such as combining variables.

#### data reduction ####

heart = hx + ekg %nin% c('normal', 'benign') ###combining variables
label(heart) = 'Heart Disease Code'

map = (2*dbp+sbp)/3
label(map) = 'Mean Arterial Pressure/10' ##comnining vairables

dd = datadist(dd, heart, map) # add combined variables to datadist(dd)


####################
# create Cox model #
####################

f = cph(S ~ rx + rcs(age,3) + wt + map + pf.coded + hg + sz + ap + bm,
        x = TRUE, y = TRUE, surv = TRUE, time.inc = 5*12)

################
# Assess Model #
################

anova(f)

#* age, pf.coded, hg, sz and bm are statistically significant predictors.
#* pf.coded, has a large coefficient - (large effect)
#* 

# model over fitting ?
df = 12
LR = 135.55


shrinkage = (LR-df)/LR
(1-shrinkage)*100

# shrinkage less than 90 - data plotting around 9% noise

AIC = LR-2*df
print(AIC)  


################################################
# visualise predicted variables based on model #
################################################

ggplot(Predict(f), sepdiscrete='vertical', nlevels=4, vnames='names')

##############################
# interactions with dosage ? #
##############################
# predict values based on model
z =  predict(f, type = 'terms')
# predicted dosage
z.dosage = z[,1]
z.other = z[,-1]

f.ia= cph(S~z.dosage*z.other, x=TRUE, y=TRUE)
print(f.ia)

anova(f.ia)
#* whist dosage and age could be interacting, p-value is often exagerated so 
#* will not include in model. This is supported by global P-value of 0.3227
#* across all predictors.

###########################
# Check model assumptions #
###########################
# are the variable proportional - Shoenfelds residuals

z = predict(f, type = 'terms')
f.short = cph(S ~ z, x = TRUE, y= TRUE)

phtest = cox.zph(f.short, transform = 'identity', terms = FALSE)
phtest

# plot residuals

par(mfrow=c(3,3), 
    oma = c(2,2,2,2) + 0.1,
    mar = c(2,2,2,2) + 0.1)

plot(phtest, var = 'rx', main = 'Shoenfelds residual plot - rx') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'age', main = 'Shoenfelds residual plot - age') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'wt', main = 'Shoenfelds residual plot - wt') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'map', main = 'Shoenfelds residual plot - map') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'pf.coded', main = 'Shoenfelds residual plot - pf.coded') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'hg', main = 'Shoenfelds residual plot - hg') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'sz', main = 'Shoenfelds residual plot - sz') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'ap', main = 'Shoenfelds residual plot - ap') +
  abline(h = 1, col = 2, lty = 2)

plot(phtest, var = 'bm', main = 'Shoenfelds residual plot - bm') +
  abline(h = 1, col = 2, lty = 2)

par(mfrow=c(1,1))
    
#* all P-values are greater than 0.05 so we can reject the null hypothesis. This
#* is also visulaised in the residual plots where predicted coefficient lies
#* within the 95% confidence interval. Therefore we can say there is no evidence 
#* to suggest that the variables arenot proportional.

####################
# Validating model #
####################

# set seed generates a predictable set of random numbers for repeatable results
set.seed(1)
# Validating the f2 model using bootstrapping of 300
v = validate(f, B=300)
print(v)

###########################
### interpreting model ####
###########################
Dxy = 0.5096
c.index = Dxy*0.5+0.5
#* Concordance index equals 0.75. meaning that the model correctly ranks the 
#* predicted outcome to observed. Ideally this would want to be > 0.8.

#* shrinkage (slope) calcualted at 0.89 meaning about 10% of noise is being fitted
#* this is just in the boundary of acceptable.

#####################
# calibrating model #
#####################

#### calibrate at time 60 months ####

# set random seed
set.seed(1)
# set time interval where survival will be measured against, Note this was set before but for clarity we are reinstating this at 60 months (5years)
f = update(f, time.inc = 60)
# calibrate data set using bootstrap of 300 permutations (b), point in time to calibrate model 60 months (u) – this needs to be the same as time.inc, maximum dimensions (m) – not sure if this is number of predictors used? 
cal = calibrate(f, B=300, u=60, maxdim = 4)
# Plot results
plot(cal, subtitles = FALSE) +
  legend(0.55,0.2, legend = c('expected', 'observed'), col=c("black", "blue"), lty = 1, cex=0.8)
#add lines to plot to visualise over estimation of survival 
abline(h=0.5, col="grey", lty = 2, cex = 1.5)
abline(v=0.5, col = 'grey', lty = 2, cex = 0.5)
# Include confidence intervals for the estimates of fraction surviving
cal = calibrate(f, cmethod = 'KM', u = 5*12, B = 150)
plot(cal, add = TRUE)


### calibrate at 30 months ###
# set random seed
set.seed(1)
# set time interval where survival will be measured against, Note this was set before but for clarity we are reinstating this at 60 months (5years)
f = update(f, time.inc = 30)
# calibrate data set using bootstrap of 300 permutations (b), point in time to calibrate model 60 months (u) – this needs to be the same as time.inc, maximum dimensions (m) – not sure if this is number of predictors used? 
cal = calibrate(f, B=100, u=30, maxdim = 4)
# Plot results
plot(cal, subtitles = FALSE) +
  legend(0.5,0.3, legend = c('expected', 'observed'), col=c("black", "blue"), lty = 1, cex=0.8)
#add lines to plot to visualise over estimation of survival 
abline(h=0.5, col="grey", lty = 2, cex = 1.5)
abline(v=0.45, col = 'grey', lty = 2, cex = 0.5)
# Include confidence intervals for the estimates of fraction surviving
cal = calibrate(f, cmethod = 'KM', u = 5*12, B = 50)
plot(cal, add = TRUE)


##################
# Interpretation #
##################
#* Looking at the median at (50% surviving at 5 years), we see that the model 
#* correctly predicts survival (~50%). 
#* The confidence intervals supports this where it can be observed that the the
#* variability of the predicted values are statistically within the variability 
#* of the observed survival.

#* It is also worth while looking at this at different survival time points. Here
#* we look at 24 months. This time the median survival fraction the model appears
#* to slightly underestimate survival (~45%).
#* I am not sure how to interret the confidence intervals.


####################
# summarising data #
####################
plot(summary(f))

#* 1mg and 5mg oestrogen treatment show better survival than 0.2mg and placebo

###############
# Conclusions #
###############

#* 1) Model consisting of 12 predictors
#* 2) A shrinkage index of 0.91 suggests model doesn't excessively over fit data
#* 3) rx, sz and bm predictors are highly significant with large coefficients
#* 4) All predictors are propotional - checked by shoenfelds residuals
#* 5) validated model using bootstrapping method, c-index of 0.75 indicates model
#* is not perfect (ideally should be >0.8), however not dreadful.
#* 6) calibration over a 60 month period showed model accurately predicted median 
#* survival.


#############
# Questions #
#############
#* 1) I haven't checked the assumption that data is linear
#* 2) How to interpret calibration at 30 month (low sample size???)
#* 3) Can you test new predictors against predicted model without jeopardising the 
#* 'honesty' of the data. i.e if not significant then doesnt need to be kept in
#*  the model


