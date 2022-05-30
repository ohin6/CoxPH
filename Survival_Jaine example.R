getHdata(prostate)


attach(prostate)

S = Surv(dtime, status != 'alive')

#full model
f.full = coxph(S ~ ., method = 'efron', data = prostate)
# remove age as a predictor
f.reduced = coxph(S ~ . -age, method = 'efron', data = prostate)
# compare whether there is loss in predictive power?
anova(f.reduced, f.full, test='LRT')

install.packages("Rtools")
require(Rtools)
install.packages("Mass")
require(MASS)

x = stepAIC(f.full)
