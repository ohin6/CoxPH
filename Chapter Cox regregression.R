require(rms)

############
# get data #
############
getHdata(kprats)

dd = datadist(kprats); options(datadist = 'dd')
units(t) = 'Day'
#################
# Plot Survival #
#################

# Create survival object
S = Surv(kprats$t, kprats$death)

# compute survival curve estimate
f = npsurv(S ~ group, type = 'fleming', data = kprats)

# Plot curve (exclude confidence intervals)
survplot(f, n.risk = TRUE, conf = 'none',
         label.curves = list(keys='lines'),
         levels.only = TRUE)

#######################
# create COX PH model #
#######################
# loop to through different methods for determining baseline survival function 
# in COX PH
for(meth in c('exact', 'breslow', 'efron')) {
  g = cph(S ~ group, method = meth, surv = TRUE, x = TRUE, y = TRUE, data = kprats)
}
print(g)

######################################################
# create other parametric Proportional hazard models #
######################################################

f.exp = psm(S ~ group, dist = 'exponential', data = kprats)
fw = psm(S ~ group, dist = 'weibull', data = kprats)
#paramateric proportial hazard form of AFT models?
phform = pphsm(fw)

# plot survival curve using parametric (f) and cox (g) proportional hazard models
co = gray(c(0, .8))
survplot(f, lty=c(1,1), lwd=c(1,3), col=co, label.curves = FALSE, conf = 'none')
survplot(g, lty=c(3,3), lwd=c(1,3), col=co, add = TRUE, label.curves = FALSE, conf.type = 'none')
legend(c(2, 160), c(.38,.54), c('Non parametric Estimates', 'Cox-Breslow Estimates'), lty = c(1,3), cex = .8, bty = 'n')
legend(c(2,160), c(.18,.34), cex = .8, c('Group 1', 'Group 2'), lwd=c(1,3), col=co, bty='n')

##################
# Interpretation #
##################

# good agreement between Cox model and the on parametric estimates 

