library(effectsize)
library(powerAnalysis)

# d, sample size, Bayes
# YP: Concerning the formula for calculation of effect size, some people use dz = t/sqrt(N) but it t_to_d 
# function actually uses dz = t / ???{df_{error}} (from the documentation of 'effectsize' package)
# I'm not sure what is the correct way to do the transformation.



p = 0.0447 
p = 0.0117
df = 22
t = qt(p/2, df) # t from p
effSize <- effectsize::t_to_d(t,df,paired = FALSE) # d from t
power.t(es = effSize$d/2, power = 0.9, sig.level = 0.01, type = "one", alternative = "left") 
