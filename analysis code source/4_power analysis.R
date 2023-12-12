library(powerSurvEpi)

ssizeEpiCont.default(power = 0.8,
                     alpha = 0.008,
                     theta = 1.2295,  # postulated hazard ratio.
                     sigma2 = 1, #variance of the covariate of interest.
                     rho2 = 0,
                     psi = 0.5 #proportion of subjects who get the disease
                     )



library(powerMediation)

pwr.r <- minEffect.VSMc(n=572, power=.8, sigma.m=1,sigma.e=1,
                        corr.xm=0,alpha = 0.05,verbose = TRUE)

pwr.r <- minEffect.VSMc(n=572, power=.8, sigma.m=1,sigma.e=1,
                        corr.xm=0.2,alpha = 0.05,verbose = TRUE)

pwr.r <- minEffect.VSMc(n=572, power=.8, sigma.m=1,sigma.e=1,
                        corr.xm=0.4,alpha = 0.05,verbose = TRUE)


# For Danias grant
pwr.r <- minEffect.VSMc(n=700, power=.8, sigma.m=1,sigma.e=1,corr.xm=.1,alpha = 0.05,verbose = TRUE)

pwr.r <- minEffect.VSMc(n=700, power=.8, sigma.m=1,sigma.e=1,corr.xm=0.2,alpha = 0.05,verbose = TRUE)

pwr.r <- minEffect.VSMc(n=700, power=.8, sigma.m=1,sigma.e=1,corr.xm=0.4,alpha = 0.05,verbose = TRUE)
