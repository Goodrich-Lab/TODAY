
g.prior.sel.logistic <- function(X=NULL, Y=NULL, U=NULL){
  
  g_prior_sel.model <- 
    "model {
  for(i in 1:N) {
    Y[i] ~ dbern(mu[i])
    logit(mu[i]) <- alpha + inprod(b[1:P], X[i,1:P]) + inprod(delta[1:Q], U[i,1:Q])
  }
  # prior on covariate effects
  for(q in 1:Q) { delta[q] ~ dnorm(0, 1.0E-06) }
  
  # prior on intercept
  alpha ~ dnorm(0, 1.0E-06)
  
  # prior on exposure effects
  beta[1:P] ~ dmnorm(mu.beta[1:P], T[1:P, 1:P])
  for(j in 1:P) {
    mu.beta[j] <- (1-gamma[j])*prop.mu.beta[j]
    b[j] <- beta[j]*gamma[j]
    gamma[j] ~ dbern(pi)
    for(k in 1:P) {
      T[j,k] <- gamma[j]*gamma[k]*XtX[j,k]/(G) + (1-gamma[j]*gamma[k])*equals(j,k)*pow(prop.sd.beta[j],-2)
    }
  }
  
  pi ~ dbeta(1,P)
  # pi ~ dbeta(P, 1)
  # semi-Bayes
  G <- w/(1-w)
  w <- .99   # w -> 0 shrink to common mean; as w -> inf toward the maximum likelihood estimate
  # Zellner and Siow prior on G
  # b0 <- 0.5*N
  # inv.G ~ dgamma(0.5, b0)
  # G <- 1/inv.G
  # w <- G/(G+1)
  # Hyper-g prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  # a <- 3
  # bw <- a/2 - 1
  # w~dbeta(1,bw)
  # G <- w/(1-w)
  # Hyper-g/n prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  # a <- 3
  # bw <- a/2 - 1
  # w~dbeta(1,bw)
  # G <- N*w/(1-w)
  # beta-prime 
  # G <- w/(1-w)
  # w ~ dbeta(bw, .25)
  # bw <- (N-P_m-1.5)/2
  # P_m <- sum(gamma[1:P])
  # g-estimation
  eta.low <- inprod(b[1:P], profiles[1,1:P])
  eta.high <- inprod(b[1:P], profiles[2,1:P])
  psi <- eta.high-eta.low 
}"

  N <- length(Y)
  P <- ncol(X)
  Q <- ncol(U)
  exposure.Names = colnames(X)
  
  # Create profiles for 1 SD
  profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P) 
  
  ### get the univariate result
  univariate.results <- t(sapply(1:P, FUN=function(p) {  # using index p facilitate write
    x <- as.matrix(X[,p])
    reg <- glm(Y~x, family=binomial)    # perform logistic regression
    s.reg <- summary(reg)                 # get the summary for the regression
    c.reg <- s.reg$coef["x",]             # select the coefficients for the exposure
    # write.table(t(c(exposure.Names[p], c.reg)), 
    #             file="ExposomeUnivariateResults.txt",
    #             append=ifelse(p==1, F, T),
    #             quote=F, 
    #             sep="\t", 
    #             col.names=ifelse(p==1, T, F), 
    #             row.names=F)
    return(c.reg)  # to avoid potential memory issues only return coefficients if small number of exposures
  }, simplify=T))
  univariate.results <- data.frame(exposure.Names, 
                                   univariate.results)
  
  ### g prior model result
  prop.mu.beta <- univariate.results$Estimate
  prop.sd.beta <- univariate.results$Std..Error
  XtX <- t(as.matrix(X))%*%as.matrix(X) 
  
  jdata <- list(N=N, Y=Y, X=X, U=U, P=P, Q=Q, XtX=XtX, profiles=profiles, prop.mu.beta=prop.mu.beta, prop.sd.beta=prop.sd.beta)
  var.s <- c("b", "beta","gamma", "psi", "pi", "G", "w")
  model.fit <- jags.model(file=textConnection(g_prior_sel.model), data=jdata, n.chains=1, n.adapt=4000, quiet=T)
  update(model.fit, n.iter=1000, progress.bar="none")
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=50000, thin=1, progress.bar="none")
  
  r <- summary(model.fit)
  var.names <- c("G", paste(exposure.Names, "b", sep="."), 
                 paste(exposure.Names, "beta", sep="."),
                 paste(exposure.Names, "gamma", sep="."),
                 "pi", "psi", "w")
  g_prior_sel.results <- data.frame(var.names, 
                                    r$statistics[,1:2],
                                    r$quantiles[,c(1,5)])
}
