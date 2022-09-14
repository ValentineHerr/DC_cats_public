NimModel <- nimbleCode({
  #detection function priors
  lam0.beta0~dnorm(0,1.0E-4)
  lam0.beta.sex~dnorm(0,1.0E-4)
  sigma.beta0~dnorm(0,1.0E-4)
  sigma.beta.sex~dnorm(0,1.0E-4)
  for(i in 1:n.trap.cov){
    lam0.beta.trap[i]~dnorm(0,1.0E-6)
  }
  r~dunif(0,25) #NB dispersion parameter
  #data augmentation prior
  psi~dunif(0,1)
  thin.beta0~dlogis(0,1) #intercept - black/grey,bicolor false, long hair false
  for(i in 1:4){
    thin.beta.coat[i]~dnorm(0,1.0E-4)
  }
  thin.beta.bi~dnorm(0,1.0E-4) #bicolor TRUE offset
  thin.beta.hair~dnorm(0,1.0E-4) #long hair TRUE offset
  
  #categorical ID covariate priors
  for(m in 1:n.cat){
    alpha[m,1:n.levels[m]] <- 1 # prior parameters
    gammaMat[m,1:n.levels[m]]~ddirch(alpha[m,1:n.levels[m]])
  }
  
  #likelihoods (except for s priors)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    for(m in 1:n.cat){
      G.true[i,m]~dcat(gammaMat[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    male.indicator[i] <- GetMale(G.true[i,3]) #convert male/female G.true to an indicator for males (female=0,male=1)
    for(j in 1:J){
      log(lam0[i,j]) <- lam0.beta0 + inprod(lam0.beta.trap[1:n.trap.cov],covs.trap[j,1:n.trap.cov]) + lam0.beta.sex*male.indicator[i]
    }
    log(sigma[i]) <- sigma.beta0 + sigma.beta.sex*male.indicator[i]
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma[i], lam0=lam0[i,1:J], z=z[i])
    p[i,1:J] <- r/(r+lam[i,1:J])
    y.true[i,1:J] ~ dNBVector(p=p[i,1:J],r=r*K1D[1:J],z=z[i]) #vectorized obs mod
    coat.indicator[i,1:4] <- GetCoat(G.true[i,1])
    bi.indicator[i] <- GetBi(G.true[i,2])
    hair.indicator[i] <- GetHair(G.true[i,4])
    logit(theta.i[i]) <- thin.beta0 + inprod(coat.indicator[i,1:4],thin.beta.coat[1:4]) + 
      bi.indicator[i]*thin.beta.bi + hair.indicator[i]*thin.beta.hair
    y.ID[i,1:J] ~ dBinomialVector(theta.i[i], y.true[i,1:J],capcounts=capcounts[i])  # Model for ID process
  }
  #calculate number of marked and unmarked inds captured and abundance
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J])
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],G.latent=G.latent[1:M,1:n.cat])
  N <- sum(z[1:M])
})# end model
