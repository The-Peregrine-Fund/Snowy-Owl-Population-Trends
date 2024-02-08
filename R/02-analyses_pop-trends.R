library ('nimble')
library ('nimbleHMC')
library ('MCMCvis')
library ('coda')
library('parallel')
load("data/data.Rdata")
set.seed(5757575)

# Function for posterior predictive checks
# to assess goodness-of-fit
plot.diag <- function(out, ratio=FALSE, lab=""){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  mx <- max(c(samps$dmape.rep, samps$dmape.obs))
  mn <- min(c(samps$dmape.rep, samps$dmape.obs))
  plot(jitter(samps$dmape.obs, amount=300), 
       jitter(samps$dmape.rep, amount=300),
       main=paste0("Mean absolute percentage error\nmodel\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(samps$dmape.rep > samps$dmape.obs),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  if (ratio==TRUE){
    # plot variance/mean ratio
    hist(samps$tvm.rep, nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=samps$tvm.obs, col="red")
    axis(1); axis(2)
  }
  return(list('Bayesian p-value'=bp1))
}

#***********************
#* Negative binomial model
#* Best fitting model
#***********************
run <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      for (j in 1:nsite){ 
        mu[j] ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
        beta[j] ~ dnorm(0, sd=5)
        r[j] ~ dexp(0.2) 
      }
      sigma.time ~ dexp(2)
      # likelihood
      for (t in 1:ntime){
        for (j in 1:nsite){
          p[t,j] <- r[j]/(r[j]+lambda[t,j])
          y[t,j] ~ dnegbin(p[t,j], r[j])
          # abundance
          log(lambda[t,j]) <- mu[j] + 
                              beta[j]*time[t] + 
                              eps.time[t] 
        } # j
        eps.time[t] ~ dnorm(0, sd=sigma.time)
      } # t

      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
        for (j in 1:nsite){
          c.obs[t,j] <- y[t,j] # observed counts
          c.exp[t,j] <- lambda[t,j] # expected counts adult breeder
          c.rep[t,j] ~ dnegbin(p[t,j],r[j]) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t,j] <- abs( ( (c.obs[t,j]) - (c.exp[t,j]) ) / (c.obs[t,j]+0.001) )
          dssm.rep[t,j] <- abs( ( (c.rep[t,j]) - (c.exp[t,j]) ) / (c.rep[t,j]+0.001) )
        } } # j # t
      dmape.obs <- sum(dssm.obs[1:ntime,1:nsite])
      dmape.rep <- sum(dssm.rep[1:ntime,1:nsite])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime,1:nsite])^2/mean(c.rep[1:ntime,1:nsite])
      tvm.obs <- sd(y[1:ntime,1:nsite])^2/mean(y[1:ntime,1:nsite])
    }
  ) # end model
  
  params <- c(    "sigma.time",
                  "mu", "beta", "r", "p",
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  y.inits <- array(NA, dim(datl$y))
  y.inits[is.na(datl$y)] <- rpois(sum(is.na(datl$y)), 2)
  inits <- function(){ list (y = y.inits,
                             sigma.time = rexp(1),
                             mu = runif(constl$nsite, -2, 2),
                             beta = runif(constl$nsite, -2, 2),
                             r = runif(constl$nsite)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  
  confhmc <- configureHMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
}

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)

params_nb <-c( "sigma.time",
                "mu", "beta", "r")

nb <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

save(out=nb, post=post, run, 
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")

# Check for convergence
MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

par(mfrow=c(1,1))
MCMCplot(object = nb, params = params_nb)
plot.diag(nb) # posterior predictive check

#***********************
#* Poisson model
#***********************
run <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      for (j in 1:nsite){ 
        mu[j] ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
        beta[j] ~ dnorm(0, sd=5)
      }
      sigma.time ~ dexp(2)
      # likelihood
      for (t in 1:ntime){
        for (j in 1:nsite){
          y[t,j] ~ dpois(lambda[t,j])
          # abundance
          log(lambda[t,j]) <- mu[j] + 
            beta[j]*time[t] + 
            eps.time[t] 
        } # j
        eps.time[t] ~ dnorm(0, sd=sigma.time)
      } # t
      
      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
        for (j in 1:nsite){
          c.obs[t,j] <- y[t,j] # observed counts
          c.exp[t,j] <- lambda[t,j] # expected counts adult breeder
          c.rep[t,j] ~ dpois(lambda[t,j]) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t,j] <- abs( ( (c.obs[t,j]) - (c.exp[t,j]) ) / (c.obs[t,j]+0.001) )
          dssm.rep[t,j] <- abs( ( (c.rep[t,j]) - (c.exp[t,j]) ) / (c.rep[t,j]+0.001) )
        } } # j # t
      dmape.obs <- sum(dssm.obs[1:ntime,1:nsite])
      dmape.rep <- sum(dssm.rep[1:ntime,1:nsite])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime,1:nsite])^2/mean(c.rep[1:ntime,1:nsite])
      tvm.obs <- sd(y[1:ntime,1:nsite])^2/mean(y[1:ntime,1:nsite])
    }
  ) # end model
  
  params <- c(    "sigma.time",
                  "mu", "beta", 
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  y.inits <- array(NA, dim(datl$y))
  y.inits[is.na(datl$y)] <- rpois(sum(is.na(datl$y)), 2)
  inits <- function(){ list (y = y.inits,
                             sigma.time = rexp(1),
                             mu = runif(constl$nsite, -2, 2),
                             beta = runif(constl$nsite, -2, 2)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  
  confhmc <- configureHMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
}

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)

params_pois <-c( "sigma.time",
               "mu", "beta")

pois <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

save(pois, post, run, 
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\pois.Rdata")

# Check for convergence
MCMCtrace(pois, params_pois, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

par(mfrow=c(1,1))
MCMCplot(object = pois, params = params_pois)
plot.diag(pois) # posterior predictive check

#***********************
#* Zero-inflated Poisson model
#***********************
run <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      for (j in 1:nsite){ 
        mu[j] ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
        beta[j] ~ dnorm(0, sd=5)
        psi[j] ~ dunif(0, 1) 
      }
      sigma.time ~ dexp(2)
      # likelihood
      for (t in 1:ntime){
        for (j in 1:nsite){
          y[t,j] ~ dpois(lambda.star[t,j])
          lambda.star[t,j] <- lambda[t,j]*z[t,j]
          z[t,j] ~ dbern(psi[j])
          # abundance
          log(lambda[t,j]) <- mu[j] + 
            beta[j]*time[t] + 
            eps.time[t] 
        } # j
        eps.time[t] ~ dnorm(0, sd=sigma.time)
      } # t
      
      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
        for (j in 1:nsite){
          c.obs[t,j] <- y[t,j] # observed counts
          c.exp[t,j] <- z[t,j]*lambda[t,j] # expected counts adult breeder
          c.rep[t,j] ~ dpois(lambda[t,j]*psi[j]) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t,j] <- abs( ( (c.obs[t,j]) - (c.exp[t,j]) ) / (c.obs[t,j]+0.001) )
          dssm.rep[t,j] <- abs( ( (c.rep[t,j]) - (c.exp[t,j]) ) / (c.rep[t,j]+0.001) )
        } } # j # t
      dmape.obs <- sum(dssm.obs[1:ntime,1:nsite])
      dmape.rep <- sum(dssm.rep[1:ntime,1:nsite])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime,1:nsite])^2/mean(c.rep[1:ntime,1:nsite])
      tvm.obs <- sd(y[1:ntime,1:nsite])^2/mean(y[1:ntime,1:nsite])
    }
  ) # end model
  
  params <- c(    "sigma.time",
                  "mu", "beta", "psi",
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  y.inits <- array(NA, dim(datl$y))
  y.inits[is.na(datl$y)] <- rpois(sum(is.na(datl$y)), 2)
  inits <- function(){ list (y = y.inits,
                             sigma.time = rexp(1),
                             mu = runif(constl$nsite, -2, 2),
                             beta = runif(constl$nsite, -2, 2),
                             psi = runif(constl$nsite)
  )}
  
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  # Use MCMC rather than HMC 
  # because HMC samplers are stuck
  confhmc <- configureMCMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
}

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)

params_zip <- c( "sigma.time", "psi",
                 "mu", "beta")

zip <- list(as.mcmc(post[[1]]), 
             as.mcmc(post[[2]]), 
             as.mcmc(post[[3]]),
             as.mcmc(post[[4]]))

save(zip, post, run, 
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\zip.Rdata")

# Check for convergence
MCMCtrace(zip, params_zip, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

par(mfrow=c(1,1))
MCMCplot(object = zip, params = params_zip)
plot.diag(zip) # posterior predictive check

#***********************
#* Negative binomial model
#* Best fitting model
#* having site and time random effects
#* in Appendix 3
#***********************
run <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      for (j in 1:nsite){ 
        mu[j] ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
        beta[j] ~ dnorm(0, sd=5)
        r[j] ~ dexp(0.2) 
      }
      sigma.time ~ dexp(2)
      sigma.site ~ dexp(2)
      # likelihood
      for (t in 1:ntime){
        for (j in 1:nsite){
          p[t,j] <- r[j]/(r[j]+lambda[t,j])
          y[t,j] ~ dnegbin(p[t,j], r[j])
          # abundance
          log(lambda[t,j]) <- mu[j] + 
            beta[j]*time[t] + 
            eps.time[t] + eps.site[j]
        } # j
        eps.time[t] ~ dnorm(0, sd=sigma.time)
      } # t
      for (j in 1:nsite){ eps.site[j] ~ dnorm(0, sigma.site) }    
      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
        for (j in 1:nsite){
          c.obs[t,j] <- y[t,j] # observed counts
          c.exp[t,j] <- lambda[t,j] # expected counts adult breeder
          c.rep[t,j] ~ dnegbin(p[t,j],r[j]) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t,j] <- abs( ( (c.obs[t,j]) - (c.exp[t,j]) ) / (c.obs[t,j]+0.001) )
          dssm.rep[t,j] <- abs( ( (c.rep[t,j]) - (c.exp[t,j]) ) / (c.rep[t,j]+0.001) )
        } } # j # t
      dmape.obs <- sum(dssm.obs[1:ntime,1:nsite])
      dmape.rep <- sum(dssm.rep[1:ntime,1:nsite])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime,1:nsite])^2/mean(c.rep[1:ntime,1:nsite])
      tvm.obs <- sd(y[1:ntime,1:nsite])^2/mean(y[1:ntime,1:nsite])
    }
  ) # end model
  
  params <- c(    "sigma.time", "sigma.site",
                  "mu", "beta", "r", "p",
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  y.inits <- array(NA, dim(datl$y))
  y.inits[is.na(datl$y)] <- rpois(sum(is.na(datl$y)), 2)
  inits <- function(){ list (y = y.inits,
                             sigma.time = rexp(1),
                             sigma.site = rexp(1),
                             mu = runif(constl$nsite, -2, 2),
                             beta = runif(constl$nsite, -2, 2),
                             r = runif(constl$nsite)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  
  confhmc <- configureHMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
}

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)
stopCluster(this_cluster)

params_nb <-c( "sigma.time", "sigma.site",
               "mu", "beta", "r")

nb <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

save(out=nb, post=post, run, 
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_site-re.Rdata")

# Check for convergence
MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

par(mfrow=c(1,1))
MCMCplot(object = nb, params = params_nb)
plot.diag(nb) # posterior predictive check
