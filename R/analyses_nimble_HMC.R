library (nimble)
library (nimbleHMC)
library (MCMCvis)
library (coda)
library (AHMbook)
library(parallel)
load("data/data.Rdata")
set.seed(5757575)

# function to plot posterior predictive checks
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
#***********************
run <- function(seed, datl, constl){
  library(nimble)
  library(coda)
  library(nimbleHMC)
  code <- nimbleCode(
    {
      # priors
      for (j in 1:nsite){ 
        mu[j] ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
        beta[j] ~ dnorm(0,sd=5)
        r[j] ~ dexp(0.2) 
      }
      
      # likelihood
      for (t in 1:ntime){
        for (j in 1:nsite){
          p[t,j] <- r[j]/(r[j]+lambda[t,j])
          y[t,j] ~ dnegbin(p[t,j], r[j])
          # abundance
          log(lambda[t,j]) <- mu[j] + beta[j]*time[t]
        } # j
        lam.tot[t] <- sum(lambda[t,1:nsite])
      } # t
      
      # derived parameters
      for (t in 2:(ntime-1)){
        gamma.tot[t] <- lam.tot[t] / (lam.tot[t-1] + 0.1)
        for (j in 1:nsite){
          gamma.site[t,j] <- lambda[t,j] / (lambda[t-1,j] + 0.1) 
        }}
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
  
  params <- c( "mu", "beta", "r", "p",
    "gamma.site", "gamma.tot",
    "lambda", "lam.tot",
    "dssm.obs", "dssm.rep",
    "dmape.obs", "dmape.rep",
    "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  y.inits <- array(NA, dim(datl$y))
  y.inits[is.na(datl$y)] <- rpois(sum(is.na(datl$y)), 2)
  inits <- function(){ list (y = y.inits,
                             mu = runif(constl$nsite, -2, 2),
                             beta = runif(constl$nsite, -2, 2),
                             r = runif(constl$nsite)
  )}
  
  #n.chains=1; n.thin=1; n.iter=2000; n.burnin=1000
  n.chains=1; n.thin=5; n.iter=10000; n.burnin=5000
  #n.chains=1; n.thin=1000; n.iter=2000000; n.burnin=1000000
  
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

params_nb <-c( "mu", "beta", "r")

nb <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

par(mfrow=c(1,1))
MCMCplot(object = nb, params = params_nb[1:2])
MCMCplot(object = nb, params = params_nb[c(3)])

MCMCplot(object = nb, params = 'mu', HPD=TRUE)
MCMCplot(object = nb, params = 'r', HPD=TRUE)
MCMCplot(object = nb, params = 'beta', HPD=TRUE)
MCMCplot(object = nb, params = 'lam.tot', HPD=TRUE, 
         ylim=c(0,200), horiz=FALSE)

plot.diag(nb)

save(out=nb, post=post,
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")

# try running without Wrangel
# data manip remove Wrangel
datl2 <- list(y=datl$y[,-5])
constl2 <- constl
constl2$nsite <- 4

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl2, 
                  const=constl2)
stopCluster(this_cluster)

nb <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

save(out2=nb, post2=post,
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_noWrangel.Rdata")
