library(nimble)
library(MCMCvis)
library (MASS)
library (pscl)
# create a function to simulate
# a time series of abundance 
# for an irruptive species 
sim.ts <- function( period = 4,
                    ntime = 30,
                    mn.psi = -1.38, # about 1/5
                    mn.abund = 10,
                    beta0 = -0.7,
                    beta1 = 2,
                    beta2 = 2,
                    plot.ts = TRUE){
          yr2 <- 1:ntime
          yr <- (yr2-(ntime/2+0.5))/(ntime/2-0.5)
          period2 <- (yr[2]-yr[1])*period
          # decreasing abundance over time
          lambda <- exp( log(mn.abund) + beta0*yr ) 
          # irruptive cycles as a zero-inflation parameter
          psi <- plogis(mn.psi +
                          beta1 * sin(2*pi*(yr/period2)) +
                          beta2 * cos(2*pi*(yr/period2)))
          # simulate abundance
          z <- rbinom(length(psi), size=1, prob=psi)
          mu <- lambda*z
          abund <- rpois(length(mu), mu)  
          if(plot.ts==TRUE){
            plot(abund, type="l", 
                 xlab="Year", ylab="Abundance")
          }
          return(list(datl = list(y = abund),
                      constl = list(ntime = ntime,
                                    yr = yr, 
                                    period=period2*5),
                       truth = list(mn.psi = mn.psi,
                            mn.abund = mn.abund,
                            beta0 = beta0,
                            beta1 = beta1,
                            beta2 = beta2)
                      ))
}

#*******************
#* Frequentist version
#* Fast AF boi
#* Negative binomial
#******************* 
set.seed(123456)
# different population trends
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lci=NA, uci=NA)
war <- rep(FALSE, 5000)
ind <- 1
beta.scenarios <- c(-0.7, -0.35, 0, 0.35, 0.7)
for (j in 1:length(beta.scenarios)){
  for (i in 1:1000){
    dat <- sim.ts(beta0=beta.scenarios[j], plot.ts=FALSE)
    dfm <- data.frame(y=dat$datl$y, yr=dat$constl$yr)
    fit <- tryCatch(glm(y ~ 1 + yr, family = negative.binomial(2), 
                        data = dfm),
                    warning= function(w) {  return(TRUE) } ) 
    if (is.logical(fit)){ war[ind] <- fit
    ind <- ind+1 
    next} else{  
      out <- summary(fit)
      uci <- out$coefficients[2,1] + 1.96*out$coefficients[2,2]
      lci <- out$coefficients[2,1] - 1.96*out$coefficients[2,2]
      estimate <- out$coefficients[2,1]
      truth <- dat$truth$beta0
      relbias <- ifelse(truth==0, 
                        ((estimate+1)-(truth+1))/(truth+1), 
                        (estimate-truth)/truth)
      coverage <- ifelse(lci>=truth | uci<=truth, 0, 1)
      df[ind, 1:6] <- c(estimate, truth, relbias, coverage, lci, uci)
      ind <- ind + 1
    }
  }}
df1 <- df[war!=TRUE,]

#* test different levels of abundance
set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lci=NA, uci=NA, abund=NA)
war <- rep(FALSE, 4000)
ind <- 1
abund.scenarios <- c(0.5, 1, 2, 3)
for (j in 1:length(abund.scenarios)){
  for (i in 1:1000){
    dat <- sim.ts(mn.abund=abund.scenarios[j], plot.ts=FALSE)
    dfm <- data.frame(y=dat$datl$y, yr=dat$constl$yr)
    fit <- tryCatch(glm(y ~ 1 + yr, family = negative.binomial(2), 
                        data = dfm),
                    warning= function(w) {  return(TRUE) } ) 
    if (is.logical(fit)){ war[ind] <- fit
    ind <- ind+1 
    next} else{      
      out <- summary(fit)
      uci <- out$coefficients[2,1] + 1.96*out$coefficients[2,2]
      lci <- out$coefficients[2,1] - 1.96*out$coefficients[2,2]
      estimate <- out$coefficients[2,1]
      truth <- dat$truth$beta0
      relbias <- ifelse(truth==0, 
                        ((estimate+1)-(truth+1))/(truth+1), 
                        (estimate-truth)/truth)
      coverage <- ifelse(lci>=truth | 
                           uci<=truth, 0, 1)
      df[ind, 1:7] <- c(estimate, truth, relbias, coverage, lci, uci, abund.scenarios[j])
      ind <- ind + 1
    }
  }}
# need to remove mods that didn't converge!
df2 <- df[war!=TRUE,]

#* test different lengths of time
set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lci=NA, uci=NA, time=NA)
war <- rep(FALSE, 3000)
ind <- 1
time.scenarios <- c(5, 10, 15)
for (j in 1:length(time.scenarios)){
  for (i in 1:1000){
    dat <- sim.ts(ntime=time.scenarios[j], plot.ts=FALSE)
    dfm <- data.frame(y=dat$datl$y, yr=dat$constl$yr)
    fit <- tryCatch(glm(y ~ 1 + yr, 
                        family = negative.binomial(2), 
                        data = dfm),
                    warning= function(w) {  return(TRUE) } ) 
    if (is.logical(fit)){ war[ind] <- fit
    ind <- ind+1 
    next} else{      
      out <- summary(fit)
      uci <- out$coefficients[2,1] + 1.96*out$coefficients[2,2]
      lci <- out$coefficients[2,1] - 1.96*out$coefficients[2,2]
      estimate <- out$coefficients[2,1]
      truth <- dat$truth$beta0
      relbias <- ifelse(truth==0, 
                        ((estimate+1)-(truth+1))/(truth+1), 
                        (estimate-truth)/truth)
      coverage <- ifelse(lci>=truth | 
                           uci<=truth, 0, 1)
      df[ind, 1:7] <- c(estimate, truth, relbias, coverage, lci, uci, time.scenarios[j])
      ind <- ind + 1
    }
  }}
# need to remove mods that didn't converge!
df3 <- df[war!=TRUE,]

# Plots
par(mfrow=c(1,3))
boxplot(df1$bias~df1$truth, 
        xlab=expression(paste("Trend (", delta[1],")")),
        ylab="Relative bias of trend")
abline(h=0, lty=2, lwd=2)
mtext("A", cex=1.5, side=3, outer=F)

boxplot(df2$bias~df2$abund, 
        xlab=expression(paste("Mean abundance (", delta[0],")")), 
        ylab="")
abline(h=0, lty=2, lwd=2)
mtext("B", cex=1.5, side=3, outer=F)

boxplot(df3$bias~df3$time, 
        xlab="Duration of monitoring",
        ylab="")
abline(h=0, lty=2, lwd=2)
mtext("C", cex=1.5, side=3, outer=F)

# plot a sample of time series
dat <- list()
for (i in 1:1000){
  dat[[i]] <- sim.ts(beta0=beta.scenarios[1], plot.ts=FALSE)[[1]]$y
}
samps <- sample(1:1000, size=6, replace = FALSE)
df <- do.call(rbind, dat[samps])
library (reshape2)
library(ggplot2)
rownames(df) <- samps
ldf <- melt(df)
  
ggplot(ldf)+
  geom_line(aes(x=Var2, y=value)) +
  facet_wrap('Var1') +
  ylab("Simulated count") +
  xlab("Time")


# tables
# different population trends
tapply(df1$bias, df1$truth, mean) |> round (2) 
tapply(df1$bias, df1$truth, sd) |> round (2)
tapply(df1$cov, df1$truth, mean) |> round (2)
#* test different levels of abundance
tapply(df2$bias, df2$abund, mean) |> round (2)
tapply(df2$bias, df2$abund, sd) |> round (2)
tapply(df2$cov, df2$abund, mean) |> round (2)
#* test different lengths of time
tapply(df3$bias, df3$time, mean) |> round (2)
tapply(df3$bias, df3$time, sd) |> round (2)
tapply(df3$cov, df3$time, mean) |> round (2)

#* Zero-inflated 
set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lci=NA, uci=NA)
war <- rep(FALSE, 5000)
ind <- 1
beta.scenarios <- c(-0.7, -0.35, 0, 0.35, 0.7)
for (j in 1:length(beta.scenarios)){
  for (i in 1:1000){
    dat <- sim.ts(beta0=beta.scenarios[j], plot.ts=FALSE)
    dfm <- data.frame(y=dat$datl$y, yr=dat$constl$yr)
    fit <- tryCatch(fit <- zeroinfl(formula = y ~ 1 + yr |
                                      1,
                                    dist    = "poisson",
                                    link    = "log",
                                    data    = dfm),
                    warning= function(w) {  return(TRUE) },
                    error= function(w) {  return(TRUE) }) 
    if (is.logical(fit)){ war[ind] <- fit
    ind <- ind+1 
    next} else{  
      out <- summary(fit)
      uci <- out$coefficients$count[2,1] + 1.96*out$coefficients$count[2,2]
      lci <- out$coefficients$count[2,1] - 1.96*out$coefficients$count[2,2]
      estimate <- out$coefficients$count[2,1]
      truth <- dat$truth$beta0
      relbias <- ifelse(truth==0, 
                        ((estimate+1)-(truth+1))/(truth+1), 
                        (estimate-truth)/truth)
      coverage <- ifelse(lci>=truth | uci<=truth, 0, 1)
      df[ind, 1:6] <- c(estimate, truth, relbias, coverage, lci, uci)
      ind <- ind + 1
    }
  }}
df1 <- df[war!=TRUE,]

#* test different levels of abundance
set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lci=NA, uci=NA, abund=NA)
war <- rep(FALSE, 4000)
ind <- 1
abund.scenarios <- c(0.5, 1, 2, 3)
for (j in 1:length(abund.scenarios)){
  for (i in 1:1000){
    dat <- sim.ts(mn.abund=abund.scenarios[j], plot.ts=FALSE)
    dfm <- data.frame(y=dat$datl$y, yr=dat$constl$yr)
    fit <- tryCatch(fit <- zeroinfl(formula = y ~ 1 + yr |
                                      1,
                                    dist    = "poisson",
                                    link    = "log",
                                    data    = dfm),
                    warning= function(w) {  return(TRUE) },
                    error= function(w) {  return(TRUE) }) 
    if (is.logical(fit)){ war[ind] <- fit
    ind <- ind+1 
    next} else{      
      out <- summary(fit)
      uci <- out$coefficients$count[2,1] + 1.96*out$coefficients$count[2,2]
      lci <- out$coefficients$count[2,1] - 1.96*out$coefficients$count[2,2]
      estimate <- out$coefficients$count[2,1]
      truth <- dat$truth$beta0
      relbias <- ifelse(truth==0, 
                        ((estimate+1)-(truth+1))/(truth+1), 
                        (estimate-truth)/truth)
      coverage <- ifelse(lci>=truth | 
                           uci<=truth, 0, 1)
      df[ind, 1:7] <- c(estimate, truth, relbias, coverage, lci, uci, abund.scenarios[j])
      ind <- ind + 1
    }
  }}
# need to remove mods that didn't converge!
df2 <- df[war!=TRUE,]

#* test different lengths of time
set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lci=NA, uci=NA, time=NA)
war <- rep(FALSE, 3000)
ind <- 1
time.scenarios <- c(5, 10, 15)
for (j in 1:length(time.scenarios)){
  for (i in 1:1000){
    dat <- sim.ts(ntime=time.scenarios[j], plot.ts=FALSE)
    dfm <- data.frame(y=dat$datl$y, yr=dat$constl$yr)
    fit <- tryCatch(fit <- zeroinfl(formula = y ~ 1 + yr |
                                      1,
                                    dist    = "poisson",
                                    link    = "log",
                                    data    = dfm),
                    warning= function(w) {  return(TRUE) },
                    error= function(w) {  return(TRUE) }) 
    if (is.logical(fit)){ war[ind] <- fit
    ind <- ind+1 
    next} else{      
      out <- summary(fit)
      uci <- out$coefficients$count[2,1] + 1.96*out$coefficients$count[2,2]
      lci <- out$coefficients$count[2,1] - 1.96*out$coefficients$count[2,2]
      estimate <- out$coefficients$count[2,1]
      truth <- dat$truth$beta0
      relbias <- ifelse(truth==0, 
                        ((estimate+1)-(truth+1))/(truth+1), 
                        (estimate-truth)/truth)
      coverage <- ifelse(lci>=truth | 
                           uci<=truth, 0, 1)
      df[ind, 1:7] <- c(estimate, truth, relbias, coverage, lci, uci, time.scenarios[j])
      ind <- ind + 1
    }
  }}
# need to remove mods that didn't converge!
df3 <- df[war!=TRUE,]

# Plots
par(mfrow=c(1,3))
boxplot(df1$bias~df1$truth, 
        xlab="Trend (delta_1)",
        ylab="Relative bias of trend")
abline(h=0, lty=2, lwd=2)
mtext("A", cex=1.5, side=3, outer=F)

boxplot(df2$bias~df2$abund, 
        xlab="Average abundance exp(delta_0)",
        ylab="Relative bias of trend")
abline(h=0, lty=2, lwd=2)
mtext("B", cex=1.5, side=3, outer=F)

boxplot(df3$bias~df3$time, 
        xlab="Duration of time series",
        ylab="Relative bias of trend")
abline(h=0, lty=2, lwd=2)
mtext("C", cex=1.5, side=3, outer=F)


# tables
tapply(df1$bias, df1$truth, mean)
tapply(df1$bias, df1$truth, sd)
tapply(df1$cov, df1$truth, mean)

tapply(df2$bias, df2$abund, mean)
tapply(df2$bias, df2$abund, sd)
tapply(df2$cov, df2$abund, mean)

tapply(df3$bias, df3$time, mean)
tapply(df3$bias, df3$time, sd)
tapply(df3$cov, df3$time, mean)




















mod2 <- nimbleCode(
  {
    # priors
      mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
      beta ~ dnorm(0,sd=5)
      r ~ dexp(0.2) 
    
    # likelihood
    for (t in 1:ntime){
        p[t] <- r/(r+lambda[t])
        y[t] ~ dnegbin(p[t], r)
        # abundance
        log(lambda[t]) <- mu + beta*yr[t]
    } # t
    
  }
) # end model

params <- c("r", "mu", "beta")

set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lhdi=NA, uhdi=NA)
ind <- 1
beta.scenarios <- c(-0.7, -0.35, 0, 0.35, 0.7)
for (j in 1:length(beta.scenarios)){
for (i in 1:100){
dat <- sim.ts(beta0=beta.scenarios[j])
ntime <- dat[[2]]$ntime
inits <- function(){ list (mu = log(runif(1,0,100)),
                           r = runif(1,0,10),
                           beta = runif(1, -5, 5),
                           lambda = runif(ntime, 0, 20))}

out <- tryCatch(post <- nimbleMCMC(
  constants = dat[["constl"]],
  data = dat[["datl"]],
  inits=inits,
  code = mod2,
  niter = 5000,
  nburnin = 2500,
  nchains = 3
) , error = function(e){ e }) 
if (any(class(out) == "error")){ 
  ind <- ind + 1
  out <- c()
  next }
else{
ps <- MCMCsummary(post, HPD=TRUE)
#MCMCtrace(post, pdf=F)
estimate <- ps["beta","mean"]
lhdi <- ps["beta","95%_HPDL"] 
uhdi <- ps["beta","95%_HPDU"]
truth <- dat$truth$beta0
#try(relbias <- (estimate-truth)/truth)
coverage <- ifelse(ps$`95%_HPDL`[1]>=truth & 
                  ps$`95%_HPDU`[1]<=truth, 0, 1)
df[ind, 1:6] <- c(estimate, truth, NA, coverage, lhdi, uhdi)
ind <- ind + 1
} # else
}}  

write.csv(df, "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\simulation-tab-trend.csv")                   
boxplot(df$bias, ylab="Relative bias of trend estimate")
abline(h=0, lty=2, lwd=2)
mean(df$cov)

#******************
#* Test whether small abundances 
#* influence the detection of trends
#******************
set.seed(123456)
df <- data.frame(est=NA, truth=NA, bias=NA, cov=NA, lhdi=NA, uhdi=NA, abund=NA)
ind <- 1
abund.scenarios <- c(0.5, 1, 2, 3)
for (j in 1:length(abund.scenarios)){
  for (i in 1:100){
    dat <- sim.ts(mn.abund=abund.scenarios[j])
    ntime <- dat[[2]]$ntime
    inits <- function(){ list (mu = log(runif(1,0,100)),
                               r = runif(1,0,10),
                               beta = runif(1, -5, 5),
                               lambda = runif(ntime, 0, 20))}
    
    out <- tryCatch(post <- nimbleMCMC(
      constants = dat[["constl"]],
      data = dat[["datl"]],
      inits=inits,
      code = mod2,
      niter = 5000,
      nburnin = 2500,
      nchains = 3
    ) , error = function(e){ e }) 
    if (any(class(out) == "error")){ 
      ind <- ind + 1
      out <- c()
      next }
    else{
      ps <- MCMCsummary(post, HPD=TRUE)
      #MCMCtrace(post, pdf=F)
      estimate <- ps["beta","mean"]
      lhdi <- ps["beta","95%_HPDL"] 
      uhdi <- ps["beta","95%_HPDU"]
      truth <- dat$truth$beta0
      relbias <- (estimate-truth)/truth
      coverage <- ifelse(ps$`95%_HPDL`[1]>=truth & 
                           ps$`95%_HPDU`[1]<=truth, 0, 1)
      df[ind, 1:7] <- c(estimate, truth, relbias, coverage, lhdi, uhdi, abund.scenarios[j])
      ind <- ind + 1
    } # else
  }}  

write.csv(df, "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\simulation-tab-abundance.csv")                   
boxplot(df$bias, ylab="Relative bias of trend estimate")
abline(h=0, lty=2, lwd=2)
mean(df$cov)


#***************
# scraps
#**************
#*#################
#* set up for model
#* run using nimble
#**************

# specify a nimble function for 
# ZIP distribution. Speeds runtimes.
dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

# specify model
mod <- nimbleCode({
  mn.psi <- logit(psi1)
  psi1 ~ dunif(0, 1)
  lmn.abund ~ dnorm(0, sd=5)
  #period2 ~ dunif(0, period*5)
  #sd.period ~ T(dnorm(0, sd=2),0, )
  beta0 ~ dunif(-5, 5)
  beta1 ~ dunif(-5, 5)
  beta2 ~ dunif(-5, 5)
  
  for (t in 1:length(ntime)){
    #period ~ dnorm(mu.period, sd=sd.period)
    log(lambda[t]) <- lmn.abund + beta0*yr[t] 
    logit(psi[t]) <- mn.psi +
      beta1 * sin(2*3.141593*(yr[t]/period)) +
      beta2 * cos(2*3.141593*(yr[t]/period))
    y[t] ~ dpois(mu[t])
    mu[t] <- lambda[t]*z[t]
    z[t] ~ dbern(psi[t])
    #y[t] <- dZIP(lambda=lambda[t], zeroProb = psi[t])
  } # t
})

params <- c("mn.psi", "psi1", "mn.abund",
            "period", "beta0", "beta1", "beta2")

dat <- sim.ts()
ntime <- dat[[2]]$ntime
inits <- function(){ list (mn.psi = runif(1,0,4),
                           psi1 = runif(1, 0.5, 1),
                           lmn.abund = log(runif(1,0,100)),
                           period = runif(1, 1, ntime),
                           #sd.period = runif(1, 0, 2),
                           beta0 = runif(1, -5, 5),
                           beta1 = runif(1, -5, 5),
                           beta2 = runif(1, -5, 5),
                           psi = runif(ntime, 0, 1),
                           lambda = runif(ntime, 0, 20),
                           z = rbinom(ntime, 1, prob=0.5)
                           
)}

post <- nimbleMCMC(
  constants = dat[["constl"]],
  data = dat[["datl"]],
  inits=inits,
  code = mod,
  niter = 10000,
  nburnin = 5000,
  nchains = 1
) 

