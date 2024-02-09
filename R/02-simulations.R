## ---- start --------
library(nimble)
library(MCMCvis)
library (MASS)
library (pscl)

## ---- simfunc --------
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

## ---- sims --------
#*******************
#* Frequentist analysis
#* using a
#* negative binomial distribution
#******************* 
set.seed(123456)
# Population trends vary
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

#* Abundance varies
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
# need to remove mods that didn't converge.
df2 <- df[war!=TRUE,]

#* Lengths of time vary
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
# need to remove mods that didn't converge
df3 <- df[war!=TRUE,]

## ---- plots --------
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

## ---- tables --------
# tables
tab1 <- data.frame(
                  # different population trends
                  varied.by= "population trends",
                  mn.abund=10,
                  trend= c(-0.7, -0.35, 0, 0.35, 0.7),
                  duration=30,
                  mean.rel.bias=tapply(df1$bias, df1$truth, mean), 
                  sd.rel.bias=tapply(df1$bias, df1$truth, sd),
                  coverage=tapply(df1$cov, df1$truth, mean)
)
tab2 <- data.frame(
                  #* test different levels of abundance
                  varied.by= "abundance",
                  mn.abund=c(0.5, 1, 2, 3),
                  trend=-0.7,
                  duration=30,
                  mean.rel.bias=tapply(df2$bias, df2$abund, mean),
                  sd.rel.bias=tapply(df2$bias, df2$abund, sd),
                  coverage=tapply(df2$cov, df2$abund, mean)
)
tab3 <- data.frame(
                  #* test different lengths of time
                  varied.by= "duration",
                  mn.abund=10,
                  trend=-0.7,
                  duration=c(5, 10, 15),
                  mean.rel.bias=tapply(df3$bias, df3$time, mean),
                  sd.rel.bias=tapply(df3$bias, df3$time, sd),
                  coverage=tapply(df3$cov, df3$time, mean)
)
tab <- rbind(tab1, tab2, tab3)
rownames(tab) <- NULL
knitr::kable(tab, digits=c(0, 2, 2, 2, 2),
             row.names=FALSE,
             caption="Table S1. Simulation scenarios and results.")
