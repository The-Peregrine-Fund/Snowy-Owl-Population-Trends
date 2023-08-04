library (MCMCvis)
library (HDInterval)
library(coda)
load("data/data.Rdata")
# No Wrangel results
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_noWrangel.Rdata")
# Results with Wrangel
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")
convert <- function(x){
  list(as.mcmc(x$samples$chain1), 
       as.mcmc(x$samples$chain2), 
       as.mcmc(x$samples$chain3))
}

# function to plot posterior predictive checks
plot.diag <- function(out, ratio=FALSE, lab="",
                      low=low, high=high){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  
  if(length(low)==1 & length(high)==1){ mn <- low; mx <- high }
  else{
    mx <- max(c(samps$dmape.rep[1,], samps$dmape.obs[1,]))
    mn <- min(c(samps$dmape.rep[1,], samps$dmape.obs[1,]))
  }
  plot(jitter(samps$dmape.obs[1,], amount=300), 
       jitter(samps$dmape.rep[1,], amount=300),
       main=paste0("Mean absolute percentage error\nmodel\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(samps$dmape.rep[1,] > samps$dmape.obs[1,]),2)
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


op <- convert(out_p)
ozip <- convert(out_zip)
ozip2 <- convert(out_zip2)
onb <- convert(out_nb)

plot.diag(onb, low=0, high=5000000)

# calculate "data" from offsets
df <- data.frame(array(NA, dim(datl$y)))
for (j in 1:constl$nsite){
  df[,j] <- datl$y[,j]/(constl$area[j]/100)
} 

MCMCsummary(onb, round=3,
            c("sigma", 
              "mu", "beta"))
all <- MCMCpstr(onb, type="chains")
mb <- apply(all$beta, 1, median)
hdi95b <- apply(all$beta, 1, HDInterval::hdi)
hdi80b <- apply(all$beta, 1, HDInterval::hdi, credMass=0.8)

MCMCplot(object = onb, params = 'lambda')

#***************
#* Plot trends
#**************
lam.md <- apply(all$lambda, c(1,2), median)
lam.67hdis <- apply(all$lambda, c(1,2), HDInterval::hdi, credMass=0.67)
lam.95hdis <- apply(all$lambda, c(1,2), HDInterval::hdi, credMass=0.95)
# add NAs for unsurveyed years
nas <- is.na(datl$y)
lam.md[nas] <- lam.67hdis[nas] <- lam.95hdis[nas] <- NA
ymx <-  max(lam.95hdis, na.rm=T)
  
plot(NA, type="n", 
     xlim=c(1988, 2020), ylim=c(0,25))
for (j in 1:constl$nsite){
  text(1988:2020, lam.md[,j], as.character(1:6)[j])
}
