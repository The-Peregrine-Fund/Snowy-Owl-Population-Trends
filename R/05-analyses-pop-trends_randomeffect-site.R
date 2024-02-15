## ---- start --------
library(readxl)
library(reshape2)
library(MCMCvis)
library(HDInterval)
library(ggplot2)
library(tidybayes)
library(bayestestR)
library(cowplot)
dat <- read.csv("data\\data.csv", na.strings="NA")
# change Bylot Island into 
# 2 separate sites
# one of 100 km^2 and 
# one of 300 km^2
dat$BylotIsland_300 <- dat$BylotIsland_400-dat$BylotIsland_100
colSums(!is.na(dat))

# Bundle data
y <- as.data.frame(dat[-c(1,2), -c(4)]) 

# create data lists for nimble
y.nim <- y[,-1]
datl <- list(y = y.nim) # exclude 8-fennoscania because no area information
constl <- list(nsite = ncol(y.nim),
               ntime = nrow(y.nim),
               time = (y$year-2004)/16
)

## ---- runmod --------
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
      mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
      beta ~ dnorm(0, sd=5)
      for (j in 1:nsite){ 
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
          log(lambda[t,j]) <- mu + 
            beta*time[t] + 
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
                             mu = runif(1, -2, 2),
                             beta = runif(1, -2, 2),
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

# approximately 15 minute runtime
this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)
stopCluster(this_cluster)

nb <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

# save(out=nb, post=post, run, 
#      file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_site-re_allsites.Rdata")

## ---- lambda --------
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_site-re_allsites.Rdata")
# Check for convergence
params_nb <-c( "sigma.time", "sigma.site",
               "mu", "beta", "r")
MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

par(mfrow=c(1,1))
MCMCplot(object = nb, params = params_nb)
coeftab95 <- MCMCsummary(nb, params_nb, 
                       HPD=TRUE,  hpd_prob = 0.95,
                       round=2, pg0 = TRUE, func = median, 
                       func_name = "Median")
coeftab85 <- MCMCsummary(nb, params_nb, 
                       HPD=TRUE,  hpd_prob = 0.85,
                       round=2, pg0 = TRUE, func = median, 
                       func_name = "Median")

mu <- MCMCpstr(nb, "mu", type="chains")[[1]]
beta <- MCMCpstr(nb, "beta", type="chains")[[1]]
pred.lam <- array(NA, dim=c(length(constl$time), ncol(mu)))
for (t in 1:length(constl$time)){
  pred.lam[t,] <- mu + 
    beta[1,]*constl$time[t] 
}

lam.post <- array (NA, dim(pred.lam), dimnames=list(Year=1988:2020, Iter=1:4000))
# calculate lambda (population growth rate)
for (t in 2:nrow(pred.lam)){
    lam.post[t,] <- exp(pred.lam[t,]) / exp(pred.lam[t-1,])
  }
lam.post.long <- melt(lam.post)

# calculate overall average pop growth rate
av.post <- apply(lam.post, 2, mean, na.rm=T)
median(av.post, na.rm=T)
HDInterval::hdi(av.post, credMass = 0.95, na.rm=T)
HDInterval::hdi(av.post, credMass = 0.80, na.rm=T)
pd(av.post, null=1, na.rm=T)

# calculate pop growth rates for each year
lam.df <- data.frame(
  year= 1988:2020,
  m= apply(lam.post, 1, median, na.rm=T),
  lci= apply(lam.post, 1, HDInterval::hdi)[1,],
  uci= apply(lam.post, 1, HDInterval::hdi)[2,],
  lci80 = apply(lam.post, 1, HDInterval::hdi, credMass=.8)[1,],
  uci80 = apply(lam.post, 1, HDInterval::hdi, credMass=.8)[2,],
  pd = apply(lam.post, 1, pd, null=1)
)

makeTransparent<-function(someColor, alpha=100){ 
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
c3 <- makeTransparent("gray40", alpha=10)
c4 <- makeTransparent(4, alpha=180)

rownames(datl$y) <- 1988:2020
names(dimnames(datl$y)) <- c("Year", "Site")
all.df <- melt(as.matrix(datl$y), value.name="Count") 
colnames(all.df) <- c("Year", "Site", "Count")
labels <- unique(all.df$Site)

sl <- all.df[!is.na(all.df$Count),]

p1 <- ggplot()  + theme_minimal() +
  coord_cartesian(clip = "off") +
  geom_line(data=all.df,aes(x=Year,y=Site),linewidth=2, color='gray60')+
  geom_point(data=sl,aes(x=Year,y=Site),shape="I",size=4)+
  xlab("Year")+ylab("Site")+scale_color_viridis_c()+
  scale_y_discrete(labels=labels)+
  labs(color='black') +
  xlim(1988, 2020) + #theme(legend.position = "top") + 
  annotate(geom = "text", x = 1990, y = 9.5, label = "A", size=6)

p2 <- ggplot() + theme_minimal() + 
  geom_line(data=lam.post.long, aes(x=Year, y=value, group=Iter), 
            color="gray40", linewidth=0.5, alpha=0.05 ) +
  geom_hline(yintercept=1, lwd=2, color="black", linetype="dashed") +
  geom_line(data=lam.df, aes(x=year, y=lci), color="deepskyblue3", linewidth=0.5 ) +
  geom_line(data=lam.df, aes(x=year, y=uci ), color="deepskyblue3", linewidth=0.5) +
  geom_line(data=lam.df, aes(x=year, y=lci80), color="deepskyblue3", linewidth=1 ) +
  geom_line(data=lam.df, aes(x=year, y=uci80), color="deepskyblue3", linewidth=1 ) +
  geom_line(data=lam.df, aes(x=year, y=m), color="deepskyblue3", linewidth=2 ) +
  xlab("Year") +
  ylab( expression(paste("Population growth (", lambda,")")) )+
  xlim(1988, 2020) + 
  annotate(geom = "text", x = 1990, y = 1.1, label = "B", size=6)+
  coord_cartesian(xlim=c(1988, 2020), ylim=c(0.85, 1.15))

ap12 <- align_plots(p1, p2, align="v", axis="l")
p12 <- plot_grid(ap12[[1]], ap12[[2]], nrow = 2, align="v")
## ---- lambdaplot2 --------
p12

# ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\pop-growth-rate_year-siteRE.tiff",
#        plot=p12,
#        device="tiff",
#        width=6,
#        height=8,
#        dpi=300)

## ---- perc --------
# plot percent change
perc.post <- array(NA, dim=dim(pred.lam), dimnames=list(Year=1988:2020, Iter=1:4000))
# calc lambda for each site relative to 1996
for (t in 1:nrow(pred.lam)){
  perc.post[t,] <- exp(pred.lam[t,]) / exp(pred.lam[9,])
}
perc.post.long <- melt(perc.post)
iucn <- (perc.post-1) * 100
iucn.post <- melt(iucn)

i.df <- data.frame(
  year = 1988:2020,
  m = apply(iucn, 1, median, na.rm=T),
  lci95 = apply(iucn, 1, HDInterval::hdi)[1,],
  uci95 = apply(iucn, 1, HDInterval::hdi)[2,],
  lci80 = apply(iucn, 1, HDInterval::hdi, credMass=0.8)[1,],
  uci80 = apply(iucn, 1, HDInterval::hdi, credMass=0.8)[2,],
  pd = apply(iucn, 1, pd)
)

knitr::kable(i.df, digits=c(0,1,1,1,1,1,2),
             row.names=FALSE,
             col.names = c("Year", "Median", "95% Lower HDI", "95% Upper HDI", 
                           "85% Lower HDI", "85% Upper HDI", "Prob. direction"),
             caption="Table S5. Percent change since 1996.")

p3 <- ggplot() + theme_minimal() + 
  geom_rect(aes(xmin=1986, xmax=2022, ymin=-20, ymax=300), color="green4", fill="green4") +
  geom_rect(aes(xmin=1986, xmax=2022, ymin=-30, ymax=-20), color="green3", fill="green3") +
  geom_rect(aes(xmin=1986, xmax=2022, ymin=-50, ymax=-30), color="yellow", fill="yellow") +
  geom_rect(aes(xmin=1986, xmax=2022, ymin=-80, ymax=-50), color="orange", fill="orange") +
  geom_rect(aes(xmin=1986, xmax=2022, ymin=-100, ymax=-80), color="red", fill="red") +
  geom_line(data=iucn.post, aes(x=Year, y=value, group=Iter), 
            color="gray40", size=0.5, alpha=0.05 ) +
  geom_hline(yintercept=1, lwd=2, color="black", linetype="dashed") +
  geom_line(data=i.df, aes(x=year, y=lci95), color="black", size=0.5 ) +
  geom_line(data=i.df, aes(x=year, y=uci95 ), color="black", size=0.5) +
  geom_line(data=i.df, aes(x=year, y=lci80), color="black", size=1 ) +
  geom_line(data=i.df, aes(x=year, y=uci80), color="black", size=1 ) +
  geom_line(data=i.df, aes(x=year, y=m), color="black", size=2 ) +
  annotate(geom = "text", x = 1995, y = 25, label = "Least concern", size=4, color="black") +
  annotate(geom = "text", x = 1995, y = -25, label = "Near threatened", size=4, color="black") +
  annotate(geom = "text", x = 1995, y = -40, label = "Vulnerable", size=4, color="black") +
  annotate(geom = "text", x = 1995, y = -70, label = "Endangered", size=4, color="black") +
  annotate(geom = "text", x = 1995, y = -90, label = "Critically endangered", size=4, color="black") +
  xlab("Year") +
  ylab( "Percent change")+
  coord_cartesian(xlim=c(1988, 2020), ylim=c(-100, 60)) +
  annotate(geom = "text", x = 1995, y = 55, label = "A", size=8)

iucn2020 <- iucn.post[iucn.post$Year==2020, ] 
i.df2020 <- i.df[i.df$year==2020,]

p4 <- ggplot() + theme_minimal() +
  geom_rect(aes(ymin=0, ymax=2, xmin=-20, xmax=300), color="green4", fill="green4") +
  geom_rect(aes(ymin=0, ymax=2, xmin=-30, xmax=-20), color="green3", fill="green3") +
  geom_rect(aes(ymin=0, ymax=2, xmin=-50, xmax=-30), color="yellow", fill="yellow") +
  geom_rect(aes(ymin=0, ymax=2, xmin=-80, xmax=-50), color="orange", fill="orange") +
  geom_rect(aes(ymin=0, ymax=2, xmin=-100, xmax=-80), color="red", fill="red") +
  stat_halfeye(data=iucn2020, aes(x = value, y=0), alpha=0.65, 
               slab_color="gray20", slab_fill="gray40",
               point_interval="median_hdi", .width = c(0.80, 0.95),
               point_size=5) +
  scale_size_continuous(range = c(7, 15)) +
  geom_vline(xintercept=0, lwd=2, color="black", linetype="dashed") +
  xlab("") + ylab("Density (scaled)\nof percent change\nover three generations") +
  annotate(geom = "text", x = 55, y = 0.25, label = "B", size=8) +
  coord_flip(xlim=c(-100, 60), ylim=c(0, 1), clip="on")

ap45 <- align_plots(p3, p4, align="h", axis="l")
p45 <- plot_grid(ap45[[1]], ap45[[2]], nrow = 1, align="h", rel_widths = c(2, 1))
## ---- perchangeplot2 --------
p45

# ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\percentchange_year-resites.tiff",
#        plot=p45,
#        device="tiff",
#        width=8,
#        height=4,
#        dpi=300)
