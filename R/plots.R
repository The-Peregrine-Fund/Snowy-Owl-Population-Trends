library ("zoo")
library ("HDInterval")
library ("MCMCvis")
library ("ggplot2")
library ("reshape2")
load("data/data.Rdata")
# No Wrangel results
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_noWrangel.Rdata")
# Results with Wrangel
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")
# helper function for transparency
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
# check for autocorrelation
par(mfrow=c(3,3))
acf(dat$Utqiagvik_213[-c(1:6)])
acf(dat$BylotIsland_100[-c(1:7)])
acf(dat$BylotIsland_300[-c(1:14,35)])
acf(dat$'Igloolik Island_114'[-c(1:25)])
acf(dat$'Karupelv Valley_75'[-c(1:2)])
acf(dat$'Hochstetter Forland_100'[!is.na(dat$`Hochstetter Forland_100`)] )
acf(dat$Fennoscandia_x[!is.na(dat$Fennoscandia_x)] )
acf(dat$Wrangel[!is.na(dat$Wrangel)] )
acf(dat$Alaska[!is.na(dat$Alaska)] )

#******************
#* Raw data plots
#******************
# calculate moving averages
mn <- array(NA, dim=dim(datl$y))
for (i in 1:ncol(datl$y)){
  mn[,i] <- rollmean(datl$y[,i], k=3, fill=NA)
}

par(mfrow=c(2,1))
# raw data
plot(1988:2020, rep(NA, length(1988:2020)), 
     xlim=c(1988,2020), 
     ylim=c(0, max(datl$y, na.rm=T)),
     xlab="Year", 
     ylab="Count",
     main="Raw data")
for (i in 1:ncol(datl$y)){
  lines(1988:2020, datl$y[,i], lty=i, lwd=2, col="red")
}

# 2-year moving average
plot(1988:2020, rep(NA, length(1988:2020)), 
     xlim=c(1988,2020), 
     ylim=c(0, max(mn, na.rm=T)),
     xlab="Year", 
     ylab="",
     main="3-year moving average")
for (i in 1:ncol(mn)){
  lines(1988:2020, mn[,i], lty=i, lwd=2)
}

#******************
#* Model plots
#******************
#* Relative abundance of each site
library (viridis)
beta <- MCMCpstr(nb, "beta", type="chains")[[1]]
mu <- MCMCpstr(nb, "mu", type="chains")[[1]]
pred.lam <- array(NA, dim=c(length(constl$time), nrow(mu),  ncol(mu)))

for (i in 1:nrow(mu)){
  for (t in 1:length(constl$time)){
    pred.lam[t,i,] <- mu[i,] + beta[i,]*constl$time[t]
  }
}
# set NAs for unsurveyed years
wna <- which(is.na(datl$y), arr.ind=TRUE)
for (i in 1:nrow(wna)){
  pred.lam[wna[i,1],wna[i,2], ] <- NA
}

beta.tab <- data.frame(
  median = apply(beta, 1, median),
  LHDI95 = apply(beta, 1, hdi, credMass=0.95)[1,],
  UHDI95 = apply(beta, 1, hdi, credMass=0.95)[2,],
  LHDI80 = apply(beta, 1, hdi, credMass=0.80)[1,],
  UHDI80 = apply(beta, 1, hdi, credMass=0.80)[2,]
)
nms <- c("Utqiagvik", "Bylot Island Core", 
         "Karupelv Valley", 
         "Fennoscandia", "Wrangel")
beta.tab$sig95 <- ifelse(beta.tab$LHDI95>0 | beta.tab$UHDI95<0, TRUE, FALSE )
beta.tab$sig80 <- ifelse(beta.tab$LHDI80>0 | beta.tab$UHDI80<0, TRUE, FALSE )
beta.tab$site <- nms
beta.tab <- beta.tab[, c(8, 1:7)]
print(beta.tab, digits=3)

lam.md <- apply(pred.lam, c(1,2), median, na.rm=F) |> exp()
lam.lci95 <- apply(pred.lam, c(1,2), hdi, na.rm=F)[1,,] |> exp()
lam.uci95 <- apply(pred.lam, c(1,2), hdi, na.rm=F)[2,,] |> exp()
lam.lci80 <- apply(pred.lam, c(1,2), hdi, na.rm=F, credMass=0.80)[1,,] |> exp()
lam.uci80 <- apply(pred.lam, c(1,2), hdi, na.rm=F, credMass=0.80)[2,,] |> exp()
lam.lci67 <- apply(pred.lam, c(1,2), hdi, na.rm=F, credMass=0.67)[1,,] |> exp()
lam.uci67 <- apply(pred.lam, c(1,2), hdi, na.rm=F, credMass=0.67)[2,,] |> exp()

cols <- viridis(nrow(mu))
cols2 <- viridis(nrow(mu))
c3 <- makeTransparent("gray40", alpha=10)

tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\trends_noOffset.tiff",
     height=4, width=6, units="in", res=300)
par(mfrow=c(1,1), 
    mar=c(0,2,0,0), 
    oma=c(5,5,1,1))
for (i in c(3)){
  plot(NA, 
       xlim=c(2010,2020), 
       ylim=c(0,  
              max(lam.uci95[,i], na.rm=T)+0.25),
       type="n",
       xlab="", ylab="", 
       xaxt="n", yaxt="n")
  # if (i==1){
  #   axis(1, at=c(1990, 2000, 2010, 2020),
  #        labels=c(NA, NA, NA, NA))
  #   axis(2, at=c(0, 10, 20, 30))
  # }
  if (i==3){
    axis(1, at=c(2010, 2015, 2020))
    axis(2, at=c(0, 1, 2, 3))
  }
  notna <- !is.na(lam.md[,i])
  l80 <- lam.lci80[ notna, i ]
  u80 <- lam.uci80[ notna, i ]
  l95 <- lam.lci95[ notna, i ]
  u95 <- lam.uci95[ notna, i ]
  years <- c(1988:2020)[ notna ]
  
  cis95 <- c(l95, rev(u95))
  cis80 <- c(l80, rev(u80))
  yrs <- c(years, rev(years))
  for (j in 1:dim(pred.lam)[[3]]){
    lines(1988:2020, exp(pred.lam[,i,j]), col = c3)
  }
  # polygon(yrs, cis95, col="gray80", border=NA)
  # polygon(yrs, cis80, col="gray60", border=NA)
  lines(1988:2020, lam.md[,i], 
        col="black", #cols[i], 
        lwd=3)
  
  mtext(side=1, "Year", outer=TRUE, 
        adj=0.55, padj=4, cex=1.3)
  #expr <- expression(Detections~per~100~km^2)
  mtext(side=2, "Relative abundance", outer=TRUE, 
        padj=-1, cex=1.3)
  title(nms[i], line = -3, cex=1.3)
}
dev.off()

#***********************
#* Plot an overall trend
#* weighted by abundance
#*********************** 
# weight by modeled abundance for each year
# weight by data, mean abundance across years
# INCLUDES WRANGEL
library(gridBase)
library(grid)
library (gridExtra)
library (cowplot)
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")
beta <- MCMCpstr(nb, "beta", type="chains")[[1]]
mu <- MCMCpstr(nb, "mu", type="chains")[[1]]
pred.lam <- array(NA, dim=c(length(constl$time), nrow(mu),  ncol(mu)))

for (i in 1:nrow(mu)){
  for (t in 1:length(constl$time)){
    pred.lam[t,i,] <- mu[i,] + beta[i,]*constl$time[t]
  }
}
# set NAs for unsurveyed years
wna <- which(is.na(datl$y), arr.ind=TRUE)
for (i in 1:nrow(wna)){
  pred.lam[wna[i,1],wna[i,2], ] <- NA
}

#wmns <- colMeans(datl$y, na.rm=T)
wmns <- as.matrix(datl$y+1)
wmns <- ifelse(is.na(wmns), 0, wmns)
zo <- ifelse(is.na(datl$y), 0, 1)
wts2 <- array(NA, dim(zo)) 

for (t in 1:nrow(zo)){
  mns <- wmns[t,]*zo[t,] # impute zeroes for unmonitored years
  wts2[t,1:5] <- unlist(mns/sum(mns, na.rm=T)) # calculate proportions
}

# calculate lambda
r.post <- lam.post <- array (NA, dim=dim(pred.lam))
for (t in 2:nrow(pred.lam)){
  for (j in 1:ncol(pred.lam)){
    lam.post[t,j,] <- exp(pred.lam[t,j,]) / exp(pred.lam[t-1,j,])
    r.post[t,j,] <- pred.lam[t,j,] - pred.lam[t-1,j,]
  }}

wm.post.r <- wm.post.lam <- array (NA, c(nrow(pred.lam),dim(pred.lam)[[3]]))
for (t in 1:nrow(pred.lam)){
  for (k in 1:dim(pred.lam)[[3]]){
    wm.post.lam[t,k] <- weighted.mean(x=lam.post[t,,k], w=wts2[t,], na.rm=T)
    wm.post.r[t,k] <- weighted.mean(x=r.post[t,,k], w=wts2[t,], na.rm=T)
  }}
dimnames(wm.post.lam) <- list(Year=1988:2020, Iter=1:4000)
lam.post <- melt(wm.post.lam)

wm.m <- apply(wm.post.r, 1, median, na.rm=T)
wm.lci <- apply(wm.post.r, 1, hdi)[1,]
wm.uci <- apply(wm.post.r, 1, hdi)[2,]
wm.lci80 <- apply(wm.post.r, 1, hdi, credMass=.8)[1,]
wm.uci80 <- apply(wm.post.r, 1, hdi, credMass=.8)[2,]

wm.df <- data.frame(
  year= 1988:2020,
  m= apply(wm.post.lam, 1, median, na.rm=T),
  lci= apply(wm.post.lam, 1, hdi)[1,],
  uci= apply(wm.post.lam, 1, hdi)[2,],
  lci80 = apply(wm.post.lam, 1, hdi, credMass=.8)[1,],
  uci80 = apply(wm.post.lam, 1, hdi, credMass=.8)[2,]
)

wm2.m <- apply(wm.post.lam, 1, median, na.rm=T)
wm2.lci <- apply(wm.post.lam, 1, hdi)[1,]
wm2.uci <- apply(wm.post.lam, 1, hdi)[2,]
wm2.lci80 <- apply(wm.post.lam, 1, hdi, credMass=.8)[1,]
wm2.uci80 <- apply(wm.post.lam, 1, hdi, credMass=.8)[2,]

c3 <- makeTransparent("gray40", alpha=10)
c4 <- makeTransparent(4, alpha=180)

rownames(wts2) <- rownames(datl$y) <- 1988:2020
colnames(wts2) <- colnames(datl$y)
names(dimnames(datl$y)) <- c("Year", "Site")
all.df <- melt(as.matrix(datl$y), value.name="Count") 
all.df$weights <- melt(as.matrix(wts2), value.name="weights")$weights
colnames(all.df) <- c("Year", "Site", "Count", "Weights")
labels <- c("Utqiagvik", "Bylot Island", "Karupelv Valley", 
            "Fennoscandia", "Wrangel Island")

sl <- all.df[!is.na(all.df$Count),]

p1 <- ggplot()  + theme_minimal() +
  coord_cartesian(clip = "off") +
  geom_line(data=all.df,aes(x=Year,y=Site,color=Weights),size=2)+
  geom_point(data=sl,aes(x=Year,y=Site),shape="I",size=4)+
  xlab("Year")+ylab("Site")+scale_color_viridis_c()+
  scale_y_discrete(labels=labels)+
  labs(color='Weight') +
  xlim(1988, 2020) + theme(legend.position = "top") + 
  annotate(geom = "text", x = 1990, y = 6, label = "A", size=6)


p2 <- ggplot() + theme_minimal() + 
      geom_line(data=lam.post, aes(x=Year, y=value, group=Iter), 
                color="gray40", size=0.5, alpha=0.05 ) +
      geom_hline(yintercept=1, lwd=2, color="black", linetype="dashed") +
      geom_line(data=wm.df, aes(x=year, y=lci), color="deepskyblue3", size=0.5 ) +
      geom_line(data=wm.df, aes(x=year, y=uci ), color="deepskyblue3", size=0.5) +
      geom_line(data=wm.df, aes(x=year, y=lci80), color="deepskyblue3", size=1 ) +
      geom_line(data=wm.df, aes(x=year, y=uci80), color="deepskyblue3", size=1 ) +
      geom_line(data=wm.df, aes(x=year, y=m), color="deepskyblue3", size=2 ) +
      xlab("Year") +
      ylab( expression(paste("Population growth (", lambda,")")) )+
      #ylab( expression("Population growth (",lambda,")") )+
      xlim(1988, 2020) + 
      annotate(geom = "text", x = 1990, y = 1.2, label = "B", size=6)

align_plots(p1, p2, align="v", axis="l")
both_aligned <- plot_grid(p1, p2, nrow = 2, align="v")
ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\pop-growth-rate_year.tiff", 
       plot=both_aligned, 
       device="tiff",
       width=6, 
       height=8,
       dpi=300)

# tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\pop-growth-rate_year.tiff",
#      height=9, width=6, units="in", res=300)
# plot(1988:2020, wm2.m,
#      xlab="Year", cex.lab=1.3,
#      ylab="Population growth rate (lambda)",
#      ylim=c(0.8, 1.2),
#      type="n", lwd=4,
#      xaxt="n", yaxt="n")
# axis(1, at=c(1990,2000,2010, 2020) )
# axis(2, at=c(0.8, 0.9, 1, 1.1, 1.2))
# for (j in 1:ncol(wm.post.lam)){
#   lines(1988:2020, wm.post.lam[,j], col = c3)
# }
# lines(1988:2020, wm2.lci, col=c4, lwd=2)
# lines(1988:2020, wm2.uci, col=c4, lwd=2)
# lines(1988:2020, wm2.lci80, col=c4, lwd=2)
# lines(1988:2020, wm2.uci80, col=c4, lwd=2)
# abline(h=1, lty=2, lwd=4)
# lines(1988:2020, wm2.m, lwd=4, col=c4)
# dev.off()

############
# Drop Wrangel
############
# No Wrangel results
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb_noWrangel.Rdata")

beta <- MCMCpstr(nb, "beta", type="chains")[[1]]
mu <- MCMCpstr(nb, "mu", type="chains")[[1]]
pred.lam <- array(NA, dim=c(length(constl$time), nrow(mu),  ncol(mu)))

for (i in 1:nrow(mu)){
  for (t in 1:length(constl$time)){
    pred.lam[t,i,] <- mu[i,] + beta[i,]*constl$time[t]
  }
}
# set NAs for unsurveyed years
wna <- which(is.na(datl$y[,-5]), arr.ind=TRUE)
for (i in 1:nrow(wna)){
  pred.lam[wna[i,1],wna[i,2], ] <- NA
}

wmns <- as.matrix(datl$y[,-5]+1)
wmns <- ifelse(is.na(wmns), 0, wmns)
zo <- ifelse(is.na(datl$y[,-5]), 0, 1)
wts2 <- array(NA, dim(zo)) 

for (t in 1:nrow(zo)){
  mns <- wmns[t,]*zo[t,] # impute zeroes for unmonitored years
  wts2[t,1:4] <- unlist(mns/sum(mns, na.rm=T)) # calculate proportions
}

# calculate lambda
r.post <- lam.post <- array (NA, dim=dim(pred.lam))
for (t in 2:nrow(pred.lam)){
  for (j in 1:ncol(pred.lam)){
    lam.post[t,j,] <- exp(pred.lam[t,j,]) / exp(pred.lam[t-1,j,])
    r.post[t,j,] <- pred.lam[t,j,] - pred.lam[t-1,j,]
  }}

wm.post.r <- wm.post.lam <- array (NA, c(nrow(pred.lam),dim(pred.lam)[[3]]))
for (t in 1:nrow(pred.lam)){
  for (k in 1:dim(pred.lam)[[3]]){
    wm.post.lam[t,k] <- weighted.mean(x=lam.post[t,,k], w=wts2[t,], na.rm=T)
    wm.post.r[t,k] <- weighted.mean(x=r.post[t,,k], w=wts2[t,], na.rm=T)
  }}
dimnames(wm.post.lam) <- list(Year=1988:2020, Iter=1:4000)
lam.post <- melt(wm.post.lam)

wm.m <- apply(wm.post.r, 1, median, na.rm=T)
wm.lci <- apply(wm.post.r, 1, hdi)[1,]
wm.uci <- apply(wm.post.r, 1, hdi)[2,]
wm.lci80 <- apply(wm.post.r, 1, hdi, credMass=.8)[1,]
wm.uci80 <- apply(wm.post.r, 1, hdi, credMass=.8)[2,]

wm2.m <- apply(wm.post.lam, 1, median, na.rm=T)
wm2.lci <- apply(wm.post.lam, 1, hdi)[1,]
wm2.uci <- apply(wm.post.lam, 1, hdi)[2,]
wm2.lci80 <- apply(wm.post.lam, 1, hdi, credMass=.8)[1,]
wm2.uci80 <- apply(wm.post.lam, 1, hdi, credMass=.8)[2,]

c3 <- makeTransparent("gray40", alpha=10)
c4 <- makeTransparent(4, alpha=180)

rownames(wts2) <- rownames(datl$y[,-5]) <- 1988:2020
colnames(wts2) <- colnames(datl$y[,-5])
names(dimnames(datl$y)) <- c("Year", "Site")
all.df <- melt(as.matrix(datl$y[,-5]), value.name="Count") 
all.df$weights <- melt(as.matrix(wts2), value.name="weights")$weights
colnames(all.df) <- c("Year", "Site", "Count", "Weights")
labels <- c("Utqiagvik", "Bylott Island", "Karupelv Valley", 
            "Northern Norway", "Wrangel")

sl <- all.df[!is.na(all.df$Count),]

p1 <- ggplot()  + theme_minimal() +
  coord_cartesian(clip = "off") +
  geom_line(data=all.df,aes(x=Year,y=Site,color=Weights),size=2)+
  geom_point(data=sl,aes(x=Year,y=Site),shape="I",size=4)+
  xlab("Year")+ylab("Site")+scale_color_viridis_c()+
  scale_y_discrete(labels=labels)+
  labs(color='Weight') +
  xlim(1988, 2020) + theme(legend.position = "top") + 
  annotate(geom = "text", x = 1990, y = 5, label = "A", size=6)


p2 <- ggplot() + theme_minimal() + 
  geom_line(data=lam.post, aes(x=Year, y=value, group=Iter), 
            color="gray40", size=0.5, alpha=0.05 ) +
  geom_line(data=wm.df, aes(x=year, y=lci), color="deepskyblue3", size=0.5 ) +
  geom_line(data=wm.df, aes(x=year, y=uci ), color="deepskyblue3", size=0.5) +
  geom_line(data=wm.df, aes(x=year, y=lci80), color="deepskyblue3", size=1 ) +
  geom_line(data=wm.df, aes(x=year, y=uci80), color="deepskyblue3", size=1 ) +
  geom_hline(yintercept=1, lwd=2, color="black", linetype="dashed") +
  geom_line(data=wm.df, aes(x=year, y=m), color="deepskyblue3", size=2 ) +
  xlab("Year") +
  ylab("Population growth rate (lambda)")+
  xlim(1988, 2020) + 
  annotate(geom = "text", x = 1990, y = 1.2, label = "B", size=6)


#both <- grid.arrange(p1, p2, nrow = 2)
align_plots(p1, p2, align="v", axis="l")
both_aligned <- plot_grid(p1, p2, nrow = 2, align="v")
ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\pop-growth-rate_year_noWrangel.tiff", 
       plot=both_aligned, 
       device="tiff",
       width=6, 
       height=8,
       dpi=300)

# tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\pop-growth-rate_year_noWrangel.tiff",
#      height=6, width=6, units="in", res=300)
# plot(1988:2020, wm2.m,
#      xlab="Year", cex.lab=1.3,
#      ylab="Population growth rate (lambda)",
#      ylim=c(0.8, 1.2),
#      type="n", lwd=4,
#      xaxt="n", yaxt="n")
# axis(1, at=c(1990,2000,2010, 2020) )
# axis(2, at=c(0.8, 0.9, 1, 1.1, 1.2))
# for (j in 1:ncol(wm.post.lam)){
#   lines(1988:2020, wm.post.lam[,j], col = c3)
# }
# lines(1988:2020, wm2.lci, col=c4, lwd=2)
# lines(1988:2020, wm2.uci, col=c4, lwd=2)
# lines(1988:2020, wm2.lci80, col=c4, lwd=2)
# lines(1988:2020, wm2.uci80, col=c4, lwd=2)
# abline(h=1, lty=2, lwd=4)
# lines(1988:2020, wm2.m, lwd=4, col=c4)
# dev.off()


# Probabity of direction
# for each year
pd.func <- function(x){
  if (median(x) < 1){
    x.tf <- x<1
    pd <- sum(x.tf)/length(x)
  } else{
    x.tf <- x>1
    pd <- sum(x.tf)/length(x)
  }
  return(pd)}

pd <- apply(wm.post.lam[-1,], 1, pd.func)
pd.tab <- data.frame( year=1989:2020, pd=pd, 
                      md=wm2.m[-1], 
                      lci80=wm2.lci80[-1], uci80=wm2.uci80[-1],
                      lci95=wm2.lci[-1], uci95=wm2.uci[-1])
print(pd.tab, digits=3)

# calculate mean over entire time period
lam <- colMeans(wm.post.lam, na.rm=T)
data.frame( md = median(lam),
            lhdi95 = hdi(lam)[1],
            uhdi95 = hdi(lam)[2],
            lhdi80 = hdi(lam, credMass=0.8)[1],
            uhdi80 = hdi(lam, credMass=0.8)[2],
            pd.all = pd.func(lam))
# Investigate time periods 
# 1989–2010 and 2011–2018 
l.early <- colMeans(wm.post.lam[1:21,], na.rm=T)
data.frame( md = median(l.early),
            lhdi95 = hdi(l.early)[1],
            uhdi95 = hdi(l.early)[2],
            lhdi80 = hdi(l.early, credMass=0.8)[1],
            uhdi80 = hdi(l.early, credMass=0.8)[2],
            pd.all = pd.func(l.early))

l.late <- colMeans(wm.post.lam[22:29,], na.rm=T)
data.frame( md = median(l.late),
            lhdi95 = hdi(l.late)[1],
            uhdi95 = hdi(l.late)[2],
            lhdi80 = hdi(l.late, credMass=0.8)[1],
            uhdi80 = hdi(l.late, credMass=0.8)[2],
            pd.all = pd.func(l.late))

pd.func2 <- function(x){
  if (median(x) < 1){
    x.tf <- x<0
    pd <- sum(x.tf)/length(x)
  } else{
    x.tf <- x>0
    pd <- sum(x.tf)/length(x)
  }
  return(pd)}
pd.func2(l.late-l.early)

#***********************
#* table with parameter estimates
#***********************
alll <- MCMCpstr(nb, params = c("mu", "beta", "sigma", "r"), 
                 type="chains")
allm <- do.call(rbind, alll)
tab <- data.frame(
  Median = apply(allm, 1, median),
  LHDI95 = apply(allm, 1, hdi, credMass=0.95)[1,],
  UHDI95 = apply(allm, 1, hdi, credMass=0.95)[2,],
  LHDI80 = apply(allm, 1, hdi, credMass=0.80)[1,],
  UHDI80 = apply(allm, 1, hdi, credMass=0.80)[2,]
)
tab$sig95 <- ifelse(tab$LHDI95>0 | tab$UHDI95<0, TRUE, FALSE )
tab$sig80 <- ifelse(tab$LHDI80>0 | tab$UHDI80<0, TRUE, FALSE )
tab$Site <- rep(nms, 3)
tab$Parameter <- c(rep("mu", 5), rep("beta", 5), 
                  rep("r", 5))
tab <- tab[, c(8, 9, 1:7)]
tab <- tab[order(tab$Site), ]
print(tab, digits=3)
tab[,c(3:7)] <- round(tab[,c(3:7)],3)
write.csv(tab, 
          file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\model_estimates.csv")
