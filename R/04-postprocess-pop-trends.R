## ---- start --------
library("zoo")
library ("HDInterval")
library ("MCMCvis")
library ("ggplot2")
library ("reshape2")
library("gridBase")
library("grid")
library ("gridExtra")
library ("cowplot")
library("bayestestR")
library ("ggdist")
library ("viridis")
options(scipen=999)
load("data/data.Rdata")
# load your negative binomial model here
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")
# Function for transparency in base R plots
makeTransparent<-function(someColor, alpha=100){ 
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

## ---- autocorr --------
# check for autocorrelation
par(mfrow=c(3,3), mar=c(4,4,3,1))
acf(dat$Utqiagvik_213[-c(1:6)], main="Utqiagvik", xlab="")
acf(dat$BylotIsland_100[-c(1:7)], main="BylotIsland_100", xlab="", ylab="")
acf(dat$BylotIsland_300[-c(1:14,35)], main="BylotIsland_300", xlab="", ylab="")
acf(dat$'Igloolik Island_114'[-c(1:25)], main="Igloolik Island", xlab="")
acf(dat$'Karupelv Valley_75'[-c(1:2)], main="Karupelv Valley", xlab="", ylab="")
acf(dat$'Hochstetter Forland_100'[!is.na(dat$`Hochstetter Forland_100`)], main="Hochstetter Forland" , xlab="", ylab="")
acf(dat$Fennoscandia_x[!is.na(dat$Fennoscandia_x)], main="Fennoscandia" )
acf(dat$Wrangel[!is.na(dat$Wrangel)], main="Wrangel" , ylab="")
acf(dat$Alaska[!is.na(dat$Alaska)] , main="Alaska", ylab="")

## ---- rawdata --------
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

## ---- parest --------
#***********************
#* table with parameter estimates
#***********************
library ("knitr")
nms <- c("Utqiagvik", "Bylot Island Core", 
"Karupelv Valley", 
"Fennoscandia", "Wrangel")
alll <- MCMCpstr(nb, params = c("mu", "beta", "r","sigma.time"), 
                 type="chains")
allm <- do.call(rbind, alll)
tab <- data.frame(
  Median = apply(allm, 1, median),
  LHDI95 = apply(allm, 1, HDInterval::hdi, credMass=0.95)[1,],
  UHDI95 = apply(allm, 1, HDInterval::hdi, credMass=0.95)[2,],
  LHDI80 = apply(allm, 1, HDInterval::hdi, credMass=0.80)[1,],
  UHDI80 = apply(allm, 1, HDInterval::hdi, credMass=0.80)[2,],
  pd = apply(allm, 1, pd)
)
tab$sig95 <- ifelse(tab$LHDI95>0 | tab$UHDI95<0, TRUE, FALSE )
tab$sig80 <- ifelse(tab$LHDI80>0 | tab$UHDI80<0, TRUE, FALSE )
tab$Site <- c(rep(nms, 3), "All sites")
tab$Parameter <- c(rep("mu", 5), rep("beta", 5), 
                   rep("r", 5), "sigma.time")
tab <- tab[, c(9, 10, 1:8)]
tab <- tab[order(tab$Site), ]
tab[,c(3:8)] <- round(tab[,c(3:8)],3)

knitr::kable(tab[,1:8], digits=2,
             caption="Table S2. Coefficient estimates from the negative binomial model.",
             row.names=FALSE,
             col.names = c("Site", "Parameter", "Median", 
                           "95% Lower HDI", "95% Upper HDI", 
                           "85% Lower HDI", "85% Upper HDI", 
                           "Prob. direction"))

#write.csv(tab, 
#          file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\model_estimates.csv")

## ---- abund --------
#******************
#* Plot abundance for each site
#******************
beta <- MCMCpstr(nb, "beta", type="chains")[[1]]
mu <- MCMCpstr(nb, "mu", type="chains")[[1]]
pred.lam <- array(NA, dim=c(length(constl$time), nrow(mu), ncol(mu)))

for (i in 1:nrow(mu)){
for (t in 1:length(constl$time)){
    pred.lam[t,i,] <- mu[i] + beta[i]*constl$time[t]
  }}

# set NAs for unsurveyed years
wna <- which(is.na(datl$y), arr.ind=TRUE)
for (i in 1:nrow(wna)){
  pred.lam[wna[i,1], wna[i,2], ] <- NA
}

beta.tab <- data.frame(
  median = apply(beta, 1, median),
  LHDI95 = apply(beta, 1, HDInterval::hdi, credMass=0.95)[1,],
  UHDI95 = apply(beta, 1, HDInterval::hdi, credMass=0.95)[2,],
  LHDI80 = apply(beta, 1, HDInterval::hdi, credMass=0.80)[1,],
  UHDI80 = apply(beta, 1, HDInterval::hdi, credMass=0.80)[2,],
  pd = apply(beta, 1, pd), 
  nms = nms
)

beta.tab$sig95 <- ifelse(beta.tab$LHDI95>0 | beta.tab$UHDI95<0, TRUE, FALSE )
beta.tab$sig80 <- ifelse(beta.tab$LHDI80>0 | beta.tab$UHDI80<0, TRUE, FALSE )
beta.tab$sigpd <- ifelse(beta.tab$pd>0.9 | beta.tab$pd<0.1, TRUE, FALSE )
beta.tab <- beta.tab[, c(7, 1:6, 8:10)]
print(beta.tab, digits=2)

lam.md <- apply(pred.lam, c(1,2), median, na.rm=F) |> exp()
lam.lci95 <- apply(pred.lam, c(1,2), HDInterval::hdi, na.rm=F)[1,,] |> exp()
lam.uci95 <- apply(pred.lam, c(1,2), HDInterval::hdi, na.rm=F)[2,,] |> exp()
lam.lci80 <- apply(pred.lam, c(1,2), HDInterval::hdi, na.rm=F, credMass=0.80)[1,,] |> exp()
lam.uci80 <- apply(pred.lam, c(1,2), HDInterval::hdi, na.rm=F, credMass=0.80)[2,,] |> exp()
lam.lci67 <- apply(pred.lam, c(1,2), HDInterval::hdi, na.rm=F, credMass=0.67)[1,,] |> exp()
lam.uci67 <- apply(pred.lam, c(1,2), HDInterval::hdi, na.rm=F, credMass=0.67)[2,,] |> exp()

cols <- viridis(nrow(mu))
cols2 <- viridis(nrow(mu))
c3 <- makeTransparent("gray40", alpha=10)

tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\trends_noOffset.tiff",
     height=4, width=6, units="in", res=300)
par(mfrow=c(1,1), 
    mar=c(0,2,0,0), 
    oma=c(5,5,1,1))
for (i in c(1)){
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

## ---- lambda --------
#***********************
#* Plot an overall trend
#* weighted by abundance
#*********************** 
# weight by modeled abundance for each year
# weight by data, mean abundance across years
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\nb.Rdata")
beta <- MCMCpstr(nb, "beta", type="chains")[[1]]
mu <- MCMCpstr(nb, "mu", type="chains")[[1]]
pred.lam <- array(NA, dim=c(length(constl$time), nrow(mu),  ncol(mu)))

for (i in 1:nrow(mu)){
  for (t in 1:length(constl$time)){
    pred.lam[t,i,] <- mu[i,] + 
                      beta[i,]*constl$time[t] 
  }}
# set NAs for unsurveyed years
wna <- which(is.na(datl$y), arr.ind=TRUE)
for (i in 1:nrow(wna)){
  pred.lam[wna[i,1],wna[i,2], ] <- NA
}

wmns <- as.matrix(datl$y+1)
wmns <- ifelse(is.na(wmns), 0, wmns)
zo <- ifelse(is.na(datl$y), 0, 1)
wts2 <- array(NA, dim(zo)) 

for (t in 1:nrow(zo)){
  mns <- wmns[t,]*zo[t,] # impute zeroes for unmonitored years
  wts2[t,1:5] <- unlist(mns/sum(mns, na.rm=T)) # calculate proportions
}

lam.post <- array (NA, dim(pred.lam))
# calculate lambda (population growth rate)
for (t in 2:nrow(pred.lam)){
  for (j in 1:ncol(pred.lam)){
    lam.post[t,j,] <- exp(pred.lam[t,j,]) / exp(pred.lam[t-1,j,])
  }}
# calculate a weighted mean (wm)
# of population growth rate
wm.post.lam <- array (NA, c(nrow(pred.lam),dim(pred.lam)[[3]]), dimnames=list(Year=1988:2020, Iter=1:4000))
for (t in 1:nrow(pred.lam)){
  for (k in 1:dim(pred.lam)[[3]]){
    wm.post.lam[t,k] <- weighted.mean(x=lam.post[t,,k], w=wts2[t,], na.rm=T)
  }}
lam.post.long <- melt(wm.post.lam)

# calculate overall average pop growth rate
av.post <- apply(wm.post.lam, 2, mean, na.rm=T)
median(av.post, na.rm=T)
HDInterval::hdi(av.post, credMass = 0.95, na.rm=T)
HDInterval::hdi(av.post, credMass = 0.80, na.rm=T)
pd(av.post, null=1, na.rm=T)

# calculate pop growth rates for each year
wm.df <- data.frame(
  year= 1988:2020,
  m= apply(wm.post.lam, 1, median, na.rm=T),
  lci= apply(wm.post.lam, 1, HDInterval::hdi)[1,],
  uci= apply(wm.post.lam, 1, HDInterval::hdi)[2,],
  lci80 = apply(wm.post.lam, 1, HDInterval::hdi, credMass=.8)[1,],
  uci80 = apply(wm.post.lam, 1, HDInterval::hdi, credMass=.8)[2,],
  pd = apply(wm.post.lam, 1, pd, null=1)
)
wm.df 

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
  geom_line(data=all.df,aes(x=Year,y=Site,color=Weights), linewidth=2)+
  geom_point(data=sl,aes(x=Year,y=Site),shape="I", size=4)+
  xlab("Year")+ylab("Site")+scale_color_viridis_c()+
  scale_y_discrete(labels=labels)+
  labs(color='Weight') +
  xlim(1988, 2020) + theme(legend.position = "top") + 
  annotate(geom = "text", x = 1990, y = 6, label = "A", size=6)

p2 <- ggplot() + theme_minimal() + 
      geom_line(data=lam.post.long, aes(x=Year, y=value, group=Iter), 
                color="gray40", linewidth=0.5, alpha=0.05 ) +
      geom_hline(yintercept=1, lwd=2, color="black", linetype="dashed") +
      geom_line(data=wm.df, aes(x=year, y=lci), color="deepskyblue3", linewidth=0.5 ) +
      geom_line(data=wm.df, aes(x=year, y=uci ), color="deepskyblue3", linewidth=0.5) +
      geom_line(data=wm.df, aes(x=year, y=lci80), color="deepskyblue3", linewidth=1 ) +
      geom_line(data=wm.df, aes(x=year, y=uci80), color="deepskyblue3", linewidth=1 ) +
      geom_line(data=wm.df, aes(x=year, y=m), color="deepskyblue3", linewidth=2 ) +
      xlab("Year") +
      ylab( expression(paste("Population growth (", lambda,")")) )+
      xlim(1988, 2020) + 
      annotate(geom = "text", x = 1990, y = 1.2, label = "B", size=6)+
      coord_cartesian(xlim=c(1988, 2020), ylim=c(0.85, 1.15))

ap12 <- align_plots(p1, p2, align="v", axis="l")
p12 <- plot_grid(ap12[[1]], ap12[[2]], nrow = 2, align="v")
## ---- lambdaplot --------
p12

# ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\pop-growth-rate_year.tiff", 
#        plot=p12, 
#        device="tiff",
#        width=6, 
#        height=8,
#        dpi=300)

## ---- perchange1 --------
#***********
#* Calculate proportion of base year
#* or percent change
#* Short generation time
#* 8-year generations
#*************
# calculate lambda from 1996
startyr <- 2020-(3*8) # 3 generations of 8 years
dnames <- list(year=1988:2020,
               site=nms,
               Iter=1:4000)

iucn <- array(NA, dim=dim(wm.post.lam), dimnames=dimnames(wm.post.lam))
# abund1996 <- apply(pred.lam[9,,], c(2), mean, na.rm=T)
iucn[9,] <- 1 #abund1996/mean(abund1996)
for (t in 10:nrow(wm.post.lam)){
iucn[t,] <- iucn[(t-1),] * wm.post.lam[t,]  
}
for (t in 8:1){
  iucn[t,] <- iucn[(t+1),] / wm.post.lam[(t+1),]  
}

# convert to percent change
iucn <- (iucn-1) * 100
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
             caption="Table S3. Percent change since 1996.")


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
  xlab("") +
  ylab( "Percent change")+
  coord_cartesian(xlim=c(1988, 2020), ylim=c(-100, 60)) +
  annotate(geom = "text", x = 1995, y = 55, label = "A", size=8) +
  theme(axis.text= element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks = element_line(color = "black"), 
        plot.title = element_text(size=22)) +
  labs(title="8-year generation time")

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
  xlab("") + ylab("") +
  annotate(geom = "text", x = 55, y = 0.25, label = "C", size=8) +
  theme(axis.text= element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks = element_line(color = "black")) +
  coord_flip(xlim=c(-100, 60), ylim=c(0, 1), clip="on") 

# calculate proportion of distribution in each category and mode
cuttab <- table(cut(iucn2020$value, breaks=c(-100,-80, -50, -30, -20, 300)))
prop <- cuttab/sum(cuttab)

df.iucn <- data.frame(IUCN.criteria=rev( c('Least concern', 'Near threatened', 'Vulnerable', 'Endangered', 'Critically endangered')),
                      Proportion.within=prop,
                      Proportion.within.and.worse= cumsum(prop)
)
df.iucn <- df.iucn[nrow(df.iucn):1,]
knitr::kable(df.iucn, digits=c(0, 0, 2, 2),
             row.names=FALSE,
             col.names = c("IUCN Category", "A2 Criteria", "Prop within", "Prop within and worse"),
             caption="Table S4. Percent change over three generations.")
# write.csv(df.iucn, 
#           "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\IUCN.csv")

# calculate the mode
dens <- density(iucn2020$value)
# mode
dens$x[dens$y==max(dens$y)]
# median
median(iucn2020$value)
# mean
mean(iucn2020$value)


#***********
#* Calculate proportion of base year
#* or percent change
#* Longer generation time
#* 10.7 yr generations
#*************
# calculate lambda from 1998
max.gen.time <- (2020-1988)/3
startyr <- 1988 # 3 generations of 8 years
dnames <- list(year=1988:2020,
               site=nms,
               Iter=1:4000)

iucn <- array(NA, dim=dim(wm.post.lam), dimnames=dimnames(wm.post.lam))
# abund1996 <- apply(pred.lam[9,,], c(2), mean, na.rm=T)
iucn[1,] <- 1 #abund1996/mean(abund1996)
for (t in 2:nrow(wm.post.lam)){
  iucn[t,] <- iucn[(t-1),] * wm.post.lam[t,]  
}

# convert to percent change
iucn <- (iucn-1) * 100
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

p5 <- ggplot() + theme_minimal() + 
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
  annotate(geom = "text", x = 1995, y = 55, label = "B", size=8) +
  theme(axis.text= element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks = element_line(color = "black"), 
        plot.title = element_text(size=22)) +
  labs(title="10.7-year generation time")

iucn2020 <- iucn.post[iucn.post$Year==2020, ] 
i.df2020 <- i.df[i.df$year==2020,]

p6 <- ggplot() + theme_minimal() +
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
  annotate(geom = "text", x = 55, y = 0.25, label = "D", size=8) +
  theme(axis.text= element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks = element_line(color = "black")) +
  coord_flip(xlim=c(-100, 60), ylim=c(0, 1), clip="on") 

ap56 <- align_plots(p5, p6, align="h", axis="l")
p56 <- plot_grid(ap56[[1]], ap56[[2]], nrow = 1, align="h", rel_widths = c(2, 1)) 

ap56 <- align_plots(p3, p4, p5, p6, align="h", axis="l")
p56 <- plot_grid(ap56[[1]], ap56[[2]],
                 ap56[[3]], ap56[[4]],
                 nrow = 2, align="h", rel_widths = c(2, 1)) 

## ---- perchangeplot2 --------
p56

# ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\figs\\percentchange_year2.tiff",
#        plot=p56,
#        device="tiff",
#        width=8,
#        height=8,
#        dpi=300)

# calculate proportion of distribution in each category and mode
cuttab <- table(cut(iucn2020$value, breaks=c(-100,-80, -50, -30, -20, 300)))
prop <- cuttab/sum(cuttab)

df.iucn <- data.frame(IUCN.criteria=rev( c('Least concern', 'Near threatened', 'Vulnerable', 'Endangered', 'Critically endangered')),
                      Proportion.within=prop,
                      Proportion.within.and.worse= cumsum(prop)
)
df.iucn <- df.iucn[nrow(df.iucn):1,]
knitr::kable(df.iucn, digits=c(0, 0, 2, 2),
             row.names=FALSE,
             col.names = c("IUCN Category", "A2 Criteria", "Prop within", "Prop within and worse"),
             caption="Table S6. Percent change over three generations.")
# write.csv(df.iucn,
#           "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\docs\\IUCN2.csv")

# calculate the mode
dens <- density(iucn2020$value)
# mode
dens$x[dens$y==max(dens$y)]
# median
median(iucn2020$value)
# mean
mean(iucn2020$value)

