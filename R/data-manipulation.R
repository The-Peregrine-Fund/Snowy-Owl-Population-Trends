library(readxl)
library(reshape2)
dat <- read_excel("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\SnowyOwl_HawkMountain\\data\\data.xlsx",
                  sheet="Modified", na="NA")
# change Bylot Island into 
# 2 separate sites
# one of 100 km^2 and 
# one of 300 km^2
dat$BylotIsland_300 <- dat$BylotIsland_400-dat$BylotIsland_100
colSums(!is.na(dat))

# Bundle data
# exclude iglook_100 and Hochstetter because <3 cycles (4 year cycles)
# remove Alaska because data are on different scale.
y <- as.data.frame(dat[-c(1,2), -c(4,5,7,9,11)]) 

# create data lists for nimble
y.nim <- y[,-c(1)]
datl <- list(y = y.nim) # exclude 8-fennoscania because no area information
constl <- list(nsite = ncol(y.nim),
               ntime = nrow(y.nim),
               area = c(213, 100, 75, NA, 45),
               time = (y$year-2004)/16
)

# Create long format data for mgcv GAMs
lf <- melt(as.matrix(y[,-1]))
colnames(lf) <- c("year.ind", "site", "count")
yr <- data.frame('year.ind'=1:33, year=1988:2020)
lf <- merge(lf, yr, by='year.ind')

save(file = "data/data.Rdata",
     datl=datl, constl=constl, 
     lf=lf, dat=dat)
