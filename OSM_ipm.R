################################################################################
# Integrated population model for Oil Sands point count, ARU, and MAPS data
#
# In this version...
#
# Last edited 21 Oct. 2024, J. Saracco
################################################################################

## load libraries, functions ---------------------------------------------------
library(data.table)
library(tidyverse)
Sys.setenv(JAGS_HOME = "C:\\Program Files\\JAGS\\JAGS-4.3.1")
library(jagsUI)
library(sf)

# function to add transparency to plot colors
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# transformations for logit-linear models 
logit <- function(x){
  y <- log(x/(1-x))
  return(y)
}
expit <- function(x){
  y <- exp(x)/(1+exp(x))
  return(y)
}

#-------------------------------------------------------------------------------
# 1. Read in and process the point, aru, and covariate data
#-------------------------------------------------------------------------------

# 1.1 - Owl Moon point count data-----------------------------------------------
# Owl Moon point count data through 2023
omei_count_data <- read.csv("data/OMEI_pointcounts_woffsets_2011-2023.csv") %>% 
  subset(select = c(Year, Method, Station.Code, Species.Code, Total.Number, lat, 
                    lon, date, dur, dis, time, offset), subset = !is.na(offset)) %>% 
  mutate(date = as.Date(parse_date_time(date, "mdy"))) 
omei_count_data <- omei_count_data %>% 
  mutate(date_time = ymd_hm(paste(date, time, sep = " "), tz = "America/Edmonton")) %>% 
  mutate(sunrise = suncalc::getSunlightTimes(data = data.frame(date = ymd(omei_count_data$date), 
                                                               lat = omei_count_data$lat,
                                                               lon = omei_count_data$lon),
                                             tz = "America/Denver")$sunrise) %>% 
  mutate(tssr = as.numeric(difftime(date_time, sunrise, units="hours"))) %>% 
  rename(sensor = Method, location_id = Station.Code, count = Total.Number, 
         species_code = Species.Code, year = Year) %>% 
  mutate(day = yday(date), task_method = "PC", project_id = "MAPS", dur = as.character(dur)) %>% 
  subset(day < 196 & day > 156 & tssr > -0.1667 & tssr < 1.5)

# 1.2 - ABMI count data (individual observations from point count and ARU data)-
# From Elly Knight: include task_method detection model: For ARU data, "1SPM" =
# 1 tag per individual per minute and "1SPT" = 1 tag per individual per 
# task/recording). 1SPM has higher counts.

abmi_count_data <- read.csv("data/OSMIPMSupplementalARUData_5.csv") %>% 
  mutate(date_time = ymd_hms(date_time), date = ymd(date), sunrise = ymd_hms(sunrise)) %>% 
  subset(day < 196 & day > 156 & tssr > -0.1667 & tssr < 1.5 & year > 2010) %>% 
  group_by(sensor, project_id, task_method, duration, location_id, lat, lon, date_time,
           species_code, date, day, year, offset) %>% 
  rename(dur = duration) %>% 
  summarise(count = sum(abundance)) %>% 
  mutate(location_id = as.character(location_id), project_id = as.character(project_id)) %>% 
  as.data.frame()

# 1.3 - covariate data -----------------------------

covariate_data <- read.csv("data/OSMIPMCovariates_5.csv")

# 1.4 - merge OMEI and ABMI data, subset ovenbirds -----------------------------
counts <- bind_rows(omei_count_data, abmi_count_data) %>% 
  group_by(location_id, project_id, year) %>% 
  mutate(rep = frank(date_time, ties.method = "dense"))


library(MODISTools) # for EVI data download

# download and process MODIS data
get_modis <- function(dat, product, band, qa, start_year, end_year, start_mmdd,
                      end_mmdd, km_buff){
  require(MODISTools)
  start_year2 <- ifelse(as.numeric(substr(end_mmdd, 1, 2)) < as.numeric(substr(start_mmdd, 1, 2)),
                        start_year + 1, start_year)
  
  x <- mt_batch_subset(dat, product, band, start = paste(start_year, start_mmdd, sep = "-"),
                       end = paste(start_year2, end_mmdd, sep = "-"),
                       km_lr = km_buff, km_ab = km_buff, internal = TRUE)
  x$value[x$value < -2000] <- NA # valid range: -2000 to 10000 (-3000 = fill value)
  if(!missing(qa)){
    q <- mt_batch_subset(dat, product, band = qa,
                         start = paste(start_year[1], start_mmdd, sep = "-"),
                         end = paste(start_year2, end_mmdd, sep = "-"),
                         km_lr = km_buff, km_ab = km_buff, internal = TRUE)
    x$value[q$value > 1] <- NA # keep only marginal and good quality pixels
    rm(q)
  }
  for (i in (start_year + 1):end_year){
    sy2 <- ifelse(start_year2 > start_year, i + 1, i)
    y <- mt_batch_subset(dat, product, band,
                         start = paste(i, start_mmdd, sep = "-"),
                         end = paste(sy2, end_mmdd, sep = "-"),
                         km_lr = km_buff, km_ab = km_buff, internal = TRUE)
    y$value[y$value < -2000] <- NA # valid range: -2000 to 10000 (-3000 = fill value)
    if(!missing(qa)){
      q <- mt_batch_subset(dat, product, band = qa,
                           start = paste(i, start_mmdd, sep = "-"),
                           end = paste(sy2, end_mmdd, sep = "-"),
                           km_lr = km_buff, km_ab = km_buff, internal = TRUE)
      y$value[q$value > 1] <- NA # keep only marginal and good quality pixels
      x <- bind_rows(x, y)
      rm(y, q, sy2)
    }
  }
  xsum <- x %>%
    group_by(latitude, longitude, site, calendar_date) %>%
    summarize(evi = mean(value, na.rm = TRUE))
  xsum$evi <- xsum$evi/10000
  xsum$calendar_date <- as.Date(xsum$calendar_date)
  xsum$Month <- as.character(format(xsum$calendar_date, format="%m"))
  xsum$Year <- as.character(format(xsum$calendar_date, format="%Y"))
  xsumSEASON <- ifelse(xsum$Month %in% c("11", "12"), as.numeric(xsum$Year) + 1, as.numeric(xsum$Year))
  
  return(list(x, xsum))
}

locs <- counts %>% subset(select = c(location_id, lat, lon, year)) %>% unique()


dat <- data.frame(site_name = locs$location_id, lat = locs$lat, lon = locs$lon)
modis_out <- get_modis(dat = dat, product = "MOD13Q1", band = "250m_16_days_NDVI",
                            qa = "250m_16_days_pixel_reliability",
                            start_year = 2011, end_year = 2023,
                            start_mmdd = "06-01", end_mmdd = "07-30",
                            km_buff = 0)

write.csv(modis_out[[1]], "data/modis_out1.csv", row.names = FALSE)
write.csv(modis_out[[2]], "data/modis_out2.csv", row.names = FALSE)

modis_sum <- modis_out[[2]] %>% group_by(latitude, longitude, site, Year) %>% 
  summarize(evi = mean(evi, na.rm = TRUE)) %>% 
  rename(lat = latitude, lon = longitude, location_id = site, year = Year) %>% 
  mutate(year = as.integer(year))



ef <- unique(subset(counts, select = c(sensor, location_id, project_id, task_method, year, lat, lon, dis, dur, date_time, sunrise, tssr, rep)))

oven_counts <- counts %>% subset(species_code %in% "OVEN")
oven_counts <- right_join(oven_counts, ef)
oven_counts$count[is.na(oven_counts$count)] <- 0
oven_counts <- as.data.frame(oven_counts)

oven_counts$dur.num <- as.numeric(oven_counts$dur)
oven_counts$dur.num[oven_counts$dur %in% "0-3-5-10min"] <- 10
oven_counts$dur.num[oven_counts$dur %in% "120s"] <- 2
oven_counts$dur.num[oven_counts$dur %in% "180s"] <- 3
oven_counts$dur.num[oven_counts$dur %in% "300s"] <- 5
oven_counts$dur.num[oven_counts$dur %in% "540s"] <- 9
oven_counts$dur.num[oven_counts$dur %in% "600s"] <- 10
oven_counts$dur.num[oven_counts$dur %in% "60s"] <- 1

x <- aggregate(count ~ task_method + dur.num, oven_counts, length)
# 1        1SPM       1     4
# 2        1SPM       2     1
# 3        1SPM       3  4689
# 4        1SPT       3   958
# 5        1SPM       5     9
# 6        1SPM       9     5
# 7        1SPM      10   231
# 8          PC      10  1002
oven_counts <- subset(oven_counts, dur.num == 3 | dur.num == 10)
oven_counts$sens_dur <- paste(oven_counts$sensor, oven_counts$dur.num, sep=".")
oven_counts$meth_dur <- paste(oven_counts$task_method, oven_counts$dur.num, sep=".")

# plot(oven_counts$offset ~ oven_counts$dur.num)
# # for nmix model
oven_counts_wide <- oven_counts %>%
  pivot_wider(id_cols = c(sensor, project_id, task_method, meth_dur, sens_dur, location_id, lat, lon, year),
              names_from = rep, names_sort = TRUE, values_from = c(count, offset)) %>%
  as.data.frame()
#-------------------------------------------------------------------------------
# 2. Read in and process MAPS CMR data
#-------------------------------------------------------------------------------

# note that in this version age codes: 1 = AHY, 2 = SY, 3 = ASY and 
# sex codes: 1 = F, 2 = M, 3 = U

ch <- foreign::read.dbf("data/OSM_CH.dbf") %>% 
  subset(!(LOSS %in% "-1") & SPEC %in% "OVEN", 
         c(STA, STATION, BAND, SPEC, AGE, SEX, Y2011:Y2023, MARKED))
ch$STATION  <- as.character(ch$STATION)
ch$STATION[ch$STA %in% "18812"] <- "BPND"
ch$STATION[ch$STA %in% "18814"] <- "VWET"
ch$STATION <- as.factor(ch$STATION)
ch[ch == "N"] <- NA # missed occasions
ch[,7:19] <- lapply(ch[,7:19], as.numeric) 
ch[,7:19] <- ch[,7:19]-1
ch$last <- ch$first <- NA
for(i in 1:nrow(ch)){
  h <- as.vector(ch[i,7:19])
  ch$first[i] <- min( (1:13)[!is.na(h) & h==1])
  ch$last[i]  <- max( (1:13)[!is.na(h) & h==1])
}

# Process age structure data
sydat <- ch %>% pivot_longer(cols = Y2011:Y2023, names_to = "year") %>% 
  mutate(year = as.numeric(factor(year))) %>% 
  subset(value == 1)
sydat$AGE[sydat$AGE == 1] <- NA
sydat$age.init <- NA
sydat$age.init[sydat$AGE == 2] <- 1
sydat$age.init[sydat$AGE == 3] <- 0
sydat$age2 <- ifelse(sydat$year == sydat$first, sydat$age.init, 0)

ch <- subset(ch, first < 13) # first encounter on last occasion does not inform CJS model

sink("ipm_zip.txt")
cat("
model{
    
# Priors------------------------------------------------------------------------

# Method intercept for detection probability. levels 1 and 2 are methods of 
# counting from ARU recordings. Level 3 is human point count
for (i in 1:nmeths){
lpct0[i] <- log(pct0[i]/(1-pct0[i]))
pct0[i] ~ dunif(0, 1)
}

pocc ~ dunif(0, 1)

for (j in 1:nlocs){
zeta[j] ~ dbern(pocc)
}
for (i in 1:ncts){ 
N[i] ~ dpois(ntot[YR[i]]*zeta[location[i]])
  for (j in 1:nrep){
# Count model-------------------------------------------------------------------
    C[i, j] ~ dbin(pct[i, j], N[i])
    logit(pct[i, j]) <- lpct0[meth[i]] 
    # Expected count 
    C.exp[i, j] <- pct[i, j] * N[i]
    # Data discrepancy 
    chi2.act[i,j] <- pow(C[i,j] - C.exp[i,j], 2)/(C.exp[i, j] +.001)
    #Simulated counts 
    C.rep[i, j] ~ dbin(pct[i, j], N[i])
    # Replicate discrepancy 
    chi2.sim[i,j] <- pow(C.rep[i,j] - C.exp[i,j], 2)/(C.exp[i, j] +.001)
  } # j
} # i

# chi-squared test statistics
chi2.fit <- sum(chi2.act[,])
chi2.fit.rep <- sum(chi2.sim[,])

# Bayesian p-value
bpvalue <- step(chi2.fit - chi2.fit.rep)
    
################################################################################
# Population process model
################################################################################
# Priors for initial population size
nrecr[1] ~ dunif(0, 10)
nsurv[1] ~ dunif(0, 10)
nimm[1] ~ dunif(0, 10)

# Population dynamics
for (t in 1:nyears){
  ntot[t] <- nsurv[t] + nrecr[t] + nimm[t]
}

tau.lomega <- pow(sigma.lomega, -2)
sigma.lomega ~ dunif(0, 10)
lomega.mn <- log(omega.mn)
omega.mn ~ dunif(0, 10)

tau.lgamma <- pow(sigma.lgamma, -2)
sigma.lgamma ~ dunif(0, 10)
lgamma.mn <- log(gamma.mn)
gamma.mn ~ dunif(0, 10)

for (t in 1:(nyears-1)){

  # Model for immigration
    log.omega[t] ~ dnorm(lomega.mn, tau.lomega)
    omega[t] <- exp(log.omega[t])

  # Model for recruitment
    log.gamma[t] ~ dnorm(lgamma.mn, tau.lgamma)
    gamma[t] <- exp(log.gamma[t])
}

for (t in 2:nyears){
  # Models for numbers of survivors, recruits, and adult immigrants
  nsurv[t] <- ntot[t-1]*phit[t-1]
  nrecr[t] <- ntot[t-1]*gamma[t-1]
  nimm[t] <- ntot[t-1]*omega[t-1]
} 

################################################################################
# MAPS CJS model for adult survival
################################################################################

# Priors for intercepts    
p0 ~ dunif(0, 1)
lp0 <- log(p0/(1-p0))
rho0 ~ dunif(0, 1)
lrho0 <- log(rho0/(1-rho0))

# Random station effects for p and rho models 
sigma.psta ~ dunif(0, 10)
sigma.rhosta ~ dunif(0,10)
tau.psta <- pow(sigma.psta, -2)
tau.rhosta <- pow(sigma.rhosta, -2)
    
for(j in 1:nsta){
  stap[j] ~ dnorm(0, tau.psta)
  starho[j] ~ dnorm(0, tau.rhosta)
}

# year variance for phi and pi models    
sigma.lphi ~ dunif(0, 10)
sigma.lpi ~ dunif(0, 10)
tau.lphi <- pow(sigma.lphi, -2)
tau.lpi <- pow(sigma.lpi, -2)

phi.mn ~ dunif(0, 1)
lphi.mn <- logit(phi.mn)
pi.mn ~ dunif(0, 1)
lpi.mn <- logit(pi.mn)
    
for (t in 1:(nyears - 1)){
  # Year-specific survival estimate
  logit.phi[t] ~ dnorm(lphi.mn, tau.lphi)
  phit[t] <- ilogit(logit.phi[t])
  logit.pi[t] ~ dnorm(lpi.mn, tau.lpi)
}

# State-space likelihood for i individuals    
for (i in 1:nind){
    
# Latent alive state at first occasion
  for(t in 1:first[i]){
    z[i,t] ~ dbern(1)
  }
  # Residency state based on first capture occasion
  R[i] ~ dbern(pi[i,first[i]])
  mu[i] <- R[i]*rho[i]
  r[i] ~ dbern(mu[i])
    
  # logit-linear models
  logit(p[i]) <- lp0 + stap[sta[i]] 
  logit(rho[i]) <- lrho0 + starho[sta[i]] 
  for (t in first[i]:(nyears-1)){
    pi[i,t] <- ilogit(logit.pi[t])
    phi[i,t] <- ilogit(logit.phi[t])
  }
  for (t in (first[i]+1):nyears){
    # State process
    mu2[i,t] <- z[i,t-1]*phi[i,t-1]*R[i]
    z[i,t] ~ dbern(mu2[i,t])
    # Observation process
    mu1[i,t] <- p[i]*z[i,t]
    y[i,t] ~ dbern(mu1[i,t])
    pred_y[i,t] <- mu1[i,t]
    ft1[i, t] <- (sqrt(y[i, t]) - sqrt(pred_y[i, t]))^2 
    y.rep[i,t] ~ dbern(mu1[i,t])
    ft1.rep[i,t] <- (sqrt(y.rep[i,t]) - sqrt(pred_y[i,t]))^2 
  } # end t
    fit_i[i] <- sum(ft1[i,(first[i]+1):nyears])          
    fit.rep_i[i] <- sum(ft1.rep[i,(first[i]+1):nyears])
} # end i
    fit <- sum(fit_i[])           # test statistic for data
    fit.rep <- sum(fit.rep_i[])   # test statistic for new predicted data
    cjs.bpvalue <- step(fit.rep - fit)   # Test whether new data set more extreme

################################################################################
# Age structure model based on MAPS SY/ASY data
################################################################################    

tau.psista <- pow(sigma.psista, -2)
sigma.psista ~ dunif(0, 10)

for (i in 1:nsta){
  sta.psi[i] ~ dnorm(0, tau.psista)
}
for (i in 1:M){
  age[i] ~ dbern(psi[i])
  logit(psi[i]) <- logit(nrecr[year[i]]/ntot[year[i]]) + sta.psi[station[i]]
}
}
        ",fill = TRUE)
sink()


# Define data for CJS model
r <- as.numeric(ch$MARKED) 
y <- as.matrix(ch[,7:19])
class(y) <- "numeric"
sta <- as.numeric(factor(ch$STA))
nsta <- max(sta)
nind <- nrow(y)
first <- ch$first
nyears <- dim(y)[2]

# Initial state values for CJS model
Zst<-y
Zst[Zst==0]<-1
Zst[is.na(Zst)]<-1
Rst<-rep(1,nrow(y))

# Define data for age-structure model
M <- nrow(sydat)
age <- sydat$age2
year <- sydat$year
station <- as.numeric(factor(sydat$STA))

C <- subset(oven_counts_wide, select = count_1:count_4)
# pc.ind <- ifelse(oven_counts_wide$sensor %in% "PC", 1, 0)

evi <- (oven_counts_wide$evi - mean(oven_counts_wide$evi, na.rm=TRUE))

ncts <- dim(C)[1]
proj <- as.numeric(factor(oven_counts_wide$project_id))
nproj <- max(proj)
meth <- as.numeric(factor(oven_counts$meth_dur))
nmeths <- max(meth)
location <- as.numeric(factor(oven_counts_wide$location_id))
nlocs <- max(location)
YR <- as.numeric(oven_counts_wide$year) - 2010


# Initial values for population size of n-mixture model
oven_counts_wide$N.st <- apply(C, 1, max, na.rm=T)
nrecr.st <- nsurv.st <- nimm.st <- aggregate(N.st ~ year, oven_counts_wide, mean)$N.st/3

jags_dat <- list(C = C, ncts = ncts, nrep = 4, meth = meth, nmeths = nmeths, pc.ind = pc.ind,
                 YR = YR, nyears = nyears, r = r, y = y, sta = sta, nsta = nsta, 
                 location=location, nlocs=nlocs,
                 nind=nind, first = first, M = M, age = age, station = station, year = year)

inits <- function(){list(N = oven_counts_wide$N.st, sigma.rhosta = runif(1), pocc = runif(1),
                         sigma.psta = runif(1), nsurv1=nsurv.st[1], nrecr1=nrecr.st[1], 
                         nimm1=nimm.st[1], pct0 = runif(nmeths, .2, .6), sigma.lpi = runif(1), 
                         method = runif(nmeths, -2, 2), sigma.lphi = runif(1), 
                         p0 = runif(1, .2, .5), rho0 = runif(1, 0.1, 0.5), 
                         phi.mn = runif(1, 0.1, 0.5), pi.mn = runif(1, 0.2, 0.8), R = Rst, 
                         z = Zst, sigma.yr = runif(1), sigma.psista = runif(1), tau.e = 1,
                         gamma.mn = runif(1, .1, 1.1), omega.mn = runif(1, .1, 1.1), 
                         sigma.lgamma = runif(1), sigma.lomega = runif(1))}

parameters <- c("nrecr", "nsurv", "nimm", "chi2.fit", "chi2.fit.rep", "bpvalue", "cjs.bpvalue",
                "phit", "phi.mn", "pimn", "p0", "gamma", "omega", "omega.mn", "pct0",
                "gamma.mn", "b1", "sigma.lgamma", "pocc", 
                "sigma.lphi", "sigma.lpi", "sigma.lomega")

# 
## MCMC settings
ni <- 120000
nc <- 3
na <- 40000
nb <- 40000
nt <- 5

## Fit model
oven_ipm_fit_newdat <- jags(data = jags_dat, inits = inits, parameters.to.save = parameters, 
                      model.file = "ipm_zip.txt", n.chains = nc, n.iter = ni, n.adapt = na,
                      n.burnin = nb, n.thin = nt, parallel = TRUE)

ntot <- oven_ipm_fit_newdat$sims.list$nrecr + oven_ipm_fit_newdat$sims.list$nsurv + 
  oven_ipm_fit_newdat$sims.list$nimm

ntot.med <- apply(ntot, 2, median)
ntot.li <- apply(ntot, 2, quantile, probs = 0.05)
ntot.ui <- apply(ntot, 2, quantile, probs = 0.95)

png("Report_figs/oven_N_newdat.png", width = 5, height = 3, units = "in", res = 600)
par(mar = c(2,4.5,1,1))
plot(x = seq(2011, 2023, 1), y= ntot.med, type = "n", las = 1, ylim = c(0, 1.2*max(ntot.ui)), 
     ylab = expression(paste(N[t]^"tot", " (birds/point)", sep = "")), xlab = "", axes=FALSE)
axis(1)
axis(2, las=1)
polygon(x = c(seq(2011, 2023, 1), rev(seq(2011, 2023, 1))), 
        y =c(ntot.li, rev(ntot.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
points(x = seq(2011, 2023, 1), y = ntot.med, type = "l", lwd=2)
dev.off()

Trend <- 100*((ntot[,13]/ntot[,1])^(1/12)-1)
median(Trend)
quantile(Trend, probs = c(0.05, 0.95))

Trend10 <- 100*((ntot[,10]/ntot[,1])^(1/9)-1)
median(Trend10)
quantile(Trend10, probs = c(0.05, 0.95))

oven_data_locs <- oven_counts_wide %>% 
  subset(select = c(sensor, project_id, location_id, lat, lon)) %>% 
  unique()
oven_data_locs$type <- factor(NA, levels = c("ARU", "PC", "PC+Banding"))
oven_data_locs$type[oven_data_locs$sensor %in% "ARU"] <- "ARU"
oven_data_locs$type[oven_data_locs$sensor %in% "PC"] <- "PC"
oven_data_locs$type[oven_data_locs$location_id %in% unique(ch$STATION)] <- "PC+Banding"

write.csv(oven_data_locs, "data/oven_data_locs.csv", row.names = FALSE)
oven_locs_sf <- st_as_sf(oven_data_locs, coords = c("lon", "lat"), crs = st_crs("+proj=longlat +datum=WGS84"))

#Write it as kml file
st_write(oven_locs_sf, "oven_locs.kml", layer="zinc", driver="KML") 

pct.med <- apply(oven_ipm_fit1$sims.list$pct0, 2, median)
pct.li <- apply(oven_ipm_fit1$sims.list$pct0, 2, quantile, probs = 0.05)
pct.ui <- apply(oven_ipm_fit1$sims.list$pct0, 2, quantile, probs = 0.95)

library(plotrix)
png("Report_figs/detec.png", height = 4, width = 6.5, res = 600, units = "in")
par(mar = c(4,4,1,1))
plotCI(y = pct.med[c(2,3,1,4)], x = c(1, 2, 3, 4), li = pct.li[c(2,3,1,4)], ui = pct.ui[c(2,3,1,4)], 
       xlim = c(0.5, 4.5), ylim = c(0, 1), axes = FALSE,
       ylab = "Detection probability", xlab = "", pch = 16, cex = 1)
axis(1, at = c(1, 2, 3, 4), labels = FALSE)
text(x = 1:4,
     y = par("usr")[3]-.05,
     labels = c("1SPM-3min", "1SPT-3min", "1SPM-10min", "PC-10min"),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     pos=2,
     ## Increase label size.
     cex = 1)
axis(2, las = 1)
box()
dev.off()

# mean demographic rates

quantile(oven_ipm_fit1$sims.list$phi.mn, probs = c(0.05, 0.5, 0.95))
quantile(oven_ipm_fit1$sims.list$omega.mn, probs = c(0.05, 0.5, 0.95))
quantile(oven_ipm_fit1$sims.list$gamma.mn, probs = c(0.05, 0.5, 0.95))

phi.med <- apply(oven_ipm_fit1$sims.list$phit, 2, median)
phi.li <- apply(oven_ipm_fit1$sims.list$phit, 2, quantile, probs = 0.05)
phi.ui <- apply(oven_ipm_fit1$sims.list$phit, 2, quantile, probs = 0.95)

gam.med <- apply(oven_ipm_fit1$sims.list$gamma, 2, median)
gam.li <- apply(oven_ipm_fit1$sims.list$gamma, 2, quantile, probs = 0.05)
gam.ui <- apply(oven_ipm_fit1$sims.list$gamma, 2, quantile, probs = 0.95)

omega.med <- apply(oven_ipm_fit1$sims.list$omega, 2, median)
omega.li <- apply(oven_ipm_fit1$sims.list$omega, 2, quantile, probs = 0.05)
omega.ui <- apply(oven_ipm_fit1$sims.list$omega, 2, quantile, probs = 0.95)

lam.med <- apply(ntot[,2:13]/ntot[,1:12], 2, median)
lam.li <- apply(ntot[,2:13]/ntot[,1:12], 2, quantile, probs = 0.05)
lam.ui <- apply(ntot[,2:13]/ntot[,1:12], 2, quantile, probs = 0.95)

# ### realized demographic rates
# phi.real.med <- apply(oven_ipm_fit1$sims.list$nsurv[,2:13]/ntot[,1:12], 2, median)
# phi.real.05 <- apply(oven_ipm_fit1$sims.list$nsurv[,2:13]/ntot[,1:12], 2, quantile, probs = 0.05)
# phi.real.25 <- apply(oven_ipm_fit1$sims.list$nsurv[,2:13]/ntot[,1:12], 2, quantile, probs = 0.25)
# phi.real.75 <- apply(oven_ipm_fit1$sims.list$nsurv[,2:13]/ntot[,1:12], 2, quantile, probs = 0.75)
# phi.real.95 <- apply(oven_ipm_fit1$sims.list$nsurv[,2:13]/ntot[,1:12], 2, quantile, probs = 0.95)
# 
# gam.real.med <- apply(oven_ipm_fit1$sims.list$nrecr[,2:13]/ntot[,1:12], 2, median)
# gam.real.05 <- apply(oven_ipm_fit1$sims.list$nrecr[,2:13]/ntot[,1:12], 2, quantile, probs = 0.05)
# gam.real.25 <- apply(oven_ipm_fit1$sims.list$nrecr[,2:13]/ntot[,1:12], 2, quantile, probs = 0.25)
# gam.real.75 <- apply(oven_ipm_fit1$sims.list$nrecr[,2:13]/ntot[,1:12], 2, quantile, probs = 0.75)
# gam.real.95 <- apply(oven_ipm_fit1$sims.list$nrecr[,2:13]/ntot[,1:12], 2, quantile, probs = 0.95)
# 
# omega.real.med <- apply(oven_ipm_fit1$sims.list$nimm[,2:13]/ntot[,1:12], 2, median)
# omega.real.05 <- apply(oven_ipm_fit1$sims.list$nimm[,2:13]/ntot[,1:12], 2, quantile, probs = 0.05)
# omega.real.25 <- apply(oven_ipm_fit1$sims.list$nimm[,2:13]/ntot[,1:12], 2, quantile, probs = 0.25)
# omega.real.75 <- apply(oven_ipm_fit1$sims.list$nimm[,2:13]/ntot[,1:12], 2, quantile, probs = 0.75)
# omega.real.95 <- apply(oven_ipm_fit1$sims.list$nimm[,2:13]/ntot[,1:12], 2, quantile, probs = 0.95)
# 
# 
# png("Report_figs/dem-realized.png", width = 6.5, height = 5, units = "in", res = 600)
# par(mar = c(2,4.5,1,1), mfrow = c(2,2))
# 
# par(mar = c(2,4.5,1,1))
# plot(x = seq(2011, 2022, 1), y= lam.real.med, type = "n", las = 1, ylim = c(0, max(lam.real.95)),
#      ylab = expression(paste("Population growth (", "N"[t+1]/"N"[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(lam.real.05, rev(lam.real.95)), border = NA, col = add.alpha("gray70", alpha=.5))
# # polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(lam.real.25, rev(lam.real.75)), border = NA, col = add.alpha("gray70", alpha=.7))
# points(x = seq(2011, 2022, 1), y = lam.real.med, type = "l", lwd=2)
# legend("topleft", "A", bty = "n")
# 
# plot(x = seq(2011, 2022, 1), y= phi.real.med, type = "n", las = 1, ylim = c(0, 1), 
#      ylab = expression(paste("Realized adult survival (", "N"[t+1]^"surv"/"N"[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(phi.real.05, rev(phi.real.95)), border = NA, col = add.alpha("gray70", alpha=.5))
# # polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(phi.real.25, rev(phi.real.75)), border = NA, col = add.alpha("gray70", alpha=.7))
# points(x = seq(2011, 2022, 1), y = phi.real.med, type = "l", lwd=2)
# legend("topleft", "B", bty = "n")
# 
# plot(x = seq(2011, 2022, 1), y= gam.real.med, type = "n", las = 1, ylim = c(0, 1.5), 
#      ylab = expression(paste("Realized recruitment (", "N"[t+1]^"recr"/"N"[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(gam.real.05, rev(gam.real.95)), border = NA, col = add.alpha("gray70", alpha=.5))
# # polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(gam.real.25, rev(gam.real.75)), border = NA, col = add.alpha("gray70", alpha=.7))
# points(x = seq(2011, 2022, 1), y = gam.real.med, type = "l", lwd=2)
# legend("topleft", "C", bty = "n")
# 
# plot(x = seq(2011, 2022, 1), y= omega.real.med, type = "n", las = 1, ylim = c(0, 1.5), 
#      ylab = expression(paste("Realized immigration (", "N"[t+1]^"imm"/"N"[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(omega.real.05, rev(omega.real.95)), border = NA, col = add.alpha("gray70", alpha=.5))
# # polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(omega.real.25, rev(omega.real.75)), border = NA, col = add.alpha("gray70", alpha=.7))
# points(x = seq(2011, 2022, 1), y = omega.real.med, type = "l", lwd=2)
# legend("topleft", "D", bty = "n")
# 
# dev.off()

# 
# 
# png("Report_figs/oven_phi.png", width = 5, height = 3, units = "in", res = 600)
# par(mar = c(2,4.5,1,1))
# plot(x = seq(2011, 2022, 1), y= phi.med, type = "n", las = 1, ylim = c(0, 1), 
#      ylab = expression(paste("Adult survival probability (", phi[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(phi.li, rev(phi.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
# points(x = seq(2011, 2022, 1), y = phi.med, type = "l", lwd=2)
# plotCI(x = seq(2011, 2022, 1), y = phi.real.med, li = phi.real.05, ui = phi.real.95, sfrac=0, pch = 16, add=TRUE)
# dev.off()
# 
# png("Report_figs/oven_gamma.png", width = 5, height = 3, units = "in", res = 600)
# par(mar = c(2,4.5,1,1))
# plot(x = seq(2011, 2022, 1), y= gam.med, type = "n", las = 1, ylim = c(0, 1), 
#      ylab = expression(paste("Recruitment rate (", gamma[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(gam.li, rev(gam.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
# points(x = seq(2011, 2022, 1), y = gam.med, type = "l", lwd=2)
# plotCI(x = seq(2011, 2022, 1), y = gam.real.med, li = gam.real.05, ui = gam.real.95, sfrac=0, pch = 16, add=TRUE)
# dev.off()
# 
# png("Report_figs/oven_omega.png", width = 5, height = 3, units = "in", res = 600)
# par(mar = c(2,4.5,1,1))
# plot(x = seq(2011, 2022, 1), y= omega.med, type = "n", las = 1, ylim = c(0, 1.5), 
#      ylab = expression(paste("Immigration rate (", omega[t], ")"), xlab = "", axes=FALSE))
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(omega.li, rev(omega.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
# points(x = seq(2011, 2022, 1), y = omega.med, type = "l", lwd=2)
# plotCI(x = seq(2011, 2022, 1), y = omega.real.med, li = omega.real.05, ui = omega.real.95, sfrac=0, pch = 16, add=TRUE)
# dev.off()
# 
# png("Report_figs/oven_lambda.png", width = 5, height = 3, units = "in", res = 600)
# par(mar = c(2,4.5,1,1))
# plot(x = seq(2011, 2022, 1), y= lam.med, type = "n", las = 1, ylim = c(0, max(lam.ui)),
#      ylab = expression(paste("Population growth (", lambda[t], ")", sep = "")), xlab = "", axes=FALSE)
# axis(1)
# axis(2, las=1)
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(lam.li, rev(lam.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
# points(x = seq(2011, 2022, 1), y = lam.med, type = "l", lwd=2)
# 
# dev.off()

# correlations between demographic rates and realized lambda

png("Report_figs/oven_phi_v_lam.png", width = 3, height = 3, units = "in", res = 600)
par(mar = c(4.5,4.5,1,1))
plotrix::plotCI(x = phi.med, y = lam.real.med, li = phi.li, ui=phi.ui, sfrac = 0, err = "x", 
                pch = 16, cex = 1, las = 1, lwd = 2,  xlim = c(0, 0.5), 
                ylim = c(0, 2), xlab = expression(paste("Adult survival probability (", phi[t], ")")),
                ylab = expression(paste("Population growth (", lambda[t], ")")), axes = FALSE)
axis(1)
axis(2, las=1)
plotrix::plotCI(x = phi.med, y = lam.real.med, li = lam.real.05, ui=lam.real.95, sfrac = 0, err = "y", 
                pch = 16, cex = 1, las = 1, lwd = 2, add=TRUE)
dev.off()



png("Report_figs/dem-expected.png", width = 6.5, height = 5, units = "in", res = 600)
par(mar = c(2,4.5,1,1), mfrow = c(2,2))

par(mar = c(2,4.5,1,1))
plot(x = seq(2011, 2022, 1), y= lam.med, type = "n", las = 1, ylim = c(0, max(lam.ui)),
     ylab = expression(paste("Population growth (", "N"[t+1]/"N"[t], ")", sep = "")), xlab = "", axes=FALSE)
axis(1)
axis(2, las=1)
polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(lam.li, rev(lam.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
# polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(lam.real.25, rev(lam.real.75)), border = NA, col = add.alpha("gray70", alpha=.7))
points(x = seq(2011, 2022, 1), y = lam.med, type = "l", lwd=2)
legend("topleft", "A", bty = "n")

plot(x = seq(2011, 2022, 1), y= phi.med, type = "n", las = 1, ylim = c(0, 1), 
     ylab = expression(paste("Adult survival probability (", phi[t], ")")), xlab = "", axes=FALSE)
axis(1)
axis(2, las=1)
polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(phi.li, rev(phi.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
points(x = seq(2011, 2022, 1), y = phi.med, type = "l", lwd=2)
legend("topleft", "B", bty = "n")

plot(x = seq(2011, 2022, 1), y= gam.med, type = "n", las = 1, ylim = c(0, 1), 
     ylab = expression(paste("Recruitment rate (", gamma[t], ")", sep = "")), xlab = "", axes=FALSE)
axis(1)
axis(2, las=1)
polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(gam.li, rev(gam.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
points(x = seq(2011, 2022, 1), y = gam.med, type = "l", lwd=2)
legend("topleft", "C", bty = "n")


plot(x = seq(2011, 2022, 1), y= omega.med, type = "n", las = 1, ylim = c(0, 1.5), 
     ylab = expression(paste("Immigration rate (", omega[t], ")"), xlab = "", axes=FALSE))
axis(1)
axis(2, las=1)
polygon(x = c(seq(2011, 2022, 1), rev(seq(2011, 2022, 1))), y =c(omega.li, rev(omega.ui)), border = NA, col = add.alpha("gray70", alpha=.5))
points(x = seq(2011, 2022, 1), y = omega.med, type = "l", lwd=2)
legend("topleft", "D", bty = "n")

dev.off()






# vital rate correlations
phi.lam.cor <- rep(NA, dim(oven_ipm_fit1$sims.list$phit)[1])
for (i in 1:length(phi.lam.cor)){
  phi.lam.cor[i] <- cor.test(oven_ipm_fit1$sims.list$phit[i,], 
                             ntot[i,2:13]/ntot[i,1:12])$estimate
}

median(phi.lam.cor)
quantile(phi.lam.cor, probs = c(0.05, 0.95))
length(phi.lam.cor[phi.lam.cor > 0])/length(phi.lam.cor)



# vital rate correlations
imm.lam.cor <- rep(NA, dim(oven_ipm_fit1$sims.list$omega)[1])
for (i in 1:length(imm.lam.cor)){
  imm.lam.cor[i] <- cor.test(oven_ipm_fit1$sims.list$omega[i,], 
                             ntot[i,2:13]/ntot[i,1:12])$estimate
}

median(imm.lam.cor)
quantile(imm.lam.cor, probs = c(0.05, 0.95))
length(imm.lam.cor[imm.lam.cor > 0])/length(imm.lam.cor)


png("oven_gam_v_lam.png", width = 3, height = 3, units = "in", res = 600)
par(mar = c(4.5,4.5,1,1))
plotrix::plotCI(x = gam.med, y = lam.real.med, li = gam.li, ui=gam.ui, sfrac = 0, err = "x", 
                pch = 16, cex = 1, las = 1, lwd = 2,  xlim = c(0, 1), 
                ylim = c(0, 2), xlab = expression(paste("Recruitment rate (", gamma[t], ")")),
                ylab = expression(paste("Population growth (", lambda[t], ")")), axes = FALSE)
axis(1)
axis(2, las=1)
plotrix::plotCI(x = gam.med, y = lam.real.med, li = lam.real.05, ui=lam.real.95, sfrac = 0, err = "y", 
                pch = 16, cex = 1, las = 1, lwd = 2, add=TRUE)
dev.off()


# vital rate correlations
gam.lam.cor <- rep(NA, dim(oven_ipm_fit1$sims.list$gamma)[1])
for (i in 1:length(gam.lam.cor)){
  gam.lam.cor[i] <- cor.test(oven_ipm_fit1$sims.list$gamma[i,], 
                             ntot[i,2:13]/ntot[i,1:12])$estimate
}

median(gam.lam.cor)
quantile(gam.lam.cor, probs = c(0.05, 0.95))
length(gam.lam.cor[gam.lam.cor > 0])/length(gam.lam.cor)


