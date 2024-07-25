library(tidyr)
library(ggplot2)
library(FactoMineR)

setwd("APPROACH/data")
######xxxxxx######xxxxxx######xxxxxx######xxxxxx######xxxxxx
setwd("/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data")
######xxxxxx######xxxxxx######xxxxxx######xxxxxx######xxxxxx
load("janclinic.Rdata")
load("hlademo.Rdata")

names(janclinic)
table(janclinic$ISMETHOD)
jansics0 <- subset(janclinic, ISMETHOD == "ICS Percentage Positive")
icsags <- unique(jansics0$ISTEST)
icsags <- icsags[order(icsags)]
icsmos <- icsags[grep("(Mos1)", icsags, fixed = T)]

ags <- goodags <- icsags[c(2,3,4,9:12,14,17,18)]
mosags <- ags[grep("(Mos1)", ags, fixed = T)]
jansics0 <- subset(jansics0, ISTEST %in% goodags & TSEQPG1N != 8)
length(unique(jansics0$PID))

########################################
# examine outliers
########################################
names(jansics0)
sd <- sd(jansics0$ISSTRESN)
mean <- mean(jansics0$ISSTRESN)
denom <- sd - mean
exp(-2)

x <- 1.5
dist <- x - mean
dens0 <- 1/(sd*sqrt(2*pi))
dense <- (1/(sd*sqrt(2*pi)))
reldense <- exp(-0.5*((x - mean)/sd)^2) 
cd8 <- subset(jansics0, ISCAT == "T cell response (CD8)")
1/nrow(cd8)

pnorm(1.25, sd = 0.28, mean = 0.0677, lower.tail = F)

########################################
# bad flow data:  0345, 0650; both visitday 183
########################################
jansics <- subset(jansics0, !(PID == "HA004-0345" & VISITDY == 183 & ISTEST == "HIV Gag pep pool (Mos1)" & ISCAT == "T cell response (CD8)"))
jansics <- subset(jansics, !(PID == "HA004-0650" & VISITDY == 183 & ISTEST == "HIV ENV gp41 pep pool 1 (Mos1)" & ISCAT == "T cell response (CD8)"))
# rm(jansics0)
# save(jansics, file = "/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output/jansics_unnorm.Rdata")
################################################################################
# initial processing removed 2 bad data points, antigens without complete data
################################################################################
######xxxxxx######xxxxxx######xxxxxx######xxxxxx######xxxxxx
load("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output/jansics_unnorm.Rdata")
######xxxxxx######xxxxxx######xxxxxx######xxxxxx######xxxxxx
setwd("APPROACH/processed data archive")
load("jansics_unnorm.Rdata")

hist(jansics$ISSTRESN, breaks = 200, xlab = "ICS raw data %positive", ylab = "Frequency", main = "Histogram, ICS %positive")
hist.data = hist(jansics$ISSTRESN, plot=F, breaks = 500)
hist.data$counts = log(hist.data$counts + 1, 10)
plot(hist.data, xlab = "ICS raw data %positive", ylab = "Log10 frequency", main = "Histogram, Log10 ICS %positive") # Supp Figure 
########################################
# negatives to 0
########################################
names(jansics)
jansics$ics <- pmax(0, as.numeric(jansics$ISORRES))
hist(jansics$ics, breaks = 200)
hist(subset(jansics, ics < 0.1)$ics, breaks = 200)
########################################
# normalize to elispot$logresponse
########################################
jansics_0values <- subset(jansics, ics == 0)
jansics_posvalues <- subset(jansics, ics > 0)
hist(jansics_posvalues$ics, breaks = 100)
load("elispot for normalization.Rdata")
quantin <- subset(elispot, logresponse > 0)$logresponse
eliquants <- quantile(quantin, probs = seq(0, 1, 1/9125), type = 7)
jansics_posvalues <- jansics_posvalues[order(jansics_posvalues$ics), ]
jansics_posvalues$ics_qscaled <- eliquants
plot(jansics_posvalues$ics_qscaled, jansics_posvalues$ics)
hist(jansics_posvalues$ics_qscaled, breaks = 200)
jansics_0values$ics_qscaled <- 0
jansics_scaled <- rbind(jansics_0values, jansics_posvalues)
jansics_scaled$ics <- jansics_scaled$ics_qscaled
hist(jansics_scaled$ics, breaks = 200, xlab = "ICS scaled data %positive", main = "Histogram, ICS data quantile scaled")

########################################
# from here use normalized ics response
########################################
# TEMP
jansicsCD4 <- subset(jansics_scaled, ISCAT == "T cell response (CD4)")
jansicsCD8 <- subset(jansics_scaled, ISCAT == "T cell response (CD8)")

######xxxxxx######xxxxxx######xxxxxx######xxxxxx######xxxxxx
save(jansicsCD8, jansicsCD4, file = "/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data/jansics qscaled.rdata")
######xxxxxx######xxxxxx######xxxxxx######xxxxxx######xxxxxx
setwd("APPROACH/processed data archive")
save(jansicsCD8, jansicsCD4, file = "jansics qscaled.rdata")

load("/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data/jansics qscaled.rdata")
for (CD in c("CD4", "CD8"))
{
  if (CD == "CD4") CDdata <- jansicsCD4
  else CDdata <- jansicsCD8
  ags <- unique(CDdata$ISTEST)
  ###########################################################
  # subset to Mosaic Ags; no done above
  ###########################################################
  # timesdata <- subset(CDdata, EPOCH %in% c("TREATMENT3", "TREATMENT4"))
  timesdata <- subset(CDdata, EPOCH %in% c("TREATMENT3", "TREATMENT4"))
  dayswide <- as.data.frame(pivot_wider(timesdata, names_from = c("ISTEST"), values_from = "ics", id_cols = c("PID", "EPOCH", "VISITDY"))) 
  responsedata <- dayswide
  subjects <- unique(responsedata$PID)
  # no baseline for ICS, so just average treatment3 and 4 (omit follow up)
  meanresponse <- matrix(nrow = length(subjects), ncol = length(ags), dimnames = list(subjects, ags))
  for (subject in subjects)
  {
    meanresponse[subject,] <- colMeans(subset(responsedata, PID == subject)[,ags], na.rm = T)
  }
  jansicsmean <- as.data.frame(meanresponse)
  
  envags <- subset(ags, substr(ags, 1, 7) == "HIV ENV")
  envpte <- envags[grep("PTE", envags)]
  envmos <- envags[grep("Mos1", envags)]
  envza <- envags[grep("ZA", envags)]
  gagags <- subset(ags, substr(ags, 1, 7) == "HIV Gag")
  polags <- subset(ags, substr(ags, 1, 7) == "HIV Pol")

  jansicsmean$mean_env_pte <- rowMeans(jansicsmean[,envpte], na.rm = T)
  jansicsmean$mean_env_mos <- rowMeans(jansicsmean[,envmos], na.rm = T)
  jansicsmean$mean_env_za <- rowMeans(jansicsmean[,envza], na.rm = T)
  env_composite <- c("mean_env_pte", "mean_env_mos", "mean_env_za")
  jansicsmean$mean_env_preavg <- rowMeans(jansicsmean[,env_composite], na.rm = T)
  jansicsmean$mean_env <- rowMeans(jansicsmean[,envags], na.rm = T)
  jansicsmean$mean_env_nopte <- rowMeans(jansicsmean[,env_composite[2:3]], na.rm = T)
  jansicsmean$mean_gag <- jansicsmean[,gagags]
  jansicsmean$mean_pol <- rowMeans(jansicsmean[,polags], na.rm = T)
  jansicsmean$env.mns.gag <- jansicsmean$mean_env - jansicsmean$mean_gag
  jansicsmean$gag.mns.env.preavg <- jansicsmean$mean_gag - jansicsmean$mean_env_preavg
  jansicsmean$gag.mns.env <- jansicsmean$mean_gag - jansicsmean$mean_env
  jansicsmean$gag.mns.env.nopte <- jansicsmean$mean_gag - jansicsmean$mean_env_nopte
  jansicsmean$env.mns.gagpol <- jansicsmean$mean_env - (jansicsmean$mean_gag + jansicsmean$mean_pol)
  jansicsmean$envgagpoll <- jansicsmean$mean_env + jansicsmean$mean_gag + jansicsmean$mean_pol
  # scale gve; "mean_all" are means over all ags and individuals
  mean_all_env <- mean(jansicsmean$mean_env, na.rm = T) 
  mean_all_gag <- mean(jansicsmean$mean_gag, na.rm = T)
  jansicsmean$mean_env_scld <- jansicsmean$mean_env*mean_all_gag/mean_all_env
  jansicsmean$gag.mns.env.scld <- jansicsmean$mean_gag - jansicsmean$mean_env_scld
  if (CD == "CD4") jansicsCD4mean <- jansicsmean
  else jansicsCD8mean <- jansicsmean
}

plot(jansicsCD8mean$mean_env_pte, jansicsCD8mean$mean_env_mos)
summary(lm(mean_env_pte ~ mean_env_mos, data = jansicsCD8mean))

# regular
setwd("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output")
save(jansicsCD4mean, jansicsCD8mean, file = "scaled jansics mean.Rdata")








