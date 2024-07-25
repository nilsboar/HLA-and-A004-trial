library(tidyr)
library(ggplot2)
library(MASS)
library(FactoMineR)
library(svd)


################################################################
# load  data
################################################################
setwd("/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data")
load("elispot.data.Rdata") 
names(elispot)
table(elispot$epoch, elispot$VISITNUM)
table(elispot$epoch, elispot$VISITDY)
# table(janclinic$EPOCH, janclinic$VISITDY)
################################################################
# use double PCA function to get means, PCA over ags
################################################################
source("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/janssen code archive/PCA function rkv.R") #
pdf(width = 11, height = 8, file = "/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output/tst elis pca.pdf")
avtime.list <- doublePCA(elispot, timepoints = unique(elispot$epoch), use.svd = F, use_meantime = T, baseline = T)
dev.off()
elispot.pca.avtimeag.coeffs <- avtime.list[[1]]
elispot.pca.avtimeag.coeffs[,2] <- -elispot.pca.avtimeag.coeffs[,2]

elismean <- avtime.list[[3]]

# get loadings for vector graph for PC3,4
elispot.pca.avtimeag.loadings <- t(avtime.list[[2]])
# loadings1 <- round(t(elispot.pca.avtimeag.loadingsq), 3)
# write.csv(loadings1, file = "/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/janssen gag mns env paper/september draft/JVI/output re reviewers queries/loadings.csv")
save(elispot.pca.avtimeag.loadings, file = "/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output/elispot loadings.Rdata")
# save(elispot.pca.avtimeag.loadingsq, file = "/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data/elispot loadings all PCS.Rdata")

################################################################
# mean ag response
################################################################

elisgve.meantime <- as.data.frame(elismean)
elisgve.meantime$mean.env <- rowMeans(elisgve.meantime[,c(2,8:11)], na.rm = T)
elisgve.meantime$mean.gag <- rowMeans(elisgve.meantime[,c(3,4,7)], na.rm = T)
elisgve.meantime$mean.pol <- rowMeans(elisgve.meantime[,c(1,5,6)], na.rm = T)
elisgve.meantime$gag.mns.env <- elisgve.meantime$mean.gag - elisgve.meantime$mean.env
hist(elisgve.meantime$gag.mns.env, breaks = 50)

save(elispot.pca.avtimeag.coeffs, elisgve.meantime, file = "/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output/elispot response.Rdata")
# save(elispot.pca.avtimeag.coeffs, elisgve.meantime, file = "/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data/elispot orig means all PCS.Rdata")


