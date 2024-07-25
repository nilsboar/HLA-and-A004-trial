####################################################################################################
# for archive 
####################################################################################################
library(ggplot2)
library(stats)
source("~/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/janssen code archive/allele assoc functions rkv.R")
####################################################################################################
# load demographic data (Hispanic ethnicity for Whites, 
# trial arm data (used here for gp140 boost as covariate); placebo never typed for HLA
####################################################################################################
setwd("/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data")
load("hlademo.Rdata")
hlademo$hispanic <- ifelse(hlademo$ETHNIC == "HISPANIC OR LATINO", 1, 0) # need PCA, who are the missings?
####################################################################################################
# load data, different assays, different analysis (PCs, means)
####################################################################################################
setwd("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output")
# need PCA data here
load("elispot response.Rdata") # elispot gve is here # includes coeffs up to 5
# load("elispot revised means.Rdata")
setwd("/Users/nelsong/Documents/analysis/HLA/HLA 2019 on/Janssen/Janssen data/data")
load("scaled jansics mean.Rdata") # correct file for published figures, per archive
# load("jansicsmean_newmeans.Rdata")
# load("jansicsmean_means_nofilter.Rdata")
elispot_hla <-  merge(hlademo, elispot.pca.avtimeag.coeffs, by.x = "PID", by.y = 0)
elispot_gve <-  merge(hlademo, elisgve.meantime, by.x = "PID", by.y = 0)
CD8gve <-  merge(hlademo, jansicsCD8mean, by.x = "PID", by.y = 0)
####################################################################################################
# select data for analysis, filenames, outcomes and class I and/or II; variables for frame
####################################################################################################
framevars <- c("assay", "outcome", "race", "group", "group count", "allele", "count", "dom freq", "coeff", "p", "FDR")
# assaydata <- elisa_hla
# assaydata <- elispot_hla
# assaydata <- elispot_gve
assaydata <- CD8gve
nperm <- 0
# suffix to name; assay, outcomes, prefixed in function
csvname <- " tst ics gag.csv" 
pdfname <- " tst ics gag.pdf" 
min_allele_ct <- 4
classes <- c("Class I")
# classes <- c("Class II")
outcomes <- c("gag.mns.env") 
outcomestext <- c("gag.mns.env") 
thisassay <- "ICS"
# thisassay <- "ELISpot"
# thisassay <- "ELISA"
setwd("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output")
regions <- c("East-Africa", "South-Africa", "USA")
Nperm <- 0
# regions <- NULL
allele.assoc <- NULL
for (race in c("WHITE", "BLACK OR AFRICAN AMERICAN"))
{
  if (race == "WHITE")
  {
    clingrp <- subset(assaydata, RACE == race & REGION1 == "USA")
    group <- "Whites"
    this_min_allele_ct <- max(ceiling(nrow(clingrp)/100), min_allele_ct)
    allele.assoc <- allele_association(clingrp, race, group, assay = thisassay, min_allele_ct = this_min_allele_ct, outcomes, classes, pdfname, csvname, tablevars = framevars, nperm = Nperm)
  } else if (race == "ASIAN")  {
    clingrp <- subset(assaydata, RACE == race & REGION1 == "Thailand")
    group <- "Asians"
    this_min_allele_ct <- max(ceiling(nrow(clingrp)/100), min_allele_ct)
    this.allele.assoc <- allele_association(clingrp, race, group, assay = thisassay, min_allele_ct = this_min_allele_ct, outcomes, classes, pdfname, csvname, tablevars = framevars, nperm = Nperm)
    if (is.null(allele.assoc)) allele.assoc <- this.allele.assoc
    else allele.assoc <- rbind(allele.assoc, this.allele.assoc)
  } else {
    clingrp <- subset(assaydata, RACE == race)
    group <- "Blacks"
    this_min_allele_ct <- max(ceiling(nrow(clingrp)/100), min_allele_ct)
    # browser()
    this.allele.assoc <- allele_association(clingrp, race, group,  assay = thisassay, min_allele_ct = this_min_allele_ct, outcomes, classes, pdfname, csvname, tablevars = framevars, nperm = Nperm)
    if (is.null(allele.assoc)) allele.assoc <- this.allele.assoc
    else allele.assoc <- rbind(allele.assoc, this.allele.assoc)
    # browser()
    # break
    for (region in regions)
    {
      group <- region
      grpclingrp <- subset(clingrp, REGION1 == region)
      this_min_allele_ct <- max(ceiling(nrow(clingrp)/100), min_allele_ct)
      this.allele.assoc <- allele_association(grpclingrp, race, group,  assay = thisassay, min_allele_ct = this_min_allele_ct, outcomes, classes, pdfname, csvname, tablevars = framevars, nperm = Nperm)
      allele.assoc <- rbind(allele.assoc, this.allele.assoc)
      # allele.counts <- rbind(allele.counts, allele.counts.rtn)
    }
  }
  write.csv(allele.assoc, file = paste(thisassay, outcomestext, Nperm, "permutations allele association.csv", sep = " "))
}


