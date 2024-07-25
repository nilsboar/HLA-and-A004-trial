library(ggplot2)
library(grid)
library(MASS)
# library(gggrid)
# code to call functions testing and plotting the association of cohort viral load with allele vaccine response, 
# or differential response, e.g. gag - env
################################################################################################
# load VL data
################################################################################################
setwd("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/janssen input data rkv")
load("onevisit blacks whites noswiss.Rdata")
load("hlademo.Rdata")
################################################################################################
# load vaccine response, HLA, and demographic data
######################################################################################load("ics CD8 pca and gve.Rdata")
setwd("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/janssen processed data rkv")
load("scaled jansics mean.Rdata")
load("elispot response.Rdata")
ics_resp <- merge(jansicsCD8mean, hlademo, by.x = 0, by.y = "PID")
elispot_resp <- merge(elisgve.meantime, hlademo, by.x = 0, by.y = "PID")
################################################################################################
# set up to run allele vs. VL function
######################################################################################load("ics CD8 pca and gve.Rdata")
setwd("/Users/nelsong/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/rkv test output")
outvars <- c("locus", "race", "region", "coeff", "p", "Rsquared", "coeff.alleles", "p.alleles", "Rsquared.alleles", "perm p", "n permute", "perm p alleles")
vl.vs.resp.mtrx <- matrix(nrow = 10, ncol = length(outvars), dimnames = list(NULL, outvars))
vl.vs.resp <- as.data.frame(vl.vs.resp.mtrx)[0,]
mincount = 4
minfreq = 0.03
npermute <- 100
permsets <- 1 # # of calls to function generating random data, each returns an npermute-long matrix
allele_plot <- T
whichresponse <- elispot_resp
whichassay <- "elispot"
# whichresponse <- ics_resp
# whichassay <- "ICS"
# whichrespvar <- c("gag.mns.env", "gag.mns.pol", "gagpol.mns.env", "gag.mns.envpol", "gagenv.mns.pol")
filenote <- "tst archive"
allele_plot <- T
# source("~/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/VL v janssen response functions.R")
regions <- c("East-Africa", "South-Africa", "USA")
allele.counts <- NULL
# regions <- NULL
source("~/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/Janssen vaccine/janssen code archive/VL v janssen response functions rkv.R")
################################################################################################
# loop to run allele vs. VL function
######################################################################################load("ics CD8 pca and gve.Rdata")
for (whichrespvar in "gag.mns.env")
{
  allele.counts <- NULL
  for (race in c("WHITE", "BLACK OR AFRICAN AMERICAN"))
  {
    if (race == "WHITE")
    {
      group <- "Whites"
      groupresponse <- subset(whichresponse, RACE == race)
      # browser()
      resp.return <- locus.resptest(vldata = onevisit_whites, jan.resp = groupresponse,  resp.var = whichrespvar, resp.mtrx = vl.vs.resp.mtrx, assay = whichassay, race = race, group = group, nperm = npermute, npermsets = permsets, mincount = mincount, minfreq = minfreq, alleleplot = allele_plot, filenamenote = filenote)
      this.vl.vs.resp <- resp.return[[1]]
      allele.counts.rtn <- resp.return[[2]]
      allele.counts.rtn$race <- "White"
      allele.counts.rtn$region <- "all"
      # write.csv(allele.counts.rtn, file = "allele counts Whites.csv")
      # write.csv(this.vl.vs.resp, file = "vlvresp Whites.csv")
      # write.csv(this.vl.vs.resp, file = paste(whichassay, whichrespvar, group, permsets*npermute, "perm.csv"))
      allele.counts <- allele.counts.rtn
      vl.vs.resp <- this.vl.vs.resp
    } else 
    {
      group <- "Blacks"
      groupresponse <- subset(whichresponse, RACE == race)
      resp.return <- locus.resptest(vldata = onevisit_blacks, jan.resp = groupresponse,  resp.var = whichrespvar, resp.mtrx = vl.vs.resp.mtrx, assay = whichassay, race = race, group = group, nperm = npermute, npermsets = permsets, mincount = mincount, minfreq = minfreq, alleleplot = allele_plot, filenamenote = filenote)
      this.vl.vs.resp <- resp.return[[1]]
      allele.counts.rtn <- resp.return[[2]]
      allele.counts.rtn$race <- "Black"
      allele.counts.rtn$region <- "all"
      # write.csv(allele.counts.rtn, file = "allele counts Blacks.csv")
      # write.csv(this.vl.vs.resp, file = paste(whichassay, whichrespvar, group, permsets*npermute, "perm.csv"))
      if (is.null(allele.counts)) 
      {
        allele.counts <- allele.counts.rtn
        vl.vs.resp <- this.vl.vs.resp 
      } else     
      {
        allele.counts <- rbind(allele.counts, allele.counts.rtn)
        vl.vs.resp <- rbind(vl.vs.resp, this.vl.vs.resp)
      }
      # break()
      for (region in regions)
      {
        groupresponse <- subset(whichresponse, RACE == race & REGION1 == region)
        group <- region
        resp.return <- locus.resptest(vldata = onevisit_blacks, jan.resp = groupresponse,  resp.var = whichrespvar, resp.mtrx = vl.vs.resp.mtrx, assay = whichassay, race = race, group = group, nperm = npermute, npermsets = permsets, mincount = mincount, minfreq = minfreq, alleleplot = allele_plot)
        this.vl.vs.resp <- resp.return[[1]]
        this.vl.vs.resp$region <- region
        allele.counts.rtn <- resp.return[[2]]
        allele.counts.rtn$race <- "Black"
        allele.counts.rtn$region <- region
        allele.counts <- rbind(allele.counts, allele.counts.rtn)
        vl.vs.resp <- rbind(vl.vs.resp, this.vl.vs.resp)
        
        # write.csv(allele.counts.rtn, file = paste("allele counts Blacks ", region, ".csv", sep = ""))
        # write.csv(this.vl.vs.resp, file = paste(whichassay, whichrespvar, group, permsets*npermute, "perm.csv"))
        
        # allele.counts <- rbind(allele.counts, allele.counts.rtn)
      }
    }
    write.csv(allele.counts, file = paste(filenote, "allele counts all", whichassay, whichrespvar, ".csv", sep = " "))
    write.csv(vl.vs.resp, file = paste(filenote, whichassay, whichrespvar, permsets*npermute, ".csv"))
  }
  
}
write.csv(vl.vs.gve, file = "VL vs gve test.csv")


