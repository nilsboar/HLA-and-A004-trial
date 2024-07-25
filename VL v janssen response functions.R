allele.resp <- function(resp.data, resp.var, min.allele.ct, min.allele.freq, locusalleles = c("A4d1", "A4d2", "B4d1", "B4d2", "C4d1", "C4d2"), callnperm = 0)
{
  # browser()
  resp.data$jan.resp <- resp.data[, resp.var]
  nsubs <- nrow(resp.data)
  min.freq <- min(min.allele.ct - 1, floor(min.allele.freq*nsubs))
  alleledata <- resp.data[,locusalleles[1]]
  for (i in 2:length(locusalleles)) alleledata <- c(alleledata, resp.data[,locusalleles[i]]) # assuming 3 loci
  allele.tbl <- table(alleledata)
  frequent <- names(which(allele.tbl > min.freq)) # 3 is 3/78 chromosomes; ~ 4%
  # browser()
  datacols <- as.matrix(resp.data[,locusalleles])
  alleles.resp <- as.data.frame(matrix(nrow = length(frequent), ncol = 3, dimnames = list(frequent, c("allele.ct", "allele.carrier.ct", "resp.score"))))
  allele.dat.list <- list()
  # names(allele.dat.list) <- frequent
  for (allele in frequent)
  {
    thisallele.dat <- resp.data
    for (i in 1:nrow(resp.data))
    {
      thisallele.dat[i, "hasallele"] <- allele %in% as.vector(datacols[i,])
    }
    thisallele.dat <- subset(thisallele.dat, hasallele) # make a list of these, for permutation
    allele.dat.list[[allele]] <- row.names(thisallele.dat)
    alleles.resp[allele, "allele.ct"] <- allele.tbl[allele]
    alleles.resp[allele, "allele.carrier.ct"] <- nrow(thisallele.dat) # use, unless there's evidence that homozygotes have stronger association
    alleles.resp[allele, "resp.score"] <- mean(thisallele.dat$jan.resp, na.rm = T)
    # if(nrow(thisallele.dat) == 0) browser()
    # browser()
  }
  row.names(alleles.resp) <- substr(row.names(alleles.resp), 1, 7)
  if (callnperm == 0) return(alleles.resp)
  perm.matrix <- matrix(nrow = length(frequent), ncol = callnperm, dimnames = list(frequent, NULL))
  for (i in 1:callnperm)
  {
    thisallele.dat <- resp.data
    thisallele.dat$permresp <- resp.data[sample.int(nrow(resp.data)), "jan.resp"]
    for (allele in frequent)
    {
      perm.matrix[allele, i] <- mean(thisallele.dat[allele.dat.list[[allele]],"permresp"], na.rm = T)
    }
    # row.names(alleles.resp) <- substr(row.names(alleles.resp), 1, 7)
  }
  # browser()
  return(as.data.frame(perm.matrix))
}

locus.resptest <- function(vldata, jan.resp, resp.var, resp.mtrx, assay, race, group, nperm = 0, npermsets = 1, mincount = 4, minfreq = 0.04, vldotplot = T, alleleplot = T, resp_powercalc = F, filenamenote = "", print_chartdata = T)
{
  # browser()
  permtotal <- npermsets*nperm
  alleles.resp <- allele.resp(jan.resp, resp.var, mincount, minfreq, callnperm = 0) # first call; no permutation
  # browser() # stop here to output allele response scores
  row.names(alleles.resp) <- substr(row.names(alleles.resp), 1, 7)
  # browser()
  vl.vs.resp <- as.data.frame(resp.mtrx)[0,]
  i <- 1
  # for (locus in c("B", "A", "C"))
  for (locus in c("B"))
  {
    locus.resp <- locus.alleles.resp <- alleles.resp[substr(row.names(alleles.resp),1,1) == locus,]
    # browser()
    row.names(locus.resp) <- locus.alleles.resp$shortallele <- substr(row.names(locus.resp), 3, 7) # vl data (onevisit) alleles have no locus prefix
    allele1 <- paste(locus, "1s", sep = "")
    allele2 <- paste(locus, "2s", sep = "")
    # browser()
    for (j in 1:nrow(locus.alleles.resp))
    {
      # CHECK THIS
      VL.withallele <- vldata[vldata[,allele1] == row.names(locus.resp)[j] | vldata[,allele2] == row.names(locus.resp)[j],]
      # browser()
      locus.alleles.resp[j, "N.VL"] <- nrow(VL.withallele)
      # browser()
      if (nchar(rownames(locus.alleles.resp)[j]) > 7) browser()
      locus.resp[j, "N.VL"] <- locus.alleles.resp[j, "N.VL"]
      locus.resp[j, "weight"] <- locus.alleles.resp[j, "weight"] <- locus.alleles.resp[j, "N.VL"] 
    }
    vldata$allelescore1 <- ifelse(is.na(locus.resp[vldata[,allele1],"resp.score"]), 0, locus.resp[vldata[,allele1],"resp.score"])
    vldata$allelescore2 <- ifelse(is.na(locus.resp[vldata[,allele2],"resp.score"]), 0, locus.resp[vldata[,allele2],"resp.score"])
    vldata$allelescore <- (vldata$allelescore1 + vldata$allelescore2)/2 # per text, is mean; 
    vldata$allelescore <- ifelse(is.na(locus.resp[vldata[,allele1],"resp.score"]) | is.na(locus.resp[vldata[,allele2],"resp.score"]), NA, vldata$allelescore) # & instead of     # regression, all VL points vs. allele score
    # browser()
    lm.vl <- lm(mVL ~ allelescore , data = vldata) # for combined blacks, can't separate groups; outcome = VL 
    ci <- confint(lm.vl)
    sumry <- summary(lm.vl)
    print(c(race, group))
    # print(sumry)
    # browser()
    vl.vs.resp[i, "locus"] <- locus
    vl.vs.resp[i, "race"] <- race
    vl.vs.resp[i, "region"] <- group
    vl.vs.resp[i, "N.VL"] <- nrow(vldata) # should be n rows with data for the locus
    vl.vs.resp[i, "N.2gve"] <- nrow(subset(vldata, !is.na(allelescore))) # should be n rows with data for the locus
    vl.vs.resp[i, "coeff"] <- sumry$coefficients[2,1]
    vl.vs.resp[i, "Rsquared"] <- sumry$r.squared
    vl.vs.resp[i, "p"] <- sumry$coefficients[2,4]
    vl.vs.resp[i, "n permute"] <- permtotal
    
    withci <- cbind(sumry$coefficients, ci)[2,]
    xlable <- paste("HLA-", locus, resp.var, sep = " ")
    locus.resp$allele <- row.names(locus.resp)
    # plot calcs moved up so we can get permutation value
    for (allele in row.names(locus.resp))
    {
      locus.resp[allele, "mean.VL"] <- mean((vldata[vldata[,allele1] == allele | vldata[,allele2] == allele,])$mVL) # mean over all carriers, dominant all_allelemodel
      # locus.resp$
    }
    mean.allelect <- mean(locus.resp$allele.carrier.ct, na.rm = T)
    mean.VLcount <- mean(locus.resp$N.VL, na.rm = T)
    # scale.factor <- mean.allelect/mean.VLcount # so bubbles for VL and for resp have similar size
    scale.factor <- 9 # rounding 8.93 for all blacks plot, need others for whites, subset (use round...)
    print(c(group, scale.factor))
    vlvresp.alleles <- lm(mean.VL ~ resp.score, weights = weight, data = locus.resp)
    sumry <- summary(vlvresp.alleles)
    # browser()
    if (resp_powercalc) 
    {
      return (resppower(locus.resp, jan.resp))
    }
    # browser()
    vl.vs.resp[i, "coeff.alleles"] <- sumry$coefficients[2,1]
    vl.vs.resp[i, "Rsquared.alleles"] <- sumry$r.squared
    vl.vs.resp[i, "p.alleles"] <- sumry$coefficients[2,4]
    intercept <- sumry$coefficients[1,1]
    slope <- sumry$coefficients[2,1]
    if(nperm > 0) # for both dot and bubble plots
    { 
      count.moresig <- 0
      count.moresig.alleles <- 0
      count.moresig.all.alleles <- 0
      permvldata <- vldata
      for (j in 1:npermsets)
      {
        # browser()
        permatrix <- allele.resp(jan.resp, resp.var, mincount, minfreq, callnperm = nperm) # returns nperm columns, each column is random resp scores for each allele
        permatrix <- permatrix[substr(row.names(permatrix),1,1) == locus,]
        permdata <- as.data.frame(permatrix)
        row.names(permdata) <- substr(row.names(permdata), 3, 7) # vl data (onevisit) alleles have no locus prefix
        for (k in 1:nperm)
        {
          if (k %% 1000 == 0) print(c(k, j))
          # browser()
          permvldata$allelescore1 <- ifelse(is.na(permdata[permvldata[,allele1],k]), 0, permdata[permvldata[,allele1],k])
          permvldata$allelescore2 <- ifelse(is.na(permdata[permvldata[,allele2],k]), 0, permdata[permvldata[,allele2],k])
          permvldata$allelescore <- permvldata$allelescore1 + permvldata$allelescore2
          # permvldata$allelescore <- ifelse(is.na(permdata[permvldata[,allele1],k]) & is.na(permdata[permvldata[,allele2],k]), NA, permvldata$allelescore) # & instead of | uses subjects with one allele with no gve score:
          permvldata$allelescore <- ifelse(is.na(permdata[permvldata[,allele1],k]) | is.na(permdata[permvldata[,allele2],k]), NA, permvldata$allelescore)
          # browser()
          perm.lm.vl <- lm(mVL ~ allelescore,  data = permvldata)
          perm.sumry  <- summary(perm.lm.vl)
          # print(c(j, perm.sumry$coefficients[2,4]))
          if (perm.sumry$coefficients[2,4] <= vl.vs.resp[i, "p"] & perm.sumry$coefficients[2,1]*vl.vs.resp[i, "coeff"] > 0) count.moresig <- count.moresig + 1 # one sided; signs must be the same
          # test for alleles
          locus.resp$perm.resp.score <- permdata[,k]
          perm.vlvresp.alleles <- lm(mean.VL ~ perm.resp.score, weights = weight, data = locus.resp) # using mean VL calculated above
          perm.sumry.alleles <- summary(perm.vlvresp.alleles)
          # print(c(j, perm.sumry.alleles$coefficients[2,4]))
          if (perm.sumry.alleles$coefficients[2,4] <= vl.vs.resp[i, "p.alleles"] & perm.sumry$coefficients[2,1]*vl.vs.resp[i, "coeff.alleles"] > 0) count.moresig.alleles <- count.moresig.alleles + 1 # one sided; signs must be the same
          ################################################################################################
          # # response to alleles themselves
          ################################################################################################
          # allele_tstdata$perm_mvl <- allele_tstdata[sample.int(n_2allele),"mVL"] 
          # perm_allelemdl <- as.formula(gsub("mVL", "perm_mvl", model))
          # perm.sumry.all_alleles <- summary(lm(perm_allelemdl, data = allele_tstdata))
          # if (perm.sumry.all_alleles$r.squared >= all_alleles_r2) count.moresig.all.alleles <- count.moresig.all.alleles + 1 
        }
      }
      vl.vs.resp[i, "perm p"] <- count.moresig/permtotal
      # vl.vs.resp[i, "perm p indiv alleles"] <- count.moresig.all.alleles/permtotal
      vl.vs.resp[i, "perm p alleles"] <- count.moresig.alleles/permtotal
      # browser()
    }
    ########################################################################################
    # plots
    ########################################################################################      
    if (vldotplot) {
      plotmain <- paste(race, group, "HLA-", locus, resp.var, "coeff =", signif(vl.vs.resp[i, "coeff"], 2), "p =", signif(vl.vs.resp[i, "p"], 2))
      pdf(paste(locus, race, group, "dotplot VL vs resp vl weights", assay, resp.var, ".pdf"), 12, 8)
      plot(vldata$allelescore, vldata$mVL, xlab = xlable, ylab = "mean logVL", main = plotmain, col = "black")
      abline(reg = lm.vl, col = "black")
      dev.off()
      if (print_chartdata) write.csv(vldata[,c("allelescore", "mVL")], file = paste(locus, race, group, "dotplot", assay, resp.var, ".csv"))
    }
    # browser()
    if (alleleplot) 
    {
      miny <- min(locus.resp$mean.VL, na.rm = T)
      maxy <- max(locus.resp$mean.VL, na.rm = T)
      # plotmain <- paste(race, group, "HLA-", locus, resp.var)
      locus.resp$plotNVL <- locus.resp$N.VL/scale.factor 
      locus.resp$largerVL <- ifelse(locus.resp$plotNVL > locus.resp$allele.carrier.ct, locus.resp$plotNVL, 0)
      locus.resp$largerJanct <- ifelse(locus.resp$plotNVL < locus.resp$allele.carrier.ct, locus.resp$allele.carrier.ct, 0)
      locus.resp$smaller <- pmin(locus.resp$plotNVL, locus.resp$allele.carrier.ct)
      plotmain <- paste(assay,race, group, "HLA-", locus, resp.var, "Rsquared =", signif(vl.vs.resp[i, "Rsquared.alleles"], 2), "p =", signif(summary(vlvresp.alleles)$coeff[2,4], 2))
      # browser()
      # three x bubbles to assign color for overlap
      pdf(paste(locus, race, group, "Allele bubble plots VL vs", resp.var, assay, filenamenote, ".pdf"), 12, 8)
      allele.bubble <- ggplot(locus.resp, aes(x=resp.score, y=mean.VL, size = largerVL)) +
        geom_point(color = "#66ccff") +
        geom_point(aes(x=resp.score, y=mean.VL, size = largerJanct), color = "#ff99ff") +
        geom_point(aes(x=resp.score, y=mean.VL, size = smaller), color = "#9977ff") +
        scale_size_continuous(range = c(1,20), breaks = c(5,10,15)) +
        geom_abline(slope = slope, intercept = intercept) +
        xlab(paste("vaccine response", resp.var)) +
        ylab("Cohort data, mean log viral load") +
        geom_text(aes(label = allele), size = 5, nudge_x = 0.002, nudge_y = 0.05) +
        ggtitle(plotmain)
      print(allele.bubble)
      dev.off()
      if (print_chartdata) write.csv(locus.resp[,c("resp.score", "mean.VL", "allele.carrier.ct")], file = paste(locus, race, group, "allele plot", assay, resp.var, ".csv"))
    }
    # browser()
    if(locus == "B") allele.counts <- locus.alleles.resp
    else allele.counts <- rbind(allele.counts, locus.alleles.resp)
    # permutation
    i <- i+1
  }
  vl.vs.resp$response <- resp.var
  vl.vs.resp$mincount <- mincount
  vvg <- list()
  # browser()
  vvg[[1]] <- vl.vs.resp
  vvg[[2]] <- allele.counts
  return(vvg)
}
