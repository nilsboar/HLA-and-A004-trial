allele_association <- function(fcn.clingrp, race, group, assay, min_allele_ct, outcomes, classes, pdfname, csvname, tablevars, nperm = 0, plot = T, boxplot = F, csvout = F, multibox = T, study = "Imbokodo")
{ 
  # browser()
  assoc <- NULL
  for (outcome in outcomes)
  {
    fcn.clingrp$depvar <- fcn.clingrp[,outcome]
    # browser()
    for (class in classes)
    {
      if (class == "Class I")
      {
        loci <- c("A", "B", "C")
      } else {
        loci <- c("DRB1", "DQA", "DQB", "DPA", "DPB")
      }
      locuscount <- 0
      allfrequent <- vector(mode = "character")
      allloci <- vector(mode = "character")
      alleleframe <- data.frame()
      firstallele <- T
      for (I in 1:length(loci))
      {
        # browser()
        locuscount <- locuscount + 1
        locus <- loci[I]
        name <- paste("HLA", locus)
        # print(c(I, locus, name))
        ##############################################################################
        # For APPROACH data
        ##############################################################################
        if (class == "Class I")
        {
          allele1 <- paste("C1", locus, "1", sep = "") # for Class I
          allele2 <- paste("C1", locus, "2", sep = "")
        } else {
          allele1 <- paste(locus, ".A1", sep = "") # for Class II
          allele2 <- paste(locus, ".A2", sep = "")
        }
        locustable <- table(c(fcn.clingrp[,allele1], fcn.clingrp[,allele2]))
        frequent <- names(which(locustable >= min_allele_ct))
        # allfrequent <- c(allfrequent, frequent)
        # allloci <- c(allloci, allele1, allele2)
        print(c(group, class, "freq allele ct", length(frequent), "min allele ct", min_allele_ct))
        locus.assoc <- as.data.frame(matrix(nrow = length(frequent), ncol = length(tablevars), dimnames = list(frequent, tablevars)))
        # browser()
        nsubs <- nrow(fcn.clingrp)
        locus.assoc$assay <- assay
        locus.assoc$outcome <- outcome
        locus.assoc$class <- class
        locus.assoc$race <- race
        locus.assoc$group <- group
        locus.assoc$'group count' <- nsubs # assumes no subsetting of nsubs
        # browser()
        print(outcome)
        print(frequent)
        for(allele in frequent)
        {
          # browser()
          allelename <- paste(locus, "*", allele, sep = "")
          fcn.clingrp$alleledom <- ifelse(fcn.clingrp[,allele1] == allele | fcn.clingrp[,allele2] == allele, 1, 0)
          sum(fcn.clingrp$alleledom)
          hasallele <- subset(fcn.clingrp, alleledom == 1)
          lacksallele <- subset(fcn.clingrp, alleledom == 0)
          # hla_assoc <- t.test(hasallele$depvar, lacksallele$depvar)
          # hla_wilcox <- wilcox.test(hasallele$depvar, lacksallele$depvar)
          # browser()
          if (group == "Blacks") hla_lm <- summary(lm(depvar ~ alleledom + gp140 + REGION1, data = fcn.clingrp))
          else if (group == "Whites") hla_lm <- summary(lm(depvar ~ alleledom + gp140 + hispanic, data = fcn.clingrp))
          else hla_lm <- summary(lm(depvar ~ alleledom + gp140, data = fcn.clingrp))
          locus.assoc[allele, "locus"] <- locus
          locus.assoc[allele, "allele"] <- allelename
          locus.assoc[allele, "count"] <- nrow(hasallele)
          locus.assoc[allele, "mean response"] <- mean(hasallele$depvar)
          locus.assoc[allele, "dom freq"] <- round(nrow(hasallele)/(nrow(hasallele) + nrow(lacksallele)), 2)
          locus.assoc[allele, "p"] <- signif(as.numeric(hla_lm$coefficients["alleledom","Pr(>|t|)"]), 2)
          locus.assoc[allele, "coeff"] <- signif(as.numeric(hla_lm$coefficients["alleledom","Estimate"]), 2)
          # set up frame for boxplot 
          haveallele <- subset(fcn.clingrp, fcn.clingrp[,allele1] == allele | fcn.clingrp[,allele2] == allele)[,c("HGAL", outcome)]
          names(haveallele)[2] <- allelename         

          if (firstallele)
          {
            allele.response <- haveallele 
            firstallele <- F 
          }
          else allele.response <- merge(allele.response, haveallele, by = 1, all = T)
          # print(c(locus, allele, locus.assoc[allele, "p"]))
          
        }
        locus.assoc <- locus.assoc[order(locus.assoc[,"p"], decreasing = F),]
        ################################################################3
        # permute above calculation
        ################################################################3
        if (nperm > 0)
        {
          # browser()
          elisrand <- fcn.clingrp
          locus.assoc[,"countmoresig"] <- 0
          locus.assoc[,"wilcoxmoresig"] <- 0
          locus.assoc[,"countmoresig.lm"] <- 0
          for (i in 1:nperm)
          {
            if(nperm == 0) break
            if (!(i%%100)) print(i)
            elisrand$rand1 <- fcn.clingrp[sample.int(nrow(fcn.clingrp)), "depvar"]
            for(allele in frequent)
            {
              # browser()
              elisrand$alleledom <- ifelse(elisrand[,allele1] == allele | elisrand[,allele2] == allele, 1, 0)
              hasallele <- subset(elisrand, alleledom == 1)
              lacksallele <- subset(elisrand, alleledom == 0)
              hla_assoc <- t.test(hasallele$rand1, lacksallele$rand1)
              hla_wilcox <- wilcox.test(hasallele$rand1, lacksallele$rand1)
              hla_lm <- summary(lm(rand1 ~ alleledom, data = elisrand))
              # browser()
              locus.assoc[allele, "countmoresig"] <- locus.assoc[allele, "countmoresig"] + as.numeric(hla_assoc$p.value < locus.assoc[allele, "p"])
              # locus.assoc[allele, "wilcoxmoresig"] <- locus.assoc[allele, "wilcoxmoresig"] + as.numeric(hla_wilcox$p.value < locus.assoc[allele, "p.wilcox"])
              locus.assoc[allele, "countmoresig.lm"] <- locus.assoc[allele, "countmoresig.lm"] + as.numeric(hla_lm$coefficients["alleledom", "Pr(>|t|)"] < locus.assoc[allele, "p"]) 
              # print(c(i, locus, allele, locus.assoc$p.value))
              # browser()
            }
          }
          locus.assoc[,"p.permute"] <- locus.assoc[, "countmoresig"]/nperm
          # locus.assoc[,"p.permute.wilcox"] <- locus.assoc[, "wilcoxmoresig"]/nperm
          locus.assoc[,"p.permute.lm"] <- (locus.assoc[, "countmoresig.lm"] + 1)/nperm
          # browser()
          for(allele in frequent) locus.assoc[allele,"p.permute.lm.bico"] <- binom.test(locus.assoc[allele, "countmoresig.lm"], nperm)$estimate
          if (min(locus.assoc["p.permute.lm.bico"] == 0)) browser()
          locus.assoc <- locus.assoc[order(locus.assoc[,"p.permute"], decreasing = F),]
        }
        locus.assoc$allele.label <- ""
        if (locuscount == 1) class.assoc <- locus.assoc
        else 
        {
          if(ncol(class.assoc) != ncol(locus.assoc))   browser()
          class.assoc <- rbind(class.assoc, locus.assoc)
        }
        # browser()
      }
      class.assoc <- class.assoc[order(class.assoc[,"p"], decreasing = F),]
      class.assoc$FDR <- signif(p.adjust(class.assoc$p, method = "BH"), 2)
      if (plot)
      {
        # outcome.name <- switch(outcome, Dim.1 = "PC1", Dim.2 = "PC2",  gag.mns.env = "gag_mns_env")
        outcome.name <- outcome
        cleangroup <- gsub("-", "", group, fixed = T)
        # browser()
        # groupname <- switch(cleangroup, Blacks = "All Blacks", Whites = "USA Whites", Asians = "Thailand Asians", EastAfrica = "East Africa Blacks",
                            # 'SouthAfrica' = "South Africa Blacks", USA = "African Americans")
        groupname <- group
        # plottitle <- paste("QQ plot", assay, group, outcome.name, unique(class.assoc$class), "> 2% allele frequency")
        # browser()
        plottitle <- paste(assay, unique(class.assoc$class), "allele association for", outcome.name, "for", groupname)
        print(plottitle)
        print(class.assoc[1:5,"p"])
        class.assoc$allele.label <- ifelse(class.assoc$p < 0.2, class.assoc$allele, "")
        topalleles <- subset(class.assoc, p < 0.2)[,"allele"]
        if (csvout)
        {
          if (nperm == 0) tablefilename <- paste(assay, outcome.name, class, group, csvname)
          else tablefilename <- paste(assay, outcome.name, class, group, "permute", nperm,  csvname)
          write.csv(class.assoc[,tablevars], file = tablefilename, row.names = F)
        }
        # browser()
        step = 1/(nrow(class.assoc))
        gap = step/2
        expected <- seq(gap, 1-gap, by = step)
        qqdata <- as.data.frame(cbind(class.assoc, expected))
        qqdata$logp <- -log10(qqdata$p)
        qqdata$logp.covars <- -log10(qqdata$p)
        # qqdata$logp.wilcox <- -log10(qqdata$p.wilcox)
        qqdata$logexpected <- -log10(qqdata$expected)
        # qqdata$correlation <- ifelse(sign(qqdata$coeff) > 0, "positive correlation", "negative correlation")
        qqdata$correlation <- ifelse(sign(qqdata$coeff) > 0, "positive", "negative")
        qqdata$ploty <- qqdata$logp.covars 
        namelength <- max(nchar(class.assoc$allele.label))
        maxval <- max(qqdata$logexpected, qqdata$ploty) + 0.55
        maxy <- maxval*(1.0 + 0.01*namelength)
        # nudge_adj <- max(1,maxval/4)
        yoffset <- 0.0 + maxval*(0.008*namelength + 0.03) 
        qqfilename <- paste(assay, outcome.name, class, group, pdfname)
        # qqfilename <- pdfname
        pdf(qqfilename, 6.5, 4.25)
        # browser()
        # ylims <- xlims <- switch(group, Blacks = c(0, 4), Whites = c(0, 2.5), Asians = c(0, 3.5))
        # ylims <- switch(group, Black = c(0, 9), White = c(0, 4.5), Asian = c(0, 3))
        qqplot <- ggplot (qqdata, aes(x = logexpected, y = ploty, color = correlation, size = count)) +
          geom_point() +
          xlim(0, maxy) +
          ylim(0, maxy) +
          scale_color_manual(values = c("positive" = "royalblue", "negative" = "red3" )) +
          # geom_text(aes(label = allele.label), nudge_y = 0.3*nudge_adj, angle = 90, size = 2.4) +
          geom_text(aes(label = allele.label), nudge_y = yoffset, angle = 90, size = 2.4) +
          # xlim(xlims) +
          # ylim(ylims) +
          geom_abline(slope = 1, intercept = 0) +
          xlab("-log10 p expected") +
          ylab("-log10 p observed") +
          ggtitle(plottitle)
        print(qqplot)
        dev.off()
        # call comp boxplot function
        response <- as.data.frame(pivot_longer(allele.response, cols = 2:ncol(allele.response), names_to = "allele", values_to = "response", values_drop_na = T))
        meanresponse <- mean(response$response, na.rm = T) # catches response var passed; other responses are in data frame
        response$response <- response$response - meanresponse
        fullresponse <- response
        if (boxplot)
        {
          response <- as.data.frame(pivot_longer(allele.response, cols = 2:ncol(allele.response), names_to = "allele", values_to = "response", values_drop_na = T))
          browser()
          meanresponse <- mean(response$response, na.rm = T) # catches response var passed; other responses are in data frame
          response$response <- response$response - meanresponse
          fullresponse <- response
          response <- subset(response, allele %in% topalleles)
          for (i in 1:nrow(response))
          {
            thisallele <- response$allele[i]
            response[i, "printnames"] <- paste(response$allele[i], " N = ", class.assoc[class.assoc$allele == thisallele, "count"])
          }
          # response$allele <- as.factor(response$allele)
          plotfilename <- paste(assay, outcome.name, class, group, "boxplot.pdf")
          plottitle <- paste("Distribution of ", outcome.name, " antigen response by allele",  sep = "")
          # KLUDGE: turn off print so I get fullresponse but not plot clutter
          pdf(plotfilename, 14, 8)
          # browser()
          responseplot <- ggplot(response, aes(x = reorder(printnames, response, FUN = median), y = response)) +
            geom_boxplot() +
          # theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.7, l = 0.5, unit = "in")) +
          # theme(axis.text.x = element_text(angle = 45, size=11, face = "bold", vjust = 1, hjust = 1)) # +
          # # xlab('Cohort') +
          theme(axis.text.x = element_text(angle = 45, size=9, face = "bold", vjust = 1, hjust = 1)) +
          xlab("Allele") +
          ylab(paste(outcome, "response"))
          # theme(axis.title = element_text(size=18)) +
          # geom_crossbar(aes(y = `logistic coeff` )) +
          # ggtitle(plottitle) +
          # theme(plot.title = element_text(size=22))
          print(responseplot)
          dev.off()
        }
      }
      rownames(class.assoc) <- class.assoc$allele
      if (class == classes[1]) outcome.assoc <- class.assoc
      else outcome.assoc <- rbind(outcome.assoc, class.assoc)
      # pc <- ifelse(outcome == "Dim.1", "PC1", "PC2")
      # if (class == "Class I")
      # {
      #   clsI.assoc <- assoc
      # } else {
      #   clsII.assoc <- assoc
      # }
      # rfilename <- paste(group, outcome, assay, "loop.Rdata")
      # save(clsI.assoc, clsII.assoc, )
      # browser()
    }
    # setting up "fullresponse" frame for multibox
    fullresponse$outcome <- outcome
    for (i in 1:nrow(fullresponse)) rownames(fullresponse)[i] <- paste(fullresponse[i, "HGAL"], fullresponse[i, "allele"], outcome, sep = "")
    if (outcome == outcomes[1])
    {
      assoc <- outcome.assoc
      allresponse <- fullresponse
    }
    else
    {
      # browser()
      assoc <- rbind(assoc, outcome.assoc)
      allresponse <- rbind(allresponse, fullresponse)
    }
  }
  return (assoc)
}
