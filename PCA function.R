require(FactoMineR)
require(svd)

janssen_PCA <- function(data, timepoints, returnfirst = F, use.svd = F, use_meantime = F, baseline = F, whichresponse = "logresponse")
{
  timesdata <- subset(data, EPOCH %in% timepoints)
  if(use_meantime)
  {
    # hardwiring mean of 3 and 4 - 1
    delvector <- vector(mode = "numeric")
    log3daywide <- as.data.frame(pivot_wider(timesdata, names_from = c("ISTEST"), values_from = whichresponse, id_cols = c("PID", "EPOCH", "VISITDY"))) # VISITDY a kludge for multiple visits same epoch
    ags <- unique(data$ISTEST)
    print(c("ags available", unique(data$ISTEST), "data columns", names(log3daywide)))
    if (baseline)
    {
      if (!("TREATMENT1" %in% timepoints))
      {
        print("code assumes baseline is TREATMENT1, exiting")
        return(NULL)
      }
      baseline <- subset(log3daywide, EPOCH == "TREATMENT1") [,ags]
      responsedata <- subset(log3daywide, EPOCH != "TREATMENT1")
      subjects <- subset(log3daywide, EPOCH == "TREATMENT1")$PID 
      meanresponse <- matrix(nrow = length(subjects), ncol = length(ags), dimnames = list(subjects, ags))
      for (subject in subjects)
      {
        meanresponse[subject,] <- colMeans(subset(responsedata, PID == subject)[,ags], na.rm = T)
      }
      meanresponse <- meanresponse - as.matrix(baseline) # response has some different meanings here
      testmat0 <- meanresponse
    } else {
      responsedata <- log3daywide
      subjects <- unique(responsedata$PID)
      meanresponse <- matrix(nrow = length(subjects), ncol = length(ags), dimnames = list(subjects, ags))
      for (subject in subjects)
      {
        meanresponse[subject,] <- colMeans(subset(responsedata, PID == subject)[,ags], na.rm = T)
      }
      testmat0 <- meanresponse
    }
  } else  
    {
    logdaywide <-log3daywide
    visitdays <- as.character(names(table(timesdata$VISITDY)))
    testmat0 <- as.matrix(logdaywide[,-1])
    testmat0 <- testmat0 + 0.1 - min(testmat0, na.rm = T)
    testmat0 <- testmat0 + rnorm(length(testmat0), sd = 0.001)
    agdays <- names(logdaywide[,-1])
    testmat <-testmat0
    badcols <- vector(mode = "character")
    goodcols <- vector(mode = "character")
    variances <- vector(mode = "numeric")
    for (i in 1:length(agdays))
    {
      ag <- agdays[i]
      mean <- mean(as.vector(testmat0[,ag]), na.rm = T)
      variance <- var(as.vector(testmat0[,ag]), na.rm = T)
      variances[i] <- variance
      if (is.na(variance)) 
      {
        badcols <- c(badcols, ag)
        next
      }
      if (variance == 0) 
      {
        badcols <- c(badcols, ag) # but shouldn't be any
        next
      }
      for (j in 1:nrow(testmat0))
      {
         if(is.na(testmat0[j,ag])) testmat[j,ag] <- 0 
        else
        {
          testmat[j,ag] <- testmat0[j,ag] - mean
          if (is.na(testmat[j,ag])) browser()
        }
      }     
    }
    browser()
    testmat <- testmat[,!(colnames(testmat) %in% badcols)]
    str(testmat)
    if(use.svd)
    {
      timesvd <- svd(testmat)
      u <- timesvd$u
      v <- as.data.frame(t(timesvd$v))
      names(v) <- colnames(testmat)
      d <- timesvd$d
      hist(as.numeric(v[1,]), breaks = 1000)
      pca.time.coords <- v 
      barplot(d^2, main = "")
    } else
    {
      pca.time.out <- PCA(testmat)
      eigens <- pca.time.out$eig
      barplot(eigens[,2], las = 2, ylab = "% of total variance")
      pca.time.coeffs <- pca.time.out[["ind"]]$coord
      pca.time.coords <- t(pca.time.out[["var"]]$coord)
      
      hist(pca.time.coeffs, breaks = 100)
      hist(pca.time.coords, breaks = 100)
      hist(pca.time.coords[1,], breaks = 200)
    }  
    subs <- unique(timesdata$PID)
    ags <- unique(timesdata$ISTEST)
    timepc1 <- matrix(nrow = length(subs), ncol = length(ags), dimnames = list(subs,ags))
    for (i in 1:ncol(pca.time.coords))
    {
      split <- strsplit(colnames(pca.time.coords)[i], "_")
      timepc1[split[[1]][2],split[[1]][1]] <- pca.time.coords[1,i]
    }
    hist(timepc1, breaks = 1000)
    if (returnfirst) return(timepc1)
    testmat0 <- timepc1
    ags <- colnames(timepc1)
  }
  testmat <- testmat0
  for (ag in ags)
  {
    # browser()
    mean <- mean(as.vector(testmat0[,ag]), na.rm = T)
    variance <- var(as.vector(testmat0[,ag]), na.rm = T)
    print(c(mean, variance))
    for (j in 1:nrow(testmat0))
    {
      if(is.na(testmat0[j,ag])) testmat[j,ag] <- 0 
      else
      {
        testmat[j,ag] <- testmat0[j,ag] - mean
        testmat[j,ag] <- testmat[j,ag]/sqrt(variance)
      }
    }     
  }
  pca.timeag.out <- PCA(testmat)
  eigens <- pca.timeag.out$eig
  barplot(eigens[,2], las = 2, ylab = "% of total variance")
  pca.timeag.coeffs <- pca.timeag.out[["ind"]]$coord
  pca.timeag.coords <- t(pca.timeag.out[["var"]]$coord)
  hist(pca.timeag.coeffs[,1], breaks = 100)
  plot(pca.timeag.coeffs[,1], pca.timeag.coeffs[,2])
  if (use_meantime) return(list(pca.timeag.coeffs, meanresponse))
  return(pca.timeag.coeffs)
}












