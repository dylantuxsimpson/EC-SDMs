plotIndex <- function(dat){
  if(is.numeric(dat$site) == F){dat$site <- as.numeric(dat$site)}
  for(i in sort(unique(dat$site))){
    dat$plotID[dat$site==i] <- as.numeric(as.factor(as.character(dat$plot[dat$site==i])))
  }
  return(dat$plotID)
}


obsIndex <- function(dat){
  if(is.numeric(dat$plotID) == F){stop("Plot must be indexed numerically in column 'plotID,' run plotIndex() first.")}
  if(is.numeric(dat$site) == F){dat$siteID <- as.numeric(dat$site)}
  for(i in sort(unique(dat$siteID))){
    for(j in sort(unique(dat$plotID[dat$siteID==i]))){
      for(t in unique(dat$year[dat$siteID==i & dat$plotID==j]))
        dat$obsID[dat$siteID==i & dat$plotID==j & dat$year==t] <- 1:length(unique(dat$obs[dat$siteID==i & dat$plotID==j & dat$year==t]))
    }
  }
  return(dat$obsID)
}

##This only works for the specific case of the prevalence model, in which only 2017 has multiple occasions
yearIndex <- function(dat){
  if(is.numeric(dat$plotID) == F){stop("Plot must be indexed numerically in column 'plotID,' run plotIndex() first.")}
  if(is.numeric(dat$obsID) == F){stop("Observations must be indexed numerically in column 'obsID,' run obsIndex() first.")}
  if(is.numeric(dat$site) == F){dat$siteID <- as.numeric(dat$site)}
  for(i in sort(unique(dat$siteID))){
    for(j in sort(unique(dat$plotID[dat$siteID==i]))){
      dat$yearID[dat$siteID==i & dat$plotID==j & dat$obsID==1] <- 1:length(dat$year[dat$siteID==i & dat$plotID==j & dat$obsID==1])
      if(length(unique(dat$obsID[dat$siteID==i & dat$plotID==j])) > 1){
        for(k in dat$obsID[dat$siteID==i & dat$plotID==j & dat$obsID>1]){
          # This bit is confusing; for each row in which obsID > 1, I take the sum of the product of the vector of year ID's of site i, plot j and a logical query that matches
          # the original, unique year ID with that of the vector obsID=1. The vector is all 0's except for the year ID that matches the original year value, so the sum is the yearID for that combination.
          dat$yearID[dat$siteID==i & dat$plotID==j & dat$obsID==k] <- sum(dat$yearID[dat$siteID==i & dat$plotID==j & dat$obsID==1] * as.numeric(dat$year[dat$siteID==i & dat$plotID==j & dat$obsID==1] == dat$year[dat$siteID==i & dat$plotID==j & dat$obsID==k]))
        }
      }
    }
  }
  return(dat$yearID)
}


datLength <- function(dat){
  if(is.numeric(dat$plotID) == F){stop("Plot must be indexed numerically in column 'plotID'")}
  if(is.numeric(dat$siteID) == F){dat$siteID <- as.numeric(dat$site)}
  nsite <- length(unique(dat$site))
  
  nyear.max <- length(unique(dat$year))
  
  nplot <- vector()
  for(i in unique(dat$siteID)){
    nplot[i] <- length(unique(dat$plotID[dat$siteID==i]))
  }
  
  # This holds the number of years each plot is surveyed.
  # This method requires years be indexed individually for each plot so that, 
  # within each plot, years count from 1:N_year,
  # which I've done here simply by how I subsetted my dataframe.
  nyear <- array(NA, dim = c(nsite, max(nplot)))
  for(i in sort(unique(dat$siteID))){
    for(j in sort(unique(dat$plotID[dat$siteID==i]))){
      nyear[i,j] <- length(unique(dat$year[dat$siteID==i & dat$plotID==j]))
    }
  }
  
  # This is an array of the number of observations during year t at plot i of site k:
  nobs <- array(NA, dim = c(nsite, max(nplot), nyear.max)) 
  for(i in sort(unique(dat$siteID))) {
    for(j in sort(unique(dat$plotID[dat$siteID==i]))) {
      for(t in sort(unique(dat$yearID[dat$siteID==i & dat$plotID==j]))){
        nobs[i,j,t] <- length(unique(dat$obsID[dat$siteID==i & dat$plotID==j & dat$yearID==t]))
      }
    }
  }
  
  return(list(
    nsite = nsite,
    nplot = nplot,
    nyear = nyear,
    nobs = nobs
  ))
}


postPlot.batch <- function(dat, mod = NULL){
  posteriorPlot <- function(dat, index, mod = mod){
    tryCatch(plotPost(dat[,index],
                      xlab = paste(mod, names(dat)[index], sep = " ")), error = function(e)
                        hist(dat[,index],
                             xlab = names(dat)[index])
    )
  }
  lapply(1:ncol(dat), function(x) posteriorPlot(dat = dat, index = x, mod = mod))
}

postPlot.batchPrint <- function(dat, mod = NULL){
  posteriorPrint <- function(dat, index, mod = mod){
    jpeg(paste0(mod, "_", names(dat)[index],".jpeg"))
    tryCatch(plotPost(dat[,index],
                      xlab = paste0(mod, "_", names(dat)[index])), error = function(e)
                        hist(dat[,index],
                             xlab = names(dat)[index])
    )
    dev.off()
  }
  lapply(1:ncol(dat), function(x) posteriorPrint(dat = dat, index = x, mod = mod))
}
