ed.cohort.bin.PFT.DBH.YR <- function(ed.path, varnames, dbh.breaks, pfts, yrs) {
  require(hdf5)

  if(FALSE) varnames=c("DDBH", "MMEAN_MORT_RATE_CO", "NPLANT")

  # Hard code these for now
  ed.res.flag <- '-E-'
  ed.prefix   <- 'analysis'


  # List of vars to extract includes the requested one, plus others needed below. 
  vars <- c(varnames, 'PFT', 'DBH', 'NPLANT', 'AREA', 'PACO_N')


  # Get file names
  pattern <- glob2rx(paste0(ed.prefix,ed.res.flag,'*h5'))
  flist <- list.files(path=ed.path, pattern=pattern, full.names=TRUE) 


  # Extract times from file names. Only actually using years in this function, but the code was copied from another application where having the full time was useful. 
    # Complicated gsub/strsplit pulls out year + month
      times <- gsub(
                "(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\1\\2",
                sapply(
                  strsplit(flist, ed.res.flag), 
                  function(x) x[2] # Select only the part of each name after res.flag
                )
              )
    # Paste paste "01" on each to indicate first day of month (needed by as.Date())
      times <- paste0(times,"01")
    # Convert to Date format
      times <- as.Date(times,"%Y%m%d")
    # Get years only
      yrs.ed <- as.numeric(format(times, "%Y"))
      
  # Subset by years specified
    ind <- yrs.ed %in% yrs
    flist <- flist[ind]
    times <- times[ind]
    yrs.ed <- yrs.ed[ind]
      
    nt <- length(times)
    nyr <- length(years)


  # Get the model outputs
  if(length(flist)==0) {
    logger.severe("No outputs found!")
  } else {
    ed.dat <- list()

    for(i in 1:nt){
      outi <- hdf5load(flist[i],load=FALSE)
        if(!is.null(vars)) outi <- outi[ names(outi) %in% vars ]
      if(length(ed.dat) == 0){
        for(j in 1:length(outi)){
          ed.dat[[j]] <- list()
          ed.dat[[j]][[1]] <- outi[[j]]        
        }
        names(ed.dat) <- names(outi)
      } else {
        t <- length(ed.dat[[1]]) + 1
        for(j in 1:length(outi)){
          k <- which(names(ed.dat) == names(outi)[j])
          if(length(k)>0){
            ed.dat[[k]][[t]] <- outi[[j]]
          } else { ## add a new outiable. ***Not checked (shouldn't come up?)
            ed.dat[[length(ed.dat)+1]] <- list()    # Add space for new outiable
            ed.dat[[length(ed.dat)]][1:(t-1)] <- NA # Give NA for all previous time points
            ed.dat[[length(ed.dat)]][t] <- outi[[j]] # Assign the value of the new outiable at this time point
            names(ed.dat)[length(ed.dat)] <- names(outi)[j]
          }
        }      
      }
    }
  }


  ndbh <- length(dbh.breaks)
  npft <- length(pfts)

  out <- list()
  for(varname in varnames) {
    out[[varname]] <- array(NA, c(nt, ndbh, npft))
  }
  
  # Aggregate over PFT and DBH bins  
  i=j=1; k=2 
  for(i in 1:nt) {
    # Get additional cohort-level variables required
    pft        <- ed.dat$PFT[[i]]
    dbh        <- ed.dat$DBH[[i]]      # cm / plant
    plant.dens <- ed.dat$NPLANT[[i]]   # plant / m2

    # Get patch areas. In general patches aren't the same area, so this is needed to area-weight when averaging up to site level. Requires minor finnagling to convert patch-level AREA to a cohort-length variable. 
    patch.area <- ed.dat$AREA[[i]]    # m2  -- one entry per patch
    pacoN      <- ed.dat$PACO_N[[i]]  # number of cohorts per patch
    patch.area <- rep(patch.area, pacoN)  # patch areas, repped out to one entry per cohort

    # Now can get number of plants per cohort, which will be used for weighting. Note that area may have been (often/always is?) a proportion of total site area, rather than an absolute measure. In which case this nplant is a tiny and meaningless number in terms of actual number of plants. But that doesn't matter for weighting purposes. 
    nplant <- plant.dens * patch.area

    # Get index of DBH bin each cohort belongs to      
    dbh.bin <- sapply(dbh, function(x) which.min(x>c(dbh.breaks,Inf))) - 1


    # For each PFT x DBH bin, average variables of interest, weighting by # of plants. Not all ED cohort variables are in per-plant units. This code would not be applicable to them without modification.
    # However, it does handle two special cases. For NPLANT, it performs no weighting, but simply sums over cohorts in the PFT x DBH bin. For MMEAN_MORT_RATE_CO, it first sums over columns representing different mortality types first, then proceeds with weighting. 
    for(j in 1:ndbh) {
      for(k in 1:npft) {
        ind <- (dbh.bin==j) & (pft == pfts[k])
        
        if(any(ind)) {
          for(varname in varnames) {
            if(varname == "NPLANT") {
              # Return the total number of plants in the bin
              out$NPLANT[i,j,k] <- sum(nplant[ind])
            } else if(varname == "MMEAN_MORT_RATE_CO") {
              # Sum over all columns 
              mort = apply(ed.dat$MMEAN_MORT_RATE_CO[[i]][ind,, drop=F], 1, sum, na.rm=T)
              out$MMEAN_MORT_RATE_CO[i,j,k] <- sum(mort * nplant[ind]) / sum(nplant[ind])
            } else {
              # For all others, just get mean weighted by nplant
              out[[varname]][i,j,k] <- sum(ed.dat[[varname]][[i]][ind] * nplant[ind]) / sum(nplant[ind])
            }
          }
        }
      }
    }
  }

  
  # Aggregate over years. There are probably faster ways.
  for(varname in varnames) {
    tmp <-  array(NA, c(nyr, ndbh, npft))
    for(i in 1:nyr) {
      tmp[i,,] <- colMeans(out[[varname]][yrs.ed==yrs[i],,, drop=FALSE], na.rm=TRUE, dims=1)
    }
    out[[varname]] = tmp
    dimnames(out[[varname]]) <- list(yrs=yrs, dbhLowBound=dbh.breaks, pft=pfts)
  }


  return(out)
}
