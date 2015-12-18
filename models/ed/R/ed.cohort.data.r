### CODE for simple diagnostic plots for ED2.1 HDF output files
## This is the monthly version. Used to be called just "diagnostics.r" for some reason.

umol2kgC <- 1.20107e-8  ## umol(CO2) -> kg(C)

ed.read <- function(path,prefix,res.flag, vars=NULL){
  require(ncdf4)

  # Some vars to always extract, because they're needed below. 
  morevars = c('NPOLYGONS_GLOBAL','NSITES_GLOBAL','AGB_PY', 'PYSI_ID', 'PYSI_N', 'SIPA_ID', 'SIPA_N', 'PACO_ID', 'PACO_N', 'AREA', 'AREA_SI')


  pattern = glob2rx(paste0(prefix,res.flag,'*h5'))
  flist = list.files(path=path, pattern=pattern, full.names=TRUE) 

  # Extract times from file names. 
    # Complicated gsub/strsplit pulls out year + month
      times = gsub(
                "(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\1\\2",
                sapply(
                  strsplit(flist, res.flag), 
                  function(x) x[2] # Select only the part of each name after res.flag
                )
              )
    # Paste paste "01" on each to indicate first day of month (needed by as.Date())
      times = paste0(times,"01")
    # Convert to Date format
      times = as.Date(times,"%Y%m%d")
      n.t   = length(times)

  # Get the data
  if(length(flist)==0) {
    stop("No files found!")
  } else {
    out = list()
    out$var = list()

    ncvar_getarray = function(nc, varname) {
      return(
        array(
          ncvar_get(nc, varname), 
          dim = nc$var[[varname]]$size
    ))}
    
  
    for(i in 1:length(flist)){
      cat(flist[i], "...\n")
      nc = nc_open(flist[i])
      ls.i = names(nc$var)
        if(!is.null(vars)) ls.i = ls.i[ ls.i %in% c(vars,morevars) ]

      if(length(out$var) == 0){
        for(j in 1:length(ls.i)){
          out$var[[j]] = list()
          out$var[[j]][[1]] = ncvar_getarray(nc,ls.i[j])
        }
        names(out$var) =  ls.i
      } else {
        t = length(out$var[[1]]) + 1
        for(j in 1:length(ls.i)){
          k = which(names(out$var) == ls.i[j])
          if(length(k)>0){
            out$var[[k]][[t]] = ncvar_getarray(nc,ls.i[j])
          } else { ## add a new outiable. ***Not checked (shouldn't come up?)
            out$var[[length(out$var)+1]] = list()    # Add space for new outiable
            out$var[[length(out$var)]][1:(t-1)] = NA # Give NA for all previous time points
            out$var[[length(out$var)]][t] = ncvar_getarray(nc,ls.i[j])# Assign the value of the new outiable at this time point
            names(out$var)[length(out$var)] = ls.i[j]
          }
        }      
      }
      
      nc_close(nc)
    }

    # Additional attributes. Make sure their names are lowercase to not interfere with ed.info below.
    # Assuming these variables are the same at all time points
      out$n.py  = as.integer(out$var$NPOLYGONS_GLOBAL[[1]])
      out$n.si  = as.integer(out$var$NSITES_GLOBAL[[1]])
      # Dim labels in metadata are wrong. Actually (ipoly,n_dbh,n_pft). To find the number of PFT, just get the 1st dimension of any variable that has 3 dimensions. 
      out$n.pft =  na.omit(unlist(sapply(out$var, function(x) dim(x[[1]])[1])))[1]


    out$time = times
    out$n.t = length(times)
    out$info = ed.info(out)

    return(out)
  }
}


ed.info = function(dat) {
  # Strip out only the variables in ALL CAPS. ED stores all variables this way, but I can add whatever additional vars I want (including this 'info' table) as long as I keep their names lowercase.
#     dat = dat[ names(dat)==toupper(names(dat)) ]
  
  nvar = length(dat$var)
  info = data.frame(matrix(NA, nrow=nvar, ncol=8))
  names(info) = c("var","dimlabel","dim1","dim2","dim3","class","longname","units")
  info$var = names(dat$var)
  info$class = sapply(dat$var, function(x) class(x[[1]]))

  for(i in 1:nvar) {
    if(length(dim(dat$var[[i]][[1]])) == 1) {
      info[i, c("dim1")] = dim(dat$var[[i]][[1]])
    } else if(length(dim(dat$var[[i]][[1]])) == 2) {
      info[i, c("dim1","dim2")] = dim(dat$var[[i]][[1]])
    } else { # array
      info[i, c("dim1","dim2","dim3")] = dim(dat$var[[i]][[1]])
    }
  #   info$n[i] = prod(na.omit(as.numeric(info[i,c("dim1","dim2","dim3")])))
    if(!is.null(attributes(dat$var[[i]][[1]])$Metadata)) {
      atts = attributes(dat$var[[i]][[1]])$Metadata
      info$longname[i] = gsub("^Long Name: |\\s+$", "", atts[grep("Long Name", atts)])
      info$units[i]    = gsub("^Units: |\\s+$", "", atts[grep("Units", atts)])
      info$dimlabel[i] = gsub("^Dimensions: |\\s+$", "", atts[grep("Dimensions", atts)])
    }
  }
  
  return(info)
}


ed.search = function(dat, pattern, exact=FALSE) {
  if(exact) {
    ind = unlist(sapply(pattern, function(x) which(dat$info$var==x)))
  } else {
    ind = unlist(sapply(pattern, function(x) grep(x,dat$info$var,ignore.case=T)))
  }
  if(length(ind)==0) warning("No matches!")
  return(dat$info[ind,])
}

