### CODE for simple diagnostic plots for ED2.1 HDF output files
## This is the monthly version. Used to be called just "diagnostics.r" for some reason.

umol2kgC <- 1.20107e-8  ## umol(CO2) -> kg(C)

ed.read <- function(path,prefix,res.flag, vars=NULL, lib='ncdf4'){
  if(lib=='hdf5') {
    require(hdf5)
  } else if(lib=='rhdf5') {
    require(rhdf5)
  } else if(lib=='h5r') {
    require(h5r)
  } else if(lib=='ncdf4') {
    require(ncdf4)
  } else {
    stop("Unrecognized lib!")
  }

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

    if(lib=='hdf5') {
      for(i in 1:length(flist)){
        outi = hdf5load(flist[i],load=FALSE)
          if(!is.null(vars)) outi = outi[ names(outi) %in% c(vars,morevars) ]
        if(length(out$var) == 0){
          for(j in 1:length(outi)){
            out$var[[j]] = list()
            out$var[[j]][[1]] = outi[[j]]        
          }
          names(out$var) = names(outi)
        } else {
          t = length(out$var[[1]]) + 1
          for(j in 1:length(outi)){
            k = which(names(out$var) == names(outi)[j])
            if(length(k)>0){
              out$var[[k]][[t]] = outi[[j]]
            } else { ## add a new outiable. ***Not checked (shouldn't come up?)
              out$var[[length(out$var)+1]] = list()    # Add space for new outiable
              out$var[[length(out$var)]][1:(t-1)] = NA # Give NA for all previous time points
              out$var[[length(out$var)]][t] = outi[[j]] # Assign the value of the new outiable at this time point
              names(out$var)[length(out$var)] = names(outi)[j]
            }
          }      
        }
      }
    } else if(lib=='rhdf5') {
      for(i in 1:length(flist)){
        cat(flist[i], "...\n")
        ls.i = h5ls(flist[i])
          if(!is.null(vars)) ls.i = ls.i[ ls.i$name %in% c(vars,morevars), ]

        if(length(out$var) == 0){
          for(j in 1:nrow(ls.i)){
            out$var[[j]] = list()
            out$var[[j]][[1]] = h5read(flist[i], ls.i$name[j])
          }
          names(out$var) =  ls.i$name
        } else {
          t = length(out$var[[1]]) + 1
          for(j in 1:nrow(ls.i)){
            k = which(names(out$var) == ls.i$name[j])
            if(length(k)>0){
              out$var[[k]][[t]] = h5read(flist[i], ls.i$name[j])
            } else { ## add a new outiable. ***Not checked (shouldn't come up?)
              out$var[[length(out$var)+1]] = list()    # Add space for new outiable
              out$var[[length(out$var)]][1:(t-1)] = NA # Give NA for all previous time points
              out$var[[length(out$var)]][t] = h5read(flist[i], ls.i$name[j]) # Assign the value of the new outiable at this time point
              names(out$var)[length(out$var)] = ls.i$name[j]
            }
          }      
        }
      }
    } else if(lib=='h5r') {
      for(i in 1:length(flist)){
        cat(flist[i], "...\n")
        h5=H5File(flist[i])
        ls.i = ls(h5)
          if(!is.null(vars)) ls.i = ls.i[ ls.i %in% c(vars,morevars) ]

        if(length(out$var) == 0){
          for(j in 1:length(ls.i)){
            out$var[[j]] = list()
            out$var[[j]][[1]] = readH5Data(getH5Dataset(h5,ls.i[j]))
          }
          names(out$var) =  ls.i
        } else {
          t = length(out$var[[1]]) + 1
          for(j in 1:length(ls.i)){
            k = which(names(out$var) == ls.i[j])
            if(length(k)>0){
              out$var[[k]][[t]] = readH5Data(getH5Dataset(h5,ls.i[j]))
            } else { ## add a new outiable. ***Not checked (shouldn't come up?)
              out$var[[length(out$var)+1]] = list()    # Add space for new outiable
              out$var[[length(out$var)]][1:(t-1)] = NA # Give NA for all previous time points
              out$var[[length(out$var)]][t] = readH5Data(getH5Dataset(h5,ls.i[j]))# Assign the value of the new outiable at this time point
              names(out$var)[length(out$var)] = ls.i[j]
            }
          }      
        }
      }
    } else if(lib=='ncdf4') {
      for(i in 1:length(flist)){
        cat(flist[i], "...\n")
        nc = nc_open(flist[i])
        ls.i = names(nc$var)
          if(!is.null(vars)) ls.i = ls.i[ ls.i %in% c(vars,morevars) ]

        if(length(out$var) == 0){
          for(j in 1:length(ls.i)){
            out$var[[j]] = list()
            out$var[[j]][[1]] = ncvar_get(nc,ls.i[j])
          }
          names(out$var) =  ls.i
        } else {
          t = length(out$var[[1]]) + 1
          for(j in 1:length(ls.i)){
            k = which(names(out$var) == ls.i[j])
            if(length(k)>0){
              out$var[[k]][[t]] = ncvar_get(nc,ls.i[j])
            } else { ## add a new outiable. ***Not checked (shouldn't come up?)
              out$var[[length(out$var)+1]] = list()    # Add space for new outiable
              out$var[[length(out$var)]][1:(t-1)] = NA # Give NA for all previous time points
              out$var[[length(out$var)]][t] = ncvar_get(nc,ls.i[j])# Assign the value of the new outiable at this time point
              names(out$var)[length(out$var)] = ls.i[j]
            }
          }      
        }
      }
    }

    # Additional attributes. Make sure their names are lowercase to not interfere with ed.info below.
    # Assuming these variables are the same at all time points
      out$n.py  = as.integer(out$var$NPOLYGONS_GLOBAL[[1]])
      out$n.si  = as.integer(out$var$NSITES_GLOBAL[[1]])
      # Dim labels in metadata are wrong. Actually (ipoly,n_dbh,n_pft). To find the number of PFT, just get the 3rd dimension of any variable that has 3 dimensions. 
      out$n.pft =  na.omit(unlist(sapply(out$var, function(x) dim(x[[1]])[3])))[1]


    out$time = times
    out$n.t = length(times)
    out$info = ed.info(out)

    return(out)
  }
}


rhdf5load = function(file) {
  h5names = h5ls(file)$name
  n.name = length(h5names)
  out = list()
  for(i in 1:n.name) {
    out[[i]] = h5read(file, h5names[i])
  }
  names(out) = h5names
  return(out)
}
  

# An alternate using the 'rhdf5' package in place of 'hdf5' (which is not on CRAN). Much slower though
# ed.read <- function(path,prefix,res.flag){
#   require(rhdf5)
#   
#   flist = dir(path,prefix,full.names=TRUE)         ## grab all files matching prefix
#   flist = flist[grep(res.flag,flist)] ## select based on resolution flag
# 
#   if(length(flist)==0) {
#     stop("No files found!")
#   } else {
#     out = list()
#     for(i in 1:length(flist)){
# #       outi = hdf5load(flist[i],load=FALSE)
#       outi = rhdf5load(flist[i])
#       if(length(out) == 0){
#         for(j in 1:length(outi)){
#           out[[j]] = list()
#           out[[j]][[1]] = outi[[j]]        
#         }
#         names(out) = names(outi)
#       } else {
#         t = length(out[[1]]) + 1
#         for(j in 1:length(outi)){
#           k = which(names(out) == names(outi)[j])
#           if(length(k)>0){
#             out[[k]][[t]] = outi[[j]]
#           }
#         }      
#       }
#     }
# 
#     return(out)
#   }
# }


ed.info = function(dat) {
  # Strip out only the variables in ALL CAPS. ED stores all variables this way, but I can add whatever additional vars I want (including this 'info' table) as long as I keep their names lowercase.
#     dat = dat[ names(dat)==toupper(names(dat)) ]
  
  nvar = length(dat$var)
  info = data.frame(matrix(NA, nrow=nvar, ncol=8))
  names(info) = c("var","dimlabel","dim1","dim2","dim3","class","longname","units")
  info$var = names(dat$var)
  info$class = sapply(dat$var, function(x) class(x[[1]]))

  for(i in 1:nvar) {
    if(info$class[i] %in% c("integer", "numeric")) {
      info$dim1[i] = length(dat$var[[i]][[1]])
    } else if(info$class[i] == "matrix") {
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

