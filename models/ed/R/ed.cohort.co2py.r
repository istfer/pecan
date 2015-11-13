# ----- ed.var.init ---------------------------------------------------------------------#
# Pull variable from edoutput data, handle as specified, output with metadata
# ---------------------------------------------------------------------------------------#
ed.var.init = function(dat, varname, type=NULL, units=NULL, longname=NULL, ...) {
  if(is.null(type)) {
    type = ed.vartype(dat, varname)
  }
  
  if(toupper(type)=='CO') {
    out = ed.co2py(dat,varname, ...)
  } else if(toupper(type)=='CO_MAT') {
    out = ed.comat2py(dat,varname, ...)
  } else if(toupper(type)=='PFT') {
    out = ed.pft2py(dat)
  }
  
  if(is.null(units)) units = dat$info$units[dat$info$var==varname]
  if(is.null(longname)) longname = dat$info$longname[dat$info$var==varname]
  
  # Set some attributes (trim possible whitespace and remove [] common in ED unit labels)
  trim = function (x) gsub("^\\s+|\\s+$", "", x)
  out$n.py     = dat$n.py
  out$n.si     = dat$n.si
  out$n.pft    = dat$n.pft
  out$time     = dat$time
  out$varname  = trim(varname)
  out$units    = trim(sub("\\[","", sub("\\]","", units)))
  out$longname = trim(longname)
  return(out)
}


# ----- ed.co2py ------------------------------------------------------------------------#
# Summarize a cohort-level variable to site + polygon level
#  Converts between area- and plant-level variables using cohort-level NPLANT variable
#  *** Site vs. polygon is basically untested since all tests had n.si = n.py = 1
# ---------------------------------------------------------------------------------------#
ed.co2py = function(dat, x, from, to=from) {
  # If var is actually the name of a variable in dat (usual case), extract it. This allows the option of passing in anything desired that has the right CO variable structure, such as the individual columns of the NCOHORT x 5 variable MMEAN_MORT_RATE_CO
  if(is.character(x)) x = dat$var[[x]]
  
  if(!all(from %in% c('area','plant') & to %in% c('area','plant'))) {
    stop("Both 'from' and 'to' must be either 'area' or 'plant'.")
  } else if(from=='plant') {
    x = ed.plant2area(x, dat$var$NPLANT)
  }
  
  out = ed.co2py.sumpft(dat,x)
  out$co = x
  out$pft = dat$var[["PFT"]]
  
  if(to=='plant') {
    # Get the site-level density of each PFT
    nplant = ed.co2py.sumpft(dat, dat$var$NPLANT)
    
    # Convert /m2 values to /plant
    out$py = out$py / nplant$py
    out$si = out$si / nplant$si
    out$co = ed.area2plant(out$co, dat$var$NPLANT)
  }
  
  return(out)
}


# ----- ed.comat2py ---------------------------------------------------------------------#
# A wrapper function that runs ed.co2py for each column in a cohort matrix var (like MORT)
# ---------------------------------------------------------------------------------------#
ed.comat2py = function(dat, varname, from, to=from) {
  x = dat$var[[varname]]

  n.col = ncol(x[[1]])
  
  x.py = array(0, dim=c(dat$n.t,dat$n.py,dat$n.pft,n.col))
  x.si = array(0, dim=c(dat$n.t,dat$n.si,dat$n.pft,n.col))
  x.co = list()
  
  for(m in 1:n.col) {
    x.m = lapply(x, function(z) z[,m])
    
    temp.m = ed.co2py(dat, x.m, from, to)
    x.py[,,,m] = temp.m$py
    x.si[,,,m] = temp.m$si
    x.co[[m]] = temp.m$co
  }
  
  return(list(py=x.py, si=x.si, co=x.co, pft=dat$var[["PFT"]]))
}


# ----- ed.co2py.sumpft -----------------------------------------------------------------#
# The function that does the work of scaling from cohort to site/polygon
#  POLY > SITE > PATCH > COHORT
#  Not tested NPOLY or NSITE > 1, but should be fine
# ---------------------------------------------------------------------------------------#
ed.co2py.sumpft = function(dat,x) {
  x.poly = array(0, dim=c(dat$n.t,dat$n.py,dat$n.pft))
  x.site = array(0, dim=c(dat$n.t,dat$n.si,dat$n.pft))
  
  i=j=k=l=1
  for(i in 1:dat$n.t){   # Month loop
    for(j in 1:dat$n.py) {    # Polygon loop
      # First and last index and count for sites within this polygon
      si0 = dat$var$PYSI_ID[[i]][j]
      nsi = dat$var$PYSI_N[[i]][j]
      siF = si0 + nsi - 1

      for(k in si0:siF){    # Loop over sites
        # First and last index and count for patches within this site
        pa0 = dat$var$SIPA_ID[[i]][k]
        npa = dat$var$SIPA_N[[i]][k]
        paF = pa0 + npa - 1

        for(l in pa0:paF){    # Loop over patches
          # First and last index and count for cohorts within this patch
          co0 = dat$var$PACO_ID[[i]][l]
          nco = dat$var$PACO_N[[i]][l]
          coF = co0 + nco - 1
          if(nco > 0){
            pft = dat$var$PFT[[i]][co0:coF]
            if(any(pft>dat$n.pft | pft<1)) warning("PFT numbers > n.pft and/or <1 !!!")

            # Weight the patch total by patch area...
            sum.pa = tapply(x[[i]][co0:coF],pft,sum,na.rm=TRUE) * dat$var$AREA[[i]][l]
            
            #...Add to running total (which at the end of this loop will be the site-level, patch-area-weighted mean). Uses kludgey trick to get column numbers from numeric PFT labels.
            x.site[i,k,as.numeric(names(sum.pa))] = x.site[i,k,as.numeric(names(sum.pa))] + sum.pa
          }
        }    # End patch loop
        
        # Weight site total by site area, and add to running polygon total
        x.poly[i,j,] = x.poly[i,j,] + x.site[i,k,]*dat$var$AREA_SI[[i]][k]

      }    # End site loop
    }    # End polygon loop
  }    # End time loop

  return(list(py=x.poly, si=x.site))
}


# ----- ed.pft2py -----------------------------------------------------------------------#
# Count up cohorts per pft at the polygon and site level
# ---------------------------------------------------------------------------------------#
ed.pft2py = function(dat) {
  x.poly = array(0, dim=c(dat$n.t,dat$n.py,dat$n.pft))
  x.site = array(0, dim=c(dat$n.t,dat$n.si,dat$n.pft))
  
  i=j=k=l=1
  for(i in 1:dat$n.t){   # Month loop
    for(j in 1:dat$n.py) {    # Polygon loop
      # First and last index and count for sites within this polygon
      si0 = dat$var$PYSI_ID[[i]][j]
      nsi = dat$var$PYSI_N[[i]][j]
      siF = si0 + nsi - 1

      for(k in si0:siF){    # Loop over sites
        # First and last index and count for patches within this site
        pa0 = dat$var$SIPA_ID[[i]][k]
        npa = dat$var$SIPA_N[[i]][k]
        paF = pa0 + npa - 1

        for(l in pa0:paF){    # Loop over patches
          # First and last index and count for cohorts within this patch
          co0 = dat$var$PACO_ID[[i]][l]
          nco = dat$var$PACO_N[[i]][l]
          coF = co0 + nco - 1
          if(nco > 0){
            pft = dat$var$PFT[[i]][co0:coF]
            if(any(pft>dat$n.pft | pft<1)) warning("PFT numbers > n.pft and/or <1 !!!")

            pft.count.l = hist(pft,1:(dat$n.pft+1),plot=F,right=F)$counts

            # Add to running total of pft counts.
            x.site[i,k, ] = x.site[i,k, ] + pft.count.l
            x.poly[i,j, ] = x.poly[i,j, ] + pft.count.l
          }
        }    # End patch loop
      }    # End site loop
    }    # End polygon loop
  }    # End time loop

  return(list(py=x.poly, si=x.site, co=dat$var[["PFT"]], pft=dat$var[["PFT"]]))
}


#----- ed.vartype -----------------------------------------------------------------------#
# Try to guess the variable type. Not really being used now, except as a convenience 
#  to quickly visualize the type of variable. The types PFT, CO, and CO_MAT are actually
#  used by ed.var.init at this point, though this fn isn't necessary since at this point 
#  I also specify the vartype explicitly in all cases.
# ---------------------------------------------------------------------------------------#
ed.vartype = function(dat, varnames) {
  n.var = length(varnames)
  vartypes = character(n.var)
  
  for(i in 1:n.var) {
    varname = varnames[i]
    q = dat$var[[varname]][[1]]
    
    # Special variables by Name
    if(varname=="PFT") { 
      vartypes[i] = "PFT"
    } else if(varname=="REPRO_PA") { 
      vartypes[i] = "PA_REPRO"
    } else if(varname=="MMEAN_RAD_PROFILE_CO") { 
      vartypes[i] = "CO_RAD"
      
    # Matrix variables
    } else if(is.matrix(q)) {
      if( nrow(q)==dat$var$NCOHORTS_GLOBAL[[1]] & ncol(q)==13 ) { 
        vartypes[i] = "CO_RUNNING_MMEAN"
      } else if( nrow(q)==dat$var$NPATCHES_GLOBAL[[1]] & ncol(q)==dat$var$NZG[[1]] ) { 
        vartypes[i] = "PA_SOIL"
      } else if( nrow(q)==dat$var$NPATCHES_GLOBAL[[1]] ) {
        vartypes[i] = "CO_MAT"
      }
    
    # Vector variables
    } else { 
      if(length(q) == dat$var$NCOHORTS_GLOBAL[[1]])   { 
        vartypes[i] = "CO"
      } else if(length(q) == dat$var$NPATCHES_GLOBAL[[1]])   { 
        vartypes[i] = "PA"
      } else if(length(q) == dat$var$NSITES_GLOBAL[[1]])     {
        vartypes[i] = "SI"
      } else if(length(q) == dat$var$NPOLYGONS_GLOBAL[[1]])  { 
        vartypes[i] = "PY"
      }
    }
  } # End for i in 1:n.var
  return(vartypes)
}


# ----- ed.cobin.calc -------------------------------------------------------------------#
# Bin a cohort-level variable
#  Background fn to do the binning on a cohort-level variable
# ---------------------------------------------------------------------------------------#
ed.cobin.calc = function(x, breaks=NULL, isfactor=FALSE) {
  # TEST: x=ed.var(dat,'DBH'); breaks=10
  # TEST: x=ed.var(dat,'PFT'); breaks=7:11; isfactor=TRUE
  if(length(breaks)==1) {
    breaks = seq(min(unlist(x$co)) - 0.001*diff(range(unlist(x$co))), 
                 max(unlist(x$co)) + 0.001*diff(range(unlist(x$co))), 
                 length=breaks+1)
  }

  out = list()
  if(isfactor) {
    if( !all(unlist(x$co) %in% breaks) )
      warning("Some data are not in the list of breaks!")
    out$cobin = lapply(x$co, function(x) match(x, breaks))
    out$nbin = length(breaks)
  } else {
    if(min(unlist(x$co))<min(breaks) | max(unlist(x$co))>max(breaks))
      warning("Some data exceed range of breaks!")
#     out$cobin = lapply(x$co, cut, breaks=breaks, labels=FALSE)
    out$cobin = lapply(x$co, findInterval, vec=breaks)
    out$nbin = length(breaks)
  }
  out$isfactor = isfactor
  out$breaks = breaks

  return(out)
}


ed.cobin = function(x, y1, y2=NULL, fn=mean,
                    breaks1=NULL, isfactor1=FALSE,
                    breaks2=NULL, isfactor2=FALSE) {
  # TEST: x=ed.var(dat,'AGB_CO'); y1=ed.var(dat,'DBH'); y2=ed.var(dat,'PFT'); fn=mean; breaks1=10; isfactor1=FALSE; breaks2=7:11; isfactor2=TRUE
  ind1 = ed.cobin.calc(y1, breaks1, isfactor1)
  n.t = length(x$time)
  if(is.null(y2)) {
    xbin = array(NA, dim=c(n.t, ind1$nbin))
    for(i in 1:n.t) 
      xbin[i,] = tapply(x$co[[i]], list(as.factor(ind1$cobin[[i]])), fn)
  } else {
    ind2 = ed.cobin.calc(y2, breaks2, isfactor2)

    xbin = array(NA, dim=c(n.t, ind1$nbin, ind2$nbin))
    for(i in 1:n.t) {
      tmp = tapply(x$co[[i]], list(as.factor(ind1$cobin[[i]]), as.factor(ind2$cobin[[i]])), fn)
      xbin[i,as.numeric(rownames(tmp)),as.numeric(colnames(tmp))] = tmp
    }
  }
  xbin[is.na(xbin)] = 0
  out = list(xbin=xbin, breaks1=ind1$breaks, breaks2=ind2$breaks)
  return(out)
}


# ----- ed.plant2area/area2plant --------------------------------------------------------#
# Convert 1/plant variable to 1/m2 variable by multiplying by nplant (plant/m2)
# ---------------------------------------------------------------------------------------#
ed.plant2area = function(x, nplant) {
  for(i in 1:length(x)) { x[[i]] = x[[i]] * nplant[[i]] }
  return(x)
}

ed.area2plant = function(x, nplant) {
  for(i in 1:length(x)) { x[[i]] = x[[i]] / nplant[[i]] }
  return(x)
}

# ----- JUNK ----------------------------------------------------------------------------#
# Works in progress not being used now
# ---------------------------------------------------------------------------------------#
# 
# ed.pa2py = function(dat, x) {
#   if(is.character(x)) x = dat$var[[x]]
# 
#   x.poly = array(0, dim=c(dat$n.t,dat$n.py))
#   x.site = array(0, dim=c(dat$n.t,dat$n.si))
#   
#   i=j=k=l=1
#   for(i in 1:dat$n.t){   # Month loop
#     for(j in 1:dat$n.py) {    # Polygon loop
#       # First and last index and count for sites within this polygon
#       si0 = dat$var$PYSI_ID[[i]][j]
#       nsi = dat$var$PYSI_N[[i]][j]
#       siF = si0 + nsi - 1
# 
#       for(k in si0:siF){    # Loop over sites
#         # First and last index and count for patches within this site
#         pa0 = dat$var$SIPA_ID[[i]][k]
#         npa = dat$var$SIPA_N[[i]][k]
#         paF = pa0 + npa - 1
# 
#         # Site value = area-weighted mean of the patches
#         x.site[i,k] = sum( x[[i]][pa0:paF]*dat$var$AREA[[i]][pa0:paF] )
# 
#         # Weight site total by site area, and add to running polygon total
#         x.poly[i,j] = x.poly[i,j] + x.site[i,k]*dat$var$AREA_SI[[i]][k]
#       }    # End site loop
#     }    # End polygon loop
#   }    # End time loop
# 
#   return(list(py=x.poly, si=x.site))
# }
# 
# ed.si2py = function(dat, x) {
#   if(is.character(x)) x = dat$var[[x]]
# 
#   x.poly = array(0, dim=c(dat$n.t,dat$n.py))
#   x.site = array(0, dim=c(dat$n.t,dat$n.si))
#   
#   i=j=k=l=1
#   for(i in 1:dat$n.t){   # Month loop
#     x.site[i,] = x[[i]]
#   
#     for(j in 1:dat$n.py) {    # Polygon loop
#       # First and last index and count for sites within this polygon
#       si0 = dat$var$PYSI_ID[[i]][j]
#       nsi = dat$var$PYSI_N[[i]][j]
#       siF = si0 + nsi - 1
# 
#       for(k in si0:siF){    # Loop over sites
#         # Weight site total by site area, and add to running polygon total
#         x.poly[i,j] = x.poly[i,j] + x.site[i,k]*dat$var$AREA_SI[[i]][k]
#       }    # End site loop
#     }    # End polygon loop
#   }    # End time loop
# 
#   return(list(py=x.poly, si=x.site))
# }
