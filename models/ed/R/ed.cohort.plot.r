


# ed.plot = function(dat,varname, ...) {
#   vartype = ed.vartype(dat,varname)
#   cat("Plotting ", varname, " (type = ", vartype, ")...\n", sep="")
#   if(vartype=="CO") {
#     ed.plot.co(dat,varname,...)
#   } else if(vartype=="SI") {
#     ed.plot.si(dat,varname,...)
#   } else if(vartype=="PFT") {
#     ed.plot.coXpft(dat,...)
#   } else if(vartype=="CO_MORT") {
#     ed.plot.mort(dat, ...)
#   } else {
#     warning("Plot function not found for ", varname, " (type = ", vartype, "). Plot skipped.", sep="")
#   }
# }

#----- COHORT-LEVEL VARIABLE AGGREGATED TO POLYGON --------------------------------------#
# For now, only plots polygon summary and assumes 1 polygon only
ed.plot.py = function(x, pft=NULL, cols=pftcols, plot.legend=TRUE, 
                      lty=1, lwd=2) {
  # TEST:  pft=NULL; cols=pftcols; plot.legend=TRUE; lty=1; lwd=2

  # For now, assume 1 py, so squash the polygon dim
    py = drop(x$py)

  if(is.null(pft)) pft = which(apply(py,2,sum) != 0)

  if(length(pft)>0) {
    plot(0, type='n', xlim=range(x$time), 
          ylim=range(py,na.rm=T)+c(0,0.1*abs(diff(range(py,na.rm=T)))), 
          xlab="time", xaxt='n', 
          ylab=paste0(x$longname, " [", x$units, "]"), main=x$varname)

    t.plot = pretty(x$time,5)
    axis(1, at=t.plot, labels=format(t.plot, "%m/%Y"), cex.axis=0.8)


    for(p in pft) {
      lines(x$time, py[,p], col=cols[p], lty=lty, lwd=lwd)
    }
    if(plot.legend) 
      legend("top", lwd=lwd, col=cols[pft],
             lty=lty, bty='n', horiz=T, legend=as.character(pft))
  }
}


#----- COHORT-LEVEL VARIABLE, NO AGGREGATION --------------------------------------------#
ed.plot.co = function(x, pft=NULL, cols=pftcols, 
                      lty=1, lwd=2, cex=0.7, alph=0.5) {
  # TEST:  x=ed.var(dat,'DBH'); pft=NULL; cols=pftcols; plot.legend=TRUE; lty=1; lwd=2

  if(is.null(pft)) pft = sort(unique(unlist(x$pft)))
    n.pft = length(pft)

  if(length(pft)>0) {
    par(mfrow=c(n.pft,1), mar=c(3,3,0,0), oma=c(1,1,1,1))
    n.t = length(x$time)
    t.plot = pretty(x$time,5)
#     t.plot = pretty(1:n.t,5)
    alph.h = as.hexmode(round(256*alph))
    t.vec = as.numeric(rep(x$time, sapply(x$co,length)))
    co.vec = unlist(x$co)
    pft.vec = unlist(x$pft)
    
    for(p in pft) {
      ind = which(pft.vec == p)
  
      plot(0, type='n', xlim=range(x$time), 
            ylim=range(co.vec[ind],na.rm=T)+c(0,0.1*abs(diff(range(co.vec[ind],na.rm=T)))), 
            xlab="", xaxt='n', 
            ylab="", main="")
        title(main=paste0(x$varname,", PFT ", p), line=-1, col.main=pftcols[p])

      points(jitter(t.vec[ind],1), co.vec[ind], cex=cex, 
        col=paste0(cols[p], alph.h))  # Add 2-digit alpha spec (range from 00 to FF)

#       boxplot(co.vec[ind] ~ t.vec[ind], add=F, col=cols[p], names=rep("",n.t),
#         ylim=range(co.vec[ind],na.rm=T)+c(0,0.1*abs(diff(range(co.vec[ind],na.rm=T)))),
#         cex.axis=0.7)
#       title(main=paste0(x$varname,", PFT ", p), line=-1, col.main=pftcols[p])

      if(p == pft[n.pft]) {
        axis(1, at=t.plot, labels=format(t.plot, "%m/%Y"), cex.axis=0.8)
        title(xlab="Time", cex.axis=0.8, line=2)
      } else {
        axis(1, labels=F)
      } 
      
      if(p == ceiling(n.pft/2)) {
        title(ylab=paste0(x$longname, " [", x$units, "]"), line=2)
      }
    } # for p in pft
  }
}


#----- COHORT-LEVEL CROSSPLOTS ----------------------------------------------------------#
# ROUGH DRAFT
ed.plot.co = function(x, y, pft=NULL, cols=pftcols, 
                      lty=1, lwd=2, cex=0.7, alph=0.5) {
  # TEST:  x=ed.var(dat,'HITE'); y=ed.var(dat,'CBR_BAR'); pft=NULL; cols=pftcols; plot.legend=TRUE; lty=1; lwd=2
  # TEST:  y=ed.var(dat,'MMEAN_MORT_RATE_CO'); y$co = y$co[[2]]

  if(!all(identical(x$time, y$time), 
          identical(sapply(x$co, length), sapply(y$co, length)),
          identical(x$pft, y$pft)))
    stop("x and y must have identical dimensions")

  if(is.null(pft)) pft = sort(unique(unlist(x$pft)))
    n.pft = length(pft)

  if(length(pft)>0) {
    plot(unlist(x$co), unlist(y$co), col=cols[unlist(x$pft)], xlab="", ylab="", main="")
    abline(h=0)
    title(xlab=paste0(x$longname, " [", x$units, "]"), line=2)
    title(ylab=paste0(y$longname, " [", y$units, "]"), line=2)
  }

}


# ----- COLOR PALETTES ------------------------------------------------------------------#
# Based on iWantHue
  # PFT
    pftcols = c("#475A71","#57B43B","#D37776","#CD50CC","#6E5824",
                "#D64A30","#57AF81","#416B35","#CD8E3A","#D24073",
                "#7D68D0","#63A7BA","#9AA53E","#CA86BD","#783D70",
                "#6E87C9","#8A3A31")

    #   plot(1:17, 1:17, pch=16, cex=5, col=pftcols)
    #   abline(h=1:17, col=pftcols, lwd=2)
    #   text(1:17, 1:17, as.character(1:17), col=grey(0.7))
    #   text(7:11, 7:11, as.character(7:11), col=1)

    #   temperate.Southern_Pine      = ED PFT 7
    #   temperate.Late_Conifer       = ED PFT 8
    #   temperate.Early_Hardwood     = ED PFT 9
    #   temperate.North_Mid_Hardwood = ED PFT 10
    #   temperate.Late_Hardwood      = ED PFT 11

  # Mortality
    mortcols = c("#9AA53E","#CA86BD","#783D70","#6E87C9","#8A3A31")
    #   plot(1:5, 1:5, pch=16, cex=5, col=mortcols)
    #   abline(h=1:5, col=mortcols, lwd=2)
    #   text(1:5, 1:5, as.character(1:5), col=1)


#----- PATCH-LEVEL VARIABLE -------------------------------------------------------------#
ed.plot.pa = function(...) {


}


#----- SITE-LEVEL VARIABLE --------------------------------------------------------------#
# As above, currently assumes just one site. Would need to figure out how multi-site vars get stored.
ed.plot.si = function(dat, varname, cols=2:7, plot.legend=TRUE, info=NULL, ...) {
  si = as.numeric(dat[[varname]])
  t = 1:length(si)  # Could specify these better, e.g. actual dates

  plot(t, si, type='l')
}


#----- MORTALITY ------------------------------------------------------------------------# 
ed.plot.mort = function(x, pft=NULL, 
                      cols=mortcols, lty=1, lwd=2) {
# cols=mortcols; lty=1; lwd=2
  par.save = par(no.readonly=T)

  # For now, assume 1 py, so squash the polygon dim
    py = drop(x$py)

  if(is.null(pft)) pft = which(apply(py,2,sum) > 0)

  n.pft = length(pft)
  
  if(n.pft>0) {
    par(mfrow=c(n.pft,1), mar=c(3,3,0,0), oma=c(1,1,1,1))
    t.plot = pretty(x$time,5)
    for(p in pft) {
      plot(0, type='n', xlim=range(x$time), ylim=c(0,1.2*max(py,na.rm=T)), xaxt='n', cex.axis=0.8)
      title(main=paste("Mortality, PFT", p), line=-1, col.main=pftcols[p])
      title(ylab="Mortality Rate", line=2)

      abline(h=1, col=grey(0.7))

      for(i in 5:1) {
        lines(x$time, py[,p,i], lwd=2, col=cols[i])
      }

      if(p==pft[1]) {
        legend("left", col=cols, lwd=lwd, cex=0.9, bty='n', horiz=F,             legend=c("age","cb","fall","cold","dist"))
      }

      if(p == pft[n.pft]) {
        axis(1, at=t.plot, labels=format(t.plot, "%m/%Y"), cex.axis=0.8)
        title(xlab="Time", cex.axis=0.8, line=2)
      } else {
        axis(1, labels=F)
      } 
    }
    par(par.save) # reset par
  } # end if(n.pft>0)
}


# Mort types -- correspond to the 5 columns of the mort variables 
#   ! 1. Ageing, PFT-dependent but otherwise constant;                                      !
#    ! 2. Negative carbon balance;                                                           !
#    ! 3. Treefall mortality;                                                                !
#    ! 4. Mortality due to cold weather.                                                     !
#    ! 5. Disturbance mortality.  This is not directly applied to the cohort population,     !
#    !    because this mortality is associated with the creation of a new patch, but it is   !
#    !    saved here for posterior analysis.                                                 !




# Plot number of cohorts by PFT
ed.plot.coXpft = function(dat, cum=FALSE, cols=2:7, plot.legend=TRUE,...) {
  pft.co = dat[["PFT"]]
  n.t = length(pft.co)
  t = 1:n.t

  pft = sort(unique(unlist(pft.co)))
    n.pft = length(pft)

  pft.sum = matrix(0, nrow=n.t, ncol=n.pft)
  for(i in 1:n.t) {
    for(j in 1:n.pft) {
      pft.sum[i,j] = sum(pft.co[[i]]==pft[j])
    }
  }
  
  if(cum) {
    pft.sum = t(apply(pft.sum, 1, cumsum))
    ylab = "# Cohorts (Cumulative)"
  } else {
    ylab = "# Cohorts"
  }
  
  plot(0, type='n', xlim=range(t), ylim=range(pft.sum)+c(0,0.1*abs(diff(range(pft.sum)))), 
        xlab="time", ylab=ylab, main="# Cohorts per PFT")
  for(j in 1:n.pft) {
    lines(t,pft.sum[,j], col=cols[(j %% length(cols)) + 1])
  }
  if(plot.legend) 
    legend("top", lwd=1, col=cols[((1:length(pft)) %% length(cols)) + 1],
           lty=1, bty='n', horiz=T, legend=as.character(pft))
}






