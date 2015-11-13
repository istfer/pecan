##-------------------------------------------------------------------------------
## Copyright (c) 2012 University of Illinois, NCSA.
## All rights reserved. This program and the accompanying materials
## are made available under the terms of the 
## University of Illinois/NCSA Open Source License
## which accompanies this distribution, and is available at
## http://opensource.ncsa.illinois.edu/license.html

##---------------------------------------------------------------------------------

##' Convert output for a single model run to NetCDF
##'
##' DEPRECATED this function will be removed in future versions, please update
##' your workflow.
##'
##' This function is a wrapper for model-specific conversion functions,
##' e.g. \code{model2netcdf.ED2}, \code{model2netcdf.BIOCRO}.
##' @title Convert model output to NetCDF 
##' @param runid 
##' @param outdir
##' @param model name of simulation model currently accepts ("ED2", "SIPNET", "BIOCRO")
##' @param lat Latitude of the site
##' @param lon Longitude of the site
##' @param start_date Start time of the simulation
##' @param end_date End time of the simulation
##' @return vector of filenames created, converts model output to netcdf as a side effect
##' @author Mike Dietze, David LeBauer
model2netcdfdep <- function(runid, outdir, model, lat, lon, start_date, end_date){
  ## load model-specific PEcAn module
  do.call(require, list(paste0("PEcAn.", model)))    


  model2nc <- paste("model2netcdf", model, sep=".")
  if(!exists(model2nc)){
    logger.warn("File conversion function model2netcdf does not exist for", model)
    return(NA)
  }
  
  do.call(model2nc, list(outdir, lat, lon, start_date, end_date))    
  
  print(paste("Output from run", runid, "has been converted to netCDF"))
  ncfiles <- list.files(path = outdir, pattern="\\.nc$", full.names=TRUE)
  if(length(ncfiles) == 0){
    logger.severe("Conversion of model files to netCDF unsuccessful")
  }
  return(ncfiles)
}

##' Convert output for a single model run to NetCDF
##'
##' DEPRECATED this function will be removed in future versions, please update
##' your workflow.
##'
##' This function is a wrapper for model-specific conversion functions,
##' e.g. \code{model2netcdf.ED2}, \code{model2netcdf.BIOCRO}.
##' @title Convert model output to NetCDF 
##' @param runid 
##' @param outdir
##' @param model name of simulation model currently accepts ("ED2", "SIPNET", "BIOCRO")
##' @param lat Latitude of the site
##' @param lon Longitude of the site
##' @param start_date Start time of the simulation
##' @param end_date End time of the simulation
##' @export
##' @return vector of filenames created, converts model output to netcdf as a side effect
##' @author Mike Dietze, David LeBauer
model2netcdf <- function(runid, outdir, model, lat, lon, start_date, end_date){
  logger.severe("model2netcdf will be removed in future versions, plase update your worklow")
}


##' Reads the output of a single model run
##'
##' Generic function to convert model output from model-specific format to 
##' a common PEcAn format. This function uses MsTMIP variables except that units of
##'  (kg m-2 d-1)  are converted to kg ha-1 y-1. Currently this function converts
##' Carbon fluxes: GPP, NPP, NEE, TotalResp, AutoResp, HeteroResp,
##' DOC_flux, Fire_flux, and Stem (Stem is specific to the BioCro model)
##' and Water fluxes: Evaporation (Evap), Transpiration(TVeg),
##' surface runoff (Qs), subsurface runoff (Qsb), and rainfall (Rainf).
##' For more details, see the
##' \href{http://nacp.ornl.gov/MsTMIP_variables.shtml}{MsTMIP variables}
##' documentation 
##' @title Read model output
##' @name read.output
##' @param runid the id distiguishing the model run. 
##' @param outdir the directory that the model's output was sent to
##' @param start.year first year of output to read (should be greater than ) 
##' @param end.year last year of output to read
##' @param variables variables to be read from model output
##' @return vector of output variable
##' @export
##' @author Michael Dietze, David LeBauer, Ryan Kelly
read.output <- function(runid, outdir, start.year=NA,
                        end.year=NA, variables = "GPP") {
  require(ncdf4)
  require(udunits2)

  ## vars in units s-1 to be converted to y-1
  cflux = c("GPP", "NPP", "NEE", "TotalResp", "AutoResp", "HeteroResp",
    "DOC_flux", "Fire_flux") # kgC m-2 d-1
  wflux = c("Evap", "TVeg", "Qs", "Qsb", "Rainf") # kgH20 m-2 d-1


  # Check if looking for an ED2 cohort variable. These must be handled differently from the usual vars stored in netcdf by model2netcdf.*.  
  ed2ind <- grep("ED2\\.", variables)
  if(length(ed2ind) > 0) {
      variables <- sub("ED2\\.", "", variables[ed2ind])
      

    # Extract the variable of interest, as well as NPLANT. Relies on a helper function meant for binning cohort level variables. The binning is superfluous here, as we're going to average over all sizes and PFTs. But that code saves us time here by collecting data from all output files, etc. It can also extract NPLANT in the process. 
    # Here, we aggregate into an arbitrarily large DBH bin (i.e., 0-Infinity cm), so the result will simply be the plant-density-weighted mean of the variable, by PFT. All that's left then is to average across PFTs, again weighting by NPLANT
    vars <- union("NPLANT", variables) # Extract NPLANT + variable of interest
    dbh.breaks <- 0 # Will represent a single DBH bin from 0 - Infinity cm
    pfts <- 1:20 # There are actually only 19 PFTs in ED currently, but this futureproofing won't hurt
    yrs <- start.year:end.year

    ed.dat <- ed.cohort.bin.PFT.DBH.YR(outdir, vars, dbh.breaks, pfts, yrs)
    
    # To conform with default usage (below), the result should be a list with one component for each variable. For each variable, results from all years are then abind()-ed together as they are read in. 
    result <- list()
    for(varname in variables) {
      # varname, yrs, pfts, and dbh.breaks should all be the same as they were for the inputs to ed.cohort.bin.PFT.DBH.YR. However, in case something went wrong (e.g. files didn't exist for some year?) check them all here.
      if(varname %in% names(ed.dat)) {
        dat = ed.dat[[varname]]
        yrs.dat = dimnames(dat)[[1]]
        dbh.breaks.dat = dimnames(dat)[[2]]
        pfts.dat = dimnames(dat)[[3]]
        
        for(i in 1:length(yrs.dat)) {
          result[[varname]] <- abind(result[[varname]], dat[yrs.dat[i],,])
        }
      } else {
        logger.warn(paste0("Variable ", varname, "not found in ED2 output."))
      }
    }

  } else {
  # Proceed as normal for a model2netcdf.* variable.
    # create list of *.nc years
    nc.years <- as.vector(unlist(strsplit(list.files(path = outdir, pattern="\\.nc$", full.names=FALSE),".nc")))
    # select only those *.nc years requested by user
    keep <- which(nc.years >= as.numeric(start.year) & nc.years <= as.numeric(end.year))
    ncfiles <- list.files(path = outdir, pattern="\\.nc$", full.names=TRUE)
    ncfiles <- ncfiles[keep]
    # throw error if no *.nc files selected/availible
    if(length(ncfiles) == 0) logger.error("no netCDF files of model output present")
  
    print(paste("Years: ",start.year," - ",end.year),sep="")
    result <- list()
    for(ncfile in ncfiles) {
      nc <- nc_open(ncfile)
      for(v in variables){
        if(v %in% c(names(nc$var),names(nc$dim))){
          newresult <- ncvar_get(nc, v)
          if(v %in% c(cflux, wflux)){
            newresult <- ud.convert(newresult, "kg m-2 s-1", "kg ha-1 yr-1")
          }
          result[[v]] <- abind(result[[v]], newresult)
        } else if (!(v %in% names(nc$var))){
          logger.warn(paste(v, "missing in", ncfile))
        }
      }
      nc_close(nc)
    }
  }  
  
  print(paste("----- Mean ", variables, " : ",
              lapply(result, mean, na.rm = TRUE)))
  print(paste("----- Median ", variables, ": ",
              lapply(result, median, na.rm = TRUE)))
  return(result)
}

##'--------------------------------------------------------------------------------------------------#
##' Converts the output of all model runs
##'
##' @title convert outputs from model specific code to 
##' @name convert.outputs
##' @param model name of simulation model currently accepts ("ED", "SIPNET", "BIOCRO")
##' @param settings settings loaded from pecan.xml
##' @param ... arguments passed to \code{\link{read.output}}, e.g. \code{variables}, \code{start.year}, \code{end.year}
##' @export
##' @author Rob Kooper
convert.outputs <- function(model, settings, ...) {
  logger.severe("This function is not longer used and will be removed in the future.")
}
####################################################################################################
