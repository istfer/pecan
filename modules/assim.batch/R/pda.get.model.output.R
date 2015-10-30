##' Get Model Output for PDA
##'
##' @title Get Model Output for PDA
##' @param all params are the identically named variables in pda.mcmc / pda.emulator
##'
##' @return A list containing model outputs extracted to correspond to each observational
##'         dataset being used for PDA. 
##'
##' @author Ryan Kelly
##' @export
pda.get.model.output <- function(settings, run.id, inputs) {

  model.out <- list()
  n.input <- length(inputs)
  for(k in 1:n.input){
    if(settings$assim.batch$inputs[[k]]$variable.id == 297) {
      ## NEE
  
      NEEm <- read.output(run.id, outdir = file.path(settings$run$host$outdir, run.id),
                          strftime(settings$run$start.date,"%Y"), 
                          strftime(settings$run$end.date,"%Y"), 
                          variables="NEE")$NEE*0.0002640674

      if(length(NEEm) == 0) {   # Probably indicates model failed entirely
        return(NA)
      }
      
      ## match model and observations
      NEEm <- rep(NEEm, each=nrow(inputs[[k]]$data)/length(NEEm))
      set <- 1:length(NEEm)  ## ***** need a more intellegent year matching!!!
        # NPPm <- rep(NPPm,each=length(NPPo)/length(NPPm))
        # set <- 1:length(NPPm) 

      model.out[[k]] <- NEEm[set]

    } else if (settings$assim.batch$inputs[[k]]$variable.id  %in% c(1000000008, 1000000009)) {
      # Get dbh breaks and pfts from the input data loaded already
      dbh.breaks <- inputs[[k]]$dbh.breaks
      pfts       <- inputs[[k]]$pfts
      
      ed.path     <- file.path(settings$modeloutdir, run.id);
        yrs <- strftime(settings$run$start.date,"%Y"):strftime(settings$run$end.date,"%Y")

      if(settings$assim.batch$inputs[[k]]$variable.id  == 1000000008) { 
        # GROWTH
        model.out[[k]] <- ed.cohort.bin.PFT.DBH.YR(ed.path, 'DDBH_DT', dbh.breaks, pfts, yrs)[[1]]
      } else if(settings$assim.batch$inputs[[k]]$variable.id  == 1000000009) { 
        # MORTALITY
        model.out[[k]] <- ed.cohort.bin.PFT.DBH.YR(ed.path, 'MMEAN_MORT_RATE_CO', dbh.breaks, pfts, yrs)[[1]]
      }
      
      
      # TEMPORARY (hopefully)... Eventually should have annual data (from Clark) to compare to. But since we don't now, just average over all years and return PFT x DBH class.
      model.out[[k]] <- apply(model.out[[k]], c(3,2), mean, na.rm=T)

    }
  }
  
  return(model.out)
}
