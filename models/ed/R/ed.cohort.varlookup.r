##' Lookup function for reading in/ processing commonly used ED variables
##'
##' @title lookup function
##' @param dat
##' @param varname
##' @export
##' @author Ryan Kelly

ed.var = function(dat, varname) {
  if(toupper(varname) == "AGB_CO") {
    out = ed.var.init(dat, "AGB_CO", type='co',
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname=NULL)                                                } else
  if(toupper(varname) == "DBH") {
    out = ed.var.init(dat, "DBH", type='co', 
                      from='plant', to='plant',
                      units="cm/plant",  
                      longname=NULL)                                                } else
  if(toupper(varname) == "BA_CO") {
    out = ed.var.init(dat, "BA_CO", type='co', 
                      from='plant', to='plant',
                      units="cm2/plant",  
                      longname=NULL)                                                } else
  if(toupper(varname) == "HITE") {
    out = ed.var.init(dat, "HITE", type='co', 
                      from='plant', to='plant',
                      units="m/plant",  
                      longname="Height")                                            } else
  if(toupper(varname) == "LAI_CO") {
    out = ed.var.init(dat, "LAI_CO", type='co', 
                      from='area', to='area',
                      units="m2/m2",  
                      longname="Leaf area index")                                   } else
  if(toupper(varname) == "DAGB_DT") {
    out = ed.var.init(dat, "DAGB_DT", type='co', 
                      from='plant', to='area',
                      units="kgC/m2/yr",  
                      longname=NULL)                                                } else
  if(toupper(varname) == "DDBH_DT") {
    out = ed.var.init(dat, "DDBH_DT", type='co', 
                      from='plant', to='plant',
                      units='cm/plant/yr',  
                      longname=NULL)                                                } else
  if(toupper(varname) == "DBA_DT") {
    out = ed.var.init(dat, "DBA_DT", type='co', 
                      from='plant', to='plant',
                      units=NULL,  
                      longname=NULL)                                                } else
  if(toupper(varname) == "MMEAN_GPP_CO") {
    out = ed.var.init(dat, "MMEAN_GPP_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Monthly mean GPP")                                  } else
  if(toupper(varname) == "MMEAN_NPP_CO") {
    out = ed.var.init(dat, "MMEAN_NPP_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Monthly mean NPP")                                  } else
  if(toupper(varname) == "CBR_BAR") {
    out = ed.var.init(dat, "CBR_BAR", type='co', 
                      from='plant', to='plant',
                      units="ratio/plant",
                      longname=NULL)                                                } else
  if(toupper(varname) == "NPLANT") {
    out = ed.var.init(dat, "NPLANT", type='co', 
                      from='area', to='area',
                      units='plant/m2',  
                      longname=NULL)                                                } else
  if(toupper(varname) == "DLNAGB_DT") {
    out = ed.var.init(dat, "DLNAGB_DT", type='co', 
                      from='plant', to='area',
                      units="log(kgC/m2/yr)",  
                      longname="Log AGB growth")                                    } else
  if(toupper(varname) == "DLNDBH_DT") {
    out = ed.var.init(dat, "DLNDBH_DT", type='co', 
                      from='plant', to='plant',
                      units="log(cm/plant/yr)",  
                      longname="Log DBH growth")                                    } else
  if(toupper(varname) == "DLNBA_DT") {
    out = ed.var.init(dat, "DLNBA_DT", type='co', 
                      from='plant', to='plant',
                      units="log(cm2/plant/yr)",  
                      longname="Log BA growth")                                     } else
  if(toupper(varname) == "BALIVE") {
    out = ed.var.init(dat, "BALIVE", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Live biomass")                                      } else
  if(toupper(varname) == "BDEAD") {
    out = ed.var.init(dat, "BDEAD", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",
                      longname="Dead biomass")                                      } else
  if(toupper(varname) == "BLEAF") {
    out = ed.var.init(dat, "BLEAF", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Leaf biomass")                                      } else
  if(toupper(varname) == "BROOT") {
    out = ed.var.init(dat, "BROOT", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Root biomass")                                      } else
  if(toupper(varname) == "BSAPWOODA") {
    out = ed.var.init(dat, "BSAPWOODA", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Aboveground sapwood biomass")                       } else
  if(toupper(varname) == "BSAPWOODB") {
    out = ed.var.init(dat, "BSAPWOODB", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Belowground sapwood biomass")                       } else
  if(toupper(varname) == "BSEEDS_CO") {
    out = ed.var.init(dat, "BSEEDS_CO", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Seed biomass")                                      } else
  if(toupper(varname) == "BSTORAGE") {
    out = ed.var.init(dat, "BSTORAGE", type='co', 
                      from='plant', to='area',
                      units="kgC/m2",  
                      longname="Storage biomass")                                   } else
  if(toupper(varname) == "MMEAN_TRANSP_CO") {
    out = ed.var.init(dat, "MMEAN_TRANSP_CO", type='co', 
                      from='area', to='area',
                      units=NULL,  
                      longname=NULL)                                                } else
  if(toupper(varname) == "MMEAN_GROWTH_RESP_CO") {
    out = ed.var.init(dat, "MMEAN_GROWTH_RESP_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Growth respiration")                                } else
  if(toupper(varname) == "MMEAN_LEAF_RESP_CO") {
    out = ed.var.init(dat, "MMEAN_LEAF_RESP_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Leaf respiration")                                  } else
  if(toupper(varname) == "MMEAN_ROOT_RESP_CO") {
    out = ed.var.init(dat, "MMEAN_ROOT_RESP_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Root respiration")                                  } else
  if(toupper(varname) == "MMEAN_STORAGE_RESP_CO") {
    out = ed.var.init(dat, "MMEAN_STORAGE_RESP_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Storage respiration")                               } else
  if(toupper(varname) == "MMEAN_CB_CO") {
    out = ed.var.init(dat, "MMEAN_CB_CO", type='co', 
                      from='plant', to='area',
                      units="kgC/m2/yr",  
                      longname="Carbon Balance")                                    } else
  if(toupper(varname) == "MMEAN_NPPCROOT_CO") {
    out = ed.var.init(dat, "MMEAN_NPPCROOT_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Coarse root NPP")                                   } else
  if(toupper(varname) == "MMEAN_NPPFROOT_CO") {
    out = ed.var.init(dat, "MMEAN_NPPFROOT_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Fine root NPP")                                     } else
  if(toupper(varname) == "MMEAN_NPPLEAF_CO") {
    out = ed.var.init(dat, "MMEAN_NPPLEAF_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Leaf NPP")                                          } else
  if(toupper(varname) == "MMEAN_NPPSAPWOOD_CO") {
    out = ed.var.init(dat, "MMEAN_NPPSAPWOOD_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Sapwood NPP")                                       } else
  if(toupper(varname) == "MMEAN_NPPWOOD_CO") {
    out = ed.var.init(dat, "MMEAN_NPPWOOD_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Heartwood NPP")                                     } else
  if(toupper(varname) == "MMEAN_NPPSEEDS_CO") {
    out = ed.var.init(dat, "MMEAN_NPPSEEDS_CO", type='co', 
                      from='plant', to='area',
                      units=NULL,  
                      longname="Seed NPP")                                          } else
  if(toupper(varname) == "MMEAN_MORT_RATE_CO") {
    out = ed.var.init(dat, "MMEAN_MORT_RATE_CO", type='co_mat',
                      from='plant', to='plant',
                      units='1/plant/yr (???)',  
                      longname="Mortality rate")                                    } else
  if(toupper(varname) == "PFT") {
    out = ed.var.init(dat, "PFT", type='pft',
                      units="#",  
                      longname="Cohort count by PFT")                               } else
  { # No Match!
    warning(paste0("Couldn't find varname ", varname, "!"))
    out = NULL
  }
  return(out)
}
