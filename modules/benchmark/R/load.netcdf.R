load.netcdf <- function(data.path, format, site, vars=NULL){
  
  # obvs.id <- "1000000384"
  # data.path <- "/fs/data1/pecan.data/input/Ameriflux_site_0-796/US-Bar.2006.nc"
  # format <- query.format.vars(obvs.id,con)  
  # site  <- query.site(settings$run$site$id, con)
  # vars=format$vars$bety_name[4]
  
  files = dir(dirname(data.path),basename(data.path),full.names = TRUE)
  i = 1
  require(ncdf4)
  
  nc <- nc_open(files[i])
  dims = names(nc$dim)
  time.var = grep(pattern = "time",dims,ignore.case = TRUE)
  time = ncvar_get(nc,dims[time.var])
  
  obs.vars <- list()
  for(varname in vars) {
    obs.vars[[varname]] <- ncvar_get(nc, varname)
  }
  
  obs.vars <- as.data.frame(do.call(cbind, obs.vars))
  
  nc_close(nc)
  return(obs.vars)

# load.netcdf <- function(path,vars,units,newunits){
#   
#   print(paste("Years: ",start.year," - ",end.year),sep="")
#   result <- list()
#   for(ncfile in ncfiles) {
#     nc <- nc_open(ncfile)
#     for(i in i:length(vars)){
#       v <- vars[i]
#       u1 <- units[i]
#       u2 <- newunits[i]
#       if(v %in% c(names(nc$var),names(nc$dim))){
#         newresult <- ncvar_get(nc, v)
#         newresult <- ud.convert(newresult, u1, u2)
#         result[[v]] <- abind(result[[v]], newresult)
#       } else if (!(v %in% names(nc$var))){
#         logger.warn(paste(v, "missing in", ncfile))
#       }
#     }
#     nc_close(nc)
#   }
#   
#   print(paste("----- Mean ", variables, " : ",
#               lapply(result, mean, na.rm = TRUE)))
#   print(paste("----- Median ", variables, ": ",
#               lapply(result, median, na.rm = TRUE)))
#   
#   
  
}