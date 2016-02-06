#--------------------------------------------------------------------------------------------------#
##' 
##' 
##' @name arate
##' @title arate 
##' @return 
##' @export
##' @author
##'

`arate` <-
function(x){1-(sum(diff(x)==0,na.rm=TRUE)/(length(x)-1))}

