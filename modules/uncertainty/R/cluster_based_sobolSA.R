##' Performs cluster-based GSA based on Roux & Buis et al 2021
##' function uses code from https://github.com/sbuis/ClusterBased_GSA/blob/master/clusterBasedGSA_on_Toycurve.Rmd
##'
##' @name cluster_based_sobolSA
##' @title Reads an ensemble time-series from PEcAn for the selected target variable
##' @param settings
##' @param settings
##' @return nothing, creates cluster-based GSA plots 
##'
##' @export
##'
##' 
cluster_based_sobolSA <- function(settings, ensemble.id = NULL, variable = "NPP", start.year = NULL, end.year = NULL){
  
  PEcAn.logger::logger.info("Performing the cluster-based GSA on ", variable, "- Please wait.")
  
  if (is.null(ensemble.id)) {
    ensemble.id <- settings$ensemble$ensemble.id
  }
  if (is.null(start.year)) {
    start.year <- settings$ensemble$start.year
  }
  if (is.null(end.year)) {
    end.year <- settings$ensemble$end.year
  }
  
  fname <- ensemble.filename(settings, "ensemble.ts", "Rdata",
                             all.var.yr = FALSE, 
                             ensemble.id = ensemble.id, 
                             variable = variable, 
                             start.year = start.year, 
                             end.year = end.year)
  
  if(!(file.exists(fname))){
    ensemble.ts <- read.ensemble.ts(settings, variable = variable)
  }else{
    load(fname)
  }
  curves <- ensemble.ts[[variable]]
    
  
  ### Clustering
  
  # the number of clusters, should be passed via settings
  if (is.null(settings$ensemble$sobolSA$nbClust)) {
    nbClust <- 3 
    PEcAn.logger::logger.info("The number of clusters in cluster-based GSA is not given. Defaulting to 3.")
  }else{
    nbClust <- as.numeric(settings$ensemble$sobolSA$nbClust)
  }
  
  tp <- ncol(curves)
  # the number of points on which the curves are discretized
  if (is.null(settings$ensemble$sobolSA$npdisc)) {
    np <- tp
    PEcAn.logger::logger.info("The number of points on which the curves are discretized is not given. Defaulting to all.")
    # or if we want something smaller by default:
    # div <- seq_len(abs(tp))
    # factors <- div[tp %% div == 0L]
    # np <- sort(factors)[ceiling(length(factors) / 1.5)] 
  }else{
    np <- as.numeric(settings$ensemble$sobolSA$npdisc)
  }
  
  
  PEcAn.logger::logger.info("Detecting", nbClust, "clusters in the simulated outputs over", np, "discretized points. Please wait.")
  
  subcurves <- curves[,seq(1,tp, length.out = np)] # can take sum / mean over?
  m <- 2 # fuzziness parameter
  clust <- fclust::FKM(subcurves,k=nbClust,m=m,conv=1e-3, maxit=50)
  u <- clust$U # Membership functions
  

  centers <- data.frame(matrix(nrow=nbClust,ncol=np)) # clusters' centers
  for (k in 1:nbClust) centers[k,] <-  apply(subcurves*(u[,k]^m),2,sum) / sum(u[,k]^m)
  
  pdfname <- paste0(sub("ts.*", "", fname), "cluster-centers.", ensemble.id, ".", variable, ".", start.year, ".", end.year ,".pdf")
  rdtname <- paste0(sub("ts.*", "", fname), "cluster-centers.", ensemble.id, ".", variable, ".", start.year, ".", end.year ,".Rdata")
  
  pdf(pdfname, width = 12, height = 9)

  # plot (a sub-set of) curves for each cluster
  # don't plot all. <100 is fine
  nc <- ifelse(nrow(curves) > 100, 100, nrow(curves)) # number of curves to plot 
  nci <- sample(1:nrow(curves), nc) 
  gr <- 0.0
  col_base <- rep(gray(gr),nbClust)
  for (cl in 1:nbClust){
    plot(seq(0,1,length.out=np),subcurves[1,],typ='n',xlab='',main=sprintf("Cluster %d",cl), ylim=range(curves))
    for (i in nci) lines(seq(0,1,length.out=np),subcurves[i,],lwd=2,col= adjustcolor(col_base[cl], alpha.f = u[i,cl]^2))
    lines(seq(0,1,length.out=np),centers[cl,],lty=1,lwd=4,ylim=c(0,1),col='red')
    grid(col='black')
    mtext(paste0(nc, " sample curves plotted."))
  }
  
  dev.off()
  
  cluster <- list()
  cluster$members <- subcurves
  cluster$centers <- centers
  save(cluster, file = rdtname)
  
  
  
  #### Compute sensitivity indices
  
  pdfname <- paste0(dirname(fname), "/cluster-basedGSA.", ensemble.id, ".", variable, ".", start.year, ".", end.year, ".pdf")
  rdtname <- paste0(dirname(fname), "/cluster-basedGSA.", ensemble.id, ".", variable, ".", start.year, ".", end.year, ".Rdata")
  
  pdf(pdfname, width = 12, height = 9)
  
  load(paste0(sub("ts.*", "", fname), "samples.", ensemble.id, ".Rdata"))
  
  # Sensitivity indices on membership functions
  Clust_SI <- vector("list",nbClust)
  paramNames <- unlist(sapply(which(!names(ens.samples) %in% c("env", "sobolSA")), function(x) names(ens.samples[[x]])))
  nbParam <- length(paramNames)
  colnames(sobolSA$X1) <- colnames(sobolSA$X2) <- colnames(sobolSA$X) <- paramNames
  for (cl in 1:nbClust) {
    Clust_SI[[cl]] <- sensitivity::tell(sobolSA, y=u[,cl])
    plot_indices(Si1=Clust_SI[[cl]]$S[,1], ST=Clust_SI[[cl]]$T[,1], 
                  lowCI_Si1=Clust_SI[[cl]]$S[,4], upCI_Si1=Clust_SI[[cl]]$S[,5], 
                  lowCI_ST=Clust_SI[[cl]]$T[,4], upCI_ST=Clust_SI[[cl]]$T[,5], 
                  graph_title=paste0("Parameters leading to cluster",cl), paramNames, nbParam)
  }
  
  # Sensitivity  indices on the difference between two membership functions 
  comb <- combn(nbClust,2) # computes the different combinations of 2 clusters among nbClust
  nbComb <- ncol(comb)
  for (icomb in 1:nbComb) {
    dClust_SI <- sensitivity::tell(sobolSA, y=u[,comb[1,icomb]]-u[,comb[2,icomb]])
    plot_indices(Si1=dClust_SI$S[,1], ST=dClust_SI$T[,1],
                 lowCI_Si1=dClust_SI$S[,4], upCI_Si1=dClust_SI$S[,5],
                 lowCI_ST=dClust_SI$T[,4], upCI_ST=dClust_SI$T[,5],
                 graph_title=paste0("Parameters influencing changes along a direction between ",comb[1,icomb],"-",comb[2,icomb]),
                 paramNames, nbParam)
  }
  
  
  # Compute Clust-GSI
  clust_GSI <- compute_GSI(sapply(Clust_SI, function(x) x$V[1,1]), 
                           sapply(Clust_SI, function(x) x$S[,1]), 
                           sapply(Clust_SI, function(x) x$T[,1]), paramNames, nbParam)
  # Compute confidence intervals on Clust-GSI
  nboot=100 # could be passed via fcn args
  clust_GSI_boot <- vector("list",nboot)
  X1 <- sobolSA$X1
  X2 <- sobolSA$X2
  n <- nrow(curves) / (nbParam + 2)  # nrow(curves) = (n * (p+2))
  for (iboot in 1:nboot){
    bt_idx <- sample(n,size=n,replace=TRUE) # resample in X1 and X2 matrices
    
    gsa_boot <- sensitivity::soboljansen(model = NULL, X1= X1[bt_idx,], X2 = X2[bt_idx,], nboot=0)
    idx_in_orig_DoE <- match(data.frame(t(gsa_boot$X)), data.frame(t(sobolSA$X))) # identify the lines of the bootstrapped DoE, gsa_boot$X, in the original DoE, gsa$X, to reuse the simulated curves
    
    Clust_SI_boot <- vector("list",nbClust)
    for (cl in 1:nbClust){
      Clust_SI_boot[[cl]] <- sensitivity::tell(sobolSA, y=u[idx_in_orig_DoE,cl])
    }
    clust_GSI_boot[[iboot]] <- compute_GSI(sapply(Clust_SI_boot, function(x) x$V[1,1]), 
                                           sapply(Clust_SI_boot, function(x) x$S[,1]), 
                                           sapply(Clust_SI_boot, function(x) x$T[,1]), paramNames, nbParam)
    
  }
  # Quantiles computation
  clust_GSI$Si1_CI95pcMin <- setNames(rep(NA,nbParam), paramNames); clust_GSI$Si1_CI95pcMax <- setNames(rep(NA,nbParam), paramNames)
  clust_GSI$ST_CI95pcMin <- setNames(rep(NA,nbParam), paramNames); clust_GSI$ST_CI95pcMax <- setNames(rep(NA,nbParam), paramNames)
  for (i in 1:nbParam){
    clust_GSI$Si1_CI95pcMin[i] = quantile(sapply(clust_GSI_boot,`[[`,"Si1")[i,],0.025)
    clust_GSI$Si1_CI95pcMax[i] = quantile(sapply(clust_GSI_boot,`[[`,"Si1")[i,],0.975)
    clust_GSI$ST_CI95pcMin[i] = quantile(sapply(clust_GSI_boot,`[[`,"ST")[i,],0.025)
    clust_GSI$ST_CI95pcMax[i] = quantile(sapply(clust_GSI_boot,`[[`,"ST")[i,],0.975)
    
  }
  

  
  plot_indices(Si1=clust_GSI$Si1, ST=clust_GSI$ST, 
               lowCI_Si1=clust_GSI$Si1_CI95pcMin, upCI_Si1=clust_GSI$Si1_CI95pcMax, 
               lowCI_ST=clust_GSI$ST_CI95pcMin, upCI_ST=clust_GSI$ST_CI95pcMax, 
               graph_title=paste0("Cluster-based  GSI for ", variable), paramNames, nbParam)
  
  dev.off()
  
  save(clust_GSI, dClust_SI, file = rdtname)
  
} 


plot_indices <- function(Si1, ST, lowCI_Si1, upCI_Si1, lowCI_ST, upCI_ST, graph_title="", paramNames, nbParam) {
  # Si1: vector of main Sobol' indices for each parameter
  # ST: vector of total Sobol' indices for each parameter
  # lowCI_Si1: vector of confidence interval lower bounds for main Sobol' indices for each parameter
  # upCI_Si1: vector of confidence interval upper bounds for main Sobol' indices for each parameter
  # lowCI_ST: vector of confidence interval lower bounds for total Sobol' indices for each parameter
  # upCI_ST: vector of confidence interval upper bounds for total Sobol' indices for each parameter
  btmmar <- max(nchar(paramNames))
  par(mar=c(round(btmmar/2),4,2,2)) # adjust as needed
  
  b <-  barplot(rbind(ST,Si1), 
                beside=TRUE, ylim=c(0,1), col=c(gray(0.5),'white'), 
                names.arg=paramNames, main=graph_title, 
                legend.text=c("TSI",expression(SI[1])), las=2)
  
  # plot confidence intervals
  for (i in 1:nbParam){
    segments(b[2*i-1],lowCI_ST[i],b[2*i-1],upCI_ST[i])
    segments(b[2*i],lowCI_Si1[i],b[2*i],upCI_Si1[i])
  }
}

compute_GSI <- function(V, Si1, ST, paramNames, nbParam) {
  # Function that computes GSI indices from a set of Sobol' indices computed for multiple variables (e.g. several time-steps of a dynamic output, membership functions for different clusters, ...)
  # V: vector of variances (one value per variable)
  # Si1: data.frame containing the Sobol' main indices per parameter (row) and variable (column)
  # ST: data.frame containing the Sobol' total indices per parameter (row) and variable (column)
  clust_GSI <- list(Si1=setNames(sapply(1:nbParam,function(x) sum(V * Si1[x,]) / sum(V)), paramNames),
                    ST=setNames(sapply(1:nbParam,function(x) sum(V * ST[x,]) / sum(V)), paramNames))
  
}
