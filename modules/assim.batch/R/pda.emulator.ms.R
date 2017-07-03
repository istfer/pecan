##' Paramater Data Assimilation using emulator on multiple sites in three modes: local, global, hierarchical
##'
##' @title Paramater Data Assimilation using emulator on multiple sites
##' @param settings = a pecan settings list
##'
##' @return settings
##'
##' @author Istem Fer
##' @export
pda.emulator.ms <- function(settings, params.id = NULL, param.names = NULL, prior.id = NULL, 
                            chain = NULL, iter = NULL, adapt = NULL, adj.min = NULL, 
                            ar.target = NULL, jvar = NULL, n.knot = NULL) {
  
  ## this bit of code is useful for defining the variables passed to this function if you are
  ## debugging
  if (FALSE) {
    params.id <- param.names <- prior.id <- chain <- iter <- NULL
    n.knot <- adapt <- adj.min <- ar.target <- jvar <- NULL
  }
  
  # check mode 
  pda.mode <- settings$assim.batch$mode
  
  # how many sites
  nsites <- length(settings)
  
  if(pda.mode == "local"){
    local <- TRUE
    global <- hierarchical <- FALSE
  }else if(pda.mode == "global"){
    global <- TRUE
    local <- hierarchical <- FALSE
  }else if(pda.mode == "hierarchical"){
    hieararchical <- TRUE
    local <- global <- FALSE
  }else{
    local <- global <- hierarchical <- TRUE
  }
  
  ## Open database connection
  if (settings$database$bety$write) {
    con <- try(db.open(settings$database$bety), silent = TRUE)
    if (is(con, "try-error")) {
      con <- NULL
    } else {
      on.exit(db.close(con))
    }
  } else {
    con <- NULL
  }
  
  bety <- src_postgres(dbname = settings$database$bety$dbname, 
                       host = settings$database$bety$host, 
                       user = settings$database$bety$user, 
                       password = settings$database$bety$password)
  
  ## Get the workflow id
  if ("workflow" %in% names(settings)) {
    workflow.id <- settings$workflow$id
  } else {
    workflow.id <- -1
  }
  
  
  ## Set model-specific functions
  do.call("library", list(paste0("PEcAn.", settings$model$type)))
  my.write.config <- paste("write.config.", settings$model$type, sep = "")
  if (!exists(my.write.config)) {
    logger.severe(paste(my.write.config, 
                        "does not exist. Please make sure that the PEcAn interface is loaded for", 
                        settings$model$type))
  }
  
  ## Create an ensemble id
  settings$assim.batch$ensemble.id <- pda.create.ensemble(settings, con, workflow.id)
  
  
  multi.settings <- settings
  
  ## history restart
  pda.restart.file <- file.path(settings$outdir,paste0("history.pda",
                                                       settings$assim.batch$ensemble.id, ".Rdata"))
  
  current.step <- "START" 
  
  
  ## -------------------------------------- Runs and build emulator ------------------------------------------ 
  
  # lists to collect emulators and run MCMC per site later
  SS.stack       <- list()
  gp.stack       <- list()
  prior.stack    <- list()
  nstack         <- list()
  
  for(s in seq_along(multi.settings)){ # site - run loop
    
    settings <- multi.settings[[s]]
    
    settings$assim.batch$inputs <- settings$run$assim.batch$inputs
    
    ## Handle settings
    settings <- pda.settings(settings)
    
    ## will be used to check if multiplicative Gaussian is requested
    any.mgauss <- sapply(settings$assim.batch$inputs, `[[`, "likelihood")
    isbias <- which(unlist(any.mgauss) == "multipGauss")
    
    ## check if scaling factors are gonna be used
    any.scaling <- sapply(settings$assim.batch$param.names, `[[`, "scaling")
    sf <- unique(unlist(any.scaling))
    
    
    ## Load priors
    #temp        <- pda.load.priors(settings, bety$con, TRUE)
    #prior.list  <- temp$prior
    load("/fs/data2/output//PEcAn_1000007695/normal_priors.Rdata")
    #settings    <- temp$settings
    pname       <- lapply(prior.list, rownames)
    n.param.all <- sapply(prior.list, nrow)
    
    
    inputs      <- load.pda.data(settings, bety)
    n.input     <- length(inputs)
    
    ## Select parameters to constrain
    prior.ind <- lapply(seq_along(settings$pfts), 
                        function(x) which(pname[[x]] %in% settings$assim.batch$param.names[[x]]))
    n.param <- sapply(prior.ind, length)
    prior.ind.orig <- lapply(seq_along(settings$pfts), 
                             function(x) which(pname[[x]] %in% settings$assim.batch$param.names[[x]] |
                                                 pname[[x]] %in% any.scaling[[x]]))
    n.param.orig <- sapply(prior.ind.orig, length)
    
    
    
    ## Set up likelihood functions
    llik.fn <- pda.define.llik.fn(settings)
    
    if(!is.null(sf)){
      sf.ind <- length(prior.list) + 1
      sf.list <- pda.generate.sf(settings$assim.batch$n.knot, sf, prior.list)
      probs.sf <- sf.list$probs
      prior.list <- sf.list$priors
    }else {
      probs.sf <- NULL
    }
    
    ## Set prior distribution functions (d___, q___, r___, and multivariate versions)
    prior.fn <- lapply(prior.list, pda.define.prior.fn)
    
    
    ## Propose parameter knots (X) for emulator design
    knots.list <- lapply(seq_along(settings$pfts), 
                         function(x) pda.generate.knots(settings$assim.batch$n.knot, sf, probs.sf,
                                                        n.param.all[x], 
                                                        prior.ind.orig[[x]], 
                                                        prior.fn[[x]], 
                                                        pname[[x]]))
    names(knots.list) <- sapply(settings$pfts,"[[",'name')
    
    knots.params <- lapply(knots.list, `[[`, "params")
    knots.probs <- lapply(knots.list, `[[`, "probs")
    
    current.step <- paste0("GENERATE KNOTS - site: ", s)
    save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
    
    ## Set up runs and write run configs for all proposed knots
    run.ids <- pda.init.run(settings, con, my.write.config, workflow.id, knots.params, 
                            n = settings$assim.batch$n.knot, 
                            run.names = paste0(settings$assim.batch$ensemble.id, 
                                               ".site.", s,
                                               ".knot.", 1:settings$assim.batch$n.knot))   
    
    current.step <- paste0("pda.init.run - site: ", s)
    save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
    
    ## start model runs
    start.model.runs(settings, settings$database$bety$write)
    
    ## Retrieve model outputs and error statistics
    model.out <- list()
    pda.errors <- list()
    
    ## read model outputs    
    for (i in seq_len(settings$assim.batch$n.knot)) {
      align.return <- pda.get.model.output(settings, run.ids[i], bety, inputs)
      model.out[[i]] <- align.return$model.out
      if(all(!is.na(model.out[[i]]))){
        inputs <- align.return$inputs
      }
    }
    
    current.step <- paste0("pda.get.model.output - site: ", s)
    save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
    
    # efficient sample size calculation
    #inputs <- pda.neff.calc(inputs)
    #for developing
    inputs[[1]]$n_eff <- 2000
    inputs[[2]]$n_eff <- 500
    
    # handle bias parameters if multiplicative Gaussian is listed in the likelihoods
    if(any(unlist(any.mgauss) == "multipGauss")) {
      # how many bias parameters per dataset requested
      nbias <- ifelse(is.null(settings$assim.batch$inputs[[isbias]]$nbias), 1,
                      as.numeric(settings$assim.batch$inputs[[isbias]]$nbias))
      bias.list <- return.bias(isbias, model.out, inputs, prior.list, nbias, run.round, settings$assim.batch$bias.path)
      bias.terms <- bias.list$bias.params
      prior.list <- bias.list$prior.list.bias
      prior.fn <- lapply(prior.list, pda.define.prior.fn)
    } else {
      bias.terms <- NULL
    }
    
    for (i in seq_len(settings$assim.batch$n.knot)) {
      if(!is.null(bias.terms)){
        all.bias <- lapply(bias.terms, function(n) n[i,])
        all.bias <- do.call("rbind", all.bias)
      } else {
        all.bias <- NULL
      }
      ## calculate error statistics and save in the DB      
      pda.errors[[i]] <- pda.calc.error(settings, con, model_out = model.out[[i]], run.id = run.ids[i], inputs, bias.terms = all.bias)
    } 
    
    current.step <- paste0("pda.calc.error - site: ", s)
    save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
    
    
    prior.all <- do.call("rbind", prior.list)
    length.pars <- 0
    prior.ind.list <- prior.ind.list.ns <- list()
    # now I need to go through all parameters for each pft, but leave out the ones that scaling factor is requested
    for(p in seq_along(settings$assim.batch$param.names)){
      param.names <- settings$assim.batch$param.names[[p]]
      prior.ind.list[[p]] <- length.pars + which(pname[[p]] %in% unlist(param.names) &
                                                   !(pname[[p]] %in% sf))
      prior.ind.list.ns[[p]] <- length.pars + which(pname[[p]] %in% unlist(param.names))
      length.pars <- length.pars + length(pname[[p]])
    }
    prior.ind.all    <- unlist(prior.ind.list)
    prior.ind.all.ns <- unlist(prior.ind.list.ns)
    # if no scaling is requested prior.ind.all == prior.ind.all.ns
    # keep this ind.all w/o bias until extracting prob values below
    
    # retrieve n
    n.of.obs <- sapply(inputs,`[[`, "n") 
    names(n.of.obs) <- sapply(model.out[[1]],names)
    
    
    # UPDATE: Use mlegp package, I can now draw from parameter space
    knots.params.all <- do.call("cbind", knots.params)
    X <- knots.params.all[, prior.ind.all, drop = FALSE]
    
    if(!is.null(sf)){
      X <- cbind(X, probs.sf)
    }
    
    # retrieve SS
    error.statistics <- list()
    SS.list <- list()
    bc <- 1
    
    # what percentage of runs is allowed to fail?
    if(!is.null(settings$assim.batch$allow.fail)){
      allow.fail <- as.numeric(settings$assim.batch$allow.fail)
    } else {
      allow.fail <- 0.5
    }
    # what is it in number of runs?
    no.of.allowed <- floor(settings$assim.batch$n.knot * allow.fail)
    
    for(inputi in seq_len(n.input)){
      error.statistics[[inputi]] <- sapply(pda.errors,`[[`, inputi)
      
      if(unlist(any.mgauss)[inputi] == "multipGauss") {
        
        # if yes, then we need to include bias term in the emulator
        #bias.probs <- bias.list$bias.probs
        #biases <- c(t(bias.probs[[bc]]))
        bias.params <- bias.list$bias.params
        biases <- c(t(bias.params[[bc]]))
        bc <- bc + 1
        
        # replicate model parameter set per bias parameter
        rep.rows <- rep(1:nrow(X), each = nbias)
        X.rep <- X[rep.rows,]
        Xnew <- cbind(X.rep, biases)
        colnames(Xnew) <- c(colnames(X.rep), paste0("bias.", names(n.of.obs)[inputi]))
        SS.list[[inputi]] <- cbind(Xnew, c(error.statistics[[inputi]]))
        
      } else {
        SS.list[[inputi]] <- cbind(X, error.statistics[[inputi]])
      } # if-block
      
      # check failed runs and remove them if you'll have a reasonable amount of param sets after removal
      # how many runs failed?
      no.of.failed <- sum(is.na(SS.list[[inputi]][, ncol(SS.list[[inputi]])]))
      
      # check if you're left with enough sets
      if(no.of.failed < no.of.allowed & (settings$assim.batch$n.knot - no.of.failed) > 1){
        SS.list[[inputi]] <- SS.list[[inputi]][!rowSums(is.na(SS.list[[inputi]])), ]
        if( no.of.failed  > 0){
          logger.info(paste0(no.of.failed, " runs failed. Emulator for ", names(n.of.obs)[inputi], " will be built with ", settings$assim.batch$n.knot - no.of.failed, " knots."))
        } 
      } else{
        logger.error(paste0("Too many runs failed, not enough parameter set to build emulator for ", names(n.of.obs)[inputi], "."))
      }
      
    } # for-loop
    
    
    SS <- SS.list
    
    
    logger.info(paste0("Using 'mlegp' package for Gaussian Process Model fitting."))
    
    ## Generate emulator on SS, return a list ##
    
    # start the clock
    ptm.start <- proc.time()
    
    # prepare for parallelization
    dcores <- parallel::detectCores() - 1
    ncores <- min(max(dcores, 1), length(SS))
    
    cl <- parallel::makeCluster(ncores, type="FORK")
    
    ## Parallel fit for GPs
    GPmodel <- parallel::parLapply(cl, SS, function(x) mlegp::mlegp(X = x[, -ncol(x), drop = FALSE], Z = x[, ncol(x), drop = FALSE], verbose = 0))
    # GPmodel <- lapply(SS, function(x) mlegp::mlegp(X = x[, -ncol(x), drop = FALSE], Z = x[, ncol(x), drop = FALSE], verbose = 0))
    
    parallel::stopCluster(cl)
    
    # Stop the clock
    ptm.finish <- proc.time() - ptm.start
    logger.info(paste0("GP fitting took ", paste0(round(ptm.finish[3])), " seconds."))
    
    SS.stack[[s]]       <- SS
    gp.stack[[s]]       <- GPmodel
    # without any bias factor and scaling factor all prior lists should be the same
    # putting them into their own sublists as a reminder only
    prior.stack[[s]]    <- prior.all 
    multi.settings[[s]] <- settings
    nstack[[s]] <- n.of.obs
    
    current.step <- paste0("fit GP - site: ", s)
    save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
    
  } # end site - run loop
  
  
  # add indice and increase n.param for scaling factor
  if(!is.null(sf)){
    prior.ind.all <- c(prior.ind.all, 
                       ((length.pars + 1): (length.pars + length(sf))))
    n.param       <- c(n.param, length(sf))
    length.pars   <- length.pars + length(sf)
  }
  
  # add indice and increase n.param for bias
  if(any(unlist(any.mgauss) == "multipGauss")){
    prior.ind.all <- c(prior.ind.all, 
                       ((length.pars + 1) : (length.pars + length(isbias))))
    prior.ind.all.ns <- c(prior.ind.all.ns, 
                          ((length.pars + 1) : (length.pars + length(isbias))))
    n.param <- c(n.param, length(isbias))
    n.param.orig <- c(n.param.orig, length(isbias))
    length.pars   <- length.pars + length(isbias)
  }
  
  
  
  if (!is.null(settings$assim.batch$mix)) {
    mix <- settings$assim.batch$mix
  } else if (sum(n.param) > 1) {
    mix <- "joint"
  } else {
    mix <- "each"
  }
  
  init.list <- list()
  jmp.list <- list()
  
  ## -------------------------------------- Local MCMC ------------------------------------------ 
  if(local){ # local - if
    
    for(s in seq_along(multi.settings)){ 
      
      settings <- multi.settings[[s]]
      
      ## Set up prior functions accordingly
      prior.fn.all <- pda.define.prior.fn(prior.stack[[s]])
      
      # define range to make sure mcmc.GP doesn't propose new values outside
      rng <- matrix(c(sapply(prior.fn.all$qprior[prior.ind.all], eval, list(p = 1e-05)),
                      sapply(prior.fn.all$qprior[prior.ind.all], eval, list(p = 0.99999))),
                    nrow = sum(n.param))
      
      resume.list <- list()
      
      for (c in seq_len(settings$assim.batch$chain)) {
        jmp.list[[c]] <- sapply(prior.fn.all$qprior, 
                                function(x) 0.1 * diff(eval(x, list(p = c(0.05, 0.95)))))[prior.ind.all]
        jmp.list[[c]] <- sqrt(jmp.list[[c]])
        
        init.x <- lapply(prior.ind.all, function(v) eval(prior.fn.all$rprior[[v]], list(n = 1)))
        names(init.x) <- rownames(prior.stack[[s]])[prior.ind.all]
        init.list[[c]] <- init.x
        resume.list[[c]] <- NA
      }
      
      current.step <- paste0("pre-MCMC - site: ", s)
      save(list = ls(all.names = TRUE), envir=environment(),file=pda.restart.file)
      
      # # start the clock
      ptm.start <- proc.time()
      
      # prepare for parallelization
      dcores <- parallel::detectCores() - 1
      ncores <- min(max(dcores, 1), settings$assim.batch$chain)
      
      logger.setOutputFile(file.path(settings$outdir, "pda.log"))
      
      cl <- parallel::makeCluster(ncores, type="FORK", outfile = file.path(settings$outdir, "pda.log"))

      ## Sample posterior from emulator
      mcmc.out <- parallel::parLapply(cl, 1:settings$assim.batch$chain, function(chain) {
        mcmc.GP(gp          = gp.stack[[s]], ## Emulator(s)
                x0          = init.list[[chain]],     ## Initial conditions
                nmcmc       = settings$assim.batch$iter,       ## Number of reps
                rng         = rng,       ## range
                format      = "lin",      ## "lin"ear vs "log" of LogLikelihood 
                mix         = mix,     ## Jump "each" dimension independently or update them "joint"ly
                jmp0        = jmp.list[[chain]],  ## Initial jump size
                ar.target   = settings$assim.batch$jump$ar.target,   ## Target acceptance rate
                priors      = prior.fn.all$dprior[prior.ind.all], ## priors
                settings    = multi.settings[[s]],
                run.block   = TRUE,  
                n.of.obs    = nstack[[s]],
                llik.fn     = llik.fn,
                resume.list = resume.list[[chain]]
        )
      })
      
      parallel::stopCluster(cl)
      
      # Stop the clock
      ptm.finish <- proc.time() - ptm.start
      logger.info(paste0("Emulator MCMC took ", paste0(round(ptm.finish[3])), " seconds for ", paste0(settings$assim.batch$iter), " iterations."))
      
      current.step <- paste0("post-MCMC - site: ", s)
      save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
      
      mcmc.samp.list <- list()
      
      for (c in seq_len(settings$assim.batch$chain)) {
        
        m <- matrix(NA, nrow =  nrow(mcmc.out[[c]]$mcmc.samp), ncol = length(prior.ind.all.ns))
        
        if(!is.null(sf)){
          sfm <- matrix(NA, nrow =  nrow(mcmc.out[[c]]$mcmc.samp), ncol = length(sf))
        }
        
        # TODO: get back to this when scaling factor is used
        # # retrieve rownames separately to get rid of var_name* structures
        # prior.all.rownames <- unlist(sapply(prior.list, rownames))
        prior.all.rownames <- rownames(prior.stack[[s]])
        
        sc <- 1
        for (i in seq_along(prior.ind.all.ns)) {
          sf.check <- prior.all.rownames[prior.ind.all.ns][i]
          idx <- grep(sf.check, rownames(prior.all)[prior.ind.all])
          if(any(grepl(sf.check, sf))){
            
            m[, i] <- eval(prior.fn.all$qprior[prior.ind.all.ns][[i]],
                           list(p = mcmc.out[[c]]$mcmc.samp[, idx]))
            if(sc <= length(sf)){
              sfm[, sc] <- mcmc.out[[c]]$mcmc.samp[, idx]
              sc <- sc + 1
            }
            
          }else{
            m[, i] <- mcmc.out[[c]]$mcmc.samp[, idx]
          }
        }
        
        colnames(m) <- prior.all.rownames[prior.ind.all.ns]
        mcmc.samp.list[[c]] <- m
        
        if(!is.null(sf)){
          colnames(sfm) <- paste0(sf, "_SF")
          sf.samp.list[[c]] <- sfm
        }
        
        resume.list[[c]] <- mcmc.out[[c]]$chain.res
      }
      
      # Separate each PFT's parameter samples (and bias term) to their own list
      mcmc.param.list <- list()
      ind <- 0
      for (i in seq_along(n.param.orig)) {
        mcmc.param.list[[i]] <- lapply(mcmc.samp.list, function(x) x[, (ind + 1):(ind + n.param.orig[i]), drop = FALSE])
        ind <- ind + n.param.orig[i]
      }
      
      # Collect non-model parameters in their own list
      if(length(mcmc.param.list) > length(settings$pfts)) { 
        # means bias parameter was at least one bias param in the emulator
        # it will be the last list in mcmc.param.list
        # there will always be at least one tau for bias
        for(c in seq_len(settings$assim.batch$chain)){
          mcmc.param.list[[length(mcmc.param.list)]][[c]] <- cbind( mcmc.param.list[[length(mcmc.param.list)]][[c]],
                                                                    mcmc.out[[c]]$mcmc.par)
        }
        
      } else if (ncol(mcmc.out[[1]]$mcmc.par) != 0){
        # means no bias param but there are still other params, e.g. Gaussian
        mcmc.param.list[[length(mcmc.param.list)+1]] <- list()
        for(c in seq_len(settings$assim.batch$chain)){
          mcmc.param.list[[length(mcmc.param.list)]][[c]] <- mcmc.out[[c]]$mcmc.par
        }
      }
      
      
      # file flag will be needed to pass as these will al go into pft folder
      settings <- pda.postprocess(settings, con, mcmc.param.list, 
                                  pname, prior.list, prior.ind.orig, sffx = paste0("_site.", settings$run$site$id))
      
      multi.settings[[s]] <- settings
      
      ## save updated settings XML
      saveXML(listToXml(settings, "pecan"), file = file.path(settings$outdir, paste0("pecan.pda", settings$assim.batch$ensemble.id, 
                                                                                     "_site.", settings$run$site$id, ".xml")))
      
    } # end - local
    
    ## close database connection
    if (!is.null(con)) {
      db.close(con)
    }
    
    
    
  }
  
  
  ## -------------------------------------- Global MCMC ------------------------------------------ 
  if(global){ # global - if
    
    s<-1
    settings <- multi.settings[[s]]
    
    ## Set up prior functions accordingly
    prior.fn.all <- pda.define.prior.fn(prior.stack[[s]])
    
    # define range to make sure mcmc.GP doesn't propose new values outside
    rng <- matrix(c(sapply(prior.fn.all$qprior[prior.ind.all], eval, list(p = 1e-05)),
                    sapply(prior.fn.all$qprior[prior.ind.all], eval, list(p = 0.99999))),
                  nrow = sum(n.param))
    
    resume.list <- list()
    
    for (c in seq_len(settings$assim.batch$chain)) {
      jmp.list[[c]] <- sapply(prior.fn.all$qprior, 
                              function(x) 0.1 * diff(eval(x, list(p = c(0.05, 0.95)))))[prior.ind.all]
      jmp.list[[c]] <- sqrt(jmp.list[[c]])
      
      init.x <- lapply(prior.ind.all, function(v) eval(prior.fn.all$rprior[[v]], list(n = 1)))
      names(init.x) <- rownames(prior.stack[[s]])[prior.ind.all]
      init.list[[c]] <- init.x
      resume.list[[c]] <- NA
    }
    
    current.step <- paste0("pre-MCMC - site: ", s)
    save(list = ls(all.names = TRUE), envir=environment(),file=pda.restart.file)
    
    gp <- unlist(gp.stack, recursive = FALSE)
    # start the clock
    ptm.start <- proc.time()
     
    # prepare for parallelization
    dcores <- parallel::detectCores() - 1
    ncores <- min(max(dcores, 1), settings$assim.batch$chain)
    # 
    logger.setOutputFile(file.path(settings$outdir, "pda.log"))
    # 
    cl <- parallel::makeCluster(ncores, type="FORK", outfile = file.path(settings$outdir, "pda.log"))
    
    ## Sample posterior from emulator
    mcmc.out <- parallel::parLapply(cl, 1:settings$assim.batch$chain, function(chain) {
      mcmc.GP(gp          = gp, ## Emulator(s)
              x0          = init.list[[chain]],     ## Initial conditions
              nmcmc       = settings$assim.batch$iter,       ## Number of reps
              rng         = rng,       ## range
              format      = "lin",      ## "lin"ear vs "log" of LogLikelihood 
              mix         = mix,     ## Jump "each" dimension independently or update them "joint"ly
              jmp0        = jmp.list[[chain]],  ## Initial jump size
              ar.target   = settings$assim.batch$jump$ar.target,   ## Target acceptance rate
              priors      = prior.fn.all$dprior[prior.ind.all], ## priors
              settings    = multi.settings[[s]], # this is just for checking llik functions downstream
              run.block   = TRUE,  
              n.of.obs    = unlist(nstack),
              llik.fn     = llik.fn,
              resume.list = resume.list[[chain]]
      )
    })
    
    parallel::stopCluster(cl)
    
    # Stop the clock
    ptm.finish <- proc.time() - ptm.start
    logger.info(paste0("Emulator MCMC took ", paste0(round(ptm.finish[3])), " seconds for ", paste0(settings$assim.batch$iter), " iterations."))
    
    current.step <- paste0("post-MCMC - site: ", s)
    save(list = ls(all.names = TRUE),envir=environment(),file=pda.restart.file)
    
    mcmc.samp.list <- list()
    
    for (c in seq_len(settings$assim.batch$chain)) {
      
      m <- matrix(NA, nrow =  nrow(mcmc.out[[c]]$mcmc.samp), ncol = length(prior.ind.all.ns))
      
      if(!is.null(sf)){
        sfm <- matrix(NA, nrow =  nrow(mcmc.out[[c]]$mcmc.samp), ncol = length(sf))
      }
      
      # TODO: get back to this when scaling factor is used
      # # retrieve rownames separately to get rid of var_name* structures
      # prior.all.rownames <- unlist(sapply(prior.list, rownames))
      prior.all.rownames <- rownames(prior.stack[[s]])
      
      sc <- 1
      for (i in seq_along(prior.ind.all.ns)) {
        sf.check <- prior.all.rownames[prior.ind.all.ns][i]
        idx <- grep(sf.check, rownames(prior.all)[prior.ind.all])
        if(any(grepl(sf.check, sf))){
          
          m[, i] <- eval(prior.fn.all$qprior[prior.ind.all.ns][[i]],
                         list(p = mcmc.out[[c]]$mcmc.samp[, idx]))
          if(sc <= length(sf)){
            sfm[, sc] <- mcmc.out[[c]]$mcmc.samp[, idx]
            sc <- sc + 1
          }
          
        }else{
          m[, i] <- mcmc.out[[c]]$mcmc.samp[, idx]
        }
      }
      
      colnames(m) <- prior.all.rownames[prior.ind.all.ns]
      mcmc.samp.list[[c]] <- m
      
      if(!is.null(sf)){
        colnames(sfm) <- paste0(sf, "_SF")
        sf.samp.list[[c]] <- sfm
      }
      
      resume.list[[c]] <- mcmc.out[[c]]$chain.res
    }
    
    # Separate each PFT's parameter samples (and bias term) to their own list
    mcmc.param.list <- list()
    ind <- 0
    for (i in seq_along(n.param.orig)) {
      mcmc.param.list[[i]] <- lapply(mcmc.samp.list, function(x) x[, (ind + 1):(ind + n.param.orig[i]), drop = FALSE])
      ind <- ind + n.param.orig[i]
    }
    
    # Collect non-model parameters in their own list
    if(length(mcmc.param.list) > length(settings$pfts)) { 
      # means bias parameter was at least one bias param in the emulator
      # it will be the last list in mcmc.param.list
      # there will always be at least one tau for bias
      for(c in seq_len(settings$assim.batch$chain)){
        mcmc.param.list[[length(mcmc.param.list)]][[c]] <- cbind( mcmc.param.list[[length(mcmc.param.list)]][[c]],
                                                                  mcmc.out[[c]]$mcmc.par)
      }
      
    } else if (ncol(mcmc.out[[1]]$mcmc.par) != 0){
      # means no bias param but there are still other params, e.g. Gaussian
      mcmc.param.list[[length(mcmc.param.list)+1]] <- list()
      for(c in seq_len(settings$assim.batch$chain)){
        mcmc.param.list[[length(mcmc.param.list)]][[c]] <- mcmc.out[[c]]$mcmc.par
      }
    }
    
    # file flag will be needed to pass as these will al go into pft folder
    settings <- pda.postprocess(settings, con, mcmc.param.list, pname, prior.list, prior.ind.orig, sffx = "_global")
    
    
    # HOW TO SAVE MULTI SETTINGS
    ## save updated settings XML
    saveXML(listToXml(settings, "pecan"), file = file.path(settings$outdir, paste0("pecan.pda", settings$assim.batch$ensemble.id, 
                                                                                   "_site.", settings$run$site$id, ".xml")))
    
    
  } # end - global
  
  
  ## -------------------------------------- Hierarchical MCMC ------------------------------------------ 
  if(hierarchical){ # hierarchical - if
    
    s<-1
    settings <- multi.settings[[s]]
    
    ## Set up prior functions accordingly
    global.prior.fn.all <- pda.define.prior.fn(prior.stack[[s]])
    

    for (c in seq_len(settings$assim.batch$chain)) {
      jmp.list[[c]] <- sapply(global.prior.fn.all$qprior, 
                              function(x) 0.1 * diff(eval(x, list(p = c(0.05, 0.95)))))[prior.ind.all]
      jmp.list[[c]] <- sqrt(jmp.list[[c]])
      
      init.x <- lapply(prior.ind.all, function(v) eval(global.prior.fn.all$rprior[[v]], list(n = 1)))
      names(init.x) <- rownames(prior.stack[[s]])[prior.ind.all]
      init.list[[c]] <- init.x
    }
    
    # develop with one chain 
    chain <- 1
    x0          = init.list[[chain]]
    nmcmc       = 1000
    jmp0        = jmp.list[[chain]]
    priors      = global.prior.fn.all$dprior[prior.ind.all]
    
          # e.g. normal priors
          # 
          # dnorm(x,3,1.25,log=TRUE)
          # dnorm(x,0.25,0.2,log=TRUE)
          # dnorm(x,18,12,log=TRUE)
          # dnorm(x,18,11,log=TRUE)
          # dnorm(x,0.13,0.12,log=TRUE)
          # dnorm(x,135,30,log=TRUE)

    ################################################################
    #
    #      mu_site    : site level parameters (nsite x nparam)
    #      tau_site   : site level precision (nsite x nparam x nparam)
    #      mu_global  : global parameters (nparam)
    #      tau_global : global precision matrix (nparam x nparam)
    #
    #
    #
    #      ### process model  
    #      mu_site ~ MVN (mu_global, tau_global)
    #      (in R semantics =)  mu_site <- mvrnorm(1, mu_global, sigma_global)
    #
    #      ### priors
    #      mu_f : prior mean vector
    #      P_f  : prior covariance matrix
    #      mu_global ~ MVN (mu_f, P_f)
    #
    #      tau_df    : Wishart degrees of freedom
    #      tau_scale : Wishart scale matrix
    #      tau_global ~ W (tau_df, tau_scale)
    #      sigma_global <- solve(tau_global)
    #
    ###############################################################

    # prior mean vector
    mu_f     <- sapply(global.prior.fn.all$qprior, function(x) eval(x, list(p = c(0.5))))[prior.ind.all] 
    P_f      <- diag((jmp0)^2) # prior covariance matrix
    P_f_inv  <- solve(P_f) # prior precision matrix
    
    # initialize mu_global
    mu_global <- mvrnorm(1, mu_f, P_f)

    # initialize tau_global, how to initialize?
    tau_df <- nsites + sum(n.param) + 1
    tau_global   <- rWishart(1, tau_df, P_f)[,,1]
    
    # initialize mu_site
    sigma_global <- solve(tau_global) 
    mu_site_curr <- mvrnorm(nsites, mu_global, sigma_global) # site mean
    

    currSS  <- sapply(seq_len(nsites), function(v) get_ss(gp.stack[[v]], mu_site_curr[v,]))
    currllp <- lapply(seq_len(nsites), function(v) pda.calc.llik.par(settings, nstack[[v]], currSS[,v]))
    
    # storage
    mu_site_samp <-  array(NA_real_, c(nmcmc, sum(n.param), nsites))
    # tau_site_samp <- array(NA_real_, c(nmcmc, nsites, sum(n.param), sum(n.param)))
    mu_global_samp  <-  matrix(NA_real_, nrow = nmcmc, ncol= sum(n.param))
    tau_global_samp <-  array(NA_real_, c(nmcmc, sum(n.param), sum(n.param)))

    accept.count <- rep(0, nsites)
    
    g <- 1
    for(g in seq_len(nmcmc)){
      
      # # adapt
      # if ((g > 2) && ((g - 1) %% settings$assim.batch$jump$adapt == 0)) {
      #   params.recent <- mu_site_samp[(g - settings$assim.batch$jump$adapt):(g - 1), , ]
      #   colnames(params.recent) <- names(x0)
      #   # accept.count <- round(jmp@arate[(g-1)/settings$assim.batch$jump$adapt]*100)
      #   jcov.list <- lapply(seq_len(nsites), function(v) pda.adjust.jumps.bs(settings, jcov.arr[,,v], accept.count[v], params.recent[,,v]))
      #   jcov.arr  <- abind(jcov.list, along=3)
      #   accept.count <- rep(0, nsites)  # Reset counter
      # }


      
      ########################################
      # update tau_global | mu_global, mu_site
      #
      # tau_global ~ W(df, Sigma)
      # 
      # tau_global   : error precision matrix
       
      # sum of pairwise deviation products
      sum_term <- matrix(0, ncol = sum(n.param), nrow = sum(n.param))
      for(i in seq_len(nsites)){
        pairwise_deviation <- as.matrix(mu_site_curr[i,] - mu_global)
        sum_term <- sum_term + pairwise_deviation %*% t(pairwise_deviation)
      }
      
      # what is V? which prior?
      V_inv <- solve(P_f) 
      tau_sigma <- solve(V_inv + sum_term)
      
      # update tau
      tau_global <- rWishart(1, df = tau_df, Sigma = tau_sigma)[,,1] # site precision
      sigma_global <- solve(tau_global) # site covariance
      

      ########################################
      # update mu_global | mu_site, tau_global
      #
      # mu_global ~ MVN(global_mu, global_Sigma)
      #
      # mu_global     : global parameters
      # global_mu     : precision weighted average between the data (mu_site) and prior mean (mu_f)
      # global_Sigma  : sum of mu_site and mu_f precision
      

      
      global_Sigma <- solve(tau_global + P_f)
      
      # global_mu dimensions need to be 1x6? because mu_global is 1x6?
      global_mu <- global_Sigma %*% ((tau_global %*% colMeans(mu_site_curr)) + P_f %*% mu_f)

      mu_global <- mvrnorm(1, global_mu, global_Sigma)
      

      # site level M-H
      ########################################
      # propose mu_site | mu_global, tau_global
      #
      mu_site_new <- mvrnorm(nsites, mu_global, sigma_global)

      # re-predict current SS
      currSS <- sapply(seq_len(nsites), function(v) get_ss(gp.stack[[v]], mu_site_curr[v,]))
      ycurr  <- sapply(seq_len(nsites), function(v) pda.calc.llik(currSS[,v], llik.fn, currllp[[v]]))
      
      # predict new SS
      newSS <- sapply(seq_len(nsites), function(v) get_ss(gp.stack[[v]], mu_site_new[v,]))
      newllp <- lapply(seq_len(nsites), function(v) pda.calc.llik.par(settings, nstack[[v]], newSS[,v]))
      ynew <- sapply(seq_len(nsites), function(v) pda.calc.llik(newSS[,v], llik.fn, newllp[[v]]))
      
      ar <- is.accepted(ycurr, ynew)
      mu_site_curr[ar, ] <- mu_site_new[ar, ]
      accept.count <- accept.count + ar
      
      
      mu_site_samp[g, , seq_len(nsites)] <- t(mu_site_curr)[,seq_len(nsites)]
      mu_global_samp[g,] <- mu_global
      tau_global_samp[g,,] <- tau_global
      
    }
    
    
    
    mcmc.samp.list <- list()
    
    for (c in seq_len(settings$assim.batch$chain)) {
      
      m <- matrix(NA, nrow =  nrow(mcmc.out[[c]]$mcmc.samp), ncol = length(prior.ind.all.ns))
      
      if(!is.null(sf)){
        sfm <- matrix(NA, nrow =  nrow(mcmc.out[[c]]$mcmc.samp), ncol = length(sf))
      }
      
      # TODO: get back to this when scaling factor is used
      # # retrieve rownames separately to get rid of var_name* structures
      # prior.all.rownames <- unlist(sapply(prior.list, rownames))
      prior.all.rownames <- rownames(prior.stack[[s]])
      
      sc <- 1
      for (i in seq_along(prior.ind.all.ns)) {
        sf.check <- prior.all.rownames[prior.ind.all.ns][i]
        idx <- grep(sf.check, rownames(prior.all)[prior.ind.all])
        if(any(grepl(sf.check, sf))){
          
          m[, i] <- eval(prior.fn.all$qprior[prior.ind.all.ns][[i]],
                         list(p = mcmc.out[[c]]$mcmc.samp[, idx]))
          if(sc <= length(sf)){
            sfm[, sc] <- mcmc.out[[c]]$mcmc.samp[, idx]
            sc <- sc + 1
          }
          
        }else{
          m[, i] <- mcmc.out[[c]]$mcmc.samp[, idx]
        }
      }
      
      colnames(m) <- prior.all.rownames[prior.ind.all.ns]
      mcmc.samp.list[[c]] <- m
      
      if(!is.null(sf)){
        colnames(sfm) <- paste0(sf, "_SF")
        sf.samp.list[[c]] <- sfm
      }
      
      resume.list[[c]] <- mcmc.out[[c]]$chain.res
    }
    
    # Separate each PFT's parameter samples (and bias term) to their own list
    mcmc.param.list <- list()
    ind <- 0
    for (i in seq_along(n.param.orig)) {
      mcmc.param.list[[i]] <- lapply(mcmc.samp.list, function(x) x[, (ind + 1):(ind + n.param.orig[i]), drop = FALSE])
      ind <- ind + n.param.orig[i]
    }
    
    # Collect non-model parameters in their own list
    if(length(mcmc.param.list) > length(settings$pfts)) { 
      # means bias parameter was at least one bias param in the emulator
      # it will be the last list in mcmc.param.list
      # there will always be at least one tau for bias
      for(c in seq_len(settings$assim.batch$chain)){
        mcmc.param.list[[length(mcmc.param.list)]][[c]] <- cbind( mcmc.param.list[[length(mcmc.param.list)]][[c]],
                                                                  mcmc.out[[c]]$mcmc.par)
      }
      
    } else if (ncol(mcmc.out[[1]]$mcmc.par) != 0){
      # means no bias param but there are still other params, e.g. Gaussian
      mcmc.param.list[[length(mcmc.param.list)+1]] <- list()
      for(c in seq_len(settings$assim.batch$chain)){
        mcmc.param.list[[length(mcmc.param.list)]][[c]] <- mcmc.out[[c]]$mcmc.par
      }
    }
    
    # file flag will be needed to pass as these will al go into pft folder
    settings <- pda.postprocess(settings, con, mcmc.param.list, pname, prior.list, prior.ind.orig, sffx = "_hierarchical")
    
    
  } # end - hierarchical
  
  
  if (FALSE) {
    gp          = gp.stack[[1]]
    x0          = init.list[[1]]
    nmcmc       = settings$assim.batch$iter
    rng         = rng
    format      = "lin"
    mix         = mix
    jmp0        = jmp.list[[1]]
    ar.target   = settings$assim.batch$jump$ar.target
    priors      = prior.fn.all$dprior[prior.ind.all]
    settings    = settings
    run.block   = TRUE  
    n.of.obs    = n.of.obs
    llik.fn     = llik.fn
    resume.list = resume.list[[1]]
  }
  
  
  return(settings)
  
}  ## end pda.emulator
