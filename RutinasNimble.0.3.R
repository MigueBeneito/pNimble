library(ntfy)
library(coda)
library(parallel)
library(nimble)
library(nimbleHMC)

pNimble<-function(code = NULL, data = NULL, constants = NULL, inits = NULL, nchains = 3, 
                  seeds = NULL, summary = TRUE, monitors = NULL, ntfyAccount = NULL,
                  email = FALSE, WAIC = FALSE, HMC = FALSE, replaceSamplers = NULL,
                  parallel = TRUE, control.model = NULL, control.compile = NULL, 
                  control.configure = NULL, control.build = NULL, ...){

  # Initial arrangements and checkings
  time.start <- Sys.time()
  if(is.null(control.model)) control.model <- list()
  if(!is.null(code)) control.model$code <- code
  if(!is.null(data)) control.model$data <- data
  if(!is.null(constants)) control.model$constants <- constants
  if(!is.null(monitors)){
    if(is.null(control.configure)){control.configure <- list()}
    control.configure$monitors <- monitors
  }
  if(is.null(seeds)) seeds <- 1:nchains
  checkArguments(inits = inits, nchains = nchains, seeds = seeds, summary = summary, monitors = monitors, 
                 ntfyAccount = ntfyAccount, email = email, WAIC = WAIC, HMC = HMC,
                 replaceSamplers = replaceSamplers, parallel = parallel, control.model = control.model, 
                 control.compile = control.compile, control.configure = control.configure, 
                 control.build = control.build)
  
  # Parallelized call to Nimble
  my.cluster <- makeCluster(nchains)
  on.exit(stopCluster(my.cluster))
  if(parallel){
    resul <- parLapply(cl = my.cluster, X = seeds, fun = runParallel,
                       inits = inits, control.model = control.model, 
                       control.compile = control.compile, control.configure = control.configure, 
                       control.build = control.build, HMC = HMC, replaceSamplers = replaceSamplers, 
                       WAIC = WAIC, ...)  
  }
  else{
    resul <- runParallel(seed = seeds[1], inits = inits,control.model = control.model, 
                         control.compile = control.compile, control.configure = control.configure, 
                         control.build = control.build, HMC = HMC, replaceSamplers = replaceSamplers, 
                         WAIC = WAIC, ...)
  }
  # Output arrangement
  names(resul) <- paste0("chain",1:length(resul))
  resul2 <- list()
  resul2$samples <- coda::as.mcmc.list(lapply(resul,function(x){coda::as.mcmc(x)}))
  if(summary){
    if(parallel){
      AllSamples <- do.call(rbind,resul2$samples)
      Var0s <- which(apply(AllSamples,2,var)==0)
      NAs <- which(is.na(apply(AllSamples,2,var)))
      if(length(Var0s)>0 | length(NAs)>0){
        aux <- lapply(resul2$samples,function(x){x[,-c(Var0s, NAs)]})
        resul2$summary <- MCMCvis::MCMCsummary(aux)
        if(length(Var0s)>0){
          newSummaries <- cbind(apply(AllSamples[,Var0s],2,mean), 
                                rep(0,length(Var0s)),
                                rep(0,length(Var0s)),
                                rep(0,length(Var0s)),
                                rep(0,length(Var0s)),
                                rep(NA,length(Var0s)),
                                rep(NA,length(Var0s)))
          rownames(newSummaries) <- colnames(AllSamples)[Var0s]
          colnames(newSummaries) <- colnames(resul2$summary)
          resul2$summary <- rbind(resul2$summary, newSummaries)
        }
        if(length(NAs)>0){
          newSummaries <- matrix(rep(NA,7*length(NAs)), ncol = 7) #
          rownames(newSummaries) <- colnames(AllSamples)[NAs]
          colnames(newSummaries) <- colnames(resul2$summary)
          resul2$summary <- rbind(resul2$summary, newSummaries)
        }
        resul2$summary <- resul2$summary[match(colnames(AllSamples),
                                               rownames(resul2$summary)),]
      }
      else{resul2$summary <- MCMCvis::MCMCsummary(resul2$samples)}
    }
    else cat("summary cannot be calculated with a single chain")
  }
  if(WAIC){
    model.nimble <- nimbleModel(code = code, data = data, constants = constants,
                                check = FALSE, calculate = FALSE, buildDerivs = FALSE)
    cmodel <- compileNimble(model.nimble)
    resul2$WAIC <- calculateWAIC(do.call(rbind,resul2$samples), cmodel)
  } 

  time.end <-Sys.time()
  total.time <- as.numeric(time.end)-as.numeric(time.start)
  notify(time = total.time, email = email, ntfyAccount = ntfyAccount, model = substitute(code))
  
  return(resul2)
}

checkArguments<-function(inits, nchains, seeds, summary, monitors, ntfyAccount, email, WAIC, HMC,
               replaceSamplers, parallel, control.model, control.compile, control.configure, control.build){
  if(!is.null(inits) && !is.function(inits)) stop("inits argument should be a function.")
  if(!is.numeric(nchains)) stop("nchains argument should be numeric.")
  if(!is.numeric(seeds)) stop("seeds argument should be a numeric vector.")
  if(length(seeds) != nchains) stop("length of seeds argument does not match nchains.")
  if(!is.logical(summary)) stop("summary argument should be logical.")
  if(!is.null(monitors) && !is.character(monitors)) stop("monitors argument should be a character vector.")
  if(!is.null(ntfyAccount) && !is.character(ntfyAccount)) stop("ntfyAcccount argument should be character.")
  if(!is.logical(email)) stop("email argument should be logical.")
  if(!is.logical(WAIC)) stop("WAIC argument should be logical.")
  if(!is.logical(HMC)) stop("HMC argument should be logical.")
  if(!is.null(replaceSamplers)){
    if(!is.list(replaceSamplers)) stop("replaceSamplers argument should be a list.")
  }
  if(!is.logical(parallel)) stop("parallel should be logical")
  if(!is.list(control.model)) stop("control.model should be a list.")
  if(!is.null(control.compile) && !is.list(control.compile)) stop("control.compile should be either NULL or a list.")
  if(!is.null(control.build) && !is.list(control.build)) stop("control.build should be either NULL or a list.")
  if(!is.null(control.configure) && !is.list(control.configure)) stop("control.configure should be either NULL or a list.")
  #XXXXXX ADDITIONAL CHECKS TO CONTROLS
}

runParallel <- function(seed, inits = NULL, control.model, control.compile, 
                        control.configure, control.build, HMC, WAIC, replaceSamplers, ...) {
  
  library(nimble)
  if(HMC) library(nimbleHMC)
  set.seed(seed)
  
  # ESTA LINEA SE DEBERIA QUITAR CUANDO LAS RUTINAS SE HAGAN LIBRERIA
  source("https://raw.githubusercontent.com/MigueBeneito/pNimble/refs/heads/main/RutinasNimble.0.3.R")
  #source("~/Trabajo/pNimble/RutinasNimble.0.2.R")
  load.leroux()
  load.SPDE()
  if(!is.null(inits)) control.model$inits <- inits()
  model.nimble <- do.call(nimbleModel, c(control.model, buildDerivs = HMC), quote = TRUE)
  model.precompiled <- do.call(compileNimble, c(model.nimble,control.compile), quote = TRUE)
  if(WAIC) control.configure$monitors <- unique(c(control.configure$monitors,
                                                  model.nimble$getParents(model.nimble$getNodeNames(dataOnly = TRUE), 
                                                                          stochOnly = TRUE)))
  if(HMC) model.configure <- do.call(configureHMC, c(model.nimble, control.configure), quote = TRUE)
  else model.configure <- do.call(configureMCMC, c(model.nimble, control.configure), quote = TRUE)

  if(!is.null(replaceSamplers)){
    for(i in 1:length(replaceSamplers[[1]])){ 
      theith<-lapply(replaceSamplers,function(y,j=i){y[[j]]})
      theith<-Filter(Negate(anyNA),theith)
      do.call(model.configure$replaceSamplers,theith)
      #model.configure$replaceSampler(target = replaceSamplers[[1]][i],
      #                               type = ifelse(!is.na(replaceSamplers[[2]][i]),replaceSamplers[[2]][i],NULL),
      #                               control = ifelse(!is.na(replaceSamplers[[3]][[i]]),replaceSamplers[[3]][[i]],NULL))
    }
    model.configure$printSamplers()
  }

  model.build <- do.call(buildMCMC, c(model.configure, control.build), quote = TRUE)
  model.compiled <- do.call(compileNimble, c(model.build, control.compile), quote =TRUE)
  model.output <- runMCMC(mcmc = model.compiled, nchains = 1, setSeed = seed, 
                          progressBar = FALSE, samplesAsCodaMCMC = TRUE, summary = FALSE, 
                          WAIC = FALSE, ...)
  return(model.output)
}

notify <- function(time, email, ntfyAccount, model){
  if(time<300){
    cat(paste0("Elapsed time: ",round(time,2)," seconds.\n"))
    theMessage <- paste0("Nimble has just finished ",format(Sys.time(), "(%b %d, %X)")," to run model ", model,
                       " with a total running time of ",round(time,2)," seconds.")
  }
  else{
    if(time<3600){
      cat(paste0("Elapsed time: ",round(time/60,2)," minutes.\n"))
      theMessage <- paste0("Nimble has just finished ",format(Sys.time(), "(%b %d, %X)")," to run model ", model,
                         " with a total running time of ",round(time/60,2)," minutes.")
    }
    else{
      cat(paste0("Elapsed time: ",round(time/3600, 2)," hours.\n"))
      theMessage <- paste0("Nimble has just finished ",format(Sys.time(), "(%b %d, %X)")," to run model ", model,
                         " with a total running time of ",round(time/3600,2)," hours.")
    }
  }
  if(is.null(ntfyAccount) & email==TRUE) cat("No email account to notify. No notification has been sent.")
  if(!is.null(ntfyAccount)){
    if(email){
      resul <- try(ntfy::ntfy_send(message = theMessage, title="Nimble has finished", topic="NimbleTools", 
                      email=ntfyAccount))
      if(inherits(resul,"try-error"))
        cat("Ntfy error. Possibly the email account is invalid or the limit for the daily number of emails has been reached. Try the option email=FALSE and use the ntfy app.")
    }
    else{
      resul <- try(ntfy::ntfy_send(message = theMessage, title="Nimble has finished", topic = ntfyAccount))
      if(inherits(resul,"try-error"))
        cat("Ntfy error. Possibly the subscribed topic in Ntfy is invalid or the limit for the daily number of push messages to the ntfy app has been reached. Check these issues.")
    }
  }
}

# Leroux distribution code for Nimble
##################################################
load.leroux<-function(){
  deregisterDistributions("dcar_leroux")
  
  dcar_leroux <- nimbleFunction(
    name = 'dcar_leroux',
    run = function(x = double(1),        # Spatial random effect (vector)
                 rho = double(0),      # Amount of spatial dependence (scalar)
                 sd.theta = double(0), # Standard deviation (scalar)
                 Lambda = double(1),   # Eigenvalues of matrix D - W
                 from.to = double(2),  # Matrix of distinct pairs of neighbors from.to[, 1] < from.to[, 2]
                 zero_mean = double(0, default = 0), 
                 log = integer(0, default = 0)) {
      #returnType(double(0))
      # Number of small areas
      NMuni <- dim(x)[1]
      # Number of distinct pairs of neighbors
      NDist <- dim(from.to)[1]
      # Required vectors
      from <- nimInteger(NDist)
      to <- nimInteger(NDist)
      x.from <- nimNumeric(NDist)
      x.to <- nimNumeric(NDist)
      for (Dist in 1:NDist) {
        from[Dist] <- ADbreak(from.to[Dist, 1])
        to[Dist] <- ADbreak(from.to[Dist, 2])
        x.from[Dist] <- x[from[Dist]]
        x.to[Dist] <- x[to[Dist]]
      }
    
      logDens <- sum(dnorm(x[1:NMuni], mean = 0, sd = sd.theta * pow(1 - rho, -1/2), log = TRUE)) -
        NMuni/2 * log(1 - rho) +  sum(log(rho * (Lambda[1:NMuni] - 1) + 1))/2 -
        pow(sd.theta, -2) * rho * sum(pow(x.from[1:NDist] - x.to[1:NDist], 2))/2
      
      if(zero_mean == 1){
        logDens <- logDens -1/2 * log(2*pi) + log(10000/(sd.theta^2)) - 5000*(mean(x[1:NMuni])/sd.theta)^2
      }
      
      returnType(double(0)) # return type
      if(log){return(logDens)} else{return(exp(logDens))}
    },
    buildDerivs = list(run = list(ignore = c('Dist','NMuni','NDist')))
  )

  rcar_leroux <- nimbleFunction(
    name = 'rcar_leroux',
    run = function(n = integer(0),
                   rho = double(0),
                   sd.theta = double(0),
                   Lambda = double(1),
                   from.to = double(2),
                   zero_mean = double(0, default = 0)) {
      returnType(double(1))
      nimStop("user-defined distribution dcar_leroux provided without random generation function.")
      x <- nimNumeric(length(Lambda))
      return(x)
    }
  )

  assign('dcar_leroux', dcar_leroux, envir = .GlobalEnv)
  assign('rcar_leroux', rcar_leroux, envir = .GlobalEnv)
}

load.SPDE <- function(){
  deregisterDistributions("dSPDE")
  
  dSPDE <- nimbleFunction(
    name = "dSPDE",
    run = function(x = double(1), tau = double(0), kappa = double(0), vecLambda = double(1), diagC = double(1), 
                   diagG1 = double(1), G1.out = double(2), diagG2 = double(1), G2.out = double(2),
                   zero_mean = double(1), log = integer(0, default = 0)) {
      sizeSPDE <- length(diagC)
      dim.G1.out <- dim(G1.out)[1]
      dim.G2.out <- dim(G2.out)[1]
      
      diagCNoD <- nimNumeric(sizeSPDE)
      diagG1NoD <- nimNumeric(sizeSPDE)
      diagG2NoD <- nimNumeric(sizeSPDE)
      vecLambdaNoD <- nimNumeric(sizeSPDE)
      for(i in 1:sizeSPDE){
        diagCNoD[i] <- ADbreak(diagC[i])
        diagG1NoD[i] <- ADbreak(diagG1[i])
        diagG2NoD[i] <- ADbreak(diagG2[i])
        vecLambdaNoD[i] <- ADbreak(vecLambda[i])
      }
      ldetC <- sum(log(diagCNoD[1:sizeSPDE]))
      
      x.idx.from.to_iG1 <- nimNumeric(dim.G1.out)
      G1.from <- nimInteger(dim.G1.out)
      G1.to <- nimInteger(dim.G1.out)
      G1.val <- nimNumeric(dim.G1.out)
      for (i in 1:dim.G1.out) {
        G1.from[i] <- ADbreak(G1.out[i,1])
        G1.to[i] <- ADbreak(G1.out[i,2])
        G1.val[i] <- ADbreak(G1.out[i,3])
        x.idx.from.to_iG1[i] <- x[G1.from[i]] * x[G1.to[i]] * G1.val[i]
      }
      
      x.idx.from.to_iG2 <- nimNumeric(dim.G2.out)
      G2.from <- nimInteger(dim.G2.out)
      G2.to <- nimInteger(dim.G2.out)
      G2.val <- nimNumeric(dim.G2.out)
      for (i in 1:dim.G2.out) {
        G2.from[i] <- ADbreak(G2.out[i,1])
        G2.to[i] <- ADbreak(G2.out[i,2])
        G2.val[i] <- ADbreak(G2.out[i,3])
        x.idx.from.to_iG2[i] <- x[G2.to[i]] * x[G2.from[i]] * G2.val[i]
      }
      
      log_pdf <- - sizeSPDE/2*log(2*pi) + # (2*pi)^(k/2)
        sizeSPDE*log(tau) + 1/2*ldetC + sum(log(kappa^2 + vecLambdaNoD[1:sizeSPDE])) - # determinant |Q| = tau^2 * |C^(1/2)|^2 * |kappa^2 * I + Lambda|^2 
        1/2*tau^2 * #(kappa^4 * 
        #   inprod(x[1:sizeSPDE]^2, diagCNoD[1:sizeSPDE]) + # exponential term of the Multivariate Gaussian
        #   2 * kappa^2 * ( inprod(x[1:sizeSPDE]^2, diagG1NoD[1:sizeSPDE]) +  2 * sum(x.idx.from.to_iG1[1:dim.G1.out])) + 
        #   inprod(x[1:sizeSPDE]^2, diagG2NoD[1:sizeSPDE]) + 2 * sum(x.idx.from.to_iG2[1:dim.G2.out])
        (inprod(x[1:sizeSPDE]^2, kappa^4 * diagCNoD[1:sizeSPDE] + 2 * kappa^2 * diagG1NoD[1:sizeSPDE] + diagG2NoD[1:sizeSPDE]) + # exponential term of the Multivariate Gaussian
           2 * 2 * kappa^2 * sum(x.idx.from.to_iG1[1:dim.G1.out]) +  2 * sum(x.idx.from.to_iG2[1:dim.G2.out]))
      if(sum(zero_mean)==0){}
      else{
        x.mean <- inprod(zero_mean[1:sizeSPDE], x[1:sizeSPDE])/sum(zero_mean[1:sizeSPDE])
        log_pdf <- log_pdf -1/2 * log(2*pi) + log(10000*sizeSPDE*tau) - 5000*sizeSPDE*tau*(x.mean^2)
      }
      returnType(double(0)) # return type
      if(log){return(log_pdf)} else{return(exp(log_pdf))}
    },
    buildDerivs = list(run = list(ignore = c('i','sizeSPDE','dim.G1.out','dim.G2.out','G1.val', 'G2.val','diagCnoD')))
    #,'diagG1NoD','diagG2NoD','vecLambdaNoD','x.mean','ldetC'))
  )
  
  rSPDE <- nimbleFunction(
    name = "rSPDE",
    run = function(n = integer(0), tau = double(0), kappa = double(0), vecLambda = double(1), diagC = double(1),
                   diagG1 = double(1), G1.out = double(2),  diagG2 = double(1),  G2.out = double(2),
                   zero_mean = double(1)
    ) {
      returnType(double(1))
      nimStop("user-defined distribution dSPDE provided without random generation function")
      x <- nimNumeric(length = length(diagC))
      return(x)
    }
  )
  
  assign('dSPDE', dSPDE, envir = .GlobalEnv)
  assign('rSPDE', rSPDE, envir = .GlobalEnv)  
}







