library(ntfy)
library(coda)
library(parallel)

pNimble<-function(code = NULL, data = NULL, constants = NULL, inits = NULL, nchains = 3, 
                  seeds = NULL, summary = TRUE, monitors =FALSE,
                  email = NULL, WAIC = FALSE, control.model = NULL,
                  control.compile = NULL, control.configure = NULL, 
                  control.build = NULL, ...){

  time.start <- Sys.time()
  
  runParallel <- function(seed, inits = NULL, control.model, control.compile, 
                          control.configure, control.build, ...) {
    
    library(nimble)
    set.seed(seed)
    
    source("RutinasNimble.R")
    load_leroux()
    
    if(!is.null(inits)) control.model$inits <- inits()
    model.nimble <- do.call(nimbleModel, control.model, quote = TRUE)
    if(WAIC) control.configure$monitors <- unique(c(control.configure$monitors,
                                                    model.nimble$getParents(model.nimble$getNodeNames(dataOnly = TRUE), 
                                                                       stochOnly = TRUE)))
    model.precompiled <- do.call(compileNimble, c(model.nimble,control.compile), quote = TRUE)
    model.configure <- do.call(configureMCMC, c(model.nimble, control.configure), quote = TRUE)
    model.build <- do.call(buildMCMC, c(model.configure, control.build), quote = TRUE)
    model.compiled <- do.call(compileNimble, c(model.build, control.compile), quote =TRUE)
    model.output <- runMCMC(mcmc = model.compiled, nchains = 1, setSeed = seed, 
                            progressBar = FALSE, samplesAsCodaMCMC = TRUE, summary = FALSE, 
                            WAIC = FALSE, ...)
    return(model.output)
  }

  if(is.null(control.model)) control.model <- list()
  if(!is.null(code)) control.model$code <- code
  if(!is.null(data)) control.model$data <- data
  if(!is.null(constants)) control.model$constants <- constants
  if(!is.null(monitors)){
    if(is.null(control.configure)){control.configure <- list()}
    control.configure$monitors=monitors
  }

  if(is.null(seeds)) seeds <- 1:nchains
  
  my_cluster <- makeCluster(nchains)
  on.exit(stopCluster(my_cluster))
  resul <- parLapply(cl = my_cluster, X = seeds, fun = runParallel,
                     inits = inits, control.model = control.model, 
                     control.compile = control.compile, control.configure = control.configure, 
                     control.build = control.build, ...
                     )

  names(resul) <- paste0("chain",1:length(resul))
  resul2 <- list()
  resul2$samples <- coda::as.mcmc.list(lapply(resul,function(x){coda::as.mcmc(x)}))
  if(summary) resul2$summary <- MCMCvis::MCMCsummary(resul2$samples)

  if(WAIC){
    model.nimble <- nimbleModel(code = code, data = data, constants = constants,
                                check = FALSE, calculate = FALSE, buildDerivs = FALSE)
    cmodel <- compileNimble(model.nimble)
    resul2$WAIC <- calculateWAIC(do.call(rbind,resul2$samples), cmodel)
  } 

  time.end <-Sys.time()
  
  total.time <- as.numeric(time.end)-as.numeric(time.start)
  
  if(total.time<300){
    cat("Elapsed time: ",round(total.time,2)," seconds.\n")
    if(!is.null(email))
      
      ntfy::ntfy_send(message = paste0("Nimble has just finished to run a model with a total running time of ",
                                       round(total.time,2)," seconds."),
                      title="Nimble has finished", topic="NimbleTools", email=email)
  }
  else{
    if(total.time<3600){
      cat("Elapsed time: ",round(total.time/60,2)," minutes.\n")
      if(!is.null(email))
        ntfy::ntfy_send(message = paste0("Nimble has just finished to run a model with a total running time of ",
                                         round(total.time/60,2)," minutes."),
                        title="Nimble has finished", topic="NimbleTools", email=email)
    }
    else{
      cat("Elapsed time: ",round(total.time/3600, 2)," hours.\n")
      if(!is.null(email))
        ntfy::ntfy_send(message = paste0("Nimble has just finished to run a model with a total running time of ",
                                         round(total.time/3600,2)," hours."),
                        title="Nimble has finished", topic="NimbleTools", email=email)
    }
  }
  
  return(resul2)
}

# Definición de distribución de Leroux para Nimble
##################################################
load_leroux<-function(){
  dcar_leroux <- nimbleFunction(
    name = 'dcar_leroux',
    run = function(x = double(1),        # Spatial random effect (vector)
                 rho = double(0),      # Amount of spatial dependence (scalar)
                 sd.theta = double(0), # Standard deviation (scalar)
                 Lambda = double(1),   # Eigenvalues of matrix D - W
                 from.to = double(2),  # Matrix of distinct pairs of neighbors from.to[, 1] < from.to[, 2]
                 log = integer(0, default = 0)) {
      returnType(double(0))
      # Number of small areas
      NMuni <- dim(x)[1]
      # Number of distinct pairs of neighbors
      NDist <- dim(from.to)[1]
      # Required vectors
      x.from <- nimNumeric(NDist)
      x.to <- nimNumeric(NDist)
      for (Dist in 1:NDist) {
        x.from[Dist] <- x[from.to[Dist, 1]]
        x.to[Dist] <- x[from.to[Dist, 2]]
      }
    
      logDens <- sum(dnorm(x[1:NMuni], mean = 0, sd = sd.theta * pow(1 - rho, -1/2), log = TRUE)) -
        NMuni/2 * log(1 - rho) + 1/2 * sum(log(rho * (Lambda[1:NMuni] - 1) + 1)) -
        1/2 * pow(sd.theta, -2) * rho * sum(pow(x.from[1:NDist] - x.to[1:NDist], 2))
      if(log) return(logDens)
      else return(exp(logDens))
    }
  )

  rcar_leroux <- nimbleFunction(
    name = 'rcar_leroux',
    run = function(n = integer(0),
                   rho = double(0),
                   sd.theta = double(0),
                   Lambda = double(1),
                   from.to = double(2)) {
      returnType(double(1))
      nimStop("user-defined distribution dcar_leroux provided without random generation function.")
      x <- nimNumeric(length(Lambda))
      return(x)
    }
  )

  assign('dcar_leroux', dcar_leroux, envir = .GlobalEnv)
  assign('rcar_leroux', rcar_leroux, envir = .GlobalEnv)
}