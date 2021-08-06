#########################################################################
# 2.SC_analysis.R
# SCR book Ch 18 code;
# created by Joanna Burgar, 07-May-2018
# updated by Joanna Burgar, 25-Jul-2021 for Leroy Gonzales (streamlined)
# script for estimating density of unmarked koala population
# using JAGS with weakly, vaguely, and strongly informed priors
#########################################################################

library(coda)
library(rjags)
library(parallel)

#############################################################
###--- WRITE JAGS MODELS
#############################################################

# uninformative prior: KUI1 = dunif(0,100)
# weakly informative prior 1: HR 4-863 ha; KWI1 = dgamma = 4,1
# weakly informative prior 2: 10-80 ha; KWI2 = dgamma = 10,4.6
# strongly informative prior 1: 40 ha; KSI1 = dgamma = 60,27

###--- Model 1 = Uninformative prior

# specify model
cat("
   model {
  sigma ~ dunif(0,100) # uninformative prior
  lam0 ~ dunif(0,100)
  psi ~ dbeta(1,1)

  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])

    for(j in 1:J) {
      distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
      lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i]
    }
  } # End of 1:M
  for(j in 1:J) {
    for(k in 1:K) {
      bigLambda[j, k] <- sum(lam[1:M, j]) * oper[j, k]
      n[j, k] ~ dpois(bigLambda[j, k])

    }
  }
  N <- sum(z[])
  D <- N/area
}
    ",file="KUI1.txt")

#############################################################
###--- Model 2 =  weakly informative prior 1 (same for all study areas)
# weakly informative HR 4-863; KWI1 = dgamma = 4,1
# specify model
cat("
    model {
    sigma ~ dgamma(4,1)
    lam0 ~ dunif(0,100)
    psi ~ dbeta(1,1)

    for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])

    for(j in 1:J) {
    distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
    lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i]
    }
    } # End of 1:M
    for(j in 1:J) {
    for(k in 1:K) {
    bigLambda[j, k] <- sum(lam[1:M, j]) * oper[j, k]
    n[j, k] ~ dpois(bigLambda[j, k])

    }
    }
    N <- sum(z[])
    D <- N/area
    }
    ",file="KWI1.txt")

#############################################################

###--- Model 3 =  weakly informative prior 2 (specific to this area)
# weakly informative 10-90; KWI2 = dgamma = 10,4.6
# specify model
cat("
    model {
    sigma ~ dgamma(10,4.6)
    lam0 ~ dunif(0,100)
    psi ~ dbeta(1,1)

    for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])

    for(j in 1:J) {
    distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
    lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i]
    }
    } # End of 1:M
    for(j in 1:J) {
    for(k in 1:K) {
    bigLambda[j, k] <- sum(lam[1:M, j]) * oper[j, k]
    n[j, k] ~ dpois(bigLambda[j, k])

    }
    }
    N <- sum(z[])
    D <- N/area
    }
    ",file="KWI2.txt")
#############################################################

###--- Model 4 =  strongly informative prior 1
# strongly informative 40; KSI1 = dgamma = 60,27
# specify model
cat("
    model {
    sigma ~ dgamma(60, 27)
    lam0 ~ dunif(0,100)
    psi ~ dbeta(1,1)

    for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])

    for(j in 1:J) {
    distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
    lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i]
    }
    } # End of 1:M
    for(j in 1:J) {
    for(k in 1:K) {
    bigLambda[j, k] <- sum(lam[1:M, j]) * oper[j, k]
    n[j, k] ~ dpois(bigLambda[j, k])

    }
    }
    N <- sum(z[])
    D <- N/area
    }
    ",file="KSI1.txt")


#####################################################################
###--- FUNCTION TO RUN 3 CHAINS IN PARALLEL USING JAGS IMPLEMENTATION
# run the model without the extra s and z parameters (i.e., for results but NOT for map)
fn.run_JAGS_model <- function(StudyArea=StudyArea, dat=dat, init=init, pars=pars, model=model, 
                              n.chains=n.chains, n.adapt=n.adapt, n.iter=n.iter){
  (start.time <- Sys.time())
  cl3 <- makeCluster(3)
  clusterExport(cl3, c("dat","init","pars", "model", "n.chains", "n.adapt", "n.iter"))
  KJAGS <- clusterEvalQ(cl3, {
    library(rjags)
    jm2 <- jags.model(paste(model,".txt", sep=""), data=dat, inits=init, n.chains=n.chains, n.adapt=n.adapt)
    jc2 <- coda.samples(jm2, pars, n.iter=n.iter)    # initial run to minimise size
    return(as.mcmc(jc2))
  })
  mc.KJAGS <- mcmc.list(KJAGS)
  
  (end.time <- Sys.time())
  mc.KJAGS.ET <- difftime(end.time, start.time, units='hours')
  print(mc.KJAGS.ET)
  
  stopCluster(cl3)  
  return(mc.KJAGS)
  
}


# run model saving s and z parameters - much larger output, only do when decided on prior
fn.run_JAGS_model_map <- function(StudyArea=StudyArea, dat=dat, init=init, pars_map=pars_map, model=model, 
                              n.chains=n.chains, n.adapt=n.adapt, n.iter=n.iter){
  (start.time <- Sys.time())
  cl3 <- makeCluster(3)
  clusterExport(cl3, c("dat","init","pars_map", "model", "n.chains", "n.adapt", "n.iter"))
  KJAGS <- clusterEvalQ(cl3, {
    library(rjags)
    jm2 <- jags.model(paste(model,".txt", sep=""), data=dat, inits=init, n.chains=n.chains, n.adapt=n.adapt)
    jc2 <- coda.samples(jm2, pars_map, n.iter=n.iter) # to create density map
    return(as.mcmc(jc2))
  })
  mc.KJAGS <- mcmc.list(KJAGS)
  
  (end.time <- Sys.time())
  mc.KJAGS.ET <- difftime(end.time, start.time, units='hours')
  print(mc.KJAGS.ET)
  
  stopCluster(cl3)  
  return(mc.KJAGS)
  
}
#####################################################################

# specify parameters to monitor
pars <- c("sigma","lam0","psi","D","N")
pars_map <- c("sigma","lam0","psi","D","N","s","z")

# specify remaining model paramters
n.chains = 1
n.adapt = 1000
n.iter = 50000

#############################################################
# Load your data [change the files paths to your data effations]
StudyArea <- "BCK"
# StudyArea <- "BELK"

Year <- 2020

load(paste("./out/",StudyArea,"_",Year,"_dat.RData",sep="")) # load output created from Koala_SC_modelling.Rmd

###--- Check appropriate data loaded
str(dat)
# specify initial values
init <-  function() {  list(sigma=rnorm(1,10), lam0=runif(1) , z=rep(1,dat$M)) }

###--- Run uninformative model
model = "KUI1"
StudyArea_KUI1 <- fn.run_JAGS_model(StudyArea = StudyArea, model = model)
save(StudyArea_KUI1, file=paste(StudyArea,"_",Year,"_",model,".RData", sep=""), compress=TRUE)

###--- Run regional weakly informative model; weakly informative HR 4-863; WIA1 = dgamma = 4,1
model = "KWI1"
StudyArea_KWI1 <- fn.run_JAGS_model(StudyArea = StudyArea, model = model)
save(StudyArea_KWI1, file=paste(StudyArea,"_",Year,"_",model,".RData", sep=""), compress=TRUE)

###--- Run study area weakly informative model; weakly informative 10-90; WI2 = dgamma = 10,4.6
model = "KWI2"
StudyArea_KWI2 <- fn.run_JAGS_model(StudyArea = StudyArea, model = model)
save(StudyArea_KWI2, file=paste(StudyArea,"_",Year,"_",model,".RData", sep=""), compress=TRUE)

###--- Run strongly informative model; weakly informative HR 4-863; WIA1 = dgamma = 4,1
model = "KSI1"
StudyArea_KSI1 <- fn.run_JAGS_model(StudyArea = StudyArea, model = model)
save(StudyArea_KSI1, file=paste(StudyArea,"_",Year,"_",model,".RData", sep=""), compress=TRUE)

#############################################################
# run model to create map
model = "KSI1" # assuming going with strongest prior for final model and density map - if not change here and in code below
StudyArea_KSI1_map <- fn.run_JAGS_model_map(StudyArea = StudyArea, model = model, pars_map = pars_map)
save(StudyArea_KSI1_map, file=paste(StudyArea,"_",Year,"_",model,"_map.RData", sep=""), compress=TRUE)
#############################################################