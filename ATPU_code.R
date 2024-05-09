library(dplyr)
library(tidyverse)
library(jagsUI)
library(MCMCvis)
library(IPMbook)
library(R2ucare)
library(tidyr)


# Load productivity data
# Number of surveyed broods
ATPU_B = read.csv("ATPU_B.csv")

# Number of fledglings (successes)
ATPU_J = read.csv("ATPU_J.csv")

# Load count data
ATPU_count=  read.csv("ATPU_count.csv")

# Load CMR data (encounter histories)
ATPU_EH = read.csv("ATPU_EH.csv")

# Load # times resighted as covariate
cov = read.csv("ATPU_cov.csv")


# Create encounter history matrix to be used in model
ATPU_df = ATPU_EH  %>% select("BirdID","Contact_Year","Age") %>% mutate(seen=1)
ch=ATPU_df %>% 
  group_by(BirdID) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(id_cols=-row,names_from = Contact_Year, values_from = seen,values_fill = 0) 

# Rearrange columns so years are in correct order and convert to matrix format
ch<-as.matrix(ch[,order(colnames(ch),decreasing=FALSE)])

# Filter birds to only include covariate information for same birds in capture histories
cov <- filter(cov, BirdID %in% ch)

# Join data together into a list
data_list <- list(cov, ch)                         

# Apply reduce function to join by BirdID
data = data_list %>% reduce(full_join, by = "BirdID", copy=TRUE)     

# Separate out EH and covariate
# For juveniles
j_data <- sapply(data[data[,"Age"]=="1",], format, trim = TRUE)
j_cov <- j_data[,3:24]
j_ch <- j_data[,25:46]

# For adults
a_data <- sapply(data[data[,"Age"]=="2",], format, trim = TRUE)
a_cov <- a_data[,3:24]
a_ch <- a_data[,25:46]


# Ensure that all encounters in cov matrix are also in capture histories 
# (some ATPU resights made it into # encounters but not into original encounter histories)
# Looks like this was an issue with original encounter history data when bird was reencountered
## on more than just MSI in a single year (some MSI reencounters not included in EH data)
ch_updated <- function(cov, ch){
  x <- ch
  for (i in 1:nrow(cov)) {
    for (j in 1:ncol(cov)) {
      if (cov[i,j] != 0) {x[i,j] <- 1}
    }
  }
  z <- as.matrix(x)
  return(x)
}

# Run function for both age groups
j_ch_updated = ch_updated (j_cov,j_ch)
a_ch_updated = ch_updated (a_cov,a_ch)


# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f.j <- apply(j_ch_updated, 1, get.first)
f.a <- apply(a_ch_updated, 1, get.first)

# Compute vector with occasion of last capture
get.last <- function(x) max(which(x!=0))
l.j <- apply(j_ch_updated, 1, get.last)
l.a <- apply(a_ch_updated, 1, get.last)


#Create matrix with number of encounters per occasion (Negative 1 (-1) at first capture, 0 before first capture and at non-captures
recode_cov <- function(cov, f){
  x <- cov
  for (i in 1:nrow(cov)) {
    for (j in 1:ncol(cov)) {
      if (j == (f[i])) {x[i,j] <- -1}
      if (j < (f[i])) {x[i,j] <- 0}
    }
  }
  z <- as.matrix(x)
  return(x)
}
# Run function for both age groups
c.j = recode_cov(j_cov,f.j)
c.a = recode_cov(a_cov,f.a)


#Create age matrix covariate for known aged individuals
# Age 4 = average age of return, using this as break for juveniles/adults
# Age covariate 1 up until third year, then 4th year and after are age covariate 2
recode_ch_j <- function(ch, f){
  z <- ch
  for (i in 1:nrow(ch)) {
    for (j in 1:ncol(ch)) {
      if (j < (f[i])){z[i,j] <- NA}
      if (j == (f[i])){z[i,j] <- 1}
      if (j == (f[i]+1)){z[i,j] <- 2}
      if (j == (f[i]+2)){z[i,j] <- 3}
      if (j >= (f[i]+3)){z[i,j] <- 4}      
    }
  }
  z <- as.matrix(z)
  return(z)
}
# Run function
age.j = recode_ch_j(j_ch_updated,f.j)

#Create age matrix covariate for unknown aged individuals
# Always age covariate 2 
recode_ch_a <- function(ch, f){
  z <- ch
  for (i in 1:nrow(ch)) {
    for (j in 1:ncol(ch)) {
      if (j < (f[i])){z[i,j] <- NA}
      if (j >= (f[i])){z[i,j] <- 4}
    }
  }
  z <- as.matrix(z)
  return(z)
}
# Run function
age.a = recode_ch_a(a_ch_updated,f.a)


#Create known alive matrix (Ones where known alive, NAs elsewhere)
known.state.cjs <- function(ch, f,l){
  z <- ch
  for (i in 1:nrow(ch)) {
    for (j in 1:ncol(ch)) {
      if (j < (f[i])){z[i,j] <- NA}
      if (j == (f[i])){z[i,j] <- 1}
      if (j == (l[i])){z[i,j] <- 1}
      if (j > (f[i]) & j<(l[i])){z[i,j] <- 1}
      
    }
  }
  z <- as.matrix(z)
  return(z)
}

# Run function for both age groups
z.dat.j = known.state.cjs(j_ch_updated,f.j,l.j)
z.dat.a = known.state.cjs(a_ch_updated,f.a,l.a)

### Known available matrix (Ones where captured only)###
# Same as raw encounter histories
a.dat.j <- j_ch_updated
a.dat.a <- a_ch_updated

# Bring the data back together
ch <- rbind(j_ch_updated,a_ch_updated)
f <- c(f.j, f.a)

c <- rbind(c.j, c.a)
c <- apply(as.matrix(c[,]), 2, as.numeric)

age <- rbind(age.j, age.a)
age <- apply(as.matrix(age[,]), 2, as.numeric)

a.dat <- rbind(a.dat.j, a.dat.a)
a.dat <- apply(as.matrix(a.dat[,]), 2, as.numeric)

z.dat <- rbind(z.dat.j, z.dat.a)
z.dat <- apply(as.matrix(z.dat[,]), 2, as.numeric)


########################################################################################
# Integrated population for Atlantic Puffins breeding on Machias Seal Island (1995-2019)
# Model description
## Age structured model (First year, second year, third year, and fourth year and after)
### Two age classes for survival (first year different than remaining years)
### Four age classes for availability for detection (gamma; to take into account temporary emigration of juveniles)
# Individual heterogeneity in detection incorporated through overdispersion parameter (theta)
# Age at first breeding = 4 years
# Assuming a pre-breeding census
# Female only model (same as number of breeding pairs)
# All vital rates assumed to be time-dependent (except gamma and omega)
## Random time effect included on survival, fecundity, and mean expected number of time detected (epsilon)
# Immigration rate plus noise parameter included to improve fit
########################################################################################

# Write JAGS model file
cat(file="ATPU_IPM.txt", "
model { 
# Priors and constraints
# Constrain survival to be age- and time-specific
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- alpha[age[i,t],t]  
      } #t
   } #i


# Making age- and time-specific survival random variable (and defining priors, standard deviation, and precision)
for (u in 1:2){
  sigma.alpha[u] ~ dgamma(1,1)
  tau.alpha[u] = 1/(sigma.alpha[u] * sigma.alpha[u])
  sigma.alpha2[u] = sigma.alpha[u] * sigma.alpha[u]
  
for (t in 1:(n.occasions-1)){
   alpha.star[u,t] ~ dnorm(0, 1)T(-5,5)      # Priors for age and time specific survival
   logit(alpha[u,t]) = mu.phi[u] + alpha.star[u,t] * sigma.alpha[u]
   }
}

# Contrain survival after first year to be same as adults
# Survival after first year of life more similar to adult survival rate
for (t in 1:(n.occasions-1)){
   alpha[3,t] = alpha[2,t]      
   alpha[4,t] = alpha[2,t]      
   }

# Defining time-specific fecundity as random variable
for (t in 1:n.occasions){
  logit(fecund[t]) = mu.fec +  eps.fec[t]   
  eps.fec[t] ~ dnorm(0,tau.fec)T(-5,5)
}


# Defining time-specific epsilon as random variable
# Constraining juveniles (first three years) to be different than adults
# Due to temporary emigration (would expect to see them less than adults)
for (t in 1:n.occasions){
  log(epsilon[1,t]) = mu.eps[1] + eps[t]
  log(epsilon[2,t]) = mu.eps[1] + eps[t]
  log(epsilon[3,t]) = mu.eps[1] + eps[t]
  log(epsilon[4,t]) = mu.eps[2] + eps[t]
  eps[t] ~ dnorm(0, tau.eps)  
}

  
# Priors for gamma (availability for detection given absense/presence at last time step)  
  for (j in 1:4){
    gamma[1,j] ~ dbeta(1, 1)  # Prior for probability of availablity at time t given absense at time t-1
    gamma[2,j] ~ dbeta(1, 1) # Prior for probability of availablity at time t given available at time t-1
  }
  

  theta ~ dunif(0,250)     # Prior for individual heterogeneity parameter
  omega ~ dbeta(3, 7)      # Prior for immigration rate plus noise
  mu.phi[1] ~ dlogis(0,1)  # Prior for mean juvenile survival
  mu.phi[2] ~ dlogis(0,1)  # Prior for mean adult survival
  mu.eps[1] ~ dgamma(1,1)     # Prior for mean epsilon (expected # times detected)
  mu.eps[2] ~ dgamma(1,1)     # Prior for mean epsilon (expected # times detected)
  mu.fec ~ dbeta(1, 1)        # Prior for mean fecundity
  sigma.eps ~ dunif(0, 5)      # Prior on standard deviation of epsilon              
  tau.eps <- pow(sigma.eps , -2)    # Precision parameter for epsilon (1/variance)
  sigma.eps2 <- pow(sigma.eps, 2)  #Temporal variance in epsilon
  sigma.fec ~ dunif(0, 5)      # Prior on standard deviation of fecundity            
  tau.fec <- pow(sigma.fec , -2)    # Precision parameter for fecundity (1/variance)
  sigma.fec2 <- pow(sigma.fec, 2)  # Temporal variance in fecundity


############################################################################
################### Population model #######################################
############################################################################

# Population count data (state-space model)
# Model for the initial stage-specific population sizes

N1  ~ dnorm(2000,0.01)I(0,)
N[1,1] <- round(N1)
N2  ~ dnorm(2000,0.01)I(0,)
N[2,1] <- round(N2)
N3  ~ dnorm(2000,0.01)I(0,)
N[3,1] <- round(N3)
N4 ~ dnorm(8000,0.01)I(0,)
N[4,1] <- round(N4) 
nimm ~ dnorm(100, 0.01)I(0,)       
Nimm[1] <- round(nimm)

# Process model over time: our model of population dynamics
for (t in 1:(n.occasions-1)){
  N[1,t+1] ~ dbin((fecund[t] / 2) * alpha[1,t], Ntot[t])
  N[2,t+1] ~ dbin(alpha[2,t], N[1,t])
  N[3,t+1] ~ dbin(alpha[2,t], N[2,t])
  N[4,t+1] ~ dbin(alpha[2,t], N[3,t]+Ntot[t])
  mpo[t+1] <- (N[1,t]+N[2,t]+N[3,t]+N[4,t]) * omega
  Nimm[t+1] ~ dpois(mpo[t+1])
}

# Observation model
for (t in 1:n.occasions){
  Ntot[t] <- N[4,t] + Nimm[t]
  C[t] ~ dpois(Ntot[t])
  }

  
############################################################################
####### Productivity model #################################################
############################################################################
# Productivity data (Binomial GLM) 
for (t in 1:n.occasions){   
    J[t] ~ dbin(fecund[t], nbrood[t])    # number young observed as fledged
    }

############################################################################
################ CMR model #################################################
############################################################################
# Zero-inflated gamma-Poisson adapted from Riecke et al. (2022)
# To account for heterogeniety in observation process and overdispersion

# CMR Model Likelihood components
for (i in 1:nind){

h[i] ~ dgamma(theta, theta)

   for (t in (f[i]+1):n.occasions){
   
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
   
      # State process
      a[i,t] ~ dbern(z[i,t] * gamma[a[i,t-1]+1, age[i,t-1]])

      # Observation process
      c[i,t] ~ dpois(a[i,t] * epsilon[age[i,t],t] * h[i])
    }
}

# Derived parameters
# Detection probability
  for (u in 1:4){
  for (t in 1:n.occasions){
    p[u,t] = 1 - exp(-(epsilon[u,t]))  
  } 
  } 

# Annual population growth rate
for (t in 1:(n.occasions-1)){
  ann.growth.rate[t] <- (Ntot[t+1]) / (Ntot[t]+0.0001)   
}

mphij <- exp(mu.phi[1])/(1+exp(mu.phi[1]))   # Mean juvenile survival probability (age 1)
mphia <- exp(mu.phi[2])/(1+exp(mu.phi[2]))   # Mean adult survival probability (age 2-4+)
mfec <- exp(mu.fec)/(1+exp(mu.fec))   # Mean adult survival productivity
mepsj <- exp(mu.eps[1])  # Mean expected # encounters per juvenile (ages 1-3)
mepsa <- exp(mu.eps[2])  # Mean expected # encounters per adult (age 4+)


} # end model

")


# Parameters monitored
parameters <- c("alpha","mu.phi","mphij","mphia","sigma.alpha",
                "sigma.alpha2","gamma","p","fecund","mu.fec","mfec",
                "sigma.fec2", "epsilon","mu.eps","mepsj","mepsa","sigma.eps2",
                "omega","Ntot","N","Nimm","ann.growth.rate","theta","alpha.star","eps","eps.fec")

# Initial values
inits <- function(){list(mu.phi = c(0,2),mu.eps = c(0,0.4), mu.fec = 0)}

# Bring together data
ATPU_IPM_data <- list(z = z.dat, 
                           a = a.dat, 
                           f = f,
                           nind = dim(ch)[1], 
                           n.occasions = dim(ch)[2],
                           J=ATPU_J, 
                           nbrood=ATPU_B,
                           C=ATPU_count,
                           c=c,
                           age=age) 
# MCMC settings
ni <- 1000000; nb <- 250000; nc <- 3; nt <- 5; na <- 50000

# Call JAGS 
ATPU_IPM_out <- jags(ATPU_IPM_data, inits=inits, parameters, "ATPU_IPM.txt", n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel = TRUE)

# Print Results
print(ATPU_IPM_out, digits = 3)

# Examine traceplots to check for convergence
MCMCtrace(ATPU_IPM_out,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)


#############################################################################################
#############################################################################################
### Run CMR model only to compare estimates from full IPM
#############################################################################################
#############################################################################################
# Write JAGS model file
cat(file="ATPU_CMR.txt", "
model { 
# Priors and constraints
# Constrain survival to be age- and time-specific
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- alpha[age[i,t],t]  
      } #t
   } #i


# Making age- and time-specific survival random variable (and defining priors, standard deviation, and precision)
for (u in 1:2){
  sigma.alpha[u] ~ dgamma(1,1)
  tau.alpha[u] = 1/(sigma.alpha[u] * sigma.alpha[u])
  
for (t in 1:(n.occasions-1)){
   alpha.star[u,t] ~ dnorm(0, 1)       # Priors for age and time specific survival
   logit(alpha[u,t]) = mu.phi[u] + alpha.star[u,t] * sigma.alpha[u]
   }
}

# Contrain survival after first year to be same as adults
# Survival after first year of life more similar to adult survival rate
for (t in 1:(n.occasions-1)){
   alpha[3,t] = alpha[2,t]      
   alpha[4,t] = alpha[2,t]      
   }

# Defining time-specific epsilon as random variable
# Constraining juveniles (first three years) to be different than adults
# Due to temporary emigration (would expect to see them less than adults)
for (t in 1:n.occasions){
  log(epsilon[1,t]) = mu.eps[1] + eps[t]
  log(epsilon[2,t]) = mu.eps[1] + eps[t]
  log(epsilon[3,t]) = mu.eps[1] + eps[t]
  log(epsilon[4,t]) = mu.eps[2] + eps[t]
  eps[t] ~ dnorm(0, tau.eps)  
}


# Priors for gamma (availability for detection given absense/presence at last time step)  
  for (j in 1:4){
    gamma[1,j] ~ dbeta(1, 1)  # Prior for probability of availablity at time t given absense at time t-1
    gamma[2,j] ~ dbeta(1, 1) # Prior for probability of availablity at time t given available at time t-1
  }
  

  theta ~ dunif(0,250)     # Prior for individual heterogeneity parameter
  mu.phi[1] ~ dlogis(0,1)  # Prior for mean juvenile survival
  mu.phi[2] ~ dlogis(0,1)  # Prior for mean adult survival
  mu.eps[1] ~ dgamma(1,1)     # Prior for mean epsilon (expected # times detected)
  mu.eps[2] ~ dgamma(1,1)     # Prior for mean epsilon (expected # times detected)
  sigma.eps ~ dunif(0, 5)      # Prior for standard deviation of epsilon              
  tau.eps <- pow(sigma.eps , -2)    # Precision parameter for epsilon (1/variance)
  sigma.eps2 <- pow(sigma.eps, 2)  # Temporal variance in epsilon

############################################################################
################ CMR model #################################################
############################################################################
# Capture-recapture data (CJS model with multinomial likelihood)
# Define the multinomial likelihood

# CMR Model Likelihood components
for (i in 1:nind){

h[i] ~ dgamma(theta, theta)

   for (t in (f[i]+1):n.occasions){
   
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
   
      # State process
      a[i,t] ~ dbern(z[i,t] * gamma[a[i,t-1]+1, age[i,t-1]])

      # Observation process
      c[i,t] ~ dpois(a[i,t] * epsilon[age[i,t],t] * h[i])
    }
}

# Derived parameters
# Detection probability
  for (u in 1:4){
  for (t in 1:n.occasions){
    p[u,t] = 1 - exp(-(epsilon[u,t]))  
  } 
  }

mphij <- exp(mu.phi[1])/(1+exp(mu.phi[1]))   # Mean juvenile survival probability (age 1)
mphia <- exp(mu.phi[2])/(1+exp(mu.phi[2]))   # Mean adult survival probability (age 2-4+)
mepsj <- exp(mu.eps[1])  # Mean expected # encounters per juvenile (ages 1-3)
mepsa <- exp(mu.eps[2])  # Mean expected # encounters per adult (age 4+)
  

} # end model
")

# Parameters monitored
parameters <- c("alpha","mu.phi","mphij","mphia","gamma","p","epsilon","mu.eps","mepsj","mepsa","sigma.eps2","theta")

# Initial values
inits <- function(){list(mu.phi = c(0,2),mu.eps = c(1,1))}


ATPU_CMR_data <- list(z = z.dat, 
                               a = a.dat, 
                               f = f,
                               nind = dim(ch)[1], 
                               n.occasions = dim(ch)[2],
                               c=c,
                               age=age) 
# MCMC settings
ni <- 1000000; nb <- 250000; nc <- 3; nt <- 5; na <- 50000

# Call JAGS (ART 2 min), check convergence and summarize posteriors
ATPU_CMR_out <- jags(ATPU_CMR_data, inits=inits, parameters, "ATPU_CMR.txt", n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel = TRUE)

# Print Results
print(ATPU_CMR_out, digits = 3)

# Examine traceplots
MCMCtrace(ATPU_CMR_out,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)

#############################################################################################
#############################################################################################
### Run prod model only to compare estimates from full IPM
#############################################################################################
#############################################################################################

# Write JAGS model file
cat(file="ATPU_prod.txt", "
model { 

# Defining time-specific fecundity as random variable
for (t in 1:n.occasions){
  logit(fecund[t]) = mu.fec +  eps.fec[t]   
  eps.fec[t] ~ dnorm(0,tau.fec)
}

  mu.fec ~ dbeta(1, 1)        # Prior for mean fecundity
  sigma.fec ~ dunif(0, 5)      # Prior for standard deviation of fecundity            
  tau.fec <- pow(sigma.fec , -2)    # Precision parameter for fecundity (1/variance)
  sigma.fec2 <- pow(sigma.fec, 2)  # Temporal variance in fecundity
############################################################################
####### Productivity model #################################################
############################################################################

# Productivity data (Binomial regression model) and GOF test (GOF adapted from Schaub et al. 2015)
for (t in 1:n.occasions){   
    J[t] ~ dbin(fecund[t], nbrood[t])    # number young observed as fledged
    }

mfec <- exp(mu.fec)/(1+exp(mu.fec))   # Mean adult survival productivity


} # end model
")


# Parameters monitored
parameters <- c("fecund","sigma.fec2","mfec")

ATPU_prod_data <- list(n.occasions = length(ATPU_J),
                           J=ATPU_J, 
                           nbrood=ATPU_B)
# MCMC settings
ni <- 1000000; nb <- 250000; nc <- 3; nt <- 5; na <- 50000

# Call JAGS  check convergence and summarize posteriors
ATPU_prod_out <- jags(ATPU_prod_data, inits=NULL, parameters, "ATPU_prod.txt", n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel = TRUE)

# Print Results
print(ATPU_prod_out, digits = 3)

#Examine traceplots
MCMCtrace(ATPU_prod_out ,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)


################################################################################
################################################################################
################################################################################
# Transient Life-table response experiments
# Adapted from example 9.3 in IPM book (Schaub and Kery 2021)
# Accounts for population structure 
## (i.e., doesn't assume population is at stable distribution like LTRE)
# Focus is the realized population growth rate (i.e., ratio of pop sizes)
###############################################################################
# First type of tLTRE
## Decompose temporal variance of realized population growth rate
## into components due to the realized variability and covariation of demographic
## rates and population structure
##############################################################################

# Number of MCMC draws
n.draws <- ATPU_IPM_out$mcmc.info$n.samples
draws <- ATPU_IPM_out$sims.list

# Define lambda
lambda <- expression(((n4 + n5) * 0.5 * phij * f + phia * n1 + phia * n2 + phia * (n3+n4+n5)) / (n1 + n2 + n3 + n4 + n5))             # Define lambda

D(lambda, "phij") # Get derivative for juvenile apparent survival

# Calculate proportional population sizes
n.years <- length(ATPU_IPM_out$mean$alpha[1,])
n1 <- draws$N[, 1, 1:n.years] / draws$Ntot[, 1:n.years]
n2 <- draws$N[, 2, 1:n.years] / draws$Ntot[, 1:n.years]
n3 <- draws$N[, 3, 1:n.years] / draws$Ntot[, 1:n.years]
n4 <- draws$N[, 4, 1:n.years] / draws$Ntot[, 1:n.years]
n5 <- draws$Nimm[, 1:n.years] / draws$Ntot[, 1:n.years]

# Extract the mean demographic rates and population sizes and store them in a list
mu <- list(phij=draws$mphij, phia=draws$mphia, f=draws$mfec,
           n1=apply(n1, 1, mean), n2=apply(n2, 1, mean), n3=apply(n3, 1, mean),n4=apply(n4, 1, mean),n5=apply(n5, 1, mean))

# Calculate sensitivities
sens <- matrix(NA, n.draws, 8)
colnames(sens) <- c("phij", "phia", "f", "n1", "n2", "n3", "n4", "n5")
sens[,"phij"] <- eval(D(lambda, "phij"), envir=mu)
sens[,"phia"] <- eval(D(lambda, "phia"), envir=mu)
sens[,"f"] <- eval(D(lambda, "f"), envir=mu)
sens[,"n1"] <- eval(D(lambda, "n1"), envir=mu)
sens[,"n2"] <- eval(D(lambda, "n2"), envir=mu)
sens[,"n3"] <- eval(D(lambda, "n3"), envir=mu)
sens[,"n4"] <- eval(D(lambda, "n4"), envir=mu)
sens[,"n5"] <- eval(D(lambda, "n5"), envir=mu)

# Define matrix to store results
# Contribution due to total variance of each demographic rate
cont_totv <- matrix(NA, nrow=n.draws, ncol=8)
# Contribution due to temporal variance of each demographic rate
cont_tempv <- matrix(NA, nrow=n.draws, ncol=8)

# Variance covariance matrix (for looking at temporal covariance and process correlation)
vcov_matrix <- array(NA, dim=c(n.draws, 8, 8))

colnames(cont_totv) <- c("phij", "phia", "f", "n1", "n2", "n3", "n4", "n5")
colnames(cont_tempv) <- c("phij", "phia", "f", "n1", "n2", "n3", "n4", "n5")

# Calculate contributions for each demographic rate and stage-structured population size at each MCMC draw
for (s in 1:n.draws){
  dp_stoch <- cbind(draws$alpha[s,1,], draws$alpha[s,2,], draws$fecund[s,1:n.years],
                    n1[s,], n2[s,], n3[s,], n4[s,], n5[s,])
  
  # Derive process variance and covariance among demographic parameters using shrunk estimates of
  #     demographic rates and proportional pop. sizes
  dp_varcov <- var(dp_stoch)
  vcov_matrix[s,,] <- var(dp_stoch)
  sensvec <- sens[s, ]
  # Calculate demographic contributions
  contmatrix <- dp_varcov * outer(sensvec, sensvec)
  cont_totv[s, ] <- rowSums(contmatrix)
  cont_tempv[s, ] <- diag(contmatrix)
}
###############################################################################
# Calculate contribution of the realized temporal variation of the demographic rates
# and of proportional population structure (n) to temporal variability of the realized
# population growth rate (vertical lines  show limits of 95% CRI)
###############################################################################
CRI <- apply(cont_totv, 2, quantile, probs=c(0.025, 0.975))

# Calculate mean
mean <- colMeans(cont_totv)

# Bind together with CRI
ATPU_tLTRE_df1 = rbind(CRI,mean) 

# Change row names
rownames(ATPU_tLTRE_df1) <- c("ll_value","ul_value","mean_value")

# Transpose dataframe
ATPU_tLTRE_df1 = as.data.frame(t(ATPU_tLTRE_df1))

# Make demographic rate a type column
ATPU_tLTRE_df1 <- ATPU_tLTRE_df1  %>%
  rownames_to_column(var = "type")

# Define the desired order of levels
desired_order <- c("phia", "phij", "f","n1","n2","n3","n4","n5")  # Replace with your specific order

# Reorder the factor levels based on the desired order
ATPU_tLTRE_df1$type <- factor(ATPU_tLTRE_df1$type, levels = desired_order)


###############################################################################
# Contribution of the realized temporal variation of the demographic rates
# and of proportional population structure (n) to temporal variability of the realized
# population growth rate (vertical lines  show limits of 95% CRI)
# Relative % contribution
###############################################################################

# Create new data frame for calculating total contributions
tot_cont = cont_totv

# Sum across rows to get total contribution for each simulation
total = rowSums(tot_cont)

# Bind together as "total" column
tot_cont = cbind(tot_cont,total)

# Divide contribution of each demographic parameter by total contribution
# to get proportion contribution
js_pc = tot_cont[,1]/tot_cont[,9]*100
as_pc = tot_cont[,2]/tot_cont[,9]*100
f_pc = tot_cont[,3]/tot_cont[,9]*100
n1_pc = tot_cont[,4]/tot_cont[,9]*100
n2_pc = tot_cont[,5]/tot_cont[,9]*100
n3_pc = tot_cont[,6]/tot_cont[,9]*100
n4_pc = tot_cont[,7]/tot_cont[,9]*100
n5_pc = tot_cont[,8]/tot_cont[,9]*100


# Get mean and 95% CRI for proportional contribution of each demographic parameter
# and create dataframe
ATPU_tLTRE_df1_2=data.frame(type=c("phij","phia","f","n1","n2","n3","n4","n5"),
                            mean=c(mean(js_pc), mean(as_pc), mean(f_pc),
                                   mean(n1_pc), mean(n2_pc), mean(n3_pc),
                                   mean(n4_pc), mean(n5_pc)),
                            ll=c(quantile(js_pc,0.025),quantile(as_pc,0.025),quantile(f_pc,0.025),
                                 quantile(n1_pc,0.025),quantile(n2_pc,0.025),quantile(n3_pc,0.025),
                                 quantile(n4_pc,0.025),quantile(n5_pc,0.025)),
                            ul=c(quantile(js_pc,0.975),quantile(as_pc,0.975),quantile(f_pc,0.975),
                                 quantile(n1_pc,0.975),quantile(n2_pc,0.975),quantile(n3_pc,0.975),
                                 quantile(n4_pc,0.975),quantile(n5_pc,0.975)))



###############################################################################
# Contribution of the realized temporal variability and of the realized
# temporal covariance among demographic rates to the temporal variability of the realized
# population growth rate (vertical lines  show limits of 95% CRI)
###############################################################################
# Mean contribution due to total variance
mean_cont_totv <- colMeans(cont_totv)

# Mean contribution due to temporal variance
mean_cont_tempv <- colMeans(cont_tempv)

# Determine contribution due to temporal covariance
cont_tcov=cont_totv - cont_tempv

# Mean contribution due to temporal covariance
mean_cont_tcov=colMeans(cont_tcov)

# Bind together 
vcov_contribution = rbind(mean_cont_tempv,mean_cont_tcov)

# Create dataframe
ATPU_tLTRE_df2 = as.data.frame(vcov_contribution)

# Change row names
rownames(ATPU_tLTRE_df2) <- c("Variance","Covariance")

# Add a column for variance or covariance
ATPU_tLTRE_df2$vcov <- row.names(ATPU_tLTRE_df2)

# Reshape the data to long format
ATPU_tLTRE_df2 <-ATPU_tLTRE_df2 %>%
  pivot_longer(cols = -vcov, names_to = "Demo_rate", values_to = "value")

# Reorder the factor levels based on the desired order (defined above)
ATPU_tLTRE_df2$Demo_rate <- factor(ATPU_tLTRE_df2$Demo_rate, levels = desired_order)

#Pull out CRI data from first df and add to new df
ATPU_tLTRE_df2$ll_data = c(ATPU_tLTRE_df1$ll_value, rep(NA, 8))
ATPU_tLTRE_df2$ul_data = c(ATPU_tLTRE_df1$ul_value, rep(NA, 8))

###############################################################################
# Contribution of the realized temporal variability and of the realized
# temporal covariance among demographic rates to the temporal variability of the realized
# population growth rate (vertical lines  show limits of 95% CRI)
# Relative % contribution
###############################################################################

# Define new column names for variance matrix
var_col_names <- c("phij_v", "phia_v", "f_v", "n1_v", "n2_v", "n3_v", "n4_v", "n5_v")

# Rename the columns
colnames(cont_tempv) <- var_col_names

# Define new column names for covvariance matrix
cvar_col_names <- c("phij_cv", "phia_cv", "f_cv", "n1_cv", "n2_cv", "n3_cv", "n4_cv", "n5_cv")

# Rename the columns
colnames(cont_tcov) <- cvar_col_names

tot_cont_vcov = cbind(cont_tempv,cont_tcov)

# Sum across rows to get total contribution for each simulation
total = rowSums(tot_cont_vcov)

# Bind together as "total" column
tot_cont_vcov = cbind(tot_cont_vcov,total)

# Divide contribution of each demographic parameter by total contribution
# to get proportion contribtion
# Variance
js_pc_v = tot_cont_vcov[,1]/tot_cont_vcov[,17]*100
as_pc_v = tot_cont_vcov[,2]/tot_cont_vcov[,17]*100
f_pc_v = tot_cont_vcov[,3]/tot_cont_vcov[,17]*100
n1_pc_v = tot_cont_vcov[,4]/tot_cont_vcov[,17]*100
n2_pc_v = tot_cont_vcov[,5]/tot_cont_vcov[,17]*100
n3_pc_v = tot_cont_vcov[,6]/tot_cont_vcov[,17]*100
n4_pc_v = tot_cont_vcov[,7]/tot_cont_vcov[,17]*100
n5_pc_v = tot_cont_vcov[,8]/tot_cont_vcov[,17]*100

#Covariance
js_pc_cv = tot_cont_vcov[,9]/tot_cont_vcov[,17]*100
as_pc_cv = tot_cont_vcov[,10]/tot_cont_vcov[,17]*100
f_pc_cv = tot_cont_vcov[,11]/tot_cont_vcov[,17]*100
n1_pc_cv = tot_cont_vcov[,12]/tot_cont_vcov[,17]*100
n2_pc_cv = tot_cont_vcov[,13]/tot_cont_vcov[,17]*100
n3_pc_cv = tot_cont_vcov[,14]/tot_cont_vcov[,17]*100
n4_pc_cv = tot_cont_vcov[,15]/tot_cont_vcov[,17]*100
n5_pc_cv = tot_cont_vcov[,16]/tot_cont_vcov[,17]*100



# Get mean and 95% CRI for proportional contribution of each demographic parameter
# and create dataframe
ATPU_tLTRE_df2_2=data.frame(type=c("phij","phia","f","n1","n2","n3","n4","n5"),
                            mean_v=c(mean(js_pc_v), mean(as_pc_v), mean(f_pc_v),
                                     mean(n1_pc_v), mean(n2_pc_v), mean(n3_pc_v),
                                     mean(n4_pc_v), mean(n5_pc_v)),
                            mean_cv = c(mean(js_pc_cv), mean(as_pc_cv), mean(f_pc_cv),
                                        mean(n1_pc_cv), mean(n2_pc_cv), mean(n3_pc_cv),
                                        mean(n4_pc_cv), mean(n5_pc_cv)),
                            ll=c(quantile(js_pc_v+js_pc_cv,0.025),quantile(as_pc_v+as_pc_cv,0.025),quantile(f_pc_v+f_pc_cv,0.025),
                                 quantile(n1_pc_v+n1_pc_cv,0.025),quantile(n2_pc_v+n2_pc_cv,0.025),quantile(n3_pc_v+n3_pc_cv,0.025),
                                 quantile(n4_pc_v+n4_pc_cv,0.025),quantile(n5_pc_v+n5_pc_cv,0.025)),
                            ul=c(quantile(js_pc_v+js_pc_cv,0.975),quantile(as_pc_v+as_pc_cv,0.975),quantile(f_pc_v+f_pc_cv,0.975),
                                 quantile(n1_pc_v+n1_pc_cv,0.975),quantile(n2_pc_v+n2_pc_cv,0.975),quantile(n3_pc_v+n3_pc_cv,0.975),
                                 quantile(n4_pc_v+n4_pc_cv,0.975),quantile(n5_pc_v+n5_pc_cv,0.975)))


# Reshape the data

ATPU_tLTRE_df2_2<- pivot_longer(
  data = ATPU_tLTRE_df2_2,
  cols = starts_with("mean_") | starts_with("ll") | starts_with("ul"),
  names_to = c(".value", "stat"),
  names_sep = "_"
) 


# Pull out posteriors for temporal covariance
# Juvenile and adult survival
juv_ad_tc = vcov_matrix[,1,2]
# Juvenile survival and fecundity
juv_fec_tc = vcov_matrix[,1,3]
# Adult survival and fecundity
ad_fec_tc = vcov_matrix[,2,3]

# Bind together in dataframe
ATPU_tLTRE_df3 = as.data.frame(cbind(juv_ad_tc,juv_fec_tc,ad_fec_tc))

# Reformat dataframe 
ATPU_tLTRE_df3 <- ATPU_tLTRE_df3   %>%
  gather(key = "type", value = "value")


# Pull individual variances
#juvenile temporal variance
juv_var = vcov_matrix[,1,1]

#adult temporal variance
ad_var = vcov_matrix[,2,2]

# Temporal variance of productivity
fec_var = vcov_matrix[,3,3]

# Bind together in dataframe
ATPU_tLTRE_df3_2 = as.data.frame(cbind(juv_var,ad_var,fec_var))

# Reformat dataframe 
ATPU_tLTRE_df3_2  <- ATPU_tLTRE_df3_2   %>%
  gather(key = "type", value = "value")

###############################################################################
# Calculate process correlation of demographic rates
# Measures correlation at specific time points
## (doesn't specifically address time or sequence but focuses 
## on the correlation between variables at a given moment)
###############################################################################
# Calculate process variance 
# Temporal covariance divided by square root of product of temporal variances

# Calculating process correlation
# Adult survival and juvenile survival
ad_juv_pv = juv_ad_tc/sqrt(juv_var*ad_var)

# Juvenile survival and fecundity
juv_fec_pv = juv_fec_tc/sqrt(juv_var*fec_var)

# Adult survival and fecundity
ad_fec_pv = ad_fec_tc/sqrt(ad_var*fec_var)


# Bind together in dataframe
ATPU_tLTRE_df4 = as.data.frame(cbind(ad_juv_pv,juv_fec_pv,ad_fec_pv))

# Reformat dataframe 
ATPU_tLTRE_df4  <- ATPU_tLTRE_df4   %>%
  gather(key = "type", value = "value")


# Calculate mean, 2.5% quantile, and 97.5% quantile for each type
ATPU_ss_pc <- tapply(ATPU_tLTRE_df4$value, ATPU_tLTRE_df4$type, function(x) {
  mean_value <- mean(x)
  ll <- quantile(x, 0.025)
  ul <- quantile(x, 0.975)
  
  return(c(mean = mean_value, ll = ll, ul = ul))
})


#####################################################################################
# Second type of tLTRE
## Understand how much differences in demographic rates and population structure of
## successive years have contributed to the difference in realized pop growth
###############################################################################
# Compute differences and means of demographic rates and population structure between successive years
###############################################################################
diff <- array(NA, dim=c(n.draws, n.years-1, 8))
dimnames(diff)[[3]] <- c("phij", "phia", "f", "n1", "n2", "n3", "n4", "n5")

# Function to compute differences over successive time steps
getDiff <- function(x) x[,2:n.years] - x[,1:(n.years-1)]

diff[,,"phij"] <- getDiff(draws$alpha[,1,])
diff[,,"phia"] <- getDiff(draws$alpha[,2,])
diff[,,"f"] <- getDiff(draws$fecund)
diff[,,"n1"] <- getDiff(draws$N[,1,])
diff[,,"n2"] <- getDiff(draws$N[,2,])
diff[,,"n3"] <- getDiff(draws$N[,3,])
diff[,,"n4"] <- getDiff(draws$N[,4,])
diff[,,"n5"] <- getDiff(draws$Nimm)

# Function to compute means over successive time steps, store them in a list
getMn <- function(x) (x[,2:n.years] + x[,1:(n.years-1)]) / 2

means <- list(phij=getMn(draws$alpha[,1,]), phia=getMn(draws$alpha[,2,]),
              f=getMn(draws$fecund),n1=getMn(draws$N[,1,]), n2=getMn(draws$N[,2,]),
              n3=getMn(draws$N[,3,]),n4=getMn(draws$N[,4,]),n5=getMn(draws$Nimm))

# Compute sensitivities
senss <- array(NA, dim=c(n.draws, n.years-1, 8))
dimnames(senss)[[3]] <- c("phij", "phia","f","n1", "n2", "n3","n4","n5")
senss[,,"phij"] <- eval(D(lambda, "phij"), envir=means)
senss[,,"phia"] <- eval(D(lambda, "phia"), envir=means)
senss[,,"f"] <- eval(D(lambda, "f"), envir=means)
senss[,,"n1"] <- eval(D(lambda, "n1"), envir=means)
senss[,,"n2"] <- eval(D(lambda, "n2"), envir=means)
senss[,,"n3"] <- eval(D(lambda, "n3"), envir=means)
senss[,,"n4"] <- eval(D(lambda, "n4"), envir=means)
senss[,,"n5"] <- eval(D(lambda, "n5"), envir=means)

# Calculate contributions of the differences 
# between successive years in the demographic parameters and population structure
conts <- diff * senss


# Sequential differences in finite growth rates
n.draws <- ATPU_IPM_out$mcmc.info$n.samples
n.years <- length(ATPU_IPM_out$mean$alpha[1,])
diff.lam <- draws$ann.growth.rate[, 2:n.years] - draws$ann.growth.rate[, 1:(n.years-1)]
differences <- cbind(1998:2017, apply(diff.lam, 2, mean))

# Mean contributions
V <- apply(conts, 3:2, mean)
V1 <- V2 <- V
V1[V1<0] <- 0
V2[V2>0] <- 0

# Create vector with beginning year
b_year = c(1998:2017)

# Create vector with differences in population growth
diffs = differences[,2]

# Bind together into a single dataframe
ATPU_tLTRE_df5 = as.data.frame(cbind(b_year,diffs))

# Covert year to factor
ATPU_tLTRE_df5$b_year <- as.factor(ATPU_tLTRE_df5$b_year)


# Get contribution of the changes in demographic rates and of population structure
# Make V a dataframe
ATPU_tLTRE_df6 = as.data.frame(V) 

# Rename columns to reflect beginning year
ATPU_tLTRE_df6 = ATPU_tLTRE_df6 %>% dplyr::rename(
  "1998"= V1,"1999"= V2,"2000"= V3,"2001"= V4,
  "2002"= V5,"2003"= V6,"2004"= V7, "2005"= V8, "2006"= V9, "2007"= V10, "2008"= V11,
  "2009"= V12, "2010"= V13, "2011"= V14, "2012"= V15, "2013"= V16, "2014"= V17, "2015"= V18,
  "2016"= V19, "2017"= V20)

# Convert the data to a long format and arrange columns
ATPU_tLTRE_df6 <- ATPU_tLTRE_df6 %>%
  rownames_to_column(var = "type") %>%
  gather(year, value, -type) %>%
  arrange(type)





