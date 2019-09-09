prevalence <- function(){
  # Declaring array dimensions:
  # Likelihood:
  for(i in 1:nsite){
    for(j in 1:nplot[i]){ # The number of plots is denoted by the vector nplot[i]
      for(t in 1:nyear[i,j]){ # The number years it is surveyed is denoted by nyear[i,j]
        for(k in 1:nobs[i,j,t]){ # The number of observations is denoted by nobs[i,j,t]
          y[i,j,k,t] ~ dbern(1-(1-rho[i,j,t])^n[i,j,k,t])
        }
        logit(rho[i,j,t]) <- beta0[i,j] 
        + inprod(beta[], X[i,j,])
        + inprod(beta.year.raw[], year[i,j,t,])
      }
      beta0[i,j] ~ dnorm(alpha[i], sigma^-2)
    }
    alpha[i] ~ dnorm(mu.raw, sigma0^-2)
  }
  
  
  # Priors:
  for(z in 1:nbetas){
    beta[z] ~ dnorm(0, 5^-2)
  }

  for(t in 1:nyear.max){
    beta.year.raw[t] ~ dnorm(0, 5^-2)
  }
  
  mu.raw ~ dnorm(logit(.018), 5^-2)
  
  sigma0 ~ dunif(0, 10)
  sigma ~ dunif(0,10)
  
  # Sum-to-zero constraint for mu and beta.year
  for(t in 1:nyear.max){
    mean.year[t] <- mu.raw + beta.year.raw[t]
  }
  
  mu <- mean(mean.year[])
  
  for(t in 1:nyear.max){
    beta.year[t] <- mean.year[t] - mu
  }
  
  
  # Derived, back-transformed parameters:
  MU <- ilogit(mu) * 100

  for(i in 1:nsite){
    ALPHA[i] <- ilogit(alpha[i]) * 100
  }

  for(i in 1:nbetas){
    BETA[i] <- ilogit(mu + beta[i]) - ilogit(mu)
  }
  
  for(t in 1:nyear.max){
    YEAR.MEAN[t] <- ilogit(mu + beta.year[t])
  }
}


