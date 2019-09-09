prev.modSel <- function(){
  # Likelihood:
  for(i in 1:nsite){
    for(j in 1:nplot[i]){ # The number of plots is denoted by the vector nplot[i]
      for(t in 1:nyear[i,j]){ # The number years it is surveyed is denoted by nyear[i,j]
        for(k in 1:nobs[i,j,t]){ # The number of observations is denoted by nobs[i,j,t]
          y[i,j,k,t] ~ dbern(1-(1-rho[i,j,t])^n[i,j,k,t])
        }
        logit(rho[i,j,t]) <- beta0[i,j] + inprod(beta.year.raw[], year[i,j,t,]) + 
       
        equals(I1, 1)*(beta[1]*prop.forest[i,j,1] + beta[2]*edge.forest[i,j,1] + beta[3]*prop.forest[i,j,1]*edge.forest[i,j,1]) +
        equals(I1, 2)*(beta[1]*prop.forest[i,j,2] + beta[2]*edge.forest[i,j,2] + beta[3]*prop.forest[i,j,2]*edge.forest[i,j,2]) +
        equals(I1, 3)*(beta[1]*prop.forest[i,j,3] + beta[2]*edge.forest[i,j,3] + beta[3]*prop.forest[i,j,3]*edge.forest[i,j,3]) +
 
        equals(I1, 4)*(-beta[1]*prop.urban[i,j,1] + beta[2]*edge.forest[i,j,1] - beta[3]*prop.urban[i,j,1]*edge.forest[i,j,1]) +
        equals(I1, 5)*(-beta[1]*prop.urban[i,j,2] + beta[2]*edge.forest[i,j,2] - beta[3]*prop.urban[i,j,2]*edge.forest[i,j,2]) +
        equals(I1, 6)*(-beta[1]*prop.urban[i,j,3] + beta[2]*edge.forest[i,j,3] - beta[3]*prop.urban[i,j,3]*edge.forest[i,j,3]) +
    
        equals(I1, 7)*(beta[1]*prop.forest[i,j,1] - beta[2]*prox.forest[i,j,1] - beta[3]*prop.forest[i,j,1]*prox.forest[i,j,1]) +
        equals(I1, 8)*(beta[1]*prop.forest[i,j,2] - beta[2]*prox.forest[i,j,2] - beta[3]*prop.forest[i,j,2]*prox.forest[i,j,2]) +
        equals(I1, 9)*(beta[1]*prop.forest[i,j,3] - beta[2]*prox.forest[i,j,3] - beta[3]*prop.forest[i,j,3]*prox.forest[i,j,3]) +
    
        equals(I1, 10)*(-beta[1]*prop.urban[i,j,1] - beta[2]*prox.forest[i,j,1] + beta[3]*prop.forest[i,j,1]*prox.forest[i,j,1]) +
        equals(I1, 11)*(-beta[1]*prop.urban[i,j,2] - beta[2]*prox.forest[i,j,2] + beta[3]*prop.forest[i,j,2]*prox.forest[i,j,2]) +
        equals(I1, 12)*(-beta[1]*prop.urban[i,j,3] - beta[2]*prox.forest[i,j,3] + beta[3]*prop.forest[i,j,3]*prox.forest[i,j,3]) +
        
        beta[4]*prop.eg[i,j,I2]
        + beta[5]*dist.water[i,j]
      }

      beta0[i,j] ~ dnorm(alpha[i], sigma^-2)
    }
    alpha[i] ~ dnorm(mu.raw, sigma0^-2)
  }
  
  # Priors:
    
  ## Precision for parameter estimates:
  prec <- 5^-1
  
  for(i in 1:5){
      beta[i] ~ dnorm(0, prec)
  }
  
  for(t in 1:nyear.max){
    beta.year.raw[t] ~ dnorm(0, prec)
  }
  
  I1 ~ dcat(rep(1,12))
  I2 ~ dcat(rep(1,3))
  
  mu.raw ~ dnorm(logit(.018), prec)
  
  sigma0 ~ dunif(0,10)
  sigma ~ dunif(0,10)
  
  
  # Sum-to-zero constraint for mu and beta.year
  for(t in 1:nyear.max){
    mean.year[t] <- mu.raw + beta.year.raw[t]
  }
  
  mu <- mean(mean.year[])
  
  for(t in 1:nyear.max){
    beta.year[t] <- mean.year[t] - mu
  }
}
