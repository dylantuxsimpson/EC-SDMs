occuMS.hier.missYear <- function(){
  # Unobserved years:
  for(i in 1:nsites){
    for(j in 1:nplots[i]){
      z[i,j,1] ~ dbern(0.8)
    }
  }
  
  for(i in 1:nsites){
    for(j in 1:nplots[i]){
      z[i,j,4] ~ dbern(z[i,j,3]*(1-eps[i,j])+(1-z[i,j,3])*gam[i,j])
    }
  }
        
  
  # Observed years:
  for(i in 1:nsites){
    for(j in 1:nplots[i]){
      for(t in yearsObs[]){
        for(k in 1:nobs){
          obs[i,j,k,t] ~ dbern(z[i,j,t]*p[i,j,t])
        }
        z[i,j,t] ~ dbern(z[i,j,t-1]*(1-eps[i,j])+(1-z[i,j,t-1])*gam[i,j])
        logit(p[i,j,t]) <- det.odds[i,j,t]
        det.odds[i,j,t] ~ dnorm(b0.p[i] + inprod(beta.p[], X.p[i,j,t,]), det.error^-2)
      }
      logit(eps[i,j]) <- b0.e[i] + inprod(beta.eps[], X.e[i,j,])
      logit(gam[i,j]) <- b0.g[i] + inprod(beta.gam[], X.g[i,j,])
    }
    # Random intercepts:
    b0.e[i] ~ dnorm(alpha.e, sigma.e^-2)
    b0.g[i] ~ dnorm(alpha.g, sigma.g^-2)
    b0.p[i] ~ dnorm(alpha.p, sigma.p^-2)
  }
  
  
  # Priors & derivations
  prec <- 5^-1

  for(i in 1:nbeta.e){
    beta.eps[i] ~ dnorm(0, prec)
    # beta.eps[i] <- ilogit(alpha.e + b.e[i]) - ilogit(alpha.e)
  }
  for(i in 1:nbeta.g){
    beta.gam[i] ~ dnorm(0, prec)
    # beta.gam[i] <- ilogit(alpha.g + b.g[i]) - ilogit(alpha.g)
  }
  
  for(i in 1:nbeta.p){
    beta.p[i] ~ dnorm(0, prec)
    # beta.p[i] <- ilogit(alpha.p + b.p[i]) - ilogit(alpha.p)
  }
  
  
  alpha.e ~ dnorm(0, prec)
  alpha.g ~ dnorm(0, prec)
  alpha.p ~ dnorm(0, prec)
  
  sigma.e ~ dunif(0,10)
  sigma.g ~ dunif(0,10)
  sigma.p ~ dunif(0,10)
  
  det.error ~ dunif(0, 10)
  
  # meanEps <- ilogit(alpha.e)
  # meanGam <- ilogit(alpha.g)
  # meanP <- ilogit(alpha.p) 
}
