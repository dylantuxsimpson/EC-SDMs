occuMS.modSel <- function(){
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
      logit(eps[i,j]) <- b0.e[i] +
        equals(ie1, 1)*(beta.eps[1]*prop.forest[i,j,1] + beta.eps[2]*edge.forest[i,j,1] + beta.eps[3]*prop.forest[i,j,1]*edge.forest[i,j,1]) +
        equals(ie1, 2)*(beta.eps[1]*prop.forest[i,j,2] + beta.eps[2]*edge.forest[i,j,2] + beta.eps[3]*prop.forest[i,j,2]*edge.forest[i,j,2]) +
        equals(ie1, 3)*(beta.eps[1]*prop.forest[i,j,3] + beta.eps[2]*edge.forest[i,j,3] + beta.eps[3]*prop.forest[i,j,3]*edge.forest[i,j,3]) +

        equals(ie1, 4)*(-beta.eps[1]*prop.urban[i,j,1] + beta.eps[2]*edge.forest[i,j,1] - beta.eps[3]*prop.urban[i,j,1]*edge.forest[i,j,1]) +
        equals(ie1, 5)*(-beta.eps[1]*prop.urban[i,j,2] + beta.eps[2]*edge.forest[i,j,2] - beta.eps[3]*prop.urban[i,j,2]*edge.forest[i,j,2]) +
        equals(ie1, 6)*(-beta.eps[1]*prop.urban[i,j,3] + beta.eps[2]*edge.forest[i,j,3] - beta.eps[3]*prop.urban[i,j,3]*edge.forest[i,j,3]) +

        equals(ie1, 7)*(beta.eps[1]*prop.forest[i,j,1] - beta.eps[2]*prox.forest[i,j,1] - beta.eps[3]*prop.forest[i,j,1]*prox.forest[i,j,1]) +
        equals(ie1, 8)*(beta.eps[1]*prop.forest[i,j,2] - beta.eps[2]*prox.forest[i,j,2] - beta.eps[3]*prop.forest[i,j,2]*prox.forest[i,j,2]) +
        equals(ie1, 9)*(beta.eps[1]*prop.forest[i,j,3] - beta.eps[2]*prox.forest[i,j,3] - beta.eps[3]*prop.forest[i,j,3]*prox.forest[i,j,3]) +

        equals(ie1, 10)*(-beta.eps[1]*prop.urban[i,j,1] - beta.eps[2]*prox.forest[i,j,1] + beta.eps[3]*prop.forest[i,j,1]*prox.forest[i,j,1]) +
        equals(ie1, 11)*(-beta.eps[1]*prop.urban[i,j,2] - beta.eps[2]*prox.forest[i,j,2] + beta.eps[3]*prop.forest[i,j,2]*prox.forest[i,j,2]) +
        equals(ie1, 12)*(-beta.eps[1]*prop.urban[i,j,3] - beta.eps[2]*prox.forest[i,j,3] + beta.eps[3]*prop.forest[i,j,3]*prox.forest[i,j,3]) +

        beta.eps[4]*prop.eg[i,j,ie2] + beta.eps[5]*dist.water[i,j]

      logit(gam[i,j]) <- b0.g[i] +
        equals(ig1, 1)*(beta.eps[1]*prop.forest[i,j,1] + beta.eps[2]*edge.forest[i,j,1] + beta.eps[3]*prop.forest[i,j,1]*edge.forest[i,j,1]) +
        equals(ig1, 2)*(beta.eps[1]*prop.forest[i,j,2] + beta.eps[2]*edge.forest[i,j,2] + beta.eps[3]*prop.forest[i,j,2]*edge.forest[i,j,2]) +
        equals(ig1, 3)*(beta.eps[1]*prop.forest[i,j,3] + beta.eps[2]*edge.forest[i,j,3] + beta.eps[3]*prop.forest[i,j,3]*edge.forest[i,j,3]) +
        
        equals(ig1, 4)*(-beta.eps[1]*prop.urban[i,j,1] + beta.eps[2]*edge.forest[i,j,1] - beta.eps[3]*prop.urban[i,j,1]*edge.forest[i,j,1]) +
        equals(ig1, 5)*(-beta.eps[1]*prop.urban[i,j,2] + beta.eps[2]*edge.forest[i,j,2] - beta.eps[3]*prop.urban[i,j,2]*edge.forest[i,j,2]) +
        equals(ig1, 6)*(-beta.eps[1]*prop.urban[i,j,3] + beta.eps[2]*edge.forest[i,j,3] - beta.eps[3]*prop.urban[i,j,3]*edge.forest[i,j,3]) +
        
        equals(ig1, 7)*(beta.eps[1]*prop.forest[i,j,1] - beta.eps[2]*prox.forest[i,j,1] - beta.eps[3]*prop.forest[i,j,1]*prox.forest[i,j,1]) +
        equals(ig1, 8)*(beta.eps[1]*prop.forest[i,j,2] - beta.eps[2]*prox.forest[i,j,2] - beta.eps[3]*prop.forest[i,j,2]*prox.forest[i,j,2]) +
        equals(ig1, 9)*(beta.eps[1]*prop.forest[i,j,3] - beta.eps[2]*prox.forest[i,j,3] - beta.eps[3]*prop.forest[i,j,3]*prox.forest[i,j,3]) +
        
        equals(ig1, 10)*(-beta.eps[1]*prop.urban[i,j,1] - beta.eps[2]*prox.forest[i,j,1] + beta.eps[3]*prop.forest[i,j,1]*prox.forest[i,j,1]) +
        equals(ig1, 11)*(-beta.eps[1]*prop.urban[i,j,2] - beta.eps[2]*prox.forest[i,j,2] + beta.eps[3]*prop.forest[i,j,2]*prox.forest[i,j,2]) +
        equals(ig1, 12)*(-beta.eps[1]*prop.urban[i,j,3] - beta.eps[2]*prox.forest[i,j,3] + beta.eps[3]*prop.forest[i,j,3]*prox.forest[i,j,3]) +

        beta.gam[4]*prop.eg[i,j,ig2] + beta.gam[5]*dist.water[i,j]
    }
    # Random intercepts:
    b0.e[i] ~ dnorm(alpha.e, sigma.e^-2)
    b0.g[i] ~ dnorm(alpha.g, sigma.g^-2)
    b0.p[i] ~ dnorm(alpha.p, sigma.p^-2)
  }
  
  
  # Priors & derivations
  
  for(i in 1:5){
    beta.eps[i] ~ dnorm(0, 5^-1)
  }
  for(i in 1:5){
    beta.gam[i] ~ dnorm(0, 5^-1)
  }
  
  for(i in 1:2){
    beta.p[i] ~ dnorm(0, 5^-1)
  }
  
  ie1 ~ dcat(rep(1,12))
  ie2 ~ dcat(rep(1,3))
  ig1 ~ dcat(rep(1,12))
  ig2 ~ dcat(rep(1,3))
  
  alpha.e ~ dnorm(0, 5^-1)
  alpha.g ~ dnorm(0, 5^-1)
  alpha.p ~ dnorm(0, 5^-1)
  
  sigma.e ~ dunif(0,10)
  sigma.g ~ dunif(0,10)
  sigma.p ~ dunif(0,10)
  
  det.error ~ dunif(0,10)
  
  meanEps <- ilogit(alpha.e)
  meanGam <- ilogit(alpha.g)
  meanP <- ilogit(alpha.p) 
}
