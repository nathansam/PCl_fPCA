#
### Parameter Clustering Bayesian functional PCA (PCl-fPCA) ###
#
# The following code perform PCl-fPCA. The method is introduced
# in "Margaritella N., Inacio V. and King R. Parameter clustering
# in Bayesian functional principal component analysis of 
# neuroscientific data.Stat Med.2020;1-18.
# https://doi.org/10.1002/sim.8768"
#
#
#
load(fPCA_data)
require(rjags)
#
#1) datalist
#
N_curv<- fPCA_info$N_curv # number of curves i
Timep <- fPCA_info$Timep  # number of time points t
Y_it <- fPCA_info$Y_it # raw observations (i*t)
E <- fPCA_info$EigFun # smooth eigenfunctions 
dim.space <- dim(E)[2] # number of eigendimensions k
clmax = 20 # upperbound J for the clusters
#
eigval <- fPCA_info$EigVal # eigenvalues (exp.var.:78%,16%,1%)
inv_eigval <- 1/ eigval # eigenvalues^-1 
sqrt_eigval <- sqrt(eigval) # square-rooted eigenvalues
upb_alph <- c(6,5) # upper-bounds for the dispersion param. alpha
#
#
dataList = list(Yc=Y_it, N_curv = N_curv, Timep = Timep,
                E = E, dim.space = dim.space, 
                inv_eigval = inv_eigval,
                sqrt_eigval = sqrt_eigval,
                clmax = clmax, upb_alph = upb_alph)

#2) Model definition:
#
modelString = "
model {
  for (i in 1:N_curv) {
    for(t in 1:Timep) {
      Yc[i,t] ~ dnorm(m[i,t],tau_eps) 
      m[i,t] <- xi[i,1:dim.space] %*% E[t,1:dim.space]
    }# close t
    for(k in 1:2){
      xi[i,k] ~ dnorm(muofclust[clmax+1-c[i,k],k], tauofclust[clmax+1-c[i,k],k])
      c[i,k] ~ dcat(pi_ord[,k])
    }# close k
# further dimensions can be added as in standard BfPCA
#    for(k in 3:dim.space){
#      xi[i,k] ~ dnorm(0, eigval[k-2])
#    }# close k
#
  }#close i
#
  for(j in 1:clmax){
    for(k in 1:2){
      muofclust[j,k] ~ dnorm(0, inv_eigval[k])
      tauofc[j,k] ~ dunif(0, sqrt_eigval[k])
      tauofclust[j,k] <- 1/(tauofc[j,k]^2)
    }# close k
  }#close j
#
 for(k in 1:2){
  p[1,k]<-v[1,k]
  for(j in 2:clmax){
   p[j,k]<-v[j,k]*(1-v[j-1,k])*(p[j-1,k]/v[j-1,k])
  }#close j
  p_sum[k]<-sum(p[,k])
  for(j in 1:clmax){
   v[j,k] ~ dbeta(1, Alpha[k])T(0.001,0.999)
      pi[j,k]<-p[j,k]/p_sum[k]
  }#close j
  pi_ord[1:clmax,k] <- sort(pi[1:clmax,k])
  Alpha[k] ~ dunif(0, upb_alph[k])
 }#close k
#
tau_eps ~ dgamma(0.001, 0.001)# 
}" 
writeLines(modelString , con="PCl_fPCA.txt")

#3) JAGS
#
jagsModel = jags.model(file="PCl_fPCA.txt", 
                       data=dataList ,n.chains=1,
                       n.adapt=1000,
                       inits=list(.RNG.name="base::Super-Duper",
                                  .RNG.seed=1671))
#
update(jagsModel , n.iter=100000)
#
PostSamples = coda.samples(jagsModel,
                           variable.names=c("xi","tau_eps","muofclust",
                                            "tauofclust","pi_ord","c",
                                            "Alpha"),
                          n.iter=100000, thin=5)
#
# The model should run in ~ 45 min. 
#
# end