#
# PCl-fPCA: Analysis of clusters
#
#1) Adjusted Rand Index:
#
require(pdfCluster)
require(lattice)
#
#true clusters  
Clust_true1 <- c(rep(1,50), rep(2,50))
Clust_true2 <- c(rep(1,20), rep(2,30), rep(3,50))
Clust_true = list(Clust_true1 = Clust_true1,
                  Clust_true2 = Clust_true2)
#
#extract cluster allocation using MAP
N_curv = fPCA_info$N_curv
dim.cl = dim.space# or # of dimensions explored for clusters
clust_label <- array(0, c(N_curv,dim.cl,
                       (dim(PostSamples[[1]])[1])))
clust_MAP <- matrix(0, N_curv, dim.cl)
maxlab = NULL
#
  for(i in 1:N_curv){
    for(k in 1:dim.cl){
      clust_label[i,k,] <- PostSamples[[1]][,sprintf("c[%s,%s]",i,k)]
      maxlab <- which.max(prop.table(table(clust_label[i,k,])))
      clust_MAP[i,k] <- as.integer(names(maxlab))
    }
  }
# it should run in few seconds
#
ARI = matrix(NA,1,dim.cl)
for(k in 1:dim.cl){
ARI[1,k] <- round(adj.rand.index(Clust_true[[k]],clust_MAP[,k]),3)
}
print(ARI)
#
#
#2) pairwise probabilities
#
# dimension 1:
d1_lab <- array(0,c(2,N_curv,(dim(PostSamples[[1]])[1])))
d1_ppm <- matrix(0,N_curv,N_curv)
#
for(i in 1:N_curv){
  for(l in 1:N_curv){
    d1_lab[1,i,] <- PostSamples[[1]][,sprintf("c[%s,1]",i)]
    d1_lab[2,l,] <- PostSamples[[1]][,sprintf("c[%s,1]",l)]
    d1_ppm[l,i] <- round((length(which(d1_lab[1,i,]==d1_lab[2,l,]))
                          /dim(PostSamples[[1]])[1]),4)
    
  }
}
rm(d1_lab)
# it should run in few seconds
# dimension 2:
d2_lab <- array(0,c(2,N_curv,(dim(PostSamples[[1]])[1])))
d2_ppm <- matrix(0,N_curv,N_curv)
#
for(i in 1:N_curv){
  for(l in 1:N_curv){
    d2_lab[1,i,] <- PostSamples[[1]][,sprintf("c[%s,2]",i)]
    d2_lab[2,l,] <- PostSamples[[1]][,sprintf("c[%s,2]",l)]
    d2_ppm[l,i] <- round((length(which(d2_lab[1,i,]==d2_lab[2,l,]))
                        /dim(PostSamples[[1]])[1]),4)
    
  }
}
rm(d2_lab)
#
# plot matrices
x11(w=20,h=10)
dim1_ppm <- levelplot(d1_ppm, col.regions = gray(0:100/100),
                    main="PPM dimension 1")
dim2_ppm <- levelplot(d2_ppm, col.regions = gray(0:100/100),
                    main="PPM dimension 2")
print(dim1_ppm, position=c(0, 0, 0.5, 1),more=TRUE)
print(dim2_ppm, position=c(0.5, 0, 1, 1))
#
#end

