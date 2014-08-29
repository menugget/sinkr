#This function reconstructs a data set using a defined set of principal components.
#arguments "pca" is the pca object from prcomp, "pcs" is a vector of principal components
#to be used for reconstruction (default includes all pcs)
prcomp.recon <- function(pca, pcs=NULL){
  if(is.null(pcs)) pcs <- seq(pca$sdev)
  recon <- as.matrix(pca$x[,pcs]) %*% t(as.matrix(pca$rotation[,pcs]))
  if(pca$scale[1] != FALSE){
  	recon <- scale(recon , center=FALSE, scale=1/pca$scale)
  }
  if(pca$center[1] != FALSE){
	recon <- scale(recon , center=-pca$center, scale=FALSE)
  }
  recon
}