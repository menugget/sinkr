#reverses the scaling of an object produced by scale()
rev.scale <- function(x, scale=TRUE, center=TRUE){
	if(scale & !is.null(attr(x, "scaled:scale"))){
		x <- scale(x, center=FALSE, scale=1/attr(x,"scaled:scale"))
		attr(x,"scaled:scale") <- NULL
	}
	if(center & !is.null(attr(x, "scaled:center"))){
		x <- scale(x, center=-1*attr(x,"scaled:center"), scale=FALSE)
		attr(x,"scaled:center") <- NULL
	}
	return(x)
}
