new.range <- function(x, new.min=0, new.max=1){

	ranx <- range(x, na.rm=TRUE)
	px <- (x-ranx[1])/(ranx[2]-ranx[1])
	res <- new.min+((new.max-new.min)*px)
	res

}

