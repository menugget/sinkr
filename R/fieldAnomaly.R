#' Calculate the daily or monthly anomaly of a field
#' 
#' The \code{fieldAnomaly} function calculates an anomaly field by substraction the daily or 
#' monthly means from each spatial location (columns) of a field.
#' 
#' @param x A vector of the class "POSIXlt" containing the dates that correspond
#' to the rows of matrix \code{y}.
#' @param y A matrix of the spatio-temporal field with rows being the temporal 
#' dimension and columns being the spatial dimension
#' @param level Character string. Anomaly to compute ("daily" or "monthly")
#' 
#' @examples
#' 
#' set.seed(1)
#' TS <- seq.Date(as.Date("1990-01-01"), as.Date("2009-12-01"), by="month")
#' GRD <- 1:10
#' SIGNAL <- sin(as.POSIXlt(TS)$mon/2)
#' NOISE <- rnorm(length(SIGNAL), sd=sd(SIGNAL)*0.2)
#' Z <- outer(SIGNAL+NOISE, GRD)
#' Z.anom <- fieldAnomaly(y=Z, x=as.POSIXlt(TS), level="monthly")
#' 
#' zran <- range(Z, Z.anom)
#' pal <- colorRampPalette(c("blue", "cyan", "grey", "yellow", "red"))
#' op <- par(no.readonly=TRUE)
#' layout(matrix(c(1,2,3,3), nrow=2, ncol=2), widths=c(4,1), heights=c(2,2))
#' par(mar=c(5,5,3,1))
#' image(TS, GRD, Z, col=pal(100), zlim=zran, main="Original")
#' image(TS, GRD, Z.anom, col=pal(100), zlim=zran, main="Anomaly")
#' par(mar=c(5,0,3,3))
#' imageScale(Z, col=pal(100), zlim=zran, axis.pos=4)
#' par(op)
#' 
#' @export
#' 
fieldAnomaly <- function(y, x, level="daily"){ 

	y <- as.matrix(y)

	if(level=="monthly"){levs=unique(x$mon)}
	if(level=="daily"){levs=unique(x$yday)}

	levs_lookup=vector("list", length(levs))
	names(levs_lookup) <- levs
	for(i in 1:length(levs)){
		if(level=="monthly"){levs_lookup[[i]] <- which(x$mon == names(levs_lookup[i]))}
		if(level=="daily"){levs_lookup[[i]] <- which(x$yday == names(levs_lookup[i]))}
	}

	for(j in seq(levs)){     #for every time level
    y[levs_lookup[[j]],] <- t(t(as.matrix(y[levs_lookup[[j]],])) - apply(as.matrix(y[levs_lookup[[j]],]), 2, mean, na.rm=TRUE))
	}

	y

}

