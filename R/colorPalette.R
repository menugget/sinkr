#' @title Color interpolation with uneven step size
#' @description Color ramp with differing number of steps between color levels. Wrapper for colorRamp
#' 
#' @param colors  colors to interpolate; must be a valid argument to \code{col2rgb()}.
#' @param bias a positive number. Higher values give more widely spaced colors at the high end.
#' @param space	a character string; interpolation in RGB or CIE Lab color spaces.
#' @param interpolate	use spline or linear interpolation.
#' @param alpha	 logical: should alpha channel (opacity) values should be returned? It is an error to give a true value if space is specified.
#' @param ...	arguments to pass to \code{\link[grDevices]{colorRamp}}.
#' 
#' @details This is a wrapper function for colorRampPalette. It allows for the
#' definition of the number of intermediate colors between the main colors.
#' Using this option one can stretch out colors that should predominate
#' the palette spectrum. Additional arguments of colorRampPalette can also
#' be added regarding the type and bias of the subsequent interpolation..
#' 
#' @keywords color
#' @export
#' @examples
#' op <- par(mfcol=c(2,1), omi=c(0.1,0.1,0.1,0.1), mai=c(1,0.2,0.2,0.2))
#' steps <- c("blue4", "cyan", "white", "yellow", "red4")
#' pal <- colorPalette(steps, space="rgb")
#' z=1:1000
#' imageScale(z, col=pal(41))
#' box()
#' steps <- c("blue4", "cyan", "white", "yellow", "red4")
#' pal <- colorPalette(steps, c(20,1,1,20), space="rgb")
#' z=1:1000
#' imageScale(z, col=pal(41))
#' box()
#' par(op)
#' 
colorPalette <- function(steps, n.steps.between=NULL, ...){
 
 if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
 if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
 
 fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
 RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
 RGB[,fill.steps] <- col2rgb(steps)
 
 for(i in which(n.steps.between>0)){
  col.start=RGB[,fill.steps[i]]
  col.end=RGB[,fill.steps[i+1]]
  for(j in seq(3)){
   vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
   RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
  }
 }
 
 new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
 pal <- colorRampPalette(new.steps, ...)
 return(pal)
}

