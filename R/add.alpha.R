#' @title add.alpha (add alpha channel - transparency - to colors)
#' @description Takes a vector of colors and adds an alpha channel 
#' at the given level of transparency
#' @export

add.alpha <- function(COLORS, ALPHA){
	if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
	RGB <- col2rgb(COLORS, alpha=TRUE)
	RGB[4,] <- round(RGB[4,]*ALPHA)
	NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
	return(NEW.COLORS)
}