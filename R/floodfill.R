#' @title Flood fill a plot region
#' @description A function to "flood fill" a region on the active plotting 
#' device. 
#' Once called, the user will be asked to click on the desired target region. 
#' The flood fill algorithm then searches neighbors in 4 directions of the 
#' target cell and checks for similar colors to the target cell. 
#' If neighboring cells are of the same color, they are added to a "queue" for 
#' further searches of their neighbors. Once a cell hace been checked, its 
#' position is added to a "done" or completed list. This algorithm is referred 
#' to as "Four-way floodfill using a queue for storage"
#' 
#' @param replCol The color to apply to the flood region. Should be given 
#' as a hexidecimal string, as is the format of the output from
#' \code{rgb()}.
#' @param res The resolution (per inch) of the temporarily exported figure (7 X 7 inches) 
#' used to create a matrix of color values for the flood fill algorithm.  
#' @examples
#' #install.packages("png")
#' library(png)
#' op <- par(mar=c(4,4,1,2), bg="white")
#' x <- 1:10
#' y <- 2*x
#' plot(x, y, t="l", col=2, xaxs="i", yaxs="i")
#' abline(10,-2)
#' abline(10,2)
#' abline(20,-4)
#' abline(20,-1)
#' floodfill(replCol=rgb(0.5,0.7,0.7,1), res=50) # Choose a 1st area
#' floodfill(replCol=rgb(t(col2rgb("pink", alpha=TRUE)), maxColorValue = 255), res=50)  # Choose a 2nd area
#' floodfill(replCol=rgb(0.9,0.7,0.3,1), res=50)  # Choose a 3rd area
#' floodfill(replCol=rgb(0,1,1,1), res=50)  # Choose a 3rd area
#' par(op)
#' 
#' # Export graphics device
#' dev.print(file="floodfill_ex.png", device = png, width=5, height=5, units="in", res=200, type="cairo")
#' 
#' @export
#'
#' 
floodfill <- function(replCol=rgb(0.5,0.7,0.7,1), res=50){
	#####################
	### Required package
	require(png)

	#####################
	### Choose target position
	print("Choose a region to flood fill")
	pos <- locator(1)
	Pars <- par()
	dev.print(file="4fill.png", device = png, width=7, height=7, units="in", res=50)
	mat <- readPNG("4fill.png")
	dim(mat)
	
	#####################
	### Make Colors matrix ('Col')
	Col <- rgb(mat[,,1], mat[,,2], mat[,,3], mat[,,4])
	Col <- array(Col, dim=dim(mat)[1:2])
	Col <- t(Col)
	Col <- Col[,dim(Col)[2]:1]
	
	#####################
	### Trim Colors matrix ('Col') to only plot area
	rows <- round(Pars$mai[2]/7*nrow(Col)):round((7-Pars$mai[4])/7*nrow(Col))
	cols <- round(Pars$mai[1]/7*ncol(Col)):round((7-Pars$mai[3])/7*ncol(Col))
	Col <- Col[rows, cols]
	
	#####################
	### Make lookup table
	rowVal <- seq(Pars$usr[1], Pars$usr[2],, nrow(Col))
	colVal <- seq(Pars$usr[3], Pars$usr[4],, ncol(Col))
	tarPos <- c(which.min((rowVal-pos$x)^2), which.min((colVal-pos$y)^2))
	tarCol <- Col[tarPos[1], tarPos[2]]
	#replCol <- rgb(0,0.5,0.5,1)
	grd <- cbind(posi=seq(Col), expand.grid(row=seq(nrow(Col)), col=seq(ncol(Col)), done=0))
	grd$Col <- c(Col)
	
	#Get neighbors
	posi <- array(grd$posi, dim=dim(Col))
	#down
	down <- posi*NaN
	down[-nrow(posi),] <- posi[-1,]
	grd$down <- c(down)
	#left
	left <- posi*NaN
	left[,-1] <- posi[,-ncol(posi)]
	grd$left <- c(left)
	#up
	up <- posi*NaN
	up[-1,] <- posi[-nrow(posi),]
	grd$up <- c(up)
	#right
	right <- posi*NaN
	right[,-ncol(posi)] <- posi[,-1]
	grd$right <- c(right)
	
	#####################
	### Search for similar colors to target cell and record "done" or checked cells
	queue <- which(grd$row==tarPos[1] & grd$col==tarPos[2])
	grd$Col[queue] <- replCol
	count <- 1
	pb <- txtProgressBar(min = count, max = nrow(grd), initial = count, style=3)
	while(length(queue) > 0){
	  queue2 <- vector(mode="list", length(queue))
		setTxtProgressBar(pb, count)
		for(i in queue){
	    incl <- which(c(grd$down[i], grd$left[i], grd$up[i], grd$right[i]) > 0) 
			check <- c(grd$down[i], grd$left[i], grd$up[i], grd$right[i])[incl]
			check <- check[grd$done[check] == 0]
	
			repl <- check[which(grd$Col[check] == tarCol)]
			grd$Col[repl] <- replCol
			queue2[[count]] <- repl
			count <- count + 1
			grd$done[i] <- 1
		}
		queue <- unlist(queue2)
		queue <- queue[grd$done[queue] == 0]
	}
	
	#####################
	### Add flood fill layer to plot region
	Fill <- array(0, dim(Col))
	Fill[which(grd$Col == replCol)] <- 1
	image(x=rowVal, y=colVal, z=Fill, add=TRUE, useRaster=TRUE, col=c(0, replCol))
	
}

