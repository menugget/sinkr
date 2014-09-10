
getcolors <- function(n){
	N <- 6
	
	X <- seq(N^2)
	Y <- seq(N)
	GRD <- expand.grid(x=X, y=Y)
	Z <- matrix(0, nrow=length(X), ncol=length(Y))

	LEV <- seq(0,1,,N) 
	R <- matrix(rep(LEV, each=N^2), nrow=length(X), ncol=length(Y))
	G <- matrix(rep(rep(LEV, each=N), N), nrow=length(X), ncol=length(Y))
	B <- matrix(rep(LEV, N^2), nrow=length(X), ncol=length(Y))


	x11(width=6, height=6)
	layout(matrix(1:2, nrow=2, ncol=1), widths=c(6), heights=c(1.5,4.5))
	op <- par(mar=c(1,3,2,1))

	image(X,Y,Z, col=NA, xlab="", ylab="", xaxt="n", yaxt="n")
	for(i in seq(nrow(GRD))){
		xs <- c(GRD$x[i]-0.5, GRD$x[i]-0.5, GRD$x[i]+0.5, GRD$x[i]+0.5)
		ys <- c(GRD$y[i]-0.5, GRD$y[i]+0.5, GRD$y[i]+0.5, GRD$y[i]-0.5)
		polygon(xs, ys, col=rgb(R[i], G[i], B[i]), border=NA)
	}
	mtext(paste("Click on", n, "colors [please]"), side=3, line=0.5)
	box()

	COLS <- NA*seq(n)
	for(i in seq(n)){
		coord <- locator(1)
		COLS[i] <- rgb(R[round(coord$x), round(coord$y)], G[round(coord$x), round(coord$y)], B[round(coord$x), round(coord$y)])
	}

	par(mar=c(1,3,0,1))
	pal <- colorRampPalette(c("black", "white"))
	image(x=1:100, y=seq(n), z=matrix(rep(1:100,n), nrow=100, ncol=n), col=pal(100), xlab="", ylab="", xaxt="n", yaxt="n")
	box()
	for(i in seq(n)){
		lines(x=c(1,100), y=c(i,i), col=COLS[i], lwd=4)
	}
	axis(2, at=seq(n))

	par(op)

	COLS
}


