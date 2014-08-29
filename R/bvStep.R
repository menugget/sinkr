#' Clarke and Ainsworth's BVSTEP routine
#' 
#' The \code{bvStep} function performs Clarke and Ainsworth's (1993) "BVSTEP" routine which is a algorithm that searches variable 
#' subsets of a matrix before calculating a dissimilarity matrix which is compared to a dissimilarity 
#' matrix derived from a fixed matrix (via a Mantel test). The test is the same as that performed by the 
#' \code{\link[sinkr]{bioEnv}} function but is more appropriate when the number of variable combinations to test 
#' is large
#' 
#' @param fix.mat The "fixed" matrix of community or environmental sample by variable values
#' @param var.mat A "variable" matrix of community or environmental sample by variable values
#' @param fix.dist.method The method of calculating dissimilarity indices bewteen samples in the fixed
#' matrix (Uses the \code{\link[vegan]{vegdist}} function from the vegan package to calculate distance matrices. See 
#' the documentation for available methods.). Defaults to Bray-Curtis dissimularity \code{"bray"}.
#' @param var.dist.method The method of calculating dissimilarity indices bewteen samples in the variable
#' matrix. Defaults to Euclidean dissimularity \code{"euclidean"}.
#' @param scale.fix Logical. Should fixed matrix be centered and scaled (Defaults to \code{FALSE}, 
#' recommended for biologic data).
#' @param scale.var Logical. Should fixed matrix be centered and scaled (Defaults to \code{TRUE}, 
#' recommended for environmental data to correct for differing units between variables).
#' @param output.best Number of best combinations to return in the results object (Default=10).
#' @param var.max Maximum number of variables to include. Defaults to all, \code{var.max=ncol(var.mat)}.
#' 
#' @details The R package "vegan" contains a version of Clarke and Ainsworth's (1993) 
#' BIOENV analysis allowing for the comparison of distance/similarity matrices between 
#' two sets of data having either samples or variables in common. 
#' The typical setup is in the exploration of environmental variables 
#' that best correlate to sample similarities of the biological community 
#' (e.g. species biomass or abundance). 
#' In this case, the similarity matrix of the community is fixed, while subsets of 
#' the environmental variables are used in the calculation of the environmental similarity matrix. 
#' A correlation coefficient (typically Spearman rank correlation coefficient) is then 
#' calculated between the two matrices and the best subset of environmental variables 
#' can then be identified and further subjected to a permutation test to determine significance.
#' Due to the inflexibility of the bioEnv() function, one has little control over how the variable similarity matrix is calculated (derived from the environmental subsets in the above example) as the method assumes the subset data to be environmental and that the resulting similarity matrix should be based on normalized "euclidean" distances. This makes sense with environmental data where one normalizes the data to remove the effect of differing units between parameters, yet in cases where the variable matrix is biological one might want more flexibility (a Bray-Curtis measure of similarity is common given its non-parametric nature). The vegan function vegdist() comes with many other possible indices that could be applied ("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial" and "chao"). For example, beyond the typical biological to environmental comparison (BIOENV setup), one can also use the routine to explore other other types of relationships; e.g.:
#' - ENVBIO: subset of biological variables that best correlate to the overall environmental pattern
#' - BIOBIO: subset of biological variables that best correlate to the overall biological pattern
#' - ENVENV: subset of environmental variables that best correlate to the overall environmental pattern

#' 
#' @examples
#' \donttest{
#' 
#' library(vegan)
#' data(varespec)
#' data(varechem)
#' 
#' # Example of a 2-round BIO-BIO search. Uses the most frequently included variables
#' # in the first round at the beginning of each restart in the second round
#' # first round
#' set.seed(1)
#' res.biobio1 <- bvStep(wisconsin(varespec), wisconsin(varespec), 
#'  fix.dist.method="bray", var.dist.method="bray",
#'  scale.fix=FALSE, scale.var=FALSE, 
#'  max.rho=0.95, min.delta.rho=0.001,
#'  random.selection=TRUE,
#'  prop.selected.var=0.3,
#'  num.restarts=50,
#'  output.best=10,
#'  var.always.include=NULL
#' )
#' res.biobio1 # Best rho equals 0.833 (10 of 44 variables)
#' 
#' #second round - always includes variables 23, 26, and 29 ("Cla.ran" "Cla.coc" "Cla.fim")
#' set.seed(1)
#' res.biobio2  <- bvStep(wisconsin(varespec), wisconsin(varespec), 
#'  fix.dist.method="bray", var.dist.method="bray",
#'  scale.fix=FALSE, scale.var=FALSE, 
#'  max.rho=0.95, min.delta.rho=0.001,
#'  random.selection=TRUE,
#'  prop.selected.var=0.3,
#'  num.restarts=50,
#'  output.best=10,
#'  var.always.include=c(23,26,29)
#' )
#' res.biobio2 # Best rho equals 0.895 (15 of 44 variables)
#' 
#' # A plot of best variables
#' MDS_res=metaMDS(wisconsin(varespec), distance = "bray", k = 2, trymax = 50)
#' bio.keep <- as.numeric(unlist(strsplit(res.biobio2$order.by.best$var.incl[1], ",")))
#' bio.fit <- envfit(MDS_res, varespec[,bio.keep], perm=999)
#' bio.fit 
#' 
#' plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2")
#' plot(bio.fit, col="gray50", cex=0.8, font=4) # display only those with p>0.1
#' text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=0.7)
#' mtext(paste("Stress =",round(MDS_res$stress, 2)), side=3, adj=1, line=0.5)
#' 
#' # Display only those with envfit p >= 0.1
#' plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2")
#' plot(bio.fit, col="gray50", p.max=0.1, cex=0.8, font=4) # p.max=0.1
#' text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=0.7)
#' mtext(paste("Stress =",round(MDS_res$stress, 2)), side=3, adj=1, line=0.5)
#' 
#' }
#' 
#' @export
#' 
bvStep <- function(fix.mat, var.mat, 
fix.dist.method="bray", var.dist.method="euclidean",
scale.fix=FALSE, scale.var=TRUE,
max.rho=0.95,
min.delta.rho=0.001,
random.selection=TRUE,
prop.selected.var=0.2,
num.restarts=10,
var.always.include=NULL,
var.exclude=NULL,
output.best=10
){

	if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
	if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
	require(vegan)

	if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
	if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}

	fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)

	#an initial removal phase
	var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
	full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method="spearman"))$estimate
	var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
	RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
	for(i in 1:dim(var.comb)[2]){
		var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
		temp <- suppressWarnings(cor.test(fix.dist, var.dist, method="spearman"))
		RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
		RES$rho[i] <- temp$estimate
	}
	delta.rho <- RES$rho - full.cor
	exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))

	if(random.selection){
		num.restarts=num.restarts
		prop.selected.var=prop.selected.var
		prob<-rep(1,ncol(var.mat))
		if(prop.selected.var< 1){
			prob[exclude]<-0
		}
		n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
	} else {
		num.restarts=1
		prop.selected.var=1		
		prob<-rep(1,ncol(var.mat))
		n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
	}

	RES_TOT <- c()
	for(i in 1:num.restarts){
		step=1
		RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
		attr(RES$step.dir, "levels") <- c("F","B")
		best.comb <- which.max(RES$rho)
		best.rho <- RES$rho[best.comb]
		delta.rho <- Inf
		selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
		while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
			#forward step
			step.dir="F"
			step=step+1
			var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
			if(RES$n.var[best.comb] == 0){
				var.comb.incl<-1:length(var.comb)
			} else {
				var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
				temp <- NA*1:length(var.comb)
				for(j in 1:length(temp)){
					temp[j] <- all(var.keep %in% var.comb[[j]]) 
				}
				var.comb.incl <- which(temp==1)
			}

			RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
			for(f in 1:length(var.comb.incl)){
				var.incl <- var.comb[[var.comb.incl[f]]]
				var.incl <- var.incl[order(var.incl)]
				var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
				temp <- suppressWarnings(cor.test(fix.dist, var.dist, method="spearman"))
				RES.f$var.incl[f] <- paste(var.incl, collapse=",")
				RES.f$rho[f] <- temp$estimate
			}

			last.F <- max(which(RES$step.dir=="F"))
			RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
			best.comb <- which.max(RES$rho)
			delta.rho <- RES$rho[best.comb] - best.rho 
			best.rho <- RES$rho[best.comb]
			
			if(best.comb == step){
				while(best.comb == step & RES$n.var[best.comb] > 1){
					#backward step
					step.dir="B"
					step <- step+1
					var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
					var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
					RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
					for(b in 1:length(var.comb)){
						var.incl <- var.comb[[b]]
						var.incl <- var.incl[order(var.incl)]
						var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
						temp <- suppressWarnings(cor.test(fix.dist, var.dist, method="spearman"))
						RES.b$var.incl[b] <- paste(var.incl, collapse=",")
						RES.b$rho[b] <- temp$estimate
					}
					RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
					best.comb <- which.max(RES$rho)
					best.rho<- RES$rho[best.comb]
				}
			} else {
				break()
			}

		}

		RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
		print(paste(round((i/num.restarts)*100,3), "% finished"))
	}
	
	RES_TOT <- unique(RES_TOT[,3:5])


	if(dim(RES_TOT)[1] > output.best){
		order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
	} else {
		order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
	}
	rownames(order.by.best)<-NULL

	order.by.i.comb <- c()
	for(i in 1:length(selected.var)){
		f1 <- which(RES_TOT$n.var==i)
		f2 <- which.max(RES_TOT$rho[f1])
		order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
	}
	rownames(order.by.i.comb)<-NULL
	
	if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
	out <- list(
		order.by.best=order.by.best,
		order.by.i.comb=order.by.i.comb,
		best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
		best.model.rho=order.by.best$rho[1],
		var.always.include=var.always.include,
		var.exclude=var.exclude
		)
	out

}

			

