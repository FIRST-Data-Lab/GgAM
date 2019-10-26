######################################################################
# Project:
# Algorithm: Group Coordinate Descent Algorithm
# Function: LM_ALASSO
# Author: Xinyi Li
# Date:
#
# Input Arguments:
#	(1) y: The response vector.
#	(2) x: The design matrix of linear part; need centered for
# 			continuous x.
#	(3) nu: Tuning parameter for mBIC/eBIC in the adaptive part of
# 			LASSO. Default is 1.
#	(4) crite: A character (vector) string specifying the crite
#			for variables selection. "BIC" to use traditional bayesian
#			information criteria; "mBIC" to use modified BIC for ultra-
#			high-dimension; "eBIC" to use extended BIC for ultra-
#			high-dimension. The default is "c("BIC","mBIC","eBIC")".
#	(5) lam: The sequence of lambda values in the adative part of
#			LASSO. Default is a grid of 20 lambda values that ranges
#			uniformly on the log scale over 0.01 to 2.
#	(6) lam1: The sequence of lambda values in adative part of LASSO.
#			Default is same with lam.
#	(7) nu1: Tuning parameter for mBIC/EBIC in the initial part of
#			LASSO. Default is the same value with nu.
#	(8) group: A vector describing the grouping of the coefficients.
#
# Output Arguments:
#	(1) theta.bic: Coefficients estimated under BIC.
#	(2) theta.mbc: Coefficients estimated under modified BIC.
#	(3) theta.ebc: Coefficients estimated under extended BIC.

######################################################################
# Library needed
#' @importFrom  grpreg grpreg
#' @importFrom mgcv uniquecombs
######################################################################
LM_ALASSO=function(y,x,group=1:ncol(x),
	family, lam=exp(seq(log(.0001),log(0.01),length.out=20)),lam1=exp(seq(log(.0001),log(0.01),length.out=20))
	                                                            ,nu=1,MIpenalty="grALasso",	MIcrite='BIC',...) {

	if (length(y)!=nrow(x)) {
		stop("Sample sizes of the response and covariates do not match.")
	}
	x.unique=uniquecombs(x)
	n=nrow(x.unique)
	np=ncol(x)

	# This routine implements LASSO,to calculate the weights in ALASSO
	if (MIpenalty=="grALasso") {
		# if (missing(lam1)) {
		# 	fit1=grpreg(x,y,group=group,MIpenalty="grLasso",
		# 		family=family[[1]])
		# 	lam1=fit1$lambda
		# } else {
			fit1=grpreg(x,y,group=group,MIpenalty="grLasso",
				family=family[[1]],lambda=lam1)
		# }
		df1=colSums(fit1$beta!=0)

		# choose lambda based on criteria
		if ("BIC" %in% MIcrite) {
			IC1=log(fit1$loss)+log(n)*df1/n
		}
		if ("mBIC" %in% MIcrite) {
			IC1=log(fit1$loss)+nu*df1*log(np)*log(n)/(2*n)
		}
		if ("EBIC" %in% MIcrite) {
			IC1=log(fit1$loss)+log(n)*df1/n+nu*df1*log(np)/n
		}
		theta1=fit1$beta[,which.min(IC1)]
		# plot(log(fit1$lambda),IC1,type="l",main="Group LASSO",
		# 	xlab=expression(paste("log(",lambda,")")),ylab=MIcrite)
		# points(log(fit1$lambda)[which.min(IC1)],min(IC1),col=2)

		# Obtain weights
		nor1=unlist(lapply(split(theta1[-1],group),
			function(x) sqrt(sum(x^2))))
		if (min(group)==0) {
			weight=ifelse(nor1[-1]==0,10^10^2.2,
				1/nor1[-1]*table(group)[-1])
		} else {
			weight=ifelse(nor1==0,10^10^2.2,1/nor1*table(group))
		}

		# This routine implements adaptive LASSO
		# if (missing(lam)) {
		# 	fit=grpreg(x,y,group=group,MIpenalty="grLasso",
		# 		family=family[[1]],group.multiplier=weight)
		# 	lam=fit$lambda
		# } else {
			fit=grpreg(x,y,group=group,MIpenalty="grLasso",
				family=family[[1]],group.multiplier=weight,lambda=lam)
		# }
	} else if (MIpenalty=="grSCAD") { # Routine for (group) SCAD
		if (missing(lam)) {
			fit=grpreg(x,y,group=group,MIpenalty="grSCAD",
				family=family[[1]])
			lam=fit$lambda
		} else {
			fit=grpreg(x,y,group=group,MIpenalty="grSCAD",
				family=family[[1]],lambda=lam)
		}
	}
	df=colSums(as.matrix(fit$beta[-1,])!=0)

	if ("BIC" %in% MIcrite) {
		IC=log(fit$loss)+log(n)*df/n
	}
	if ("mBIC" %in% MIcrite) {
		IC=log(fit$loss)+nu*df*log(np)*log(n)/(2*n)
	}
	if ("EBIC" %in% MIcrite) {
		IC=log(fit$loss)+log(n)*df/n+nu*df*log(np)/n
	}
	theta=fit$beta[,which.min(IC)]
	# plot(log(fit$lambda),IC,type="l",main=MIpenalty,
	# 	xlab=expression(paste("log(",lambda,")")),ylab=MIcrite)
	# points(log(fit$lambda)[which.min(IC)],min(IC),col=2)

	list(theta.hat=theta,group=group,MIpenalty=MIpenalty,family=family,
		lambda=lam,nu=nu,lambda1=lam1)
}
