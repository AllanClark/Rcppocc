\name{PGocc}
\alias{PGocc}

\title{PGocc - single-season Bayesian occupancy model using Polya-Gamma random variables.}

\description{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian occupancy model that uses Polya-Gamma random variables.
}

\usage{
PGocc(formula, data_inputs, ndraws, alpha_m, beta_m,
sigma_inv_alpha_p, sigma_inv_beta_p, percent_burn_in, store_z)
}
\arguments{
  \item{formula }{ - Double right-hand side formula describing covariates of detection and occupancy in that order.}
  \item{data_inputs }{ - A named list that contains the elements W, X, and y.

  W - A named list of data.frames of covariates that vary within sites. i.e. The dataframes are of dimension \code{n} (number of sites surveyed) by \code{J} where each row is associated with a site and each column represents a site visit. \code{NA} values should be entered at locations that were not surveyed.

e.g. Suppose \code{W} contained three dataframes \code{W1}, \code{W2} and \code{W3}.

\code{W$W1[1, ]} = the covariate values for site 1 for all of the visits. Note that some of the entries might be 'NA' meaning that no visit took place at those occasions.

 X - A named data.frame that varies at site level.

 y - An \code{n} by \code{J} matrix of the detection, non-detection data, where \code{n} is the number of sites, \code{J} is the maximum number of sampling periods per site.

}
  \item{ndraws }{ - The number of MCMC samples undertaken. This is set as the number of burn-in iterations plus the number of post burn-in samples. See 'percent_burn_in' below.}
  \item{alpha_m }{ - The prior mean alpha vector. Supply a matrix of the correct dimension.}
  \item{beta_m }{ - The prior mean beta vector. Supply a matrix of the correct dimension.}
  \item{sigma_inv_alpha_p }{ - The alpha precision matrix (i.e. the inverse of the covariance matrix).}
  \item{sigma_inv_beta_p }{ - The beta precision matrix.}
  \item{percent_burn_in }{ - The percentage of \code{ndraws} used as burn-in samples. The argument entered should be between 0.0 and 1.0. No check however is done to assess whether the user has entered the correct argument here. }
  \item{store_z }{ - A logical. \code{store_z=TRUE} indicates that the \code{z} and \code{psi} samples are stored. Note that these are not stored by default. i.e. \code{store_z=FALSE} is the default value for this argument.}
}
\details{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian occupancy model that uses Polya-Gamma random variables.
}
\value{
\item{ }{  - The MCMC samples are returned as a named list with entries \code{alpha} and \code{beta} (\code{real.occ} and \code{psi} are optional).}
\item{}{  - \code{alpha}, \code{beta} are matrices where the number of rows equals to the number of regression effects and the number of columns equals to number of retained MCMC samples. The dimension of \code{z} and \code{psi} are \code{number of surveyed sites} by \code{number of retained MCMC samples}.}
\item{ }{  - \code{alpha} = the regression effects of the detection process;}
\item{ }{  - \code{beta} = the regression effects of the occupancy process;}
\item{ }{  - \code{real.occ} = the estimate of the true occupancy variable;}
\item{ }{  - \code{psi} = conditional occupancy probabilities.}
}
\examples{
%% ~~objects to See Also as \code{\link{help}}, ~~~
data("SpatSimData")
#--------------------------------------------------------------------------

Xmat <- SimTable$Xmat[1:805,] #only use surveyed sites (805)
colnames(Xmat) #PC1, PC2

Wmat <- list( W1=SimTable$Wmat )
names(Wmat) #W1

Ysim <- SimTable$Ysim

data_list <- list(W=Wmat, X=Xmat, y=Ysim)
names(data_list$Y) # V1
#--------------------------------------------------------------------------

#specification of the prior distributions
beta_m<-matrix(0, nrow=3, ncol=1) #prior mean for beta
sigma_inv_beta<-diag(3)/1000 #prior covariance for beta

alpha_m <- matrix(0, nrow=2, ncol=1) #prior mean for alpha
sigma_inv_alpha <- diag(2)/1000 #prior covariance for alpha

# ~ occupancy ~ detection
#ndraws should be increased
occ <- PGocc(formula= V1  ~ PC1 + PC2 ~ W1,
            data_inputs = data_list,
            ndraws=1000,
            alpha_m=alpha_m, beta_m=beta_m,
            sigma_inv_alpha_p=sigma_inv_alpha,
            sigma_inv_beta_p=sigma_inv_beta,
            percent_burn_in=2/7)
#--------------------------------------------------------------------------

#not run
#require(coda)
#alpha_mcmc <- as.mcmc( t(spat_logit$alpha) ) #note the transpose
#summary(alpha_mcmc)
#effectiveSize(alpha_mcmc) #the effective sample size
#
#par(mfrow=c(2,2))
#traceplot(alpha_mcmc)
#densplot(alpha_mcmc)
#autocorr.plot(alpha_mcmc)
#run three chains to do convergence testing!
#--------------------------------------------------------------------------

attributes(occ)
dim(occ$alpha) #note the dimensions
dim(occ$beta) #note the dimensions

#some posterior stats
#--------------------------------------------------------------------------
apply(occ$alpha,1,mean)
apply(occ$alpha,1,sd)

apply(occ$beta,1,mean)
apply(occ$beta,1,sd)
#--------------------------------------------------------------------------

#Some trace plots
par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(occ$alpha[1,], main="", xlab=expression(alpha[0]), type="l", ylab="")
plot(occ$alpha[2,], main="", xlab=expression(alpha[1]), type="l", ylab="")

plot(occ$beta[1,], main="", xlab=expression(beta[0]), type="l", col="red", ylab="")
plot(occ$beta[2,], main="", xlab=expression(beta[1]), type="l", col="red", ylab="")
#--------------------------------------------------------------------------

#plot of the posterior densities
par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(density(occ$alpha[1,]), main="", xlab=expression(alpha[0]))
plot(density(occ$alpha[1,]), main="", xlab=expression(alpha[1]))

plot(density(occ$beta[1,]), main="", xlab=expression(beta[0]), col="red")
plot(density(occ$beta[1,]), main="", xlab=expression(beta[1]), col="red")
#--------------------------------------------------------------------------
}
