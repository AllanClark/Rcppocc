\name{PGocc2}
\alias{PGocc2}

\title{PGocc2 - single-season Bayesian occupancy model using Polya-Gamma random variables. The arguments of this function are different to those of PGocc and uses a similar style to that used in the 'stocc' package. THIS PAGE HAS TO BE UPDATED!}

\description{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian occupancy model that uses Polya-Gamma random variables.

The arguments used are MOSTLY the same as those for the 'spatial.occupancy' function in the 'stocc' package although the 'initial.values' argument is not included at present.
}

\usage{
PGocc2(detection.model, occupancy.model, so.data, prior, control, store_z)
}
\arguments{
  \item{formula }{ - Double right-hand side formula describing covariates of detection and occupancy in that order.}
  \item{data_inputs }{ - A named list that contains the elements W, X, and y.

  W - A named list of data.frames of covariates that vary within sites. i.e. The dataframes are of dimension \code{n} (number of sites surveyed) by \code{J} where each row is associated with a site and each column represents a site visit. 'NA' values should be entered at locations that were not surveyed.

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
  \item{percent_burn_in }{ - The percentage of 'ndraws' used as burn-in samples. The argument entered should be between 0.0 and 1.0. No check however is done to assess whether the user has entered the correct argument here. }
  \item{store_z }{ - A logical. 'store_z=TRUE' indicates that the 'z' and 'psi' samples are stored. Note that these are not stored by default. i.e. 'store_z=FALSE' is the default value for this argument.}
}
\details{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian occupancy model that uses Polya-Gamma random variables.
}
\value{
\item{ }{  - The MCMC samples are returned as a named list with entries 'alpha', 'beta', 'z' and 'psi'. 'alpha', 'beta' are matrices where the number of rows equals to the number of regression effects and the number of columns equals to number of retained MCMC samples. The dimension of 'z' and 'psi' are 'number of surveyed sites' by 'number of retained MCMC samples'.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%% ~~objects to See Also as \code{\link{help}}, ~~~
data("SpatSimData")
#--------------------------------------------------------------------------

#Below we use simulated data to fit a RSR occupancy model
#occupancy dataframe at the subset of sites where surveys were performed
Num_surveys <- SimTable$Num_surveys
surveyIndex <- 1:Num_surveys #surveyIndex are the ids of the surveyed locations
nsites <- length(surveyIndex) #only include the surveyed sites
siteIndex <- 1:nsites #ids for all sites considered
xycords <- SimTable$xycords[surveyIndex,]
Xmat <- SimTable$Xmat[surveyIndex,] #only include the surveyed sites

Site.Data <- as.data.frame( cbind( siteIndex, xycords, Xmat ) )
colnames(Site.Data)<- c("SiteName", "Longitude", "Latitude", "PC1", "PC2")
head(Site.Data)

#detection dataframe - only at surveyed locations
PAdata <- na.omit(c(t(Ysim))) #presence-absence data
Nvisits_surveyed <- SimTable$Nvisits_surveyed
SiteNamereps <- rep(surveyIndex, Nvisits_surveyed) #the site names

Wmat <- SimTable$Wmat

nspp <- na.omit(c(t(Wmat))) #detection covariate
Longitudereps <- rep(xycords[surveyIndex,1], Nvisits_surveyed)
Latitudereps <- rep(xycords[surveyIndex,2], Nvisits_surveyed)

#at surveyed locations
Visit.data <- data.frame(SiteNamereps, Longitudereps, Latitudereps, nspp, PAdata)
colnames(Visit.data) <- c("SiteName", "Longitude", "Latitude", "nspp", "PAdata")

Names <- list(visit = list(site = "SiteName", obs = "PAdata"),
              site = list(site = "SiteName", coords = c("Longitude", "Latitude")) )

#uncomment when running the example
require(stocc) #loaded the stocc package!
Make.so.data <- make.so.data(visit.data = Visit.data, site.data = Site.Data, names = Names)
#--------------------------------------------------------------------------

#Set the prior distributions used
nalphas<-2
nbetas<-3
beta_m<-matrix( rep(0,nbetas), ncol=1)
sigma_inv_beta_p<-diag(nbetas)/1000 #prior inverse covariance for beta
alpha_m<-matrix( rep(0,nalphas), ncol=1)
sigma_inv_alpha_p<-diag(nalphas)/1000 #prior inverse covariance for alpha
#--------------------------------------------------------------------------

#mcmc settings
#These should be increased and are simply set to make it fairly fast to run!
niter <- 1000
nthin <- 1
#--------------------------------------------------------------------------

occ <- PGocc2(detection.model= ~ nspp,
                           occupancy.model= ~ PC1 + PC2,
                           so.data = Make.so.data,
                           prior = list(Q.d=sigma_inv_alpha_p,
                                        mu.d = alpha_m,
                                        Q.o=sigma_inv_beta_p,
                                        mu.o = beta_m),
                           control = list(ndraws =niter,
                                          percent_burn_in = 3/7))
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

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
