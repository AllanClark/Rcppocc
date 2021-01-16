\name{PGocc2}
\alias{PGocc2}

\title{PGocc2 - single-season Bayesian occupancy model using Polya-Gamma random variables. The arguments of this function are different to those of PGocc and uses a similar style to that used in the \code{stocc} package. }

\description{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian occupancy model that uses Polya-Gamma random variables.

The arguments used are MOSTLY the same as those for the \code{spatial.occupancy} function in the \code{stocc} package although the \code{initial.values} argument is not included at present.
}

\usage{
PGocc2(detection.model, occupancy.model, so.data, prior, control, store_z)
}
\arguments{
  \item{detection.model}{A formula object describing the detection portion of
  the occupancy model. The variables described by the detection model are
  located in the \code{visit} data frame of the \code{so.data}.}

  \item{occupancy.model}{A formula object describing the fixed effects portion
  of the spatial occupancy process. The variables described by the occupancy
  model are located in the \code{site} data frame of an \code{so.data} object.}

  \item{so.data}{An \code{so.data} object containing the observed occupancies, detection covariates, site covariates, and site coordinates. This is created via the \code{\link{make.so.data}}.}

  \item{prior}{A named list that provides the parameter values for the prior distributions. At the current time the elements of the list MUST contain
  \code{a.tau}, \code{b.tau}, \code{tau}, \code{Q.d}, \code{mu.d}, \code{Q.o} and \code{mu.o}.

  \code{a.tau} and \code{b.tau} are the parameters for the gamma prior on the spatial
  process parameter in the occupancy model;

  \code{tau} is a starting value for the spatial precision parameter;

  \code{Q.d} and \code{mu.d} are the prior precision matrix and mean vector for the detection regression effects while \code{Q.o} and \code{mu.o} are similar inputs for the occupancy process.

  All of these arguments must be supplied by the user at present.

  At present, no checks are made to ensure that the correct inputs are supplied! BE CAREFUL!}

  \item{control}{A named list with the control parameters for the MCMC.

  The elements of the list must include: (1) \code{ndraws} is the total number of iterations
  UNDERTAKEN for the MCMC sample, and (2) \code{percent_burn_in} is the percentage of 'ndraws' used as burn-in samples.

  The argument entered should be between 0.0 and 1.0. No check however is done to assess whether the user has entered the correct argument here. Note: The user has to thin the sample afterwards (if required).

  No thinning is undertaken.}
  \item{store_z }{ - A logical. \code{store_z=TRUE} indicates that the \code{z} and \code{psi} samples are stored. Note that these are not stored by default. i.e. \code{store_z=FALSE} is the default value for this argument.}
}
\details{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian occupancy model that uses Polya-Gamma random variables.

The arguments used are MOSTLY the same as those for the \code{spatial.occupancy} function in the \code{stocc} package although the \code{initial.values} argument is not included at present.
}
\value{
\item{ }{  - The MCMC samples are returned as a named list with entries \code{alpha} and \code{beta} (\code{real.occ} and \code{psi} are optional).}
\item{}{  - \code{alpha}, \code{beta} are matrices where the number of rows equals to the number of regression effects and the number of columns equals to number of retained MCMC samples. The dimension of \code{real.occ} and \code{psi} are \code{number of surveyed sites} by \code{number of retained MCMC samples}.}
\item{ }{  - \code{alpha} = the regression effects of the detection process;}
\item{ }{  - \code{beta} = the regression effects of the occupancy process;}
\item{ }{  - \code{real.occ} = the estimate of the true occupancy variable;}
\item{ }{  - \code{psi} = conditional occupancy probabilities.}
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
Ysim <- SimTable$Ysim
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
