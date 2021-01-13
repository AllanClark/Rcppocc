\name{occSPATlogit}
\alias{occSPATlogit}

\title{occSPATlogit - single-season Bayesian spatial occupancy model using Polya-Gamma random variables.}

\description{
The function undertakes MCMC sampling of the parameters of a single-season Bayesian Restricted spatial regression model that uses Polya-Gamma random variables.

The arguments used are MOSTLY the same as those for the \code{spatial.occupancy} function in the \code{stocc} package although the \code{initial.values} argument is not included at present.
}

\usage{
occSPATlogit(detection.model, occupancy.model, spatial.model, so.data, prior, control)
}
\arguments{
\item{detection.model}{A formula object describing the detection portion of
the occupancy model. The variables described by the detection model are
located in the \code{visit} data frame of the \code{so.data}.}

\item{occupancy.model}{A formula object describing the fixed effects portion
of the spatial occupancy process. The variables described by the occupancy
model are located in the \code{site} data frame of an \code{so.data} object.}

\item{spatial.model}{A named list object describing the spatial component of
the occupancy process.

Currently the only possible model is the restricted spatial regression model. Thus, \code{spatial.model=list(threshold=, moran.cut=)} is the only form that is accepted at present. The \code{threshold} component is used the create neighborhoods in the RSR model.  All sites within distance \code{threshold} of site i are considered neighbors of site.
i.

The \code{moran.cut} component is the number of eigenvalues of the Moran matrix to be used in the restricted spatial regression model. It is recommended to set this value equal to 0.1*(the number of sites in the analysis). (This argument is thus different to the \code{stocc} usage).}

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

The
elements of the list must include: (1) \code{ndraws} is the total number of iterations
UNDERTAKEN for the MCMC sample, and (2) \code{percent_burn_in} is the percentage of 'ndraws' used as burn-in samples.

The argument entered should be between 0.0 and 1.0. No check however is done to assess whether the user has entered the correct argument here. Note: The user has to thin the sample afterwards (if required).

No thinning is undertaken.}
}

\details{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\value{
\item{ }{ - The MCMC samples are returned as a named list. The elements of the list are \code{alpha}, \code{beta}, \code{theta}, \code{tau_mat}, \code{real.occ} and \code{psi}.}
\item{ }{  - \code{alpha} = the regression effects of the detection process;}
\item{ }{  - \code{beta} = the regression effects of the occupancy process;}
\item{ }{  - \code{theta} = the regression effects of the spatial component;}
\item{ }{  - \code{tau_mat} = the spatial precision parameter;}
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
data("SpatSimData")
#--------------------------------------------------------------------------

#Below we use simulated data to fit a RSR occupancy model
xycords <- SimTable$xycords
psi <- SimTable$psi
z <- SimTable$z
alpha <- SimTable$alpha

Nvisits_surveyed <- SimTable$Nvisits_surveyed

Xmat <- SimTable$Xmat
Minv <- SimTable$Minv
Q <- SimTable$Q
Kmat <- SimTable$Kmat
Wmat <- SimTable$Wmat

pij <- SimTable$pij
Ysim <- SimTable$Ysim

nsites <- SimTable$nsites
Num_surveys <- SimTable$Num_surveys
N_maxvisits <- SimTable$N_maxvisits
#--------------------------------------------------------------------------


#occupancy dataframe at the subset of sites where surveys were performed
surveyIndex <- 1:Num_surveys #surveyIndex are the ids of the surveyed locations
siteIndex <- 1:nsites #ids for all sites considered

Site.Data <- as.data.frame( cbind( siteIndex,xycords, Xmat ) )
colnames(Site.Data)<- c("SiteName", "Longitude", "Latitude", "PC1", "PC2")
head(Site.Data)

#detection dataframe - only at surveyed locations
PAdata <- na.omit(c(t(Ysim))) #presence-absence data
SiteNamereps <- rep(surveyIndex, Nvisits_surveyed) #the site names
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
numSpatre<-188
#--------------------------------------------------------------------------

#mcmc settings
#These should be increased and are simply set to make it fairly fast to run!
niter <- 500
nthin <- 1
#--------------------------------------------------------------------------

t1_logit <- proc.time()

spat_logit <- occSPATlogit(detection.model= ~ nspp, occupancy.model= ~ PC1 + PC2,
                                 spatial.model = list(threshold = 0.36,
                                 moran.cut = numSpatre, rho=1),
                                 so.data = Make.so.data,
                                 prior = list(a.tau = 0.5, b.tau = 0.0005,
                                              tau=1,
                                              Q.d=sigma_inv_alpha_p,
                                              mu.d = alpha_m,
                                              Q.o=sigma_inv_beta_p,
                                              mu.o = beta_m),
                                 control = list(ndraws =niter,
                                 percent_burn_in = 0.5))

t2_logit <- proc.time()
RcppTimer <- t2_logit-t1_logit
RcppTimer #in seconds
#--------------------------------------------------------------------------

names(spat_logit)
dim(spat_logit$alpha) #note the dimensions! i.e. row matrices
#--------------------------------------------------------------------------

#not run
#require(coda)
#alpha_mcmc <- as.mcmc( t(spat_logit$alpha) )
#summary(alpha_mcmc)
#effectiveSize(alpha_mcmc) #the effective sample size
#
#par(mfrow=c(2,2))
#traceplot(alpha_mcmc)
#densplot(alpha_mcmc)
#autocorr.plot(alpha_mcmc)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
