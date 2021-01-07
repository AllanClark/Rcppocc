\name{occSPATlogitConsBinom}
\alias{occSPATlogitConsBinom}

\title{occSPATlogitConsBinom - Implementation of the consensus model with \code{k}
subsets using a Gaussian approximation to a Polya Gamma distribuition.}

\description{
The function implements of the consensus model with k subsets using using a Gaussian approximation to a Polya Gamma distribuition.
}

\usage{
occSPATlogitConsBinom(detection.model, occupancy.model, spatial.model, so.data, prior, control, stat, k)
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

The \code{moran.cut} component is the number of eigenvalues of the Moran matrix to be used in the restricted spatial regression model. It is recommended to set this value equal to 0.1*(the number of sites in the analysis). (This argument is thus different to the 'stocc' usage).}

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
\item{stat}{The detection covariate information is aggregated across the number of surveys at a site. \code{stat} allows the user to choose the method of aggregation for the detection covariates. The default is set to \code{1}, which represents the mean, and the other option is 2, which represents the median of the covariates at the site.}

\item{k}{The number of independent subsets to split the sites in. Default is set to 3. Errors might occur if you set this value too large since as the number of subsets increase the number of data points in the subset data sets decreases!}
}
\details{In order to attempt to improve the efficiency of \eqn{occSPATlogitBinom}, we use the Consensus Monte Carlo method of Huang and Gelman (2005). This entails extracting \eqn{k} subsets of a users data and then applying \eqn{occSPATlogitBinom} to each of the subsets after adjusting the prior distribution of the parameters of the model apprpriately. The final posterior distribution of the regression effects results in a Gaussian distribution where the mean and covariance matrices are an appropriately weighted average of the posterior distributions obtained from each of the subsets.

Here the detection process is modelled as a Binomial random variable. Specifically we assume that \ifelse{html}{\out{y<sub>i</sub>|z<sub>i</sub>~Binom(K<sub>i</sub>,z<sub>i</sub>*p<sub>i</sub>)}}{} where \ifelse{html}{\out{y<sub>i</sub>}}{} and \ifelse{html}{\out{p<sub>i</sub>}}{} are ther number of surveys to site \eqn{i} and \ifelse{html}{\out{p<sub>i</sub>}}{} is the conditional detection probability for site \eqn{i}.
}
\value{
\item{ }{  - MCMC samples similar to \eqn{occSPATlogitBinom} are outputted, one for each subset of the data. A list of length \eqn{k} is returned. Within each list item another list is returned with the following named items: \eqn{locc}, \eqn{alpha}, \eqn{beta}, \eqn{theta}, \eqn{tau} and \eqn{psi}. \eqn{locc} is the sample mean occupancy probability for a subset.}
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
niter <- 100
nthin <- 1
#--------------------------------------------------------------------------

t1_logit <- proc.time()

spat_logitConsBinomPG <- occSPATlogitConsBinom(detection.model= ~ nspp, occupancy.model= ~ PC1 + PC2,
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
                                    percent_burn_in = 0.5),
                                    stat=1,
                                    k=3)

t2_logit <- proc.time()
RcppTimer <- t2_logit-t1_logit
RcppTimer #in seconds
#--------------------------------------------------------------------------

#the results for the first subset
names(spat_logitConsBinomPG[[1]])
dim(spat_logitConsBinomPG$alpha) #note the dimensions! i.e. row matrices
#--------------------------------------------------------------------------
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
