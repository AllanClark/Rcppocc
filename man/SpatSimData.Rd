\name{SpatSimData}
\alias{SpatSimData}

\title{SpatSimData - Simulated single-season occupancy spatial data.}

\description{
A RData file that contains simulated single-season occupancy spatial data.
}

\usage{
data(SpatSimData)
}

\details{
The .RData file contains a number of objects all related to a simulated spatial occupancy data set. The data relates to 1881 grid-cells that span South Africa; 805 of which were surveyed a number of times (at least 3 times and at most 50 times) while the remaining were unsurveyed. Site specific covariate data exists (2 covariates) at all grid-cell while detection specific data (1 covariate) exists at 805 locations. The first 805 locations refers to the surveyed locations while the rest are unsurveyed locations.
}
\value{
\item{SimTable}{ - A list that contains a number of objects related to a spatial
occupancy data set.}
\item{SimTable$xycords}{ - the latitude and longitude of 1881 grid-cells in South Africa.}
\item{SimTable$psi}{ - true occupancy probabilities associated with the 1881 grid-cells.}
\item{SimTable$z}{ - the latitude and longitude of 1881 grid-cells in South Africa.}
\item{SimTable$alpha}{ - Detection regression effects (intercept and one covariate effect).}
\item{SimTable$Nvisits_surveyed}{ - The number of survey visits to 805 of the grid-cells.}
\item{SimTable$Q}{ - the spatial precision matrix. A 1881 by 1881 matrix.}
\item{SimTable$Minv}{ - t(K)\%*\%Q\%*\%K}
\item{SimTable$Xmat}{ - survey specific covariates. Two covariates are included.}
\item{SimTable$Kmat}{ - a 1881 by 188 matrix. Referred to as K above. Used to fit a RSR spatial occupancy model. The columns relate to the number of eigenvalues of the Moran matrix retained.}
\item{SimTable$Wmat}{ - a 805 by 50 matrix that contains the detection covariate (only one covariate is included in the data set)}
\item{SimTable$pij}{ - the simulated detection probabilities.}
\item{SimTable$Ysim}{ - the simulated detection/non-detection data. This should be a n (number of surveys) by J (number of visits) matrix. Unequal survey-effort is allowed. Compare the number of surveyes at each site to the row-entries of this matrix.}
\item{SimTable$nsites}{ - the number of grid cells in the data set.}
\item{SimTable$Num_surveys}{ - the number of surveyed sites.}
\item{SimTable$N_maxvisits}{ - the maximum number of visits to a grid-cell.}
}
\examples{
#Below we use simulated data to fit a RSR occupancy model

data(SpatSimData)
#read the help file for the data ?SpatSimData

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
#-----------------------------------------------------

#occupancy dataframe at the subset of sites where surveys were performed
surveyIndex <- 1:Num_surveys #surveyIndex are the ids of the surveyed locations
siteIndex <- 1:nsites #ids for all sites considered

Site.Data <- as.data.frame( cbind( siteIndex, xycords,Xmat ) )
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
#-----------------------------------------------------

#Set the prior distributions used
nalphas<-2
nbetas<-3
beta_m<-matrix( rep(0,nbetas), ncol=1)
sigma_inv_beta_p<-diag(nbetas)/1000 #prior inverse covariance for beta
alpha_m<-matrix( rep(0,nalphas), ncol=1)
sigma_inv_alpha_p<-diag(nalphas)/1000 #prior inverse covariance for alpha
numSpatre<-188
#-----------------------------------------------------

#mcmc settings
#These should be increased and are simply set to make it fairly fast to run!
niter <- 500
nthin <- 1
#-----------------------------------------------------

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
                                 control = list(ndraws =niter, percent_burn_in = 0.5))

t2_logit <- proc.time()
RcppTimer <- t2_logit-t1_logit
RcppTimer
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
