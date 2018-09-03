\name{vb_Designs}
\alias{vb_Designs}

\title{vb_Designs}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
vb_Designs(W, X, y)
}
\arguments{
  \item{W }{ - A named list of data.frames of covariates that vary within sites. i.e. The dataframes are of dimension \code{n} (number of sites surveyed) by \code{J} where each row is associated with a site and each column represents a site visit. 'NA' values should be entered at locations that were not surveyed.

e.g. Suppose \code{W} contained three dataframes \code{W1}, \code{W2} and \code{W3};

\code{W$W1[1, ]} = the covariate values for site 1 for all of the visits. Note that some of the entries might be 'NA' meaning that no visit took place at those occasions.}
  \item{X }{ - A named data.frame that varies at site level.}
  \item{y }{ - An \code{n} by \code{J} matrix of the detection, non-detection data, where \code{n} is the number of sites, \code{J} is the maximum number of sampling periods per site.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

Xmat <- SimTable$Xmat[1:805,] #only use surveyed sites (805)
colnames(Xmat) #PC1, PC2

Wmat <- list( W1=SimTable$Wmat )
names(Wmat) #W1

Ysim <- SimTable$Ysim
#--------------------------------------------------------------------------

design_list <- vb_Designs(W=Wmat, X=Xmat, y=Ysim)
names(design_list$Y) # V1

#specification of the prior distributions
beta_m<-matrix(0, nrow=3, ncol=1) #prior mean for beta
sigma_inv_beta<-diag(3)/1000 #prior covariance for beta

alpha_m <- matrix(0, nrow=2, ncol=1) #prior mean for alpha
sigma_inv_alpha <- diag(2)/1000 #prior covariance for alpha
#--------------------------------------------------------------------------

# ~ occupancy covariates ~ detection covariates
#fit a single-season model using Polya Gamma sampling
occ <- PGocc4(formula= V1  ~ PC1 + PC2 ~ W1,
            design_mats = design_list,
            ndraws=1000,
            alpha_m=alpha_m, beta_m=beta_m,
            sigma_inv_alpha_p=sigma_inv_alpha,
            sigma_inv_beta_p=sigma_inv_beta,
            percent_burn_in=2/7)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
