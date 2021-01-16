\name{vb_Designs}
\alias{vb_Designs}

\title{vb_Designs}

\description{
Creates the required 'response' and 'regressor matrices' for \code{PGocc}.
}

\usage{
vb_Designs(W, X, y)
}
\arguments{
  \item{W }{ - A named list of data.frames of covariates that vary within sites. i.e. The dataframes are of dimension \code{n} (number of sites surveyed) by \code{J} where each row is associated with a site and each column represents a site visit. \code{NA} values should be entered at locations that were not surveyed.

e.g. Suppose \code{W} contained three dataframes \code{W1}, \code{W2} and \code{W3};

\code{W$W1[1, ]} = the covariate values for site 1 for all of the visits. Note that some of the entries might be \code{NA} meaning that no visit took place at those occasions.}
  \item{X }{ - A named data.frame that varies at site level.}
  \item{y }{ - An \code{n} by \code{J} matrix of the detection, non-detection data, where \code{n} is the number of sites, \code{J} is the maximum number of sampling periods per site.}
}
\details{
Creates the required 'response' and 'regressor matrices' for \code{PGocc}. Note however that it is suggested that you use \code{PGocc2} instead of \code{PGocc}.
}
\examples{
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

design_mats <- vb_Designs(W=Wmat, X=Xmat, y=Ysim)
}
