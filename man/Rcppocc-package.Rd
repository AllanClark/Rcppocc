\name{Rcppocc-package}
\alias{Rcppocc-package}
\alias{Rcppocc}
\docType{package}

\title{Rcppocc}

\description{Fits various Bayesian occupancy models (non-spatial and spatial) assuming all regression effects uses a logit link function.}

\details{The following package allows the user to build various Bayesian
  occupancy models. We specifically assume that the regression
  effects all use logit link functions. Gibbs sampling algorithms are developed
  to undertake all MCMC sampling.

  At present the following models can be fit:
\itemize{
  \item{A single season Bayesian occupancy (SSO) model using \link{PGocc}, \link{PGocc2} and \link{dRUMocc}. It is suggested that you use \link{PGocc2} if you are familiar with the \link{stocc} package.}
  \item{A SSO Bayesian spatial model (specifically a restricted spatial regression model) using \link{occSPATlogit}.}
  \item{A SSO Bayesian spatial occupancy model using a Gaussian approximation to a Polya-Gamma distribuition. This model is termed the Binomial Gaussian model. See \link{occSPATlogitBinom}.}
  \item{A SSO Bayesian spatial occupancy (Binomial) model using Polya-Gamma latent variables. See \link{occSPATlogitBinomPG}.}
  \item{The implementation of a consensus model with \code{k} subsets using a Gaussian approximation to a Polya Gamma distribution. See \link{occSPATlogitConsBinom}.}
  \item{The implementation of a consensus model with \code{k} subsets using the Binomial Polya-Gamma model. See \link{occSPATlogitConsPG}.}
  }
}
\author{A.E. Clark (allan.clark@uct.ac.za)}

\references{
\itemize{
\item{Broms, K.M., Johnson, D.S., Altwegg, R. and Conquest, L.L., (2014). Spatial occupancy models applied to atlas data show Southern Ground Hornbills strongly depend on protected areas. Ecological Applications, 24(2), pp.363-374.}
\item{Clark, A.E. and Altwegg, R., (2019). Efficient Bayesian analysis of occupancy models with logit link functions. Ecology and evolution, 9(2), pp.756-768.}
\item Dorazio, R. M. and Rodriguez, D. T. (2012), A Gibbs sampler for Bayesian analysis of site-occupancy data. Methods in Ecology and Evolution, 3: 1093-1098. \url{doi: 10.1111/j.2041-210X.2012.00237.x}
\item Johnson, D. S., Conn, P. B., Hooten, M., Ray, J., and Pond, B. (2013). Spatial occupancy
models for large data sets. Ecology, \url{http://dx.doi.org/10.1890/12-0564.1.}
\item MacKenzie, D.I., Nichols, J.D., Lachman, G.B., Droege, S., Royle, J.A., Langtimm, C.A. (2002). Estimating site occupancy rates when detection probabilities are less than one. Ecology, 83(8): 2248-2255}
}
