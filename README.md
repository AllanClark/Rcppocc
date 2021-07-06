'Rcppocc' is a package that allows the user to build various Bayesian occupancy models. We specifically assume that the regression effects all use logit link functions. Gibbs sampling algorithms are developed to undertake all MCMC sampling.

At present the following models can be fit:
1. A single season Bayesian occupancy (SSO) model using 'PGocc', 'PGocc2' and 'dRUMocc'. It is suggested that you use 'PGocc2' if you are familiar with the 'stocc' package.
2. A SSO Bayesian spatial model (specifically a restricted spatial regression model) using 'occSPATlogit'.
3. A SSO Bayesian spatial occupancy model using a Gaussian approximation to a Polya-Gamma distribuition. This model is termed the Binomial Gaussian model. See 'occSPATlogitBinom'.
4. A SSO Bayesian spatial occupancy (Binomial) model using Polya-Gamma latent variables. See 'occSPATlogitBinomPG'.
5. The implementation of a consensus model with 'k' subsets using a Gaussian approximation to a Polya Gamma distribution. See 'occSPATlogitConsBinom'.
6. The implementation of a consensus model with 'k' subsets using the Binomial Polya-Gamma model. See 'occSPATlogitConsPG'.


The following package can be installed as follows: 
--------------------------------------------------

1. Install the 'devtools' package in R or RStudio. 
   (i.e. install.packages("devtools"))
2. Now type the following R lines of code in your R workspace: 

require(devtools)  
install_github("AllanClark/Rcppocc")

'Rcppocc' should now install. If not, see the additional notes.

Additional notes:
-----------------
1. Ensure that the latest version of the 'stocc' package is also installed. 
2. You might also have to update your version of R so that it is compatible with 'stocc'. 
3. Install the latest version of 'Rtools'. 
4. If you are using old versions of R and using a Mac, have a look at 'https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/' since there might be compilation errors specific to your laptop/computer.
5. If you are a Mac user you might also have to install 'clang'. See 'https://github.com/rmacoslib/r-macos-clang'.

If the above installation does not work email your compilation errors and the specs of the your system to allan.clark@uct.ac.za.

A quick note on exiting 'Rcppocc':
----------------------------------

If you do not want to use the package in your current project, detach the package
as follows:

detach("package:Rcppocc", unload=TRUE)
