'Rcppocc' is a package that can be used to fit single-season Bayesian occupancy 
models as well as a restricted spatial regression (RSR) occupancy model. 
It is assumed that the regression effects are modelled using logit link 
functions.

The following package can be installed as follows: 
--------------------------------------------------

1. Install the 'devtools' package in R or RStudio. 
   (i.e. install.packages("devtools"))
2. Now type the following R lines of code in your R workspace: 

require(devtools)  
install_github("AllanClark/Rcppocc")

'Rcppocc' should now install.

Ensure that the latest version of the 'stocc' package is also installed. You might also have to update your version of R so that it is compatible with 'stocc'. 


A quick note on exiting 'Rcppocc':
----------------------------------

If you do not want to use the package in your current project, detach the package
as follows:

detach("package:Rcppocc", unload=TRUE)

