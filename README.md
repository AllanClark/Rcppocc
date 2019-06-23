'Rcppocc' is a package that can be used to fit single-season Bayesian occupancy 
models as well as a restricted spatial regression (RSRS) occupancy model. 
models as well as a restricted spatial regression (RSR) occupancy model. 
It is assumed that the regression effects are modelled using logit link 
functions.

The following package can be installed as follows: 
--------------------------------------------------

1. Install the 'devtools' package in R or RStudio. 
   (i.e. install.packages("devtools"))
2. Now type the following R lines of code in your R workspace: 

<<<<<<< HEAD
require(devtools)
install_github("AllanClark/RcppoccDev")
=======
require(devtools)  
install_github("AllanClark/Rcppocc")
>>>>>>> 453554764c2699626c9aa7564242b99a31d0a893

'Rcppocc' should now install.

Ensure that the latest version of the 'stocc' package is also installed. You might also have to update your version of R so that it is compatible with 'stocc'. 


A quick note on exiting 'RcppoccDev':
----------------------------------

If you do not want to use the package in your current project, detach the package
as follows:

detach("package:RcppoccDev", unload=TRUE)

