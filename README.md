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

require(devtools)  
install_github("AllanClark/Rcppocc")

'Rcppocc' should now install. If not, see the additional notes.

Additional notes:
-----------------
1. Ensure that the latest version of the 'stocc' package is also installed. 
2. You might also have to update your version of R so that it is compatible with 'stocc'. 
3. Install the latest version of 'Rtools'. 
4. If you ae using old versions of R and using a Mac, have a look at 'https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/' since there might be compilation errors specific to your laptop/computer.
5. If you are a Mac user you might also have to install 'clang'. See 'https://github.com/rmacoslib/r-macos-clang'.

If the above installation does not work email your compilation errors and the specs of the your system to allan.clark@uct.ac.za.

A quick note on exiting 'Rcppocc':
----------------------------------

If you do not want to use the package in your current project, detach the package
as follows:

detach("package:Rcppocc", unload=TRUE)
