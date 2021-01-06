#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>

#include "RNG.h"
#include "PolyaGamma.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <exception>

#define pi           3.14159265358979323846  /* pi */
using namespace arma;
using namespace Rcpp;
using namespace R;
using namespace sugar;
using namespace std;

#include <RcppParallel.h>
using namespace RcppParallel;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * output is row vector
   */
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);  //1 by ncols matrix
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat mvrnormArma2(int n, arma::vec mu, arma::mat sigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * output is column vector
   */
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(ncols, n);  //ncols by 1 matrix

  //the output is ncols by 1
  return arma::repmat(mu, 1, n) + (arma::chol(sigma)).t()*Y;
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat mvrnormArma3(int n, arma::vec A, arma::mat invSigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * output is column vector
   */
  int ncols = invSigma.n_cols;
  arma::mat Y = arma::randn(ncols, n);  //ncols by 1 matrix
  arma::vec mu = solve(invSigma, A);

  return mu + solve(arma::chol(invSigma),Y);
}

/*
 // [[Rcpp::depends("RcppArmadillo")]]
 // [[Rcpp::export]]
 arma::mat randn_p(int ncols, int n){

 randn(ncols, n)

#ifdef _OPENMP
#pragma omp parallel for
#endif
 for(int j=0;j<p;++j){
 C(j,i)=dot(x.col(j),xi);
 }

 }
 */

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma4(int n, arma::vec A, arma::mat invSigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * output is column vector
   */
  int ncols = invSigma.n_cols;
  //arma::mat xc = arma::randn(ncols, n);
  //Rcpp::Rcout << "ncols = " << ncols << std::endl;
  //Rcpp::Rcout << "nrow arma::randn(ncols, n) = " << xc.n_rows << std::endl;
  //Rcpp::Rcout << "ncol arma::randn(ncols, n) = " << xc.n_cols << std::endl;
  //Rcpp::Rcout << "random = " << xc << std::endl;

  /*
   Rcpp::Rcout << "thread number is  = " << omp_get_num_threads() << std::endl;

   omp_set_num_threads(20);

#pragma omp parallel for
   for (int i = 0; i < 5; i++) {
   //printf("thread is %d\n", omp_get_thread_num());
   Rcpp::Rcout << "thread number is  = " << omp_get_num_threads() << "\n" << std::endl;
   }*/

  //return solve(invSigma, A) + solve(arma::chol(invSigma) , arma::randn(ncols, n));

  return solve(invSigma + arma::chol(invSigma) , A + arma::randn(ncols, n));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma5(int n, arma::mat invSigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * A is a precision matrix
   * output is column vector
   */
  wall_clock timer;

  int ncols = invSigma.n_cols;
  //arma::mat Y = arma::randn(ncols, n);  //ncols by 1 matrix
  //arma::vec mu = solve(invSigma, A);

  //the output is ncols by 1
  //return solve(invSigma, A) + solve(arma::chol(invSigma), arma::randn(ncols, n));

  arma::mat z = arma::randn(ncols, n);
  //Rcpp::Rcout << "'z' = \n" << z << std::endl;
  //Rcpp::Rcout << "'arma::chol(invSigma)' = \n" << arma::chol(invSigma) << std::endl;

  timer.tic();
  //Rcpp::Rcout << "'y (solve)' = \n " << solve(arma::chol(invSigma), z) << std::endl;
  arma::mat s1 = solve(arma::chol(invSigma), z);
  Rcpp::Rcout << "'timer' = " << timer.toc() << std::endl;

  timer.tic();
  //Rcpp::Rcout << "'y (tri)' = \n" << solve(trimatu( arma::chol(invSigma) ), z) << std::endl;
  arma::mat s2 = solve(trimatu( arma::chol(invSigma) ), z);
  Rcpp::Rcout << "'timer' = " << timer.toc() << std::endl;

  return solve(arma::chol(invSigma), z);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma6(int n, arma::vec A, arma::mat invSigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * output is column vector
   */
  //int ncols = invSigma.n_cols;
  //int temp = invSigma.n_cols; //delete later
  return solve( invSigma , arma::mvnrnd(A, invSigma) );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void mvrnormArma4_timers(int n, arma::vec A, arma::mat invSigma) {
  /* Sampling from a Gaussian (multivariate) distribution
   * output is column vector
   */
  int ncols = invSigma.n_cols;
  wall_clock timer;
  double td, td6, td7;

  timer.tic();
  arma::mat mat1 = solve(invSigma, A);
  td = timer.toc();
  Rcpp::Rcout << "'mat1' = " << td << std::endl;

  timer.tic();
  arma::mat mat2 = arma::chol(invSigma);
  td = timer.toc();
  Rcpp::Rcout << "'mat2' = " << td << std::endl;

  timer.tic();
  arma::mat mat3 = arma::randn(ncols, n);
  td = timer.toc();
  Rcpp::Rcout << "'mat3' = " << td << std::endl;

  timer.tic();
  arma::mat mat4 = solve(arma::chol(invSigma) , arma::randn(ncols, n));
  td = timer.toc();
  Rcpp::Rcout << "'mat4' = " << td << std::endl;

  timer.tic();
  arma::mat mat5 = solve(invSigma, A) + solve(arma::chol(invSigma) , arma::randn(ncols, n));
  td = timer.toc();
  Rcpp::Rcout << "'mat5' = " << td << std::endl;

  timer.tic();
  arma::mat mat6 = solve( invSigma , arma::mvnrnd(A, invSigma) );
  td6 = timer.toc();
  Rcpp::Rcout << "'mat6 - new MV sampler' = " << td6 << std::endl;

  timer.tic();
  arma::mat mat7 = mvrnormArma4( 1 , A, invSigma );
  td7 = timer.toc();
  Rcpp::Rcout << "'mat7 - old MV sampler' = " << td7 << std::endl;
  Rcpp::Rcout << "'Efficiency' = " << td7/td6 << std::endl;

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::cube setCubeslices(arma::mat K) {
  /* Set cube slices for K'K
   */

  int n = K.n_rows;
  int p = K.n_cols;
  cube S(p, p, n);

  arma::mat A(p, p );
  A.fill(0);

  for (int i=0; i<n; i++){
    S.slice(i) = K.row(i).t()*K.row(i) ;
    //A += S.slice(i);
    //Rcpp::Rcout << "i=" << i << "; 'S.slice(i)' = \n" << S.slice(i) << std::endl;
  }
  //Rcpp::Rcout << "'A' = \n" << A << std::endl;

  return S;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat multiplybyconstants(arma::mat S, arma::vec Dv) {
  /*double cc = 2.5;
   arma::mat tt = cc*S;*/

  arma::mat tt = Dv(0)*S;
  return tt;
}


// [[Rcpp::depends("RcppArmadillo")]]
double posterior_r3(double x, arma::vec wr, arma::vec sr){

  /*sample the components of the linear gaussian distribution
   * x = z_i_u - log(lambda_beta_i) or y_ij_u - log(lambda_alpha_i)
   * wr and sr are the constants required to approximate the
   * Logistic distribution
   */

  NumericVector w(3);

  for (int comp_j=0; comp_j<3; ++comp_j){
    w(comp_j) = wr(comp_j)*exp( -0.5*std::pow(x/sr(comp_j), 2) )/ (sr(comp_j) * sqrt(2 * pi));
  }
  w = w/sum(w);

  //writing my own sampling function here
  arma::mat u;
  int r;
  u.randu(1);

  if ( u(0)< w(0) ){
    r = 1;
  }
  else if ( u(0)< (w(0)+w(1)) )
  {
    r = 2;
  }
  else
  {
    r = 3;
  }
  return r;
}

NumericVector arma2vec(arma::vec x) {
  //converts from arma::vec to NumericVector
  return Rcpp::NumericVector(x.begin(), x.end());
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::vec invlogit(arma::mat lincomb){
  //inverse logit
  return 1.0/( 1.0 + exp( -lincomb ) );
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::vec ln_invlogit(arma::mat lincomb){
  //log inverse logit
  return -log( 1 + exp( -lincomb ) );
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::vec log_not_prob(arma::mat lincomb){
  //log(1-prob(i))
  return -log( 1 + exp( lincomb ) );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccDA4(arma::mat X, arma::mat Y, arma::mat W_vb, NumericVector siteids,
                 int ndraws,
                 arma::vec ysum, arma::vec z,
                 arma::vec nvisits,
                 arma::mat alpha_m, arma::mat beta_m,
                 arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p){

  /*-----------------------------------------------------------------------
   * The following function uses the dRUM formulation  to fit a single season
   * site occupancy model using Gibbs sampling
   * X = design matrix associated with occupancy process
   * y = observed data
   * W_vb = the design matrix associated with detection process
   * ndraws = number of MCMC iterations
   * Called by the R function, dRUMocc
   *-----------------------------------------------------------------------*/

  //vector of constants related to the Logistic approximation
  arma::vec wr = {0.2522, 0.58523, 0.16257};
  arma::vec s2r = {1.2131, 2.9955, 7.5458};
  //arma::vec sr = {1.101408, 1.730751, 2.746962};

  arma::vec s2r_inv = 1/s2r;
  arma::vec sr = sqrt(s2r);

  arma::mat X_t = X.t(); //the transpose of the design matrix
  int n = X.n_rows; //the number of sites
  int N = W_vb.n_rows; //the total number of visits

  NumericVector siteindex(n); //siteindex = 1, 2, ..., n
  for (int idown =0; idown<n; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0); //identify the indices associated with sites where species was seen
  uvec not_observed_index = find(ysum==0); //identify the indices associated with sites where species was not seen

  /*-----------------------------------------------------------------------*/

  //#1. set initial values for z_i
  //z<-apply(y, 1, max) //this is read in now
  arma::vec psi = z; //set the starting value of psi to z
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's

  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the alpha matrix
  arma::mat sigma_inv_alpha = sigma_inv_alpha_p; //initialize the alpha matrix
  arma::mat sigma_inv_beta = sigma_inv_beta_p; //initialize the alpha matrix
  /*-----------------------------------------------------------------------*/

  //declare some matrices and vectors to use below
  //for 2. below
  arma::mat u;
  arma::mat z_i_u;
  /*-----------------------------------------------------------------------*/

  //for 3. below
  arma::mat log_lambda_beta; arma::mat lambda_beta;
  uvec r_i(n);
  /*-----------------------------------------------------------------------*/

  //for 4. below
  arma::mat s2_beta_inv;
  arma::mat constant_mat1;
  arma::mat beta_cov;
  arma::mat beta_mu;
  /*-----------------------------------------------------------------------*/

  //for 5. below
  uvec z_equals1_rows;//identify all indices with z==1
  uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  //int counter = 0;

  //construct a matrix with the start and end for each site
  arma::vec starts(n);
  starts(0) = 1; //might be better to make it start at 0!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, n-1 ) = ends.rows(0,n-2)+1;
  /*-----------------------------------------------------------------------*/

  //for 6. below
  arma::mat log_lambda_alpha;
  arma::mat lambda_alpha;
  /*-----------------------------------------------------------------------*/

  //for 7. below
  int n_iter;
  arma::mat y_ij_u;
  arma::mat u_Y;
  /*-----------------------------------------------------------------------*/

  //for 7.1 below
  uvec r_ij(N); //not all of these elements are used n_iter<<N
  arma::mat s2_alpha_inv;
  arma::mat constant_mat2;
  arma::mat alpha_cov;
  arma::mat alpha_mu;
  /*-----------------------------------------------------------------------*/

  //for 8. below
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  /*-----------------------------------------------------------------------*/

  //The posterior samples
  //arma::mat post(alpha.n_rows + beta.n_rows, floor(ndraws/2)); //place the samples in this matrix
  arma::mat post_alpha(alpha.n_rows , floor(ndraws/2)); //place the samples in this matrix
  arma::mat post_beta(beta.n_rows, floor(ndraws/2)); //place the samples in this matrix

  //now do the sampling here
  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    log_lambda_beta = X*beta;
    lambda_beta =  exp(log_lambda_beta); //#exp(X%*%beta)

    //2. posterior samples for z_i_u
    z_i_u = zeros<arma::mat>(n, 1);
    u.randu(n);

    /*for (int idown=0; idown<n; idown++){
     z_i_u(idown) = log( lambda_beta(idown)*u(idown) + z(idown) ) - log( 1-u(idown) + lambda_beta(idown)*(1-z(idown)) );
    }*/
    z_i_u = log( lambda_beta%u + z ) - log( 1-u + lambda_beta%(1-z) );

    /*3. posterior sample for r_i
     * See if this can be done more efficiently
     */
    for (int idown=0; idown<n; idown++){
      r_i(idown) = posterior_r3( z_i_u(idown) - log_lambda_beta(idown) , wr, sr) ;
    }

    //4. posterior samples for beta
    s2_beta_inv = diagmat(s2r_inv.rows(r_i-1)); //creates a diagonal matrix!

    constant_mat1 = X_t*s2_beta_inv;
    beta_cov = inv_sympd(sigma_inv_beta + constant_mat1*X); //replaced inv by inv_sympdp
    beta_mu =  beta_cov*constant_mat1*z_i_u;

    //#posterior sample for the beta vector
    beta = mvrnormArma2(1, beta_mu, beta_cov);
    /*-----------------------------------------------------------------------*/

    //5. construct the W matrix associated with z=1

    //Rcpp::Rcout << "z = " << z << std::endl;
    z_equals1_rows = find(z==1);  //row number as specified by c++. if an element say 0 ==> row 0 of z is 1.

    //convert arma::vec z to NumericVector zNM
    NumericVector z_NM = wrap(z);
    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    IntegerVector ind_z = match(siteids, z_pos);
    arma::vec convert_z_vec= as<arma::vec>(ind_z);
    W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    Y_iter = Y.rows( find(convert_z_vec>0) ); //an alternate method of obtaining W_iter and Y_iter

    /*z_equals1_allrows = row indices of W_vb that is associated with z==1
     * these row entries relates to the row indices as per R
     * e.g. is the first entry in z_equals1_allrows = 4 then it means that the
     * c++ row index will actually be 3. NB. (in future this should rather
     * be changed to be consistent with c++)
     * Look at :
     * W_iter = W_vb.rows(z_equals1_allrows.rows(0, counter-1)-1) ;
     * Y_iter = Y.rows(z_equals1_allrows.rows(0, counter-1)-1);
     * That is why I remove 1 from the row indices.
     */

    /*
     counter = 0;
     for (int idown=0; idown< z_equals1_rows.n_rows; idown++){
     for (int jacross=0; jacross< nvisits(z_equals1_rows(idown)); jacross++ ){
     counter = ++counter; //so first entry is 1 (idown=0; jacros=0)
     z_equals1_allrows(counter-1) = starts(z_equals1_rows(idown)) + jacross;
     }
     }

     //Rcpp::Rcout << "nvisits=\n" << nvisits.t() << std::endl;
     //Rcpp::Rcout << "z_equals1_allrows=\n" << z_equals1_allrows.t() << std::endl;

     W_iter = W_vb.rows(z_equals1_allrows.rows(0, counter-1)-1) ;
     Y_iter = Y.rows(z_equals1_allrows.rows(0, counter-1)-1);

     Rcout << find(convert_z_vec>0).t()<< endl;
     Rcout << (z_equals1_allrows.rows(0, counter-1)-1).t()<< endl;*/

    /*-----------------------------------------------------------------------*/

    //6.code up the posterior samples for y_ij(u) and r_ij as well as beta

    log_lambda_alpha = W_iter*alpha; //#log(lambda_alpha)
    lambda_alpha = exp(log_lambda_alpha); //#exp(W.iter%*%alpha)

    //7. posterior samples for z_i_u
    n_iter = Y_iter.n_rows;
    y_ij_u = zeros<arma::mat>(n_iter, 1);
    u_Y.randu(n_iter);
    /*for (int idown=0; idown<n_iter; idown++){
     y_ij_u(idown) = log( lambda_alpha(idown)*u_Y(idown) + Y_iter(idown) ) - log( 1-u_Y(idown) + lambda_alpha(idown)*(1-Y_iter(idown)) );
    }//The loop and the line below produces the same speed!
     */
    //y_ij_u(span(0, n_iter-1),0) = log( lambda_alpha(span(0, n_iter-1),0)%u_Y(span(0, n_iter-1),0) + Y_iter(span(0, n_iter-1),0) ) - log( 1-u_Y(span(0, n_iter-1),0) + lambda_alpha(span(0, n_iter-1),0)%(1-Y_iter(span(0, n_iter-1),0)) );
    y_ij_u.head_rows(n_iter) = log( lambda_alpha.head_rows(n_iter)%u_Y.head_rows(n_iter) + Y_iter.head_rows(n_iter) ) - log( 1-u_Y.head_rows(n_iter) + lambda_alpha.head_rows(n_iter)%(1-Y_iter.head_rows(n_iter)) );
    /*-----------------------------------------------------------------------*/

    //7.1
    /* posterior sample for r_ij
     * See if this can be done more efficiently
     */
    for (int idown=0; idown<n_iter; idown++){
      r_ij(idown) = posterior_r3( y_ij_u(idown) - log_lambda_alpha(idown) , wr, sr) ;
    }

    // posterior samples for alpha
    s2_alpha_inv = diagmat(s2r_inv.rows(r_ij.rows(0, n_iter-1)-1)); //creates a diagonal matrix!

    constant_mat2 = W_iter.t()*s2_alpha_inv;
    alpha_cov = inv_sympd(sigma_inv_alpha + constant_mat2*W_iter);
    alpha_mu =  alpha_cov*constant_mat2*y_ij_u;

    //#posterior sample for the alpha vector
    alpha = mvrnormArma2(1, alpha_mu, alpha_cov);
    /*-----------------------------------------------------------------------*/

    //8. update the p and the psi matrices
    p = 1/(1+exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //Rcpp::Rcout << "p = " << p << std::endl;
    //update the psi matrix for those sites the species has not been seen
    psi.elem(not_observed_index) = 1/(1+exp(-X.rows(not_observed_index)*beta)); //plogis(X[not.observed.index,]%*%beta)
    //Rcpp::Rcout << "psi = " << psi << std::endl;

    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */

    for (int idown=0; idown< not_observed_index.n_rows; idown++){
      /*Rcpp::Rcout << "p.rows(1,3) = \n" << p.rows(1,3) << std::endl;
       Rcpp::Rcout << "prod(p.rows(1,3)) = \n" << prod(p.rows(1,3)) << std::endl;*/

      prod_p_start = starts.row(not_observed_index(idown))-1;
      prod_p_end = ends.row(not_observed_index(idown))-1;

      /*Rcpp::Rcout << "p.rows(...) = \n" << prod(1-p.rows(prod_p_start(0), prod_p_end(0)))  << std::endl;
       Rcpp::Rcout << "Cant do this: prod(p.rows(1,3))(0)" << std::endl;
       1/(1+ (1-psi[x])/(psi[x]*prod(1-p[start.end.index[x,1]:start.end.index[x,2]])) )*/

      prob = 1/( 1+ ( 1-psi.row(not_observed_index(idown)) )/( psi.row(not_observed_index(idown))*prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );
      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) =zdraw(0);
    }//finish sampling the z matrix

    //Rcpp::Rcout << "z(not_observed_index) = \n" << z(not_observed_index).t() << std::endl;//correct

    if (isamples>= floor(ndraws/2)){
      //post.col(isamples-floor(ndraws/2)) = join_cols(alpha, beta);
      post_alpha.col(isamples-floor(ndraws/2)) = alpha;
      post_beta.col(isamples-floor(ndraws/2)) = beta;
    }
  }
  return List::create(_["alpha"]=post_alpha, _["beta"]=post_beta);
}

double rpg4(double scale) {
  //draws 1 PG(1, scale) random variables. here scale is double
  RNG r;
  PolyaGamma pg;

  double result;
  result= pg.draw(1, scale, r);

  return result;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec rpg5(arma::mat scale) {
  /*C++-only interface to PolyaGamma class
   draws random PG variates from arma::mat
   shape = 1
   scale is a arma::mat

   Code adapted from the BayesLogit-master github repository
   YOU NEED THE FOLLOWING FILES IN THE FOLDER: PolyaGamma.h,
   RcppExports.cpp, RNG.cpp, RNG.h, RRNG.cpp, RRNG.h

   */

  int d = scale.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = rpg4(scale(i));
  }

  return result;
}

/*
 Attempt to sample PG variables using RcppParallel
 THIS SHOULD BE RELOOKED!!! NOT IMPLEMENTED AT PRESENT

 struct Rpg : public Worker
 {
 // source matrix
 const RMatrix<double> input;

 // destination matrix
 RMatrix<double> output;

 // initialize with source and destination
 Rpg(const NumericMatrix input, NumericMatrix output)
 : input(input), output(output) {}

 // take the square root of the range of elements requested
 void operator()(std::size_t begin, std::size_t end) {
 std::transform(input.begin() + begin,
 input.begin() + end,
 output.begin() + begin,
 rpg4);
 }
 };

 NumericMatrix parallelMatrixRpg(NumericMatrix x) {

 // allocate the output matrix
 NumericMatrix output(x.nrow(), x.ncol());

 // SquareRoot functor (pass input and output matrixes)
 Rpg rpg_pass(x, output);

 // call parallelFor to do the work
 parallelFor(0, x.length(), rpg_pass);

 // return the output matrix
 return output;
 }*/


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccPG3(arma::mat X, arma::mat Y, arma::mat W_vb, NumericVector siteids,
                 int ndraws,
                 arma::vec ysum, arma::vec z,
                 arma::vec nvisits,
                 arma::mat alpha_m, arma::mat beta_m,
                 arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                 double percent_burn_in){

  /*-----------------------------------------------------------------------
   * The following function uses the POLYA-GAMMA formulation to fit a single season
   * site occupancy model using Gibbs sampling
   * X = design matrix associated with occupancy process
   * y = observed data
   * W_vb = the design matrix associated with detection process
   * ndraws = number of MCMC iterations
   * Called by the R function, PGocc4
   *-----------------------------------------------------------------------*/

  arma::mat X_t = X.t(); //the transpose of the design matrix
  int n = X.n_rows; //the number of sites
  int N = W_vb.n_rows; //the total number of visits

  NumericVector siteindex(n); //siteindex = 1, 2, ..., n
  for (int idown =0; idown<n; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0); //identify the indices associated with sites where species was seen
  uvec not_observed_index = find(ysum==0); //identify the indices associated with sites where species was not seen

  //construct a matrix with the start and end for each site
  arma::vec starts(n);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, n-1 ) = ends.rows(0,n-2)+1;
  //-----------------------------------------------------------------------

  //set initial values for z_i
  arma::vec psi = z; //set the starting value of psi to z
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the beta matrix

  //arma::mat sigma_inv_alpha = sigma_inv_alpha_p; //initialize the alpha precision matrix
  //arma::mat sigma_inv_beta = sigma_inv_beta_p; //initialize the beta precision matrix
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat X_beta; //X*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat diag_pg_beta;
  //arma::mat constant_mat1;
  arma::mat beta_cov;
  arma::mat beta_mu;
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  uvec z_equals1_rows;//identify all indices with z==1
  uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat diag_pg_alpha;
  arma::mat constant_mat2;
  arma::mat alpha_cov;
  arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //arma::mat post(alpha.n_rows + beta.n_rows, floor(ndraws/2)); //place the samples in this matrix
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept); //place the samples in this matrix
  arma::mat post_beta(beta.n_rows, num_samples_kept); //place the samples in this matrix

  //now do the sampling here
  arma::mat tempa;
  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    X_beta = X*beta; //arma::mat
    pg_beta = rpg5(X_beta); //arma::vec
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(X_beta)) ); //bound to be slower!
    //Rcpp::Rcout << "is this possible" << as<arma::vec>( parallelMatrixRpg(wrap(X_beta)) ) << std::endl;

    //posterior samples for beta
    //diag_pg_beta = diagmat( pg_beta ); //creates a diagonal matrix!
    //constant_mat1 = X_t*diag_pg_beta; //arma::mat
    //beta_cov = inv_sympd( sigma_inv_beta_p + X_t*diag_pg_beta*X ); //arma::mat
    beta_cov = inv_sympd( sigma_inv_beta_p + X_t*diagmat( pg_beta )*X ); //arma::mat
    beta_mu =  beta_cov*(X.t()*( z - 0.5 ) + sigma_inv_beta_p*beta_m); // arma::mat

    //#posterior sample for the beta vector
    beta = mvrnormArma2(1, beta_mu, beta_cov); //arma::mat
    //-----------------------------------------------------------------------

    //construct the W matrix associated with z=1
    z_equals1_rows = find(z==1);  //row number as specified by c++. if an element say 0 ==> row 0 of z is 1.

    //convert arma::vec z to NumericVector zNM
    NumericVector z_NM = wrap(z);
    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    IntegerVector ind_z = match(siteids, z_pos);
    arma::vec convert_z_vec= as<arma::vec>(ind_z);
    W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    Y_iter = Y.rows( find(convert_z_vec>0) ); //an alternate method of obtaining W_iter and Y_iter
    //-----------------------------------------------------------------------

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat

    pg_alpha = rpg5(W_alpha); //arma::vec
    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) );

    //posterior samples for alpha
    alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    alpha_mu =  alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + sigma_inv_alpha_p*alpha_m);
    alpha = mvrnormArma2(1, alpha_mu, alpha_cov);
    //-----------------------------------------------------------------------

    //update the p, psi and z
    p = 1/(1+exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    psi.elem(not_observed_index) = 1/(1+exp(-X.rows(not_observed_index)*beta)); //plogis(X[not.observed.index,]%*%beta)

    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    for (int idown=0; idown< not_observed_index.n_rows; idown++){
      prod_p_start = starts.row(not_observed_index(idown))-1;
      prod_p_end = ends.row(not_observed_index(idown))-1;

      /*Rcpp::Rcout << "p.rows(...) = \n" << prod(1-p.rows(prod_p_start(0), prod_p_end(0)))  << std::endl;
       Rcpp::Rcout << "Cant do this: prod(p.rows(1,3))(0)" << std::endl;
       1/(1+ (1-psi[x])/(psi[x]*prod(1-p[start.end.index[x,1]:start.end.index[x,2]])) )*/

      prob = 1/( 1+ ( 1-psi.row(not_observed_index(idown)) )/( psi.row(not_observed_index(idown))*arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );
      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }//finish sampling the z matrix

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;
    if (isamples_counter>= 0){
      //post.col(isamples-floor(ndraws/2)) = join_cols(alpha, beta);
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
    }
    //-----------------------------------------------------------------------
  }//end of sampling

  return List::create(_["alpha"]=post_alpha, _["beta"]=post_beta);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccPG3_z(arma::mat X, arma::mat Y, arma::mat W_vb, NumericVector siteids,
                   int ndraws,
                   arma::vec ysum, arma::vec z,
                   arma::vec nvisits,
                   arma::mat alpha_m, arma::mat beta_m,
                   arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                   double percent_burn_in){

  /*-----------------------------------------------------------------------
   * The following function uses the POLYA-GAMMA formulation to fit a single season
   * site occupancy model using Gibbs sampling
   * X = design matrix associated with occupancy process
   * y = observed data
   * W_vb = the design matrix associated with detection process
   * ndraws = number of MCMC iterations
   * Called by the R function, PGocc4
   *-----------------------------------------------------------------------*/

  arma::mat X_t = X.t(); //the transpose of the design matrix
  int n = X.n_rows; //the number of sites
  int N = W_vb.n_rows; //the total number of visits

  NumericVector siteindex(n); //siteindex = 1, 2, ..., n
  for (int idown =0; idown<n; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0); //identify the indices associated with sites where species was seen
  uvec not_observed_index = find(ysum==0); //identify the indices associated with sites where species was not seen

  //construct a matrix with the start and end for each site
  arma::vec starts(n);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, n-1 ) = ends.rows(0,n-2)+1;
  //-----------------------------------------------------------------------

  //set initial values for z_i
  arma::vec psi = z; //set the starting value of psi to z
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the beta matrix

  //arma::mat sigma_inv_alpha = sigma_inv_alpha_p; //initialize the alpha precision matrix
  //arma::mat sigma_inv_beta = sigma_inv_beta_p; //initialize the beta precision matrix
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat X_beta; //X*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat diag_pg_beta;
  //arma::mat constant_mat1;
  arma::mat beta_cov;
  arma::mat beta_mu;
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  uvec z_equals1_rows;//identify all indices with z==1
  uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat diag_pg_alpha;
  arma::mat constant_mat2;
  arma::mat alpha_cov;
  arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //arma::mat post(alpha.n_rows + beta.n_rows, floor(ndraws/2)); //place the samples in this matrix
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept); //place the samples in this matrix
  arma::mat post_beta(beta.n_rows, num_samples_kept); //place the samples in this matrix
  arma::mat post_z(z.n_rows, num_samples_kept); //place the samples in this matrix
  arma::mat post_psi(z.n_rows, num_samples_kept); //place the samples in this matrix

  //Rcout << "\n Before sampling" << std::endl;

  //now do the sampling here
  arma::mat tempa;
  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    X_beta = X*beta; //arma::mat
    pg_beta = rpg5(X_beta); //arma::vec
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(X_beta)) ); //bound to be slower!
    //Rcpp::Rcout << "is this possible" << as<arma::vec>( parallelMatrixRpg(wrap(X_beta)) ) << std::endl;

    //posterior samples for beta
    //diag_pg_beta = diagmat( pg_beta ); //creates a diagonal matrix!
    //constant_mat1 = X_t*diag_pg_beta; //arma::mat
    //beta_cov = inv_sympd( sigma_inv_beta_p + X_t*diag_pg_beta*X ); //arma::mat
    beta_cov = inv_sympd( sigma_inv_beta_p + X_t*diagmat( pg_beta )*X ); //arma::mat
    beta_mu =  beta_cov*(X.t()*( z - 0.5 ) + sigma_inv_beta_p*beta_m); // arma::mat

    //#posterior sample for the beta vector
    beta = mvrnormArma2(1, beta_mu, beta_cov); //arma::mat
    //Rcout << "\n beta" << std::endl;
    //-----------------------------------------------------------------------

    //construct the W matrix associated with z=1
    z_equals1_rows = find(z==1);  //row number as specified by c++. if an element say 0 ==> row 0 of z is 1.

    //convert arma::vec z to NumericVector zNM
    NumericVector z_NM = wrap(z);
    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    IntegerVector ind_z = match(siteids, z_pos);
    arma::vec convert_z_vec= as<arma::vec>(ind_z);
    W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    Y_iter = Y.rows( find(convert_z_vec>0) ); //an alternate method of obtaining W_iter and Y_iter
    //-----------------------------------------------------------------------
    //Rcout << "\n Before alpha" << std::endl;

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat

    pg_alpha = rpg5(W_alpha); //arma::vec
    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) );
    //Rcout << "\n pg_alpha" << std::endl;
    //posterior samples for alpha
    //Rcout << "\n sigma_inv_alpha_p" << size(sigma_inv_alpha_p) << std::endl;
    //Rcout << "\n W_iter" << size(W_iter) << std::endl;
    //Rcout << "\n W_iter.t()*diagmat( pg_alpha )*W_iter" << size(W_iter.t()*diagmat( pg_alpha )*W_iter) << std::endl;
    alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    //Rcout << "\n alpha_cov" << std::endl;
    alpha_mu =  alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + sigma_inv_alpha_p*alpha_m);
    //Rcout << "\n alpha_mu" << std::endl;
    alpha = mvrnormArma2(1, alpha_mu, alpha_cov);
    //  Rcout << "\n alpha" << std::endl;
    //-----------------------------------------------------------------------

    //update the p, psi and z
    p = 1/(1+exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    //psi.elem(not_observed_index) = 1/(1+exp(-X.rows(not_observed_index)*beta)); //plogis(X[not.observed.index,]%*%beta)
    //changed here since we calculate the residuals as defoned by Wright et al 2019
    psi = 1/(1+exp(-X*beta));

    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    for (int idown=0; idown< not_observed_index.n_rows; idown++){
      prod_p_start = starts.row(not_observed_index(idown))-1;
      prod_p_end = ends.row(not_observed_index(idown))-1;

      /*Rcpp::Rcout << "p.rows(...) = \n" << prod(1-p.rows(prod_p_start(0), prod_p_end(0)))  << std::endl;
       Rcpp::Rcout << "Cant do this: prod(p.rows(1,3))(0)" << std::endl;
       1/(1+ (1-psi[x])/(psi[x]*prod(1-p[start.end.index[x,1]:start.end.index[x,2]])) )*/

      prob = 1/( 1+ ( 1-psi.row(not_observed_index(idown)) )/( psi.row(not_observed_index(idown))*arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );
      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }//finish sampling the z matrix

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;
    if (isamples_counter>= 0){
      //post.col(isamples-floor(ndraws/2)) = join_cols(alpha, beta);
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
      post_z.col(isamples_counter) = z;
      post_psi.col(isamples_counter) = psi;
    }
    //-----------------------------------------------------------------------
  }//end of sampling

  return List::create(_["alpha"]=post_alpha,
                      _["beta"]=post_beta,
                      _["z"]=post_z,
                      _["psi"]=post_psi);
}

// [[Rcpp::depends("RcppArmadillo")]]
double rgammadouble(int a, double b, double c)
{   //from http://iamrandom.com/rgamma-rgammadouble
  Rcpp::NumericVector x = rgamma(a,b,1.0/c);
  return x(0);
}

// [[Rcpp::depends("RcppArmadillo")]]
double rnormdouble(double b, double c)
{   //from http://iamrandom.com/rnorm-rnormdouble
  //b = mean; c = sd
  Rcpp::NumericVector x = rnorm(1,b,c);
  return x(0);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat quadform(arma::mat X, arma::mat Xt, arma::vec dmat){
  //Calculate the quadratic form X.t()*S*X
  //where S is diagonal with positive elements

  for (int i =0; i < X.n_rows; i++){
    //X.row(i) = X.row(i)*dmat(i);
    X.row(i)*=sqrt(dmat(i));
    Xt.col(i)*=sqrt(dmat(i));
  }
  return Xt*X;
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat diagtimesX(arma::mat X, arma::vec dmat){
  //Calculate S*X
  //where S is diagonal with positive elements
  //X is a matrix

  for (int i =0; i < X.n_rows; i++){
    //X.row(i) = X.row(i)*dmat(i);
    X.row(i)*=dmat(i);
  }
  return X;
}

// [[Rcpp::export]]
arma::vec seq_int(long int a, long int b){
  /*
   *From https://github.com/coatless/r-to-armadillo/blob/master/src/seq.cpp
   *create a sequence of numbers from a to b
   */

  long int d = std::abs(b-a)+1;

  return arma::linspace(a, b, d);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat quadform2(arma::mat X, arma::vec dmat) {
  //Calculation of a quadratic form X'*D*X where D is diagonal

  X.each_col() %= sqrt(dmat);
  return X.t()*X;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
int rbinom2(arma::vec prob) {
  return rbinom(1,1, prob(0))(0);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
int rbinom3(double prob) {

  double u = randu();

  if (u <= prob){
    return 1;
  }else{
    return 0;
  }
}






// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat matrix_multiplication3(arma::mat x, arma::mat y){
  return x*y;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat matrix_multiplication4(arma::mat x){
  return x.t()*x;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat transpose_sq(arma::mat x){
  /* Taken from the Rfast package
   * transpose a square matrix
   */

  const int p=x.n_cols,n=x.n_rows;

  for(int i=1;i<p;++i){
    for(int u=0;u<i;++u){
      swap(x(u,i),x(i,u));
    }
  }
  return x;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat Transpose2(arma::mat x){
  //.t() is faster than transpose_sq and transpose_notsq
  return x.t();
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat diagtimesX2(arma::mat X, arma::vec dmat){
  //Calculate S*X
  //where S is diagonal with positive elements
  //X is a matrix (but here its a column vector)

  return X % dmat;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccSPAT(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits,
                  arma::mat K, arma::mat Minv, //arma::mat Mt,
                  double n_obs,
                  NumericVector siteids, arma::vec unsurveyed_ind,
                  double tau_0, double a_tau, double b_tau,
                  arma::mat alpha_m, arma::mat beta_m,
                  arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                  int ndraws, double percent_burn_in){

  /*
   * Called by the R function, occSPATlogit
   *
   * Some Rcpp notes
   * Some R translations put in brackets at times.
   *
   * 1. Indices
   * arma::mat x;
   * x.row(0) is the first row a matrix. (x[1,])
   * x.row(1) is the second row a matrix. (x[2,])
   * Indices start at 0!
   *
   * 2. Multiple rows
   * x.rows(1,3) is x[2:4, ] using the R notation.
   */

  //some constants or constnat matrices
  int const1 = n_obs -1;
  arma::mat const2 = sigma_inv_beta_p*beta_m;
  arma::mat const3 = sigma_inv_alpha_p*alpha_m;

  // define some matrices
  arma::mat Xs = X.rows(0, const1);
  arma::mat Xs_t = Xs.t(); //the transpose of the design matrix
  //int n = X.n_rows; //the number of sites (survey + unsurvey)
  int N = W_vb.n_rows; //the total number of visits

  arma::mat Ks = K.rows(0, const1);
  arma::mat Ks_t = Ks.t(); //the transpose of the spatial design matrix
  int r = K.n_cols; //the number of column of K; the number of spatial random effects added

  NumericVector siteindex(n_obs); //siteindex = 1, 2, ..., number of surveyed sites
  for (int idown =0; idown<n_obs; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0);
  //identify the indices associated with sites where species was seen (from surveyd sites)

  uvec not_observed_index = find(ysum==0);
  //identify the indices associated with sites where species was not seen (from surveyd sites)

  int n_not_observed_surveyed_sites = size(find(ysum.rows(0, const1)==0))(0);
  //the number of sites with 0 detections based on surveyed sites

  uvec unsurveyed_index = find(unsurveyed_ind>0) + n_obs;
  int n_unsurveyed = size(unsurveyed_index)(0);

  //construct a matrix with the start and end for each site
  arma::vec starts(n_obs);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, const1 ) = ends.rows(0, const1-1)+1;
  //-----------------------------------------------------------------------

  //Declare parameter vectors
  arma::vec psi = z; //set the starting value of psi to z
  //arma::vec psi_occ;
  //arma::vec common_term;
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the alpha matrix
  mat theta = zeros<mat>(r,1); //initialise the theta matrix to 0's; i.e. no spatial effects initiall
  double tau=tau_0; //initialize tau =1; the spatial precision scalar
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat Xs_beta; //Xs*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat beta_cov;
  //arma::mat beta_mu;

  //Posterior sampling of theta
  arma::mat Ks_theta; //Ks*theta
  arma::mat theta_cov;
  //arma::mat theta_mu;
  //-----------------------------------------------------------------------

  //Posterior sampling of tau
  mat inv_scale;
  mat tau_mat(1,1);
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  //uvec z_equals1_rows;//identify all indices with z==1
  //uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  //int counter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat alpha_cov;
  //arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //The outputs
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept);
  arma::mat post_beta(beta.n_rows, num_samples_kept);
  arma::mat post_theta(theta.n_rows, num_samples_kept);
  arma::mat post_tau(1, num_samples_kept);
  arma::mat post_z(z.n_rows, num_samples_kept);
  arma::mat post_psi(z.n_rows, num_samples_kept);
  //-----------------------------------------------------------------------

  //now do the sampling here
  Xs_beta = Xs*beta; //arma::mat
  Ks_theta = Ks*theta; //arma::mat

  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(Xs_beta + Ks_theta)) );
    pg_beta = rpg5(Xs_beta + Ks_theta);
    //-----------------------------------------------------------------------

    //posterior samples for beta
    beta_cov = inv_sympd( sigma_inv_beta_p + Xs_t*diagmat( pg_beta )*Xs ); //arma::mat
    beta = mvnrnd(beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov, 1);
    //beta = mvrnormArma6(1, beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov);
    Xs_beta = Xs*beta; //arma::mat
    //-----------------------------------------------------------------------

    //posterior samples for theta
    theta_cov = inv_sympd( tau*Minv + Ks_t*diagmat( pg_beta )*Ks );
    theta = mvnrnd(theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov, 1);
    //theta = mvrnormArma6(1, theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov);
    Ks_theta = Ks*theta; //arma::mat
    //-----------------------------------------------------------------------

    //posterior samples for tau
    /*inv_scale =  theta.t()*Minv*theta *0.5 + b_tau ;
     tau = randg(1, distr_param(a_tau + 0.5*r, 1.0/inv_scale(0)))(0);
     tau_mat(0,0) = tau;*/

    inv_scale =  theta.t()*Minv*theta*0.5 + b_tau ; //not inverse scale. rgammadouble uses 1/inv_scale = 1/( theta.t()*Minv*theta *0.5 + i2 );
    tau = rgammadouble(1, a_tau + 0.5*r, inv_scale(0) );
    tau_mat(0,0) = tau;
    //-----------------------------------------------------------------------

    //convert arma::vec z to NumericVector zNM (but only for surveyed locations)
    NumericVector z_NM = wrap(z.rows(0, const1));
    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    IntegerVector ind_z = match(siteids, z_pos);
    arma::vec convert_z_vec= as<arma::vec>(ind_z);
    W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    Y_iter = Y.rows( find(convert_z_vec>0) );
    //-----------------------------------------------------------------------

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat
    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) ); //arma::vec
    pg_alpha = rpg5(W_alpha); //arma::vec

    //posterior samples for alpha
    alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    alpha = mvnrnd(alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov, 1);
    //alpha = mvrnormArma6(1, alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov);
    //-----------------------------------------------------------------------

    //update the p, psi and z
    p = 1/(1 + exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    //psi.elem(not_observed_index) = 1/(1 + exp(-X.rows(not_observed_index)*beta - K.rows(not_observed_index)*theta ));
    psi = 1/(1 + exp(-X*beta - K*theta ));

    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    //for surveyed sites with zero detections
    //this can be made more efficient!
    for (int idown=0; idown< n_not_observed_surveyed_sites; idown++){
      prod_p_start = starts.row(not_observed_index(idown))-1;
      prod_p_end = ends.row(not_observed_index(idown))-1;

      prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/( psi.row(not_observed_index(idown))*arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );

      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }

    //for unsurveyed sites
    for (int idown=0; idown< n_unsurveyed; idown++){
      prob = psi.row(n_obs + idown);
      zdraw = rbinom(1,1, prob(0));
      z(n_obs + idown) = zdraw(0);
    } //finish sampling the z matrix

    //idea from stocc package
    /*
     common_term = exp(ln_invlogit( X*beta + K*theta ) + Mt*log_not_prob(W_vb*alpha));
     psi_occ = common_term/( common_term + exp( log_not_prob(X*beta) ) );

     //sample z
     for (int idown=0; idown< n; idown++){
     prob = psi_occ.row(idown);
     zdraw = rbinom(1,1, prob(0));
     z(idown) = zdraw(0);
     }//finish sampling the z matrix
     */

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;

    if (isamples_counter>=0){
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
      post_theta.col(isamples_counter) = theta;
      post_tau.col(isamples_counter) = tau_mat;
      post_z.col(isamples_counter) = z;
      post_psi.col(isamples_counter) = psi;
    }
  }

  return List::create(_["alpha"]=post_alpha,
                      _["beta"]=post_beta,
                      _["theta"]=post_theta,
                      _["tau"]=post_tau,
                      _["real.occ"]= post_z,
                      _["psi"]=post_psi );

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccSPAT3(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits,
                   arma::mat K, arma::mat Minv, //arma::mat Mt,
                   double n_obs,
                   NumericVector siteids, arma::vec unsurveyed_ind,
                   double tau_0, double a_tau, double b_tau,
                   arma::mat alpha_m, arma::mat beta_m,
                   arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                   int ndraws, double percent_burn_in){

  /*
   * Updated 20 April 2020
   * Called by the R function, occSPATlogit
   *
   * Some Rcpp notes
   * Some R translations put in brackets at times.
   *
   * 1. Indices
   * arma::mat x;
   * x.row(0) is the first row a matrix. (x[1,])
   * x.row(1) is the second row a matrix. (x[2,])
   * Indices start at 0!
   *
   * 2. Multiple rows
   * x.rows(1,3) is x[2:4, ] using the R notation.
   */

  //some constants or constnat matrices
  int const1 = n_obs -1;
  arma::mat const2 = sigma_inv_beta_p*beta_m;
  arma::mat const3 = sigma_inv_alpha_p*alpha_m;

  // define some matrices
  arma::mat Xs = X.rows(0, const1);
  arma::mat Xs_t = Xs.t(); //the transpose of the design matrix
  //int n = X.n_rows; //the number of sites (survey + unsurvey)
  int N = W_vb.n_rows; //the total number of visits

  arma::mat Ks = K.rows(0, const1);
  arma::mat Ks_t = Ks.t(); //the transpose of the spatial design matrix
  int r = K.n_cols; //the number of column of K; the number of spatial random effects added

  NumericVector siteindex(n_obs); //siteindex = 1, 2, ..., number of surveyed sites
  for (int idown =0; idown<n_obs; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0);
  //identify the indices associated with sites where species was seen (from surveyd sites)

  uvec not_observed_index = find(ysum==0);
  //identify the indices associated with sites where species was not seen (from surveyd sites)

  int n_not_observed_surveyed_sites = size(find(ysum.rows(0, const1)==0))(0);
  //the number of sites with 0 detections based on surveyed sites

  uvec unsurveyed_index = find(unsurveyed_ind>0) + n_obs;
  int n_unsurveyed = size(unsurveyed_index)(0);

  //construct a matrix with the start and end for each site
  arma::vec starts(n_obs);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, const1 ) = ends.rows(0, const1-1)+1;
  //-----------------------------------------------------------------------

  //Declare parameter vectors
  arma::vec psi = z; //set the starting value of psi to z
  //arma::vec psi_occ;
  //arma::vec common_term;
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the alpha matrix
  mat theta = zeros<mat>(r,1); //initialise the theta matrix to 0's; i.e. no spatial effects initiall
  double tau=tau_0; //initialize tau =1; the spatial precision scalar
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat Xs_beta; //Xs*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat beta_cov;
  //arma::mat beta_mu;

  //Posterior sampling of theta
  arma::mat Ks_theta; //Ks*theta
  arma::mat theta_cov;
  //arma::mat theta_mu;
  //-----------------------------------------------------------------------

  //Posterior sampling of tau
  mat inv_scale;
  mat tau_mat(1,1);
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  //uvec z_equals1_rows;//identify all indices with z==1
  //uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  //int counter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat alpha_cov;
  //arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //The outputs
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept);
  arma::mat post_beta(beta.n_rows, num_samples_kept);
  arma::mat post_theta(theta.n_rows, num_samples_kept);
  arma::mat post_tau(1, num_samples_kept);
  arma::mat post_z(z.n_rows, num_samples_kept);
  arma::mat post_psi(z.n_rows, num_samples_kept);
  //-----------------------------------------------------------------------

  //now do the sampling here
  Xs_beta = Xs*beta; //arma::mat
  Ks_theta = Ks*theta; //arma::mat

  wall_clock timer;

  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(Xs_beta + Ks_theta)) );
    pg_beta = rpg5(Xs_beta + Ks_theta);
    //-----------------------------------------------------------------------

    //posterior samples for beta
    //beta_cov = inv_sympd( sigma_inv_beta_p + Xs_t*diagmat( pg_beta )*Xs ); //arma::mat
    //beta = mvnrnd(beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov, 1);
    //timer.tic();
    beta_cov = inv_sympd( sigma_inv_beta_p + quadform2(Xs, pg_beta) ); //arma::mat
    beta = mvnrnd(beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - Ks_theta % pg_beta ) + const2), beta_cov, 1);
    //Rcpp::Rcout << "Beta current timer" << timer.toc() << std::endl;

    /*timer.tic();
     beta = mvrnormArma4(1, Xs_t*( z.rows(0,n_obs-1) - 0.5 - Ks_theta % pg_beta ) + const2 , sigma_inv_beta_p + Xs_t*diagmat( pg_beta )*Xs);
     Rcpp::Rcout << "Beta new timer" << timer.toc() << std::endl;
     Rcpp::Rcout << " \n" << std::endl;
     */

    //beta = mvrnormArma2(1, beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov);
    Xs_beta = Xs*beta; //arma::mat
    //-----------------------------------------------------------------------

    //posterior samples for theta
    //theta_cov = inv_sympd( tau*Minv + Ks_t*diagmat( pg_beta )*Ks );
    //theta = mvnrnd(theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov, 1);
    //timer.tic();
    theta_cov = inv_sympd( tau*Minv + quadform2(Ks, pg_beta) );
    theta = mvnrnd(theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - Xs_beta % pg_beta  ), theta_cov, 1);
    //Rcpp::Rcout << "Theta current timer = " << timer.toc() << std::endl;

    /*timer.tic();
     theta = mvrnormArma4(1, Ks_t*( z.rows(0, const1) - 0.5 - Xs_beta % pg_beta  ) , tau*Minv + quadform2(Ks, pg_beta) );
     Rcpp::Rcout << "Theta new timer = " << timer.toc() << std::endl;
     Rcpp::Rcout << " \n" << std::endl;
     */

    //theta = mvrnormArma2(1, theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov);
    Ks_theta = Ks*theta; //arma::mat
    //-----------------------------------------------------------------------

    //posterior samples for tau
    /*inv_scale =  theta.t()*Minv*theta *0.5 + b_tau ;
     tau = randg(1, distr_param(a_tau + 0.5*r, 1.0/inv_scale(0)))(0);
     tau_mat(0,0) = tau;*/

    inv_scale =  theta.t()*Minv*theta*0.5 + b_tau ; //not inverse scale. rgammadouble uses 1/inv_scale = 1/( theta.t()*Minv*theta *0.5 + i2 );
    tau = rgammadouble(1, a_tau + 0.5*r, inv_scale(0) );
    tau_mat(0,0) = tau;
    //-----------------------------------------------------------------------

    //convert arma::vec z to NumericVector zNM (but only for surveyed locations)
    NumericVector z_NM = wrap(z.rows(0, const1));
    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    IntegerVector ind_z = match(siteids, z_pos);
    arma::vec convert_z_vec= as<arma::vec>(ind_z);
    W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    Y_iter = Y.rows( find(convert_z_vec>0) );
    //-----------------------------------------------------------------------

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat
    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) ); //arma::vec
    pg_alpha = rpg5(W_alpha); //arma::vec

    //posterior samples for alpha
    //alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    //timer.tic();
    alpha_cov = inv_sympd( sigma_inv_alpha_p + quadform2(W_iter, pg_alpha) );
    alpha = mvnrnd(alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov, 1);
    //Rcpp::Rcout << " Alpha current timer = " << timer.toc() << std::endl;

    /*timer.tic();
     alpha = mvrnormArma4(1, W_iter.t()*( Y_iter - 0.5 ) + const3 , sigma_inv_alpha_p + quadform2(W_iter, pg_alpha) );
     Rcpp::Rcout << " Alpha new timer = " << timer.toc() << std::endl;
     Rcpp::Rcout << " \n" << std::endl;
     */

    //alpha = mvrnormArma2(1, alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov);
    //-----------------------------------------------------------------------

    //update the p, psi and z
    p = 1/(1 + exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    //psi.elem(not_observed_index) = 1/(1 + exp(-X.rows(not_observed_index)*beta - K.rows(not_observed_index)*theta ));
    psi = 1/(1 + exp(-X*beta - K*theta ));

    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    //for surveyed sites with zero detections
    //this can be made more efficient!
    for (int idown=0; idown< n_not_observed_surveyed_sites; idown++){
      prod_p_start = starts.row(not_observed_index(idown))-1;
      prod_p_end = ends.row(not_observed_index(idown))-1;

      prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/( psi.row(not_observed_index(idown))*arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );

      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }

    //for unsurveyed sites
    for (int idown=0; idown< n_unsurveyed; idown++){
      prob = psi.row(n_obs + idown);
      zdraw = rbinom(1,1, prob(0));
      z(n_obs + idown) = zdraw(0);
    } //finish sampling the z matrix

    //idea from stocc package
    /*
     common_term = exp(ln_invlogit( X*beta + K*theta ) + Mt*log_not_prob(W_vb*alpha));
     psi_occ = common_term/( common_term + exp( log_not_prob(X*beta) ) );

     //sample z
     for (int idown=0; idown< n; idown++){
     prob = psi_occ.row(idown);
     zdraw = rbinom(1,1, prob(0));
     z(idown) = zdraw(0);
     }//finish sampling the z matrix
     */

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;

    if (isamples_counter>=0){
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
      post_theta.col(isamples_counter) = theta;
      post_tau.col(isamples_counter) = tau_mat;
      post_z.col(isamples_counter) = z;
      post_psi.col(isamples_counter) = psi;
    }
  }

  return List::create(_["alpha"]=post_alpha,
                      _["beta"]=post_beta,
                      _["theta"]=post_theta,
                      _["tau"]=post_tau,
                      _["real.occ"]= post_z,
                      _["psi"]=post_psi );

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccSPATsplit(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits,
                       arma::mat K, arma::mat Minv, //arma::mat Mt,
                       double n_obs,
                       NumericVector siteids, arma::vec unsurveyed_ind,
                       double tau_0, double a_tau, double b_tau,
                       arma::mat alpha_m, arma::mat beta_m,
                       arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                       int ndraws, double percent_burn_in){

  /*
   * Called by the R function, occSPATlogit
   *
   * Some Rcpp notes
   * Some R translations put in brackets at times.
   *
   * 1. Indices
   * arma::mat x;
   * x.row(0) is the first row a matrix. (x[1,])
   * x.row(1) is the second row a matrix. (x[2,])
   * Indices start at 0!
   *
   * 2. Multiple rows
   * x.rows(1,3) is x[2:4, ] using the R notation.
   */

  //some constants or constnat matrices
  int const1 = n_obs -1;
  arma::mat const2 = sigma_inv_beta_p*beta_m;
  arma::mat const3 = sigma_inv_alpha_p*alpha_m;

  // define some matrices
  arma::mat Xs = X.rows(0, const1);
  arma::mat Xs_t = Xs.t(); //the transpose of the design matrix
  //int n = X.n_rows; //the number of sites (survey + unsurvey)
  int N = W_vb.n_rows; //the total number of visits

  arma::mat Ks = K.rows(0, const1);
  arma::mat Ks_t = Ks.t(); //the transpose of the spatial design matrix
  int r = K.n_cols; //the number of column of K; the number of spatial random effects added

  NumericVector siteindex(n_obs); //siteindex = 1, 2, ..., number of surveyed sites
  for (int idown =0; idown<n_obs; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0);
  //identify the indices associated with sites where species was seen (from surveyd sites)

  uvec not_observed_index = find(ysum==0);
  //identify the indices associated with sites where species was not seen (from surveyd sites)

  int n_not_observed_surveyed_sites = size(find(ysum.rows(0, const1)==0))(0);
  //the number of sites with 0 detections based on surveyed sites

  uvec unsurveyed_index = find(unsurveyed_ind>0) + n_obs;
  int n_unsurveyed = size(unsurveyed_index)(0);

  //construct a matrix with the start and end for each site
  arma::vec starts(n_obs);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, const1 ) = ends.rows(0, const1-1)+1;
  //-----------------------------------------------------------------------

  //Declare parameter vectors
  arma::vec psi = z; //set the starting value of psi to z
  //arma::vec psi_occ;
  //arma::vec common_term;
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the alpha matrix
  mat theta = zeros<mat>(r,1); //initialise the theta matrix to 0's; i.e. no spatial effects initiall
  double tau=tau_0; //initialize tau =1; the spatial precision scalar
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat Xs_beta; //Xs*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat beta_cov;
  //arma::mat beta_mu;

  //Posterior sampling of theta
  arma::mat Ks_theta; //Ks*theta
  arma::mat theta_cov;
  //arma::mat theta_mu;
  //-----------------------------------------------------------------------

  //Posterior sampling of tau
  mat inv_scale;
  mat tau_mat(1,1);
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  //uvec z_equals1_rows;//identify all indices with z==1
  //uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  //int counter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat alpha_cov;
  //arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //The outputs
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept);
  arma::mat post_beta(beta.n_rows, num_samples_kept);
  arma::mat post_theta(theta.n_rows, num_samples_kept);
  arma::mat post_tau(1, num_samples_kept);
  arma::mat post_z(z.n_rows, num_samples_kept);
  arma::mat post_psi(z.n_rows, num_samples_kept);
  //-----------------------------------------------------------------------

  //now do the sampling here
  Xs_beta = Xs*beta; //arma::mat
  Ks_theta = Ks*theta; //arma::mat

  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(Xs_beta + Ks_theta)) );
    pg_beta = rpg5(Xs_beta + Ks_theta);
    //-----------------------------------------------------------------------

    //posterior samples for beta
    beta_cov = inv_sympd( sigma_inv_beta_p + Xs_t*diagmat( pg_beta )*Xs ); //arma::mat
    beta = mvnrnd(beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov, 1);
    //beta = mvrnormArma2(1, beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov);
    Xs_beta = Xs*beta; //arma::mat
    //-----------------------------------------------------------------------

    //posterior samples for theta
    theta_cov = inv_sympd( tau*Minv + Ks_t*diagmat( pg_beta )*Ks );
    theta = mvnrnd(theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov, 1);
    //theta = mvrnormArma2(1, theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov);
    Ks_theta = Ks*theta; //arma::mat
    //-----------------------------------------------------------------------

    //posterior samples for tau
    /*inv_scale =  theta.t()*Minv*theta *0.5 + b_tau ;
     tau = randg(1, distr_param(a_tau + 0.5*r, 1.0/inv_scale(0)))(0);
     tau_mat(0,0) = tau;*/

    inv_scale =  theta.t()*Minv*theta*0.5 + b_tau ; //not inverse scale. rgammadouble uses 1/inv_scale = 1/( theta.t()*Minv*theta *0.5 + i2 );
    tau = rgammadouble(1, a_tau + 0.5*r, inv_scale(0) );
    tau_mat(0,0) = tau;
    //-----------------------------------------------------------------------

    //convert arma::vec z to NumericVector zNM (but only for surveyed locations)
    NumericVector z_NM = wrap(z.rows(0, const1));
    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    IntegerVector ind_z = match(siteids, z_pos);
    arma::vec convert_z_vec= as<arma::vec>(ind_z);
    W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    Y_iter = Y.rows( find(convert_z_vec>0) );
    //-----------------------------------------------------------------------

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat
    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) ); //arma::vec
    pg_alpha = rpg5(W_alpha); //arma::vec

    //posterior samples for alpha
    alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    alpha = mvnrnd(alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov, 1);
    //alpha = mvrnormArma2(1, alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov);
    //-----------------------------------------------------------------------

    //update the p, psi and z
    p = 1/(1 + exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    //psi.elem(not_observed_index) = 1/(1 + exp(-X.rows(not_observed_index)*beta - K.rows(not_observed_index)*theta ));
    psi = 1/(1 + exp(-X*beta - K*theta ));

    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    //for surveyed sites with zero detections
    //this can be made more efficient!
    for (int idown=0; idown< n_not_observed_surveyed_sites; idown++){
      prod_p_start = starts.row(not_observed_index(idown))-1;
      prod_p_end = ends.row(not_observed_index(idown))-1;

      prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/( psi.row(not_observed_index(idown))*arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );

      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }

    //for unsurveyed sites
    for (int idown=0; idown< n_unsurveyed; idown++){
      prob = psi.row(n_obs + idown);
      zdraw = rbinom(1,1, prob(0));
      z(n_obs + idown) = zdraw(0);
    } //finish sampling the z matrix

    //idea from stocc package
    /*
     common_term = exp(ln_invlogit( X*beta + K*theta ) + Mt*log_not_prob(W_vb*alpha));
     psi_occ = common_term/( common_term + exp( log_not_prob(X*beta) ) );

     //sample z
     for (int idown=0; idown< n; idown++){
     prob = psi_occ.row(idown);
     zdraw = rbinom(1,1, prob(0));
     z(idown) = zdraw(0);
     }//finish sampling the z matrix
     */

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;

    if (isamples_counter>=0){
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
      post_theta.col(isamples_counter) = theta;
      post_tau.col(isamples_counter) = tau_mat;
      post_z.col(isamples_counter) = z;
      post_psi.col(isamples_counter) = psi;
    }
  }

  return List::create(_["alpha"]=post_alpha,
                      _["beta"]=post_beta,
                      _["theta"]=post_theta,
                      _["tau"]=post_tau,
                      _["real.occ"]= post_z,
                      _["psi"]=post_psi );

}





/*
 * Binomial Model Implementations
 *
 */

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double rpg6(int n, double scale) {
  //draws 1 PG(n, scale) random variables. here scale is double
  RNG r;
  PolyaGamma pg;

  double result;
  result= pg.draw(n, scale, r);

  return result;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec rpg7(arma::vec n, arma::mat scale) {
  /*C++-only interface to PolyaGamma class
   draws random PG variates from arma::mat
   shape = n
   scale is a arma::mat

   Code adapted from the BayesLogit-master github repository
   YOU NEED THE FOLLOWING FILES IN THE FOLDER: PolyaGamma.h,
   RcppExports.cpp, RNG.cpp, RNG.h, RRNG.cpp, RRNG.h

   */

  int d = scale.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = rpg6(n(i), scale(i));
  }

  return result;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double rnorm6(double a, double c) {
  double mu;
  double x2;
  double var;
  if (abs(c) < 0.00001) {
    mu = a/4;
    var = a/24;
  }

  else {
    mu = (a/(2*c))*(std::tanh(c/2));
    double nom =  a*(-(2+a)*(pow(c,2))+ a*pow(c,2)*std::cosh(c)+2*c*std::sinh(c));
    double denom = 8*pow(c,4)*pow(std::cosh(c/2),2);
    x2 = nom/denom;
    var= x2-(pow(mu,2));
  }

  double sd = sqrt(var);

  double result = arma::randn<double>()*sd + mu;

  if (result <= 0.001){
    result = rpg6(a,c);
  }

  return result;


}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double rnorm_1(double mu, double sd){

  double result = arma::randn<double>()*sd + mu;
  return(result);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec rnorm7(arma::vec n, arma::mat scale){
  int d = scale.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = rnorm6(n(i),scale(i));
  }

  return result;

}





/*
 * New spatial occupancy model using the binomial model
 * Approximating Polya-Gamma by a Normal distribution
 */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccSPATBINOM(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits,
                       arma::mat K, arma::mat Minv, //arma::mat Mt,
                       double n_obs,
                       NumericVector siteids, arma::vec unsurveyed_ind,
                       double tau_0, double a_tau, double b_tau,
                       arma::mat alpha_m, arma::mat beta_m,
                       arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                       int ndraws, double percent_burn_in){

  /*
   * Called by the R function, occSPATlogit
   *
   * Some Rcpp notes
   * Some R translations put in brackets at times.
   *
   * 1. Indices
   * arma::mat x;
   * x.row(0) is the first row a matrix. (x[1,])
   * x.row(1) is the second row a matrix. (x[2,])
   * Indices start at 0!
   *
   * 2. Multiple rows
   * x.rows(1,3) is x[2:4, ] using the R notation.
   */


  //Rcpp::Rcout << "Enters cpp function" << std::endl;

  //some constants or constnat matrices
  int const1 = n_obs -1;
  arma::mat const2 = sigma_inv_beta_p*beta_m;
  arma::mat const3 = sigma_inv_alpha_p*alpha_m;

  // define some matrices
  arma::mat Xs = X.rows(0, const1);
  arma::mat Xs_t = Xs.t(); //the transpose of the design matrix
  //int n = X.n_rows; //the number of sites (survey + unsurvey)
  int N = W_vb.n_rows; //the total number of visits

  arma::mat Ks = K.rows(0, const1);
  arma::mat Ks_t = Ks.t(); //the transpose of the spatial design matrix
  int r = K.n_cols; //the number of column of K; the number of spatial random effects added

  NumericVector siteindex(n_obs); //siteindex = 1, 2, ..., number of surveyed sites
  for (int idown =0; idown<n_obs; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0);
  //identify the indices associated with sites where species was seen (from surveyd sites)

  uvec not_observed_index = find(ysum==0);
  //identify the indices associated with sites where species was not seen (from surveyd sites)

  int n_not_observed_surveyed_sites = size(find(ysum.rows(0, const1)==0))(0);
  //the number of sites with 0 detections based on surveyed sites

  uvec unsurveyed_index = find(unsurveyed_ind>0) + n_obs;
  int n_unsurveyed = size(unsurveyed_index)(0);

  //construct a matrix with the start and end for each site
  arma::vec starts(n_obs);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, const1 ) = ends.rows(0, const1-1)+1;
  //-----------------------------------------------------------------------

  //Declare parameter vectors
  arma::vec psi = z; //set the starting value of psi to z
  //arma::vec psi_occ;
  //arma::vec common_term;
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the alpha matrix
  mat theta = zeros<mat>(r,1); //initialise the theta matrix to 0's; i.e. no spatial effects initiall
  double tau=tau_0; //initialize tau =1; the spatial precision scalar
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat Xs_beta; //Xs*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat beta_cov;
  //arma::mat beta_mu;

  //Posterior sampling of theta
  arma::mat Ks_theta; //Ks*theta
  arma::mat theta_cov;
  //arma::mat theta_mu;
  //-----------------------------------------------------------------------

  //Posterior sampling of tau
  mat inv_scale;
  mat tau_mat(1,1);
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  //uvec z_equals1_rows;//identify all indices with z==1
  //uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  arma::mat ysum_iter;
  //int counter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat alpha_cov;
  //arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //The outputs
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept);
  arma::mat post_beta(beta.n_rows, num_samples_kept);
  arma::mat post_theta(theta.n_rows, num_samples_kept);
  arma::mat post_tau(1, num_samples_kept);
  arma::mat post_z(z.n_rows, num_samples_kept);
  arma::mat post_psi(z.n_rows, num_samples_kept);
  //-----------------------------------------------------------------------

  //now do the sampling here
  Xs_beta = Xs*beta; //arma::mat
  Ks_theta = Ks*theta; //arma::mat


  //Rcpp::Rcout << "Starts sampling" << std::endl;


  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(Xs_beta + Ks_theta)) );
    pg_beta = rpg5(Xs_beta + Ks_theta);
    //-----------------------------------------------------------------------

    //posterior samples for beta
    beta_cov = inv_sympd( sigma_inv_beta_p + Xs_t*diagmat( pg_beta )*Xs ); //arma::mat
    beta = mvnrnd(beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov, 1);
    //beta = mvrnormArma2(1, beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov);
    Xs_beta = Xs*beta; //arma::mat
    //-----------------------------------------------------------------------


    //Rcpp::Rcout << "Posterior samples for beta" << std::endl;


    //posterior samples for theta
    theta_cov = inv_sympd( tau*Minv + Ks_t*diagmat( pg_beta )*Ks );
    theta = mvnrnd(theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov, 1);
    //theta = mvrnormArma2(1, theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov);
    Ks_theta = Ks*theta; //arma::mat
    //-----------------------------------------------------------------------



    //Rcpp::Rcout << "Posterior samples for theta" << std::endl;

    //posterior samples for tau
    /*inv_scale =  theta.t()*Minv*theta *0.5 + b_tau ;
     tau = randg(1, distr_param(a_tau + 0.5*r, 1.0/inv_scale(0)))(0);
     tau_mat(0,0) = tau;*/

    inv_scale =  theta.t()*Minv*theta*0.5 + b_tau ; //not inverse scale. rgammadouble uses 1/inv_scale = 1/( theta.t()*Minv*theta *0.5 + i2 );
    tau = rgammadouble(1, a_tau + 0.5*r, inv_scale(0) );
    tau_mat(0,0) = tau;
    //-----------------------------------------------------------------------


    //Rcpp::Rcout << "Posterior samples for tau" << std::endl;


    //convert arma::vec z to NumericVector zNM (but only for surveyed locations)
    NumericVector z_NM = wrap(z.rows(0, const1));
    //Rcpp::Rcout << "1 arma::vec z to NumericVector zNM" << std::endl;

    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    //Rcpp::Rcout << "2 arma::vec z to NumericVector zNM" << std::endl;

    IntegerVector ind_z = match(siteids, z_pos);
    //Rcpp::Rcout << "3 arma::vec z to NumericVector zNM" << std::endl;

    arma::vec convert_z_vec= as<arma::vec>(ind_z);

    // Added this
    arma::vec convert_z_vec1 = as<arma::vec>(z_pos);

    //Rcpp::Rcout << "4 arma::vec z to NumericVector zNM" << std::endl;


    // W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    W_iter = W_vb.rows(find(convert_z_vec1>0));
    //Rcpp::Rcout << "5 arma::vec z to NumericVector zNM" << std::endl;

    //Look here
    Y_iter = Y.rows( find(convert_z_vec>0) );
    //Rcpp::Rcout << "y iterations" << std::endl;
    ysum_iter = ysum.rows(find(convert_z_vec1>0));
    arma::mat nvisits_iter;
    nvisits_iter = nvisits.rows(find(convert_z_vec1>0));

    //-----------------------------------------------------------------------

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat
    //Rcpp::Rcout << W_iter << std::endl;

    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) ); //arma::vec

    pg_alpha = rnorm7(nvisits_iter, W_alpha);


    //posterior samples for alpha



    alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    alpha = mvnrnd(alpha_cov*( W_iter.t()*( ysum_iter - nvisits_iter/2 ) + const3), alpha_cov, 1);
    //alpha = mvrnormArma2(1, alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov);
    //-----------------------------------------------------------------------


    //update the p, psi and z
    p = 1/(1 + exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    //psi.elem(not_observed_index) = 1/(1 + exp(-X.rows(not_observed_index)*beta - K.rows(not_observed_index)*theta ));
    psi = 1/(1 + exp(-X*beta - K*theta ));

    //Rcpp::Rcout << "p, psi updated" << std::endl;


    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    //for surveyed sites with zero detections
    //this can be made more efficient!
    for (int idown=0; idown< n_not_observed_surveyed_sites; idown++){
      //prod_p_start = starts.row(not_observed_index(idown))-1;
      //prod_p_end = ends.row(not_observed_index(idown))-1;

      // what is going on here // need tp check !!!
      // prob of occupying given that it is P(not d|o) + p(not occupied)
      prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/
        ( psi.row(not_observed_index(idown))*
          pow((1-p.row(not_observed_index(idown))),nvisits(idown) ) ));

      //prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/
      //( psi.row(not_observed_index(idown))*
      //arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );


      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }

    // Rcpp::Rcout << "z updated not obs surveyed" << std::endl;
    //
    // Rcpp::Rcout << "NOT SAMPLED FOR UNSURVEYED" << std::endl;

    //for unsurveyed sites

    for (int idown=0; idown< n_unsurveyed; idown++){
      prob = psi.row(n_obs + idown);
      zdraw = rbinom(1,1, prob(0));
      z(n_obs + idown) = zdraw(0);
    } //finish sampling the z matrix
    //Rcpp::Rcout << "z updated for unsurveyed" << std::endl;

    //idea from stocc package
    /*
     common_term = exp(ln_invlogit( X*beta + K*theta ) + Mt*log_not_prob(W_vb*alpha));
     psi_occ = common_term/( common_term + exp( log_not_prob(X*beta) ) );

     //sample z
     for (int idown=0; idown< n; idown++){
     prob = psi_occ.row(idown);
     zdraw = rbinom(1,1, prob(0));
     z(idown) = zdraw(0);
     }//finish sampling the z matrix
     */

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;

    if (isamples_counter>=0){
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
      post_theta.col(isamples_counter) = theta;
      post_tau.col(isamples_counter) = tau_mat;
      post_z.col(isamples_counter) = z;
      post_psi.col(isamples_counter) = psi;
    }
  }

  return List::create(_["alpha"]=post_alpha,
                      _["beta"]=post_beta,
                      _["theta"]=post_theta,
                      _["tau"]=post_tau,
                      _["real.occ"]= post_z,
                      _["psi"]=post_psi,
                      _["K"]=K);

}





/*
 * New spatial occupancy model using the binomial model
 */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List logitoccSPATBINOMPG(arma::mat X, arma::mat W_vb, arma::mat Y, arma::mat z, arma::vec ysum, arma::vec nvisits,
                         arma::mat K, arma::mat Minv, //arma::mat Mt,
                         double n_obs,
                         NumericVector siteids, arma::vec unsurveyed_ind,
                         double tau_0, double a_tau, double b_tau,
                         arma::mat alpha_m, arma::mat beta_m,
                         arma::mat sigma_inv_alpha_p, arma::mat sigma_inv_beta_p,
                         int ndraws, double percent_burn_in){

  /*
   * Called by the R function, occSPATlogit
   *
   * Some Rcpp notes
   * Some R translations put in brackets at times.
   *
   * 1. Indices
   * arma::mat x;
   * x.row(0) is the first row a matrix. (x[1,])
   * x.row(1) is the second row a matrix. (x[2,])
   * Indices start at 0!
   *
   * 2. Multiple rows
   * x.rows(1,3) is x[2:4, ] using the R notation.
   */


  //Rcpp::Rcout << "Enters cpp function" << std::endl;

  //some constants or constnat matrices
  int const1 = n_obs -1;
  arma::mat const2 = sigma_inv_beta_p*beta_m;
  arma::mat const3 = sigma_inv_alpha_p*alpha_m;

  // define some matrices
  arma::mat Xs = X.rows(0, const1);
  arma::mat Xs_t = Xs.t(); //the transpose of the design matrix
  //int n = X.n_rows; //the number of sites (survey + unsurvey)
  int N = W_vb.n_rows; //the total number of visits

  arma::mat Ks = K.rows(0, const1);
  arma::mat Ks_t = Ks.t(); //the transpose of the spatial design matrix
  int r = K.n_cols; //the number of column of K; the number of spatial random effects added

  NumericVector siteindex(n_obs); //siteindex = 1, 2, ..., number of surveyed sites
  for (int idown =0; idown<n_obs; idown++){siteindex(idown) = idown+1;}

  uvec observed_index = find(ysum>0);
  //identify the indices associated with sites where species was seen (from surveyd sites)

  uvec not_observed_index = find(ysum==0);
  //identify the indices associated with sites where species was not seen (from surveyd sites)

  int n_not_observed_surveyed_sites = size(find(ysum.rows(0, const1)==0))(0);
  //the number of sites with 0 detections based on surveyed sites

  uvec unsurveyed_index = find(unsurveyed_ind>0) + n_obs;
  int n_unsurveyed = size(unsurveyed_index)(0);

  //construct a matrix with the start and end for each site
  arma::vec starts(n_obs);
  starts(0) = 1; //might be better to make it start at 0 at a later stage!!!
  arma::vec ends = cumsum(nvisits);
  starts.rows( 1, const1 ) = ends.rows(0, const1-1)+1;
  //-----------------------------------------------------------------------

  //Declare parameter vectors
  arma::vec psi = z; //set the starting value of psi to z
  //arma::vec psi_occ;
  //arma::vec common_term;
  arma::mat p = zeros<arma::mat>(N, 1); //initialize p to contain 0's
  arma::mat alpha = alpha_m; //initialize the alpha matrix
  arma::mat beta = beta_m; //initialize the alpha matrix
  mat theta = zeros<mat>(r,1); //initialise the theta matrix to 0's; i.e. no spatial effects initiall
  double tau=tau_0; //initialize tau =1; the spatial precision scalar
  //-----------------------------------------------------------------------

  //Posterior sampling of beta
  arma::mat Xs_beta; //Xs*beta
  arma::vec pg_beta; //Polya-Gamma variates used to sample beta
  arma::mat beta_cov;
  //arma::mat beta_mu;

  //Posterior sampling of theta
  arma::mat Ks_theta; //Ks*theta
  arma::mat theta_cov;
  //arma::mat theta_mu;
  //-----------------------------------------------------------------------

  //Posterior sampling of tau
  mat inv_scale;
  mat tau_mat(1,1);
  //-----------------------------------------------------------------------

  //construct the W matrix associated with z=1
  //uvec z_equals1_rows;//identify all indices with z==1
  //uvec z_equals1_allrows(sum(nvisits)); //used to identify all rows of W_vb where z==1
  arma::mat W_iter;
  arma::mat Y_iter;
  arma::mat ysum_iter;
  //int counter;
  //-----------------------------------------------------------------------

  //Posterior sampling of alpha
  arma::mat W_alpha;
  arma::vec pg_alpha; //Polya-Gamma variates used to sample alpha
  arma::mat alpha_cov;
  //arma::mat alpha_mu;
  //-----------------------------------------------------------------------

  //update the p, psi and z
  NumericVector zdraw(1);
  arma::vec prob(1);
  arma::vec prod_p_start; arma::vec prod_p_end;
  //-----------------------------------------------------------------------

  //The outputs
  int isamples_counter;
  int num_burnin = floor(ndraws*percent_burn_in);
  int num_samples_kept = ndraws - num_burnin;
  arma::mat post_alpha(alpha.n_rows , num_samples_kept);
  arma::mat post_beta(beta.n_rows, num_samples_kept);
  arma::mat post_theta(theta.n_rows, num_samples_kept);
  arma::mat post_tau(1, num_samples_kept);
  arma::mat post_z(z.n_rows, num_samples_kept);
  arma::mat post_psi(z.n_rows, num_samples_kept);
  //-----------------------------------------------------------------------

  //now do the sampling here
  Xs_beta = Xs*beta; //arma::mat
  Ks_theta = Ks*theta; //arma::mat


  //Rcpp::Rcout << "Starts sampling" << std::endl;


  for (int isamples=0; isamples<ndraws; isamples++){

    //add in an interuptor. i.e. escape if the user cancels operations
    //checks every 1000 iterations
    if (isamples % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //Generate the Poly-Gamma variables for beta vector
    //pg_beta = as<arma::vec>( parallelMatrixRpg(wrap(Xs_beta + Ks_theta)) );
    pg_beta = rpg5(Xs_beta + Ks_theta);
    //-----------------------------------------------------------------------

    //posterior samples for beta
    beta_cov = inv_sympd( sigma_inv_beta_p + Xs_t*diagmat( pg_beta )*Xs ); //arma::mat
    beta = mvnrnd(beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov, 1);
    //beta = mvrnormArma2(1, beta_cov*( Xs_t*( z.rows(0,n_obs-1) - 0.5 - diagmat( pg_beta )*( Ks_theta) ) + const2), beta_cov);
    Xs_beta = Xs*beta; //arma::mat
    //-----------------------------------------------------------------------


    //Rcpp::Rcout << "Posterior samples for beta" << std::endl;


    //posterior samples for theta
    theta_cov = inv_sympd( tau*Minv + Ks_t*diagmat( pg_beta )*Ks );
    theta = mvnrnd(theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov, 1);
    //theta = mvrnormArma2(1, theta_cov*Ks_t*( z.rows(0, const1) - 0.5 - diagmat( pg_beta )*( Xs_beta  ) ), theta_cov);
    Ks_theta = Ks*theta; //arma::mat
    //-----------------------------------------------------------------------



    //Rcpp::Rcout << "Posterior samples for theta" << std::endl;

    //posterior samples for tau
    /*inv_scale =  theta.t()*Minv*theta *0.5 + b_tau ;
     tau = randg(1, distr_param(a_tau + 0.5*r, 1.0/inv_scale(0)))(0);
     tau_mat(0,0) = tau;*/

    inv_scale =  theta.t()*Minv*theta*0.5 + b_tau ; //not inverse scale. rgammadouble uses 1/inv_scale = 1/( theta.t()*Minv*theta *0.5 + i2 );
    tau = rgammadouble(1, a_tau + 0.5*r, inv_scale(0) );
    tau_mat(0,0) = tau;
    //-----------------------------------------------------------------------


    //Rcpp::Rcout << "Posterior samples for tau" << std::endl;


    //convert arma::vec z to NumericVector zNM (but only for surveyed locations)
    NumericVector z_NM = wrap(z.rows(0, const1));
    //Rcpp::Rcout << "1 arma::vec z to NumericVector zNM" << std::endl;

    NumericVector z_pos = siteindex[z_NM>0]; //find out where z>0
    //Rcpp::Rcout << "2 arma::vec z to NumericVector zNM" << std::endl;

    IntegerVector ind_z = match(siteids, z_pos);
    //Rcpp::Rcout << "3 arma::vec z to NumericVector zNM" << std::endl;

    arma::vec convert_z_vec= as<arma::vec>(ind_z);

    // Added this
    arma::vec convert_z_vec1 = as<arma::vec>(z_pos);

    //Rcpp::Rcout << "4 arma::vec z to NumericVector zNM" << std::endl;


    // W_iter = W_vb.rows( find(convert_z_vec>0) ) ;
    W_iter = W_vb.rows(find(convert_z_vec1>0));
    //Rcpp::Rcout << "5 arma::vec z to NumericVector zNM" << std::endl;

    //Look here
    Y_iter = Y.rows( find(convert_z_vec>0) );
    //Rcpp::Rcout << "y iterations" << std::endl;
    ysum_iter = ysum.rows(find(convert_z_vec1>0));
    arma::mat nvisits_iter;
    nvisits_iter = nvisits.rows(find(convert_z_vec1>0));

    //-----------------------------------------------------------------------

    //Generate the Poly-Gamma variables for alpha vector
    W_alpha = W_iter*alpha; // arma::mat
    //Rcpp::Rcout << W_iter << std::endl;

    //pg_alpha = as<arma::vec>( parallelMatrixRpg(wrap(W_alpha)) ); //arma::vec





    pg_alpha = rpg7(nvisits_iter, W_alpha); //arma::vec
    //pg_alpha = rnorm7(nvisits_iter, W_alpha);

    //posterior samples for alpha


    alpha_cov = inv_sympd( sigma_inv_alpha_p + W_iter.t()*diagmat( pg_alpha )*W_iter );
    alpha = mvnrnd(alpha_cov*( W_iter.t()*( ysum_iter - nvisits_iter/2 ) + const3), alpha_cov, 1);
    //alpha = mvrnormArma2(1, alpha_cov*( W_iter.t()*( Y_iter - 0.5 ) + const3), alpha_cov);
    //-----------------------------------------------------------------------


    //update the p, psi and z
    p = 1/(1 + exp(-W_vb*alpha)); //plogis(W_vb*alpha)

    //update the psi matrix for those sites the species has not been seen
    //psi.elem(not_observed_index) = 1/(1 + exp(-X.rows(not_observed_index)*beta - K.rows(not_observed_index)*theta ));
    psi = 1/(1 + exp(-X*beta - K*theta ));

    //Rcpp::Rcout << "p, psi updated" << std::endl;


    /*sample z
     * remember rbinom from Rcpp requires the prob to be a 'double'
     * Also p.rows(double, double). One cannot use vectors.
     */
    //for surveyed sites with zero detections
    //this can be made more efficient!
    for (int idown=0; idown< n_not_observed_surveyed_sites; idown++){
      //prod_p_start = starts.row(not_observed_index(idown))-1;
      //prod_p_end = ends.row(not_observed_index(idown))-1;

      // what is going on here // need tp check !!!
      // prob of occupying given that it is P(not d|o) + p(not occupied)
      prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/
        ( psi.row(not_observed_index(idown))*
          pow((1-p.row(not_observed_index(idown))),nvisits(idown) ) ));

      //prob = 1/( 1 + ( 1-psi.row(not_observed_index(idown)) )/
      //( psi.row(not_observed_index(idown))*
      //arma::prod(1-p.rows(prod_p_start(0), prod_p_end(0))) ) );


      zdraw = rbinom(1,1, prob(0));
      z(not_observed_index(idown)) = zdraw(0);
    }

    // Rcpp::Rcout << "z updated not obs surveyed" << std::endl;
    //
    // Rcpp::Rcout << "NOT SAMPLED FOR UNSURVEYED" << std::endl;

    //for unsurveyed sites

    for (int idown=0; idown< n_unsurveyed; idown++){
      prob = psi.row(n_obs + idown);
      zdraw = rbinom(1,1, prob(0));
      z(n_obs + idown) = zdraw(0);
    } //finish sampling the z matrix
    //Rcpp::Rcout << "z updated for unsurveyed" << std::endl;

    //idea from stocc package
    /*
     common_term = exp(ln_invlogit( X*beta + K*theta ) + Mt*log_not_prob(W_vb*alpha));
     psi_occ = common_term/( common_term + exp( log_not_prob(X*beta) ) );

     //sample z
     for (int idown=0; idown< n; idown++){
     prob = psi_occ.row(idown);
     zdraw = rbinom(1,1, prob(0));
     z(idown) = zdraw(0);
     }//finish sampling the z matrix
     */

    //-----------------------------------------------------------------------

    //store the samples
    isamples_counter = isamples - num_burnin;

    if (isamples_counter>=0){
      post_alpha.col(isamples_counter) = alpha;
      post_beta.col(isamples_counter) = beta;
      post_theta.col(isamples_counter) = theta;
      post_tau.col(isamples_counter) = tau_mat;
      post_z.col(isamples_counter) = z;
      post_psi.col(isamples_counter) = psi;
    }
  }

  return List::create(_["alpha"]=post_alpha,
                      _["beta"]=post_beta,
                      _["theta"]=post_theta,
                      _["tau"]=post_tau,
                      _["real.occ"]= post_z,
                      _["psi"]=post_psi,
                      _["K"]=K);

}

