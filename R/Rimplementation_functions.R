Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

temp2 <- function(){}

vb_Designs_check<-function(formula, Names){
     #perform some checks on the design matrices
     #check that the names in the formula call are in the dataframes provided
     detVars<-all.vars(formula) #return names contained in formula

     if ((sum(detVars[-1]%in% Names)==length(detVars[-1]))!= 1){
          #if the condition '(sum(detVars[-1]%in% Names)==length(detVars[-1]))'
          #is false print the error below and terminate the function
          stop(print("\n CHECK YOUR FORMULA CALL. \n MISMATCH BETWEEN CALL AND
                     DESIGN MATRICES. \n ."))
     }
}

vb_ReqDesigns<-function(formula, design_mats){

     vb_Designs_check(formula, design_mats$Names)
     #create the W matrix
     W <- model.matrix(as.formula(paste("~",formula[3],sep="")), data=design_mats$W)

     #create the X matrix
     f_xmat <- paste(formula[[2]])
     X <- model.matrix(as.formula(paste(f_xmat[1],f_xmat[3],sep="")),
                     data=design_mats$X)

     list(W=W, X=X)
}

vb_Designs<-function(W, X, y){
     #create the required 'response' and 'regressor matrices'
     #using all of the X and W data
     #the output is stored as a named list
     #create the Y matrix that will be used

     #Y is a n*J times 1 matrix
     Y <- matrix(na.omit(matrix(t(y), ncol=1)))

     pres_abs <- apply(y,1,max,na.rm=T)

     #create the W matrix
     W.temp <- NULL
     nv <- length(W)
     for (i in 1:nv){ W.temp <- cbind(W.temp, W[[i]]) }

     nvisits <- apply(W[[1]],1,function(x){length(na.omit(x))})
     n <- length(nvisits)

     W.out<-NULL
     for (i in 1:n){ W.out <- rbind(W.out, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )}
     colnames(W.out) <- names(W)

     siteids <- rep(1:n, nvisits)

     list(Y=as.data.frame(Y), X=as.data.frame(X), W=as.data.frame(W.out),
          Names=c( colnames(X), colnames(W.out)), nvisits=nvisits,
          pres_abs=pres_abs, siteids=siteids, y=y)
}

dRUMocc <- function(formula, design_mats,
                    ndraws=1,
                    alpha_m, beta_m,
                    sigma_inv_alpha_p, sigma_inv_beta_p){
    #The R function that uses the dRUM algorithm to fit a Bayesian
    #single season model
    #This function does not allow the user to specify informative
    #priors for alpha and beta
    #Not much functionality added to this function since PGocc4 is preferred

    req_design_mats  <-  vb_ReqDesigns(formula, design_mats)
    W_vb <- req_design_mats$W
    Y <- matrix(design_mats$Y$V1)
    X <- req_design_mats$X
    nvisits <- design_mats$nvisits #nvisits
    y <- design_mats$y
    siteids <- design_mats$siteids

    ysum <- apply(y,1,sum, na.rm=TRUE) #a vector!
    z <- design_mats$pres_abs #the inital z vector

    logitoccDA4(X, Y, W_vb, as.matrix(siteids, ncol=1), ndraws, ysum, z,
                nvisits,
                alpha_m, beta_m,
                sigma_inv_alpha_p, sigma_inv_beta_p)
}

#This is the function used for all non-spatial SSO analysis
PGocc4 <- function(formula, design_mats, ndraws=1,
                 alpha_m, beta_m,
                 sigma_inv_alpha_p, sigma_inv_beta_p,
                 percent_burn_in){
     #The R function that uses the Polya-Gamma algorithm to fit a Bayesian
     #single season model
     #This function does allow the user to specify informative priors for
     #alpha and beta
     #samples are outputted as a list

     req_design_mats <- vb_ReqDesigns(formula, design_mats)
     W_vb <- req_design_mats$W
     Y <- matrix(design_mats$Y$V1)
     X <- req_design_mats$X
     nvisits <- design_mats$nvisits #nvisits

     y <- design_mats$y
     siteids <- design_mats$siteids

     ysum <- apply(y,1,sum, na.rm=TRUE) #a vector!
     z <- design_mats$pres_abs #the inital z vector

     logitoccPG3(X, Y, W_vb, as.matrix(siteids, ncol=1), ndraws, ysum, z,
                 nvisits,
                 alpha_m, beta_m,
                 sigma_inv_alpha_p, sigma_inv_beta_p,
                 percent_burn_in)
}

#This is the function used for all non-spatial SSO analysis
PGocc <- function(formula, data_inputs,
                  ndraws=1,
                  alpha_m, beta_m,
                  sigma_inv_alpha_p, sigma_inv_beta_p,
                  percent_burn_in,
                  store_z=FALSE){
  #The R function that uses the Polya-Gamma algorithm to fit a Bayesian
  #single season model
  #This function does allow the user to specify informative priors for
  #alpha and beta
  #samples are outputted as a list
  #'data_inputs' is a list that contains W, X and y
  #W = A named list of data.frames of covariates that vary within sites. i.e.
  #The dataframes are of dimension n (number of sites surveyed) by J where each
  #row is associated with a site and each column represents a site visit. 'NA'
  #values should be entered at locations that were not surveyed.
  #X = A named data.frame that varies at site level.
  #y = An n by J matrix of the detection, non-detection data, where n is the
  #number of sites, J is the maximum number of sampling periods per site.

  design_mats <- vb_Designs(W=data_inputs$W, X=data_inputs$X, y=data_inputs$y)
  req_design_mats <- vb_ReqDesigns(formula, design_mats)
  W_vb <- req_design_mats$W
  Y <- matrix(design_mats$Y$V1)
  X <- req_design_mats$X
  nvisits <- design_mats$nvisits #nvisits

  y <- design_mats$y
  siteids <- design_mats$siteids

  ysum <- apply(y,1,sum, na.rm=TRUE) #a vector!
  z <- design_mats$pres_abs #the inital z vector

  if (store_z==FALSE){
    logitoccPG3(X, Y, W_vb, as.matrix(siteids, ncol=1), ndraws, ysum, z,
              nvisits,
              alpha_m, beta_m,
              sigma_inv_alpha_p, sigma_inv_beta_p,
              percent_burn_in)
  }else{
    logitoccPG3_z(X, Y, W_vb, as.matrix(siteids, ncol=1), ndraws, ysum, z,
                nvisits,
                alpha_m, beta_m,
                sigma_inv_alpha_p, sigma_inv_beta_p,
                percent_burn_in)
  }


}

#Code for running the Bayesian analysis using jagsUI
#a simple single season occupancy model
writeLines("
      model{

      #likelihood part
      for (i in 1:n) #loop over each site
      {
          #loop over the number of visits per site
          #make general later
          #ie the number of surveys should be generalized
          #beta vector and alpha vector should be generalized

          #state occupancy
          zb[i] ~ dbern(pz[i])
          logit(pz[i]) <- beta0 + beta1*X[i]

          for (j in 1:J) #loop over the number of visits per site
          {
               #observation process
               Y[i,j] ~ dbern(py[i,j])
               py[i,j] <- zb[i]*pd[i,j]
               logit(pd[i,j])<- alpha0 + alpha1*W[i,j,2]
          }#end loop over surveys
      }#end loop over sites

      alpha0 ~ dnorm(0, 0.001)
      alpha1 ~ dnorm(0, 0.001)

      beta0 ~ dnorm(0, 0.001)
      beta1 ~ dnorm(0, 0.001)

      #occupied <-sum(zb[])
}##model
", con = "single_season_occupancy_model.txt")
#-------------------------------------------------------------------------------

#rewriting the code so that it has a similar form to the stocc package!
#this function is used to fit spatial occupancy models
occSPATlogit <- function(detection.model, occupancy.model, spatial.model,
                       so.data,
                       prior,
                       control){
     #Date = 30 Aug 2018

     cat("\n ---------------------------------------------------------------------------------------")
     cat("\n You are fitting a Bayesian spatial occupancy model using Restricted Spatial Regression.")
     cat("\n Be patient while the sampling continues.")
     cat("\n ---------------------------------------------------------------------------------------\n")

     #Design matrix formulation - taken from stocc package
     site <- so.data$site
     visit <- so.data$visit
     site$site.idx <- factor(site[, attr(site, "site")])
     visit$site.idx <- factor(visit[, attr(visit, "site")],
                              levels = levels(site$site.idx))

     xy <- as.matrix(site[, attr(site, "coords")])
     X <- as.matrix(model.matrix(occupancy.model, site)) #occupancy design matrix
     W_vb <- as.matrix( model.matrix(detection.model, visit) ) #detection design matrix
     Y <- matrix( visit[, attr(visit, "obs")], ncol=1 )

     #starting values for the occupancy status
     #all unsurveyed sites are initially set to 0
     z <- as.matrix(ifelse(table(visit$site.idx, visit[, attr(visit, "obs")])[, "1"] > 0, 1, 0), ncol=1)
     nvisits <- as.numeric(table(visit$site.idx)) #the number of surveys to sites
     n.obs <- length(nvisits[ which(nvisits >0) ]) #only keep the sites that were surveyed
     n.site <- nrow(X) #the number of sites
     nvisits <- nvisits[nvisits>0] #the number of surveys to surveyed sites
     unsurveyed_index <- ((n.obs+1):n.site) #index for the unsurveyed locations

     siteids <- matrix(rep(1:n.obs, nvisits), ncol=1)
     #create index of siteids. ie 1,1,1,2,2,3,... 3 visits to site 1, 2 to site 2 etc
     ysum <- tapply(so.data$visit$PAdata, so.data$visit$SiteName, sum)
     #the number of detections at each surveyed site
     #--------------------------------------------------------------------------------

     #objects associated with the prior specification
     #so checking should be done here if one wants this to be "user proof"
     alpha_m <- matrix(prior$mu.d, ncol=1)
     beta_m <-  matrix(prior$mu.o, ncol=1)
     sigma_inv_alpha_p <- as.matrix(prior$Q.d)
     sigma_inv_beta_p <- as.matrix(prior$Q.o)
     #--------------------------------------------------------------------------------

     #mcmc control objects
     #so checking should be done here if one wants this to be "user proof"
     ndraws <- control$ndraws
     percent_burn_in <- control$percent_burn_in
     #--------------------------------------------------------------------------------

     cat("\n ---------------------------------------------------------------------------------------")
     cat("\n Creating (R)estricted (S)patial (R)egression matrices ...\n")
     cat(" ---------------------------------------------------------------------------------------\n")
     #format taken from stocc package

     #this is a slow bit!
     Q <- stocc::icar.Q(xy=xy, threshold=spatial.model$threshold, rho = 1)
     A <- diag(diag(Q)) - Q
     n <- NCOL(Q)
     P <- diag(n) - X %*% solve(crossprod(X), t(X))
     Op <- (nrow(A)/sum(A)) * (P %*% (A %*% P))
     numSpatre <- spatial.model$moran.cut
     e <- rARPACK::eigs(Op, numSpatre) # e<-eigen(Op, symmetric = TRUE)

     K <- e$vectors[,1:numSpatre] #we take the first q eigenvalues
     Minv <- crossprod(K,Q)%*%K #Using the entire spatial lattice
     #--------------------------------------------------------------------------------

     #so checking should be done here if one wants this to be "user proof"
     tau_0 <- prior$tau
     a.tau <- prior$a.tau
     b.tau <- prior$b.tau
     #--------------------------------------------------------------------------------
     #send things to Rcpp

     cat("\n ---------------------------------------------------------------------------------------")
     cat("\n Doing sampling! ...\n")
     cat(" ---------------------------------------------------------------------------------------\n")

     logitoccSPAT(X, W_vb, Y, z, ysum, nvisits,
                  K, Minv,
                  n.obs,
                  siteids, unsurveyed_index,
                  tau_0, a.tau, b.tau,
                  alpha_m, beta_m,
                  sigma_inv_alpha_p, sigma_inv_beta_p,
                  ndraws, percent_burn_in)
}

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
