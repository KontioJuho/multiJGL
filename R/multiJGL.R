

#' @title Linear and nonlinear multiclass network estimation with joint regularization over categorical groups

#' @param node.covariates An nxp dimensional matrix of p covariates measured over n samples.
#' #Input data is given as a single data matrix, and the creation of observational classes is done automatically by the grouping factor.
#' @param grouping.factor A grouping factor for creating observational classes.
#'
#' @param penalty.lin Specify "fused" or "group" penalty type for the linear JGL algorithm.
#' @param penalty.nonlin Specify "fused" or "group" penalty type for the nonlinear JGL algorithm.
#' #See the explanation for the fused and group penalties in the JGL package
#' #The original JGL CRAN repository: https://CRAN.R-project.org/package=JGL
#' #The following penalty parameters are given in pairs to --
#' #separately assign the amount of regularizations for linear and nonlinear parts
#' @param lin_lambda1 The l1-penalty parameter for the linear JGL to regulate within group network densities
#' @param lin_lambda2 The l1-penalty parameter for the nonlinear JGL.
#'
#' @param nonlin_lambda1 The fusion penalty parameter for the linear JGL.
#' @param nonlin_lambda2 The fusion penalty parameter for the nonlinear JGL.
#' @param tol.linear Convergence criterion for the linear part (see the JGL package for details).
#' @param tol.nonlinear Convergence criterion for the nonlinear part.
#' #Subsampling procedure over pseudo-observations if the number of observations is large already in the original sets.
#' @param subsample.pseudo_obs Should the subsampling procedure be used over the pseudo-observations.
#' @param omit.rate An integer: Omit rate for the subsampling pcocedure between 2L and 5L
#' @param ... Additional parameter for the nonlinear JGL.
#' @import whitening
#' @import JGL
#' @import crayon
#' @importFrom CVglasso CVglasso
#' @importFrom matrixStats rowDiffs
#' @importFrom stats cov
#' @importFrom utils combn
#' @export
#'
#' @examples  print("net <- multiJGL(node.covariates, grouping.factor)")

multiJGL <- function(node.covariates = node.covariates,
                           grouping.factor = grouping.factor,
                           penalty.lin = "fused",
                           penalty.nonlin = "fused",
                           lin_lambda1 = 0.1, lin_lambda2 = 0.025,
                           nonlin_lambda1 = 0.1,nonlin_lambda2 = 0.025,
                           tol.linear = 1e-05,
                           tol.nonlinear = 1e-05,
                           subsample.pseudo_obs = FALSE,
                           omit.rate = 2L,
                           ...){

  #This algorithm is built upon the linear JGL R ´
  #CRAN repository: https://CRAN.R-project.org/package=JGL

  #Check if the input dataset contains missing values
  if(anyNA(node.covariates)) stop("NAs in node.covariates data not allowed")

  #Return a warning message if the number of observations in some group is smaller than 10
  if(is.element(TRUE, as.vector(table(grouping.factor)) < 10))
    warning("Sample size in one class is smaller than 10")

  if(is.element(TRUE, as.vector(table(grouping.factor)) < 5))
    stop("Sample size in one class is smaller than 5")


  #Check if the grouping factor contains missing values
  if(anyNA(grouping.factor))
    stop("NAs in grouping.factor data not allowed")
  #Assign the number of observational classes objects
  num.of.classes <- length(unique(grouping.factor))

  obs.classes <- vector(mode = "list", length = num.of.classes)

  #Initialize other objects for the algorithm
  whitened_data_matrices <- vector(mode = "list", length = num.of.classes)
  pseudo_obs <- vector(mode = "list", length = num.of.classes)


  for(i in 1:num.of.classes){
    obs.classes[[i]] <- node.covariates[which(as.numeric(as.factor(grouping.factor)) == i),]
    Ctest <- CVglasso(X = obs.classes[[i]], S = cov(obs.classes[[i]]),

                      #Here a relatively small lambda 0.001 should be used to provide nonsingular solution.
                      nlam = 10, lam.min.ratio = 0.0,
                      lam = 0.001, diagonal = FALSE, path = FALSE,
                      #Modify and use the following part if whitening step is based on AIC-based lambda
                      tol = 1e-04,crit.cv = "AIC",
                      maxit = 10000)

    S <- Ctest$Sigma
    if(!(is.null(S))) {
      cat(crayon::green$bold("Done: Inversion of sample covariance matrix\n")) }
    #The whitening step is based on the ZCA-cor procedure
    W.ZCAcor = whiteningMatrix(S, method="ZCA-cor")
    #Compute the whitened data matrix
    whitened_data_matrices[[i]] = tcrossprod(as.matrix(obs.classes[[i]]), W.ZCAcor)
    whitened_domain <- scale(whitened_data_matrices[[i]])
    #Sample index ordering is required for generating the pseudo-observations
    whitened_domain <- whitened_domain[order(whitened_domain[,1]),]
    tr_whitened_domain <- t(whitened_domain)

    #Enumerate index pairs of all possible pairs of observations
    cols <- combn( ncol(tr_whitened_domain), 2)
    #Calculate the kernel specific differences among observations
    tripleSums <- apply( cols, 2, function(z) rowDiffs(tr_whitened_domain[,z]))
    whitened_domain <- t(tripleSums)
    pseudo_obs[[i]] <- scale(sqrt(abs(whitened_domain)))
    if(!(is.null(pseudo_obs[[i]]))) {
      cat(crayon::green$bold("Done: Pseudo-observations generated succesfully")) }
    colnames(pseudo_obs[[i]]) <- colnames(node.covariates) #These are the final pseudo-observations

  }

  #Estimate the linear and nonlinear network structures with the original JGL-algorithm

  lin_fgl.results = JGL(Y=obs.classes, penalty=penalty.lin, lambda1=lin_lambda1,
                        lambda2=lin_lambda2,
                        return.whole.theta = TRUE,
                        tol = tol.linear)

  #This part subsamples the pseudo-observation set
  if(subsample.pseudo_obs == TRUE){

    if(!(is.integer(omit.rate))) stop("omit.rate must be an integer (2L, 3L, 4L or 5L)")

    if(omit.rate < 2 | omit.rate > 5) stop("omit.rate must be between 2 and 5")
    sub.pseudosample <-  lapply(pseudo_obs, function(class){

      return(class[c(TRUE,rep(FALSE, omit.rate-1)),])})

    nonlin_fgl.results = JGL(Y=sub.pseudosample, penalty=penalty.nonlin,
                             lambda1=nonlin_lambda1,
                             lambda2=nonlin_lambda2,
                             return.whole.theta = TRUE,
                             tol = tol.nonlinear,...)
  } else {

    nonlin_fgl.results = JGL(Y=pseudo_obs, penalty=penalty.nonlin,
                             lambda1=nonlin_lambda1,
                             lambda2=nonlin_lambda2,
                             return.whole.theta = TRUE,
                             tol = tol.nonlinear,...)
  }


  #The linear and nonlinear networks are stored separately into the same list
  results <- list(linear_networks = lin_fgl.results,
                  nonlinear_networks = nonlin_fgl.results)

}