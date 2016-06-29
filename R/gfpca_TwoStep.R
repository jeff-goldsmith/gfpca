#' gfpca_TwoStep
#' 
#' Implements a two-step approach to generalized functional principal
#' components analysis for sparsely observed binary curves
#' 
#' 
#' @param data A dataframe containing observed data. Should have column names
#' \code{index} for observation times, \code{value} for observed responses,
#' and \code{id} for curve indicators.
#' @param npc prespecified value for the number of principal components (if
#' given, this overrides \code{pve}).
#' @param pve proportion of variance explained; used to choose the number of
#' principal components.
#' @param output_index Grid on which estimates should be computed. Defaults to
#' \code{NULL} and returns estimates on the timepoints in the observed dataset
#' @param type Type of estimate for the FPCs; either \code{approx} or
#' \code{naive}
#' @param basis Basis used as FPCs; in simulations, allows the use of the true
#' FPCs
#' @param nbasis Number of basis functions used in spline expansions
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de}
#' @seealso \code{\link{gfpca_TwoStep}}, \code{\link{gfpca_Bayes}}.
#' @references Gertheiss, J., Goldsmith, J., and Staicu, A.-M. (2016). A note
#' on modeling sparse exponential-family functional response curves.
#' \emph{Under Review}.
#' @examples
#' 
#' \dontrun{
#' library(mvtnorm)
#' library(boot)
#' library(refund.shiny)
#' 
#' ## set simulation design elements
#' 
#' bf = 10                           ## number of bspline fns used in smoothing the cov
#' D = 101                           ## size of grid for observations
#' Kp = 2                            ## number of true FPC basis functions
#' grid = seq(0, 1, length = D)
#' 
#' ## sample size and sparsity
#' I <- 300
#' mobs <- 7:10
#' 
#' ## mean structure
#' mu <- 8*(grid - 0.4)^2 - 3
#' 
#' ## Eigenfunctions /Eigenvalues for cov:
#' psi.true = matrix(NA, 2, D)
#' psi.true[1,] = sqrt(2)*cos(2*pi*grid)
#' psi.true[2,] = sqrt(2)*sin(2*pi*grid)
#' 
#' lambda.true = c(1, 0.5)
#' 
#' ## generate data
#' 
#' set.seed(1)
#' 
#' ## pca effects: xi_i1 phi1(t)+ xi_i2 phi2(t)
#' c.true = rmvnorm(I, mean = rep(0, Kp), sigma = diag(lambda.true))
#' Zi = c.true %*% psi.true
#' 
#' Wi = matrix(rep(mu, I), nrow=I, byrow=T) + Zi
#' pi.true = inv.logit(Wi)  # inverse logit is defined by g(x)=exp(x)/(1+exp(x))
#' Yi.obs = matrix(NA, I, D)
#' for(i in 1:I){
#'   for(j in 1:D){
#'     Yi.obs[i,j] = rbinom(1, 1, pi.true[i,j])
#'   }
#' }
#' 
#' ## "sparsify" data
#' for (i in 1:I)
#' {
#'   mobsi <- sample(mobs, 1)
#'   obsi <- sample(1:D, mobsi)
#'   Yi.obs[i,-obsi] <- NA
#' }
#' 
#' Y.vec = as.vector(t(Yi.obs))
#' subject <- rep(1:I, rep(D,I))
#' t.vec = rep(grid, I)
#' 
#' data.sparse = data.frame(
#'   index = t.vec,
#'   value = Y.vec,
#'   id = subject
#' )
#' 
#' data.sparse = data.sparse[!is.na(data.sparse$value),]
#' 
#' ## fit models
#' 
#' # 2-step with eigenfunctions from FPCA
#' fit.2step = gfpca_TwoStep(data = data.sparse, type="naive")
#' plot(mu)
#' lines(fit.2step$mu, col=2)
#' 
#' # 2-step with eigenfunctions estimated following Hall et al. (2008)
#' fit.2step = gfpca_TwoStep(data = data.sparse, type="approx")
#' plot(mu)
#' lines(fit.2step$mu, col=2)
#' plot_shiny(fit.2step)
#' 
#' }
#' 
#' @export gfpca_TwoStep
#' @importFrom gamm4 gamm4
#' @import Rcpp
gfpca_TwoStep <-  function(data, npc = NULL, pve = .9, output_index = NULL, 
                           type=c("approx", "naive", "fixed"),
                           basis=NULL, nbasis=10){
  
  
  # some data checks
  if (is.null(output_index)) { output_index = sort(unique(data['index'][[1]])) }
  
  type <- match.arg(type)
  type <- switch(type, approx = "approx", naive = "naive", fixed = "fixed")
  
  # data
  Y.vec <- data['value'][[1]]
  input_index <- data['index'][[1]]
  id.vec <- data['id'][[1]]
  
  D <- length(output_index)
  I <- length(unique(id.vec))
  Y.obs <- matrix(NA, nrow = I, ncol = D)
  
  for (i in 1:I) {
    Yi <- Y.vec[id.vec == i]
    ti <- input_index[id.vec == i]
    indexi <- sapply(ti, function(t) which(output_index == t))
    Y.obs[i,indexi] <- Yi
  }
  
  dta = data.frame(
    value = as.vector(t(Y.obs)),
    id = rep(1:I, rep(D,I)),
    index = rep(output_index, I)
  )
  
  # obtain efunctions and evalues using selected method
  if (type == "naive") {
    # a simple way to get estimates of the eigenfunctions
    # Y.pca <- fpca.sc(Yi.obs, pve=pve, npc=npc)
    Y.pca <- fpca.sc(Y.obs, pve = pve, npc = npc)
    npc <- ncol(Y.pca$efunctions)

    efunctions = Y.pca$efunctions[,1:npc]
    evalues = Y.pca$evalues[1:npc]
  } else if (type == "approx") {
    # use HMY approach to estimate the eigenfunctions
    hmy_cov <- covHall(data = data, u = output_index, bf = 10, pve = pve, eps = 0.01, nu = 1)
    
    # obtain spectral decomposition of the covariance of X
    eigen_HMY = eigen(hmy_cov)
    fit.lambda = eigen_HMY$values
    fit.phi = eigen_HMY$vectors
    
    # remove negative eigenvalues
    wp <- which(fit.lambda > 0)
    fit.lambda_pos = fit.lambda[wp]
    fit.phi <- fit.phi[,wp]
    
    if (is.null(npc)) {
      # truncate using the cumulative percentage of explained variance
      npc <- which((cumsum(fit.lambda_pos)/sum(fit.lambda_pos)) > pve)[1]
    }
    
    efunctions = fit.phi[,1:npc]
    evalues = fit.lambda[1:npc]
  } else {
    if (is.null(basis))
      stop("basis has to be supplied if type = fixed")
    if (nrow(basis) != D)
      stop("nrow(basis) has to match grid")
    
    if (!is.null(npc) && (ncol(basis) != npc))
      stop("ncol(basis) has to match npc")
    else
      npc <- ncol(basis)
  
    efunctions = basis[,1:npc]
    evalues = NULL
  }
  
  # construct elements for fitting using mixed model
  for (i in 1:npc) {
    dta <- cbind(dta, rep(efunctions[,i], I))
  }
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  
  random.structure = paste(paste0("psi", 1:npc), collapse = "+")
  random.formula = formula(paste("~(0+", random.structure, "|| id)"))
  
  # fit using mixed model
  outre <- gamm4(value ~ s(index, k = nbasis), family = "binomial", data = dta,
                 random = random.formula)
  
  scores = as.matrix(coef(outre$mer)$id[,npc:1])
  Z.gamm.fpca <- scores %*% t(efunctions[,1:npc])
  mu <- as.vector(predict.gam(outre$gam, newdata = data.frame(index = output_index)))
  Yhat <- matrix(rep(mu, I), nrow = I, byrow = TRUE) + Z.gamm.fpca
  
  ## format output
  Y = data
  index = output_index
  family = "binomial"
  Yhat = data.frame(
    value = as.vector(t(Yhat)),
    index = rep(output_index, I),
    id = rep(1:I, each = D)
  )
  
  ret.objects = c("Yhat", "Y", "scores", "mu", "efunctions", "evalues", "npc", "index", "family")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "fpca"
  return(ret)
  
}

