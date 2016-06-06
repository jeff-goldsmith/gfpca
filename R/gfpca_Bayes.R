#' gfpca_Bayes
#' 
#' Implements a Bayesian approach to generalized functional principal
#' components analysis for sparsely observed binary curves
#' 
#' 
#' @param data A dataframe containing observed data. Should have column names
#' \code{.index} for observation times, \code{.value} for observed responses,
#' and \code{.id} for curve indicators.
#' @param npc Number of FPC basis functions to estimate
#' @param grid Grid on which estimates should be computed. Defaults to
#' \code{NULL} and returns estimates on the timepoints in the observed dataset
#' @param nbasis Number of basis functions used in spline expansions
#' @param iter Number of sampler iterations
#' @param warmup Number of iterations discarded as warmup
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de}
#' @references Gertheiss, J., Goldsmith, J., and Staicu, A.-M. (2016). A note
#' on modeling sparse exponential-family functional response curves.
#' \emph{Under Review}.
#' @examples
#' 
#' \dontrun{
#' library(mvtnorm)
#' library(boot)
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
#'   .index = t.vec,
#'   .value = Y.vec,
#'   .id = subject
#' )
#' 
#' data.sparse = data.sparse[!is.na(data.sparse$.value),]
#' 
#' ## fit Bayesian models
#' fit.Bayes = gfpca_Bayes(data = data.sparse)
#' plot(mu)
#' lines(fit.Bayes$mu, col=2)
#' 
#' }
#' 
#' @export gfpca_Bayes
#' @import rstan
#' @importFrom splines bs
gfpca_Bayes <- function(data, npc=3, grid = NULL, nbasis=10, iter=1000, warmup=400){
  
  ## implement some data checks
  
  if(is.null(grid)){ grid = sort(unique(data['.index'][[1]])) }
  
  I = length(unique(data['.id'][[1]]))
  D = length(grid)
  
  X.des = matrix(1, nrow = I, ncol = 1)
  p = 1
  
  BS.pen = bs(grid, df=nbasis, intercept=TRUE, degree=3)
  BS = bs(data['.index'][[1]], df=nbasis, intercept=TRUE, degree=3)
  
  alpha = .1
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(BS.pen) %*% t(diff0) %*% diff0 %*% BS.pen
  P2 = t(BS.pen) %*% t(diff2) %*% diff2 %*% BS.pen
  P.mat = alpha * P0 + (1-alpha) * P2
  
  
  ## format sparse data for stan fitting
  Y.vec.obs = data['.value'][[1]]
  subject.obs = data['.id'][[1]]
  n.total = length(Y.vec.obs)
  
  # stanfile = file.path("exec", "gfpca.stan")
  # stanfile = system.file(stanfile, package = "gfpca")
  # stanfit <- rstan::stanc_builder(file = stanfile)
  # stanfit = stan_model(stanfile, model_name = stanfit$model_name)
  # stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
  #                           model_cppcode = stanfit$cppcode)
  # stanmodels = do.call(
  #   methods::new, args = c(
  #   stanfit[-(1:3)], 
  #   Class = "stanmodel", 
  #   mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))

  ## fit model using STAN
  stanmodels = make_stanmodel("gfpca")
  stanfit = stanmodels$gfpca
  
  dat = list(Y = Y.vec.obs, X = X.des, BS = BS,
             subjId = subject.obs,
             N = n.total, I = I, D = D, p = p, Kt = nbasis, Kp = npc, 
             PenMat = P.mat)

  GenFPCA.fit = sampling(stanfit,
                         data=dat, iter = iter, 
                         warmup = warmup,
                         control = list(adapt_delta = .65), 
                         chains = 1, verbose = FALSE)
  
  ## post process to obtain point estimates
  beta.post = extract(GenFPCA.fit, "beta")$beta
  beta_psi.post = extract(GenFPCA.fit, "beta_psi")$beta_psi
  c.post = extract(GenFPCA.fit, "c")$c
  
  betaHat.post = array(NA, dim = c(p, D, dim(c.post)[1]))
  for(i in 1:dim(c.post)[1]) {
    betaHat.post[,,i] = (beta.post[i,,] %*% t(BS.pen))
  }
  
  y.post = z.post = array(NA, dim = c(I, D, dim(c.post)[1]))
  for(i in 1:dim(c.post)[1]) {
    y.post[,,i] = X.des %*% (beta.post[i,,] %*% 
                               t(BS.pen)) + c.post[i,,] %*% 
      (beta_psi.post[i,,] %*% t(BS.pen))
    z.post[,,i] = c.post[i,,] %*% (beta_psi.post[i,,] %*% t(BS.pen))
  }
  Zstan = apply(z.post, c(1,2), mean)
  W.bayes = apply(y.post, c(1,2), mean)
  
  ret = list(apply(betaHat.post, 2, mean), Zstan, W.bayes)
  names(ret) = c("mu", "z", "yhat")
  ret

}
  
