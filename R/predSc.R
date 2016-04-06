#' predSc
#' 
#' Internal function for prediction of scores using the Hall et al. approach.
#' 
#' @param ev vector of eigenvalues
#' @param psi matrix of eigenfunctions
#' @param Yi.obs matrix of observed data
#' @param mu mean function
#' @param gs dispersion parameter (gamma in Hall et al.). This is gm in gfpcaMar. In
#' the simulation studies this parameter is chosen such the MSE is smallest,
#' which gives somewhat over-optimistic results for the marginal approach.
#' 
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de} and 
#' Ana-Maria Staicu \email{astaicu@@ncsu.edu}
#' 
predSc <- function(ev, psi, Yi.obs, mu, gs){
  n <- nrow(Yi.obs)
  Scores <- matrix(NA, n, length(ev))
  for (i in 1:n)
  {
    ti <- which(!is.na(Yi.obs[i,]))
    yi <- Yi.obs[i,ti]
    di <- (yi - plogis(mu[ti]))/deriv.inv.logit(mu[ti])
    Sigma <- psi%*%diag(ev)%*%t(psi)
    dS <- (gs^2)*(plogis(mu)*(1-plogis(mu)))/(deriv.inv.logit(mu))^2
    diag(Sigma) <- diag(Sigma) + dS
    xi <- ev*(t(psi[ti,])%*%solve(Sigma[ti,ti])%*%di)
    Scores[i,] <- xi
  }
  return(Scores)
}
