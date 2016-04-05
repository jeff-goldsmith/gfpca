#' predSc
#' 
#' Internal function used to compute predicted scores
#' 
#' @param ev
#' @param psi
#' @param Yi.obs
#' @param mu
#' @param gs
#' 
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de}
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
