#' covHall
#' 
#' Internal function used to compute the latent covariance
#' 
#' @param data
#' @param u
#' @param bf
#' @param pve
#' @param spectral_decomp
#' @param eps
#' @param nu
#' @param mu.fit
#' 
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de}
#' 
covHall <- function(data, u, bf=10, pve=.9, spectral_decomp=T, eps=0.01, nu=1,
                    mu.fit=NULL){
  
  Y.vec <- data['.value'][[1]]
  t.vec <- data['.index'][[1]]
  id.vec <- data['.id'][[1]]
  
  grid <- unique(t.vec)
  D <- length(grid)
  I <- length(unique(id.vec))
  Y_miss <- matrix(NA, nrow=I, ncol=D)
  
  if (length(mu.fit)==0)
  {
    out <- gam(Y.vec ~ s(t.vec), family = "binomial")
    mu.fit <- as.vector(predict.gam(out, newdata=data.frame(t.vec = u)))
  }
  
  for(i in 1:I){
    Yi <- Y.vec[id.vec==i]
    ti <- t.vec[id.vec==i]
    indexi <- sapply(ti, function(t) which(grid==t))
    Y_miss[i,indexi] <- Yi
  }
  
  Yi_2 <- sapply(c(1:I), function(k)  Y_miss[k,] %x%   Y_miss[k,])
  Yi_2mom <- matrix(rowMeans(Yi_2, na.rm=T), ncol=D)
  diag(Yi_2mom) <- NA
  
  row.vec <- rep(grid, each = D)       # set up row variable for bivariate smoothing
  col.vec <- rep(grid, D)              # set up column variable for bivariate smoothing
  
  Yi_2mom_sm <- matrix(predict.gam(gam(as.vector(Yi_2mom)~te(row.vec,col.vec,k=bf)) ,
                                   data.frame(row.vec =rep(u, each=length(u)),col.vec=rep(u, length(u)) )), ncol=length(u))
  Yi_2mom_sm <- (Yi_2mom_sm +t(Yi_2mom_sm))/2
  
  Y_miss_mean <- colMeans(Y_miss, na.rm=T)
  Y.mean_sm <- as.vector(predict(gam(Y_miss_mean~s(grid), method="REML"),
                                 newdata=data.frame(grid=u)))
  
  Yi.cov_sm <- Yi_2mom_sm - matrix(Y.mean_sm, ncol=1)%*%matrix(Y.mean_sm, nrow=1)
  
  Zi.cov_sm <- diag(1/deriv.inv.logit(nu*mu.fit)) %*% Yi.cov_sm %*%diag(1/deriv.inv.logit(nu*mu.fit))
  
  ddd <- diag(Zi.cov_sm)
  diag(Zi.cov_sm) <- ifelse(ddd<eps, eps, ddd)
  
  Zi.cov_sm
}
