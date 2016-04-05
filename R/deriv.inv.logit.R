#' deriv.inv.logit
#' 
#' Internal function used to compute the derivative of the logit function
#' 
#' @param x Value at which derivative is computed
#' 
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de}
deriv.inv.logit = function(x){exp(x)/((1+exp(x))^2)}