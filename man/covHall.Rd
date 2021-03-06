% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covHall.R
\name{covHall}
\alias{covHall}
\title{covHall}
\usage{
covHall(data, u, bf = 10, pve = 0.9, eps = 0.01, nu = 1, mu = NULL)
}
\arguments{
\item{data}{tha data as described in \code{\link{gfpca_Mar}},
\code{\link{gfpca_TwoStep}}, and \code{\link{gfpca_Bayes}}}

\item{u}{grid for evaluation}

\item{bf}{number of basis functions for smoothing}

\item{pve}{percentage of explained variance}

\item{eps}{small constant used for the diagonal if diagonal elements are
smaller}

\item{nu}{additional shrinkage parameter for the estimated mean function
(should usually be 1)}

\item{mu}{fitted mean function (if available)}
}
\description{
Internal function for estimation of the covariance matrix of the latent
process using the Hall et al. approach
}
\author{
Jan Gertheiss \email{jan.gertheiss@agr.uni-goettingen.de} and
Ana-Maria Staicu \email{astaicu@ncsu.edu}
}

