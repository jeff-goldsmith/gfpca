.onLoad <- function(libname, pkgname) { # nocov start
	modules <- paste0("stan_fit4", names(stanmodels), "_mod")
	for (m in modules) loadModule(m, what = TRUE)
} 

.onAttach <- function(...) {
	gfpcaLib <- dirname(system.file(package = "gfpca"))
	pkgdesc <- suppressWarnings(utils::packageDescription("gfpca", lib.loc = gfpcaLib))
	if (length(pkgdesc) > 1) { # FALSE if called from devtools::document()
		builddate <- gsub(';.*$', '', pkgdesc$Packaged)
		packageStartupMessage(paste("gfpca (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
	}
}