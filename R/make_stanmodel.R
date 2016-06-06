#' @title Make Stan Model from .Stan file
#' @description Simple wrapper to make stanmodel
#' @param mod_name Name of the .stan file
#' @param ... arguments passed to \code{\link{stan_model}}
#' @return Object of class \code{stanmodel}
#' @import methods
make_stanmodel = function(mod_name = "gfpca", ...) {
  exec_folder = system.file("exec", 
                            package = "gfpca")
  stanfile = file.path(exec_folder, 
                       paste0(mod_name, ".stan"))
  stanmodels = stan_model(file = stanfile, 
                          model_name = paste0("model_", mod_name),
                          ...)
  
  return(stanmodels)
}
