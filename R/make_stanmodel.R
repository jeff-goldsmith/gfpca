#' @title Make Stan Model from .Stan file
#' @description Simple wrapper to make stanmodel
#' @param mod_name Name of the .stan file
#' @return Object of class \code{stanmodel}
#' @import methods
make_stanmodel = function(mod_name = "gfpca") {
  exec_folder = system.file("exec", 
                            package = "gfpca")
  stanfile = file.path(exec_folder, 
                       paste0(mod_name, ".stan"))
  stopifnot(file.exists(stanfile))
  chunk_fol = system.file("chunks", package = "gfpca") 

    
  stanfit <- rstan::stanc_builder(stanfile, isystem = chunk_fol)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
                            model_cppcode = stanfit$cppcode)
  stanmodels = do.call(methods::new, 
                       args = c(stanfit[-(1:3)], 
                                Class = "stanmodel", 
                                mk_cppmodule = function(x) {
                                  get(paste0("model_", mod_name))
                                  
                                })
  )
  # names(stanmodels) <- mod_name
  
  return(stanmodels)
}
