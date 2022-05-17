
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("meta\033[1m\033[31mKIN\033[39m\033[22m",
                              utils::packageVersion("metaKIN")))

  if(!"metaSEM" %in% .packages()){
  suppressPackageStartupMessages(attachNamespace("metaSEM"))
    suppressPackageStartupMessages(attachNamespace("OpenMx"))
  }

  message("metaKIN is under active development.\nReport bugs to: https://github.com/conig/metaKIN/issues")

}


#' meta_list plot method
#' @param x model to print
#' @param ... additional arguments passed to ninjaForest.
#' @export
plot.meta_list = function(x, ...) {
  forest_plot(x, ..., envir = parent.frame())
}
