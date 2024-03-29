
.onLoad <- function(...) {
  s3_register("papaja::apa_table", "meta_list")
  s3_register("papaja::apa_table", "metaKIN_table")
  invisible()
}

s3_register <- function (generic, class, method = NULL)
{ # exported from vctrs with thanks
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]
  caller <- parent.frame()
  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      asNamespace(environmentName(top))
    }
    else {
      caller
    }
  }
  get_method <- function(method, env) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    }
    else {
      method
    }
  }
  register <- function(...) {
    envir <- asNamespace(package)
    method_fn <- get_method(method)
    stopifnot(is.function(method_fn))
    if (exists(generic, envir)) {
      registerS3method(generic, class, method_fn, envir = envir)
    }
    else if (identical(Sys.getenv("NOT_CRAN"), "true")) {
      warning(sprintf("Can't find generic `%s` in package %s to register S3 method.",
                      generic, package))
    }
  }
  setHook(packageEvent(package, "onLoad"), register)
  if (isNamespaceLoaded(package)) {
    register()
  }
  invisible()
}

apa_table.meta_list <-
  function(x,
           caption = NULL,
           note = NULL,
           added_stub_head = NULL,
           col_spanners = NULL,
           midrules = NULL,
           placement = "tbp",
           landscape = FALSE,
           font_size = NULL,
           escape = FALSE,
           span_text_columns = TRUE,
           format.args = NULL,...) {

  x <- format_nicely(x, ...)
  names(x) = gsub("\\%", "\\\\%", names(x))
  names(x) = gsub("_", "\\_", names(x))
  x$Moderation <- gsub("_", "\\_", x$Moderation, fixed = TRUE)

  papaja::apa_table(
    x,
    caption = caption,
    note = note,
    added_stub_head = added_stub_head,
    col_spanners = col_spanners,
    midrules = midrules,
    placement = placement,
    landscape = landscape,
    font_size = font_size,
    escape = escape,
    span_text_columns = span_text_columns,
    format.args = format.args
  )
}

#' apa_table.metaKIN_table
#'
#' format_nicely method for apa_table
#' @param x passed to papaja::apa_table
#' @param caption passed to papaja::apa_table
#' @param note passed to papaja::apa_table
#' @param stub_indents passed to papaja::apa_table
#' @param added_stub_head passed to papaja::apa_table
#' @param col_spanners passed to papaja::apa_table
#' @param midrules passed to papaja::apa_table
#' @param placement passed to papaja::apa_table
#' @param landscape passed to papaja::apa_table
#' @param font_size passed to papaja::apa_table
#' @param escape passed to papaja::apa_table
#' @param span_text_columns passed to papaja::apa_table
#' @param format.args passed to papaja::apa_table
#' @param ... passed to papaja::apa_table

apa_table.metaKIN_table <-
  function(x,
           caption = NULL,
           note = NULL,
           stub_indents = NULL,
           added_stub_head = NULL,
           col_spanners = NULL,
           midrules = NULL,
           placement = "tbp",
           landscape = FALSE,
           font_size = NULL,
           escape = FALSE,
           span_text_columns = TRUE,
           format.args = NULL,
           ...) {

  class(x) <- "data.frame"

  if(is.null(stub_indents)){
    stub_indents <- list(attr(x, "indent"))
  }

  papaja::apa_table(
    x,
    caption = caption,
    note = note,
    stub_indents = stub_indents,
    added_stub_head = added_stub_head,
    col_spanners = col_spanners,
    midrules = midrules,
    placement = placement,
    landscape = landscape,
    font_size = font_size,
    escape = escape,
    span_text_columns = span_text_columns,
    format.args = format.args,
    ...
  )
}
