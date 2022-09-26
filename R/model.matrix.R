#' meta_matrix
#'
#' Creates a matrix for use in the meta-analysis
#' @param formula formula yi + vi + cluster ~ predictor1 + predictor_n
#' @param data data object
#' @param intercept bool, TRUE includes an intercept
#' @param warn bool. If TRUE, you will be warned about variables with no variance
#' @export

meta_matrix <- function(formula, data, intercept = FALSE, warn = TRUE){
  model_frame <- eval(call("model.frame", formula = str2lang(formula),
                      data = data,
                      na.action = na.pass))

  matrix.invariant <- eval(call("model.frame", formula = str2lang(formula),
                      data = data,
                      na.action = stats::na.omit))

  if(length(stats::na.omit(unique(unlist(model_frame[,-1])))) <= 1){
    warning("No variance in predictor matrix", formula)
    return(NULL)
  }

  matrx <- model.matrix2(eval(parse(text = formula)), data = model_frame)
  matrx.invariant <- model.matrix2(eval(parse(text = formula)), data = matrix.invariant)

  if(!intercept & any(grepl("(Intercept)",colnames(matrx)))){
    matrx <- subset(matrx, select = -`(Intercept)`)
    matrx.invariant <- subset(matrx.invariant, select = -`(Intercept)`)
  }

  var0 <- rep(TRUE, ncol(matrx)) # initiate variable
  names(var0) <- colnames(matrx)

  # Determine predictors with variance.
  includes_intercept <- FALSE
  for(v in names(var0)){
    if(!v %in% colnames(matrx.invariant)){
      next
    }
    # If there is not already a variable which is all 1s. Mark as intercept and skip
    if(all(matrx.invariant[,v, drop = TRUE] == 1) & !includes_intercept){
      includes_intecept <- TRUE
      var0[v] <- FALSE
      next
    }

    if(stats::sd(matrx.invariant[,v, drop = TRUE]) != 0){
      var0[v] <- FALSE
    }

  }

  var0[names(var0) == "(Intercept)"] <- FALSE

  sum_var0 <- sum(var0)

  if (sum(sum_var0) != 0) {
    if (warn) {
      warning(
        sum(var0),
        " columns in the predictor matrix had no variation and were removed:\n",
        paste(names(which(var0 == 1)), collapse = ", "),
        call. = FALSE
      )
    }
    matrx <- matrx[,which(!var0)]
  }
  matrx
}

#' is_binary
#'
#' This function tests whether a matrix column is binary
#' @param x a matrix

is_binary = function(x){
  x = stats::na.omit(x)
  apply(as.matrix(x),2,function(y) { all(y %in% 0:1) })
}

#' get_kn
#'
#' get numbers of studies and effect sizes
#' @param model model to extract kn from
#' @param param_names the names of slopes

get_kn <- function(model, param_names) {
  dat <- model$data
  dat <- dat[, grepl("^x|(cluster)", colnames(dat))]
  colnames(dat)[2:length(colnames(dat))]  <- param_names

  if (is.null(model$call$incercept.constraints)) {
    dat$Intercept <- 1
    dat <- dat[, c("cluster", "Intercept", param_names)]
  }
  out <- data.frame(moderation = names(dat)[-1])

  long_dat <-
    data.table(tidyr::pivot_longer(dat,-cluster, names_to = "moderation"))[value == 1]
  long_dat <-
    long_dat[, .(k = length(unique(cluster)), n = length(value)), by = "moderation"]

  out <- dplyr::left_join(out, long_dat, by = "moderation")
  out[is.na(out)] <- 0

  bin <- is_binary(dat[-1])
  if ("Intercept" %in% names(bin)) {
    bin["Intercept"] <- FALSE
  }

  out[!bin, c("k", "n")] <- NA

  data.frame(out)
}

#' try_even_harder
#'
#' reruns models with problems
#' @param model the model.
#' @param extraTries the number of extra tries


try_even_harder = function(model, extraTries = 25) {
  if (!summary(model)$Mx.status %in% c(0, 1)) {
    suppressMessages(model <- metaSEM::rerun(model,
                                             extraTries = extraTries))
  }

  if (!summary(model)$Mx.status %in% c(0, 1)) {
    suppressMessages(model <-
                       #all the rerun messages in bulk just get in the way.
                       metaSEM::rerun(
                         model,
                         extraTries = extraTries,
                         finetuneGradient = FALSE
                       ))
  }

  if (any(is.na(summary(model)$coefficients["Std.Error"]))) {

    pre_summary <- summary(model)
    #pre_tau <- summary(model)$coefficients[c("Tau2_2", "Tau2_3"),]

    suppressMessages(model2 <-
                       #all the rerun messages in bulk just get in the way.
                       metaSEM::rerun(model,
                                      extraTries = extraTries,
                                      autofixtau2 = TRUE))

    post_summary <- summary(model2)
    taus <-
      function(x)
        x$coefficients[rownames(x$coefficients) %in% c("Tau2_2", "Tau2_3"), "Estimate", drop = FALSE]
    pre_taus <- taus(pre_summary)
    post_taus <- taus(post_summary)

    missing_tau <- which(! rownames(pre_taus) %in% rownames(post_taus))

    tau_change <-
      merge(pre_taus, post_taus, by = "row.names", all = TRUE)
    tau_change[is.na(tau_change)] <- 0
    tau_change$diff <-
      abs(tau_change$Estimate.x - tau_change$Estimate.y)

    if (any(tau_change$diff > 0.001) | length(missing_tau) > 1) {
      return(model)
    } else{
      attr(model2, "fixed_tau") <- rownames(pre_taus)[missing_tau]
      attr(model2, "fixed_tau.originalvalue") <- pre_taus[missing_tau, "Estimate"]
      return(model2)
    }

  }

  model
}


#' mlm
#'
#' Meta regression for meta3 objects
#' @param m meta3 object
#' @param formula formula passed to model.matrix
#' @param model.name String. Name your model (optional)
#' @param .envir the environment to run in
#' @param na.adjust a bool. Should the baseline model be adjusted to use the same data as a moderated model when missing covariate data is present
#' @param store_var a function which will be returned on data object by cluster. Useful for storing variable information such as number of participants
#' @export

mlm <- function(m, formula, model.name = NULL, .envir = parent.frame(), na.adjust = TRUE, store_var = NULL){
  formula = as.character(formula)[as.character(formula) != "~"]
  formula <- gsub("^~","",formula)

    formula <-
      paste0(m$call$y, " + ", m$call$v, " + as.numeric(as.factor(",m$call$cluster, ")) ~ ", formula)

  # check meta3 used
  if(!methods::is(m, "meta3")) stop("Only meta3 models are supported. To perform a 2-level meta-analysis use the meta3 argument: RE3.constraints = 0")

  # check intercept constraints are not zero
  if (!is.null(m$call$intercept.constraints)) {
      stop("Intercepts can only be modified to zero. To do so, use the formula interface e.g. ~ predictor - 1")
  }

  call = m$call

  # .internalData <- eval(call$data, envir = .envir)

  matrx_call <- call("meta_matrix", formula = formula, data = call$data)

  matrx <- eval(
    call(
      "meta_matrix",
      formula = formula,
      data = call$data,
      intercept = TRUE,
      warn = FALSE
    ),
    envir = .envir
  )
  if(is.null(matrx)) return(NA)

  includes_intercept = any(grepl("\\(Intercept\\)",colnames(matrx)))

  if(!includes_intercept){
    call$intercept.constraints <- 0
  }

  call$x = matrx_call
  call$model.name = model.name

  m_out <- eval(call, envir = .envir)



  status <- m_out$mx.fit$output$status[[1]]
  SEs <- summary(m_out)$coefficients["Std.Error"]
  missing_SEs <- any(is.na(SEs))

  if((!status %in% c(0,1)) | missing_SEs){
    # if status is not OK, or some SEs are NAs
    m_out <- try_even_harder(m_out)
  }

  param_names <- colnames(matrx)[!colnames(matrx) %in% "(Intercept)"]
  KN <- get_kn(m_out, param_names)

  row_mismatch <- nrow(m_out$data) < nrow(m$data) & na.adjust
  tau_fixed <- attr(m_out, "fixed_tau")

  # Adjust baseline for comparison if needed ------------------------
  if(row_mismatch | !is.null(tau_fixed)){

    if(row_mismatch){
    cli::cli_alert_info("anova(): baseline model has been adjusted due to missing data in covariance matrix")
    }
    adjusted_dat <- m_out$data
    new_baseline_call <- as.list(m$call)
    new_baseline_call[c("y","v","cluster", "data")] <- sapply(c("y","v","cluster", "adjusted_dat"), as.name)

    baseline <- eval(as.call(new_baseline_call))


    if (!is.null(tau_fixed)) {
      # Refit baseline with constraied Tau value to match comparator
      # Use basleine's tau value to maintain model characteristics
      baseline_coefs <- summary(baseline)$coefficients
      if (tau_fixed == "Tau2_3") {
        new_baseline_call$RE3.constraints <- baseline_coefs["Tau2_3","Estimate"]
        m_out$call$RE3.constraints <- attributes(m_out)$fixed_tau.originalvalue
      }
      if (tau_fixed == "Tau2_2") {
        new_baseline_call$RE2.constraints <- baseline_coefs["Tau2_2","Estimate"]
        m_out$call$RE2.constraints <- attributes(m_out)$fixed_tau.originalvalue
      }
      baseline <- eval(as.call(new_baseline_call))
    }

  } else {
    baseline <- m
  }

  # Add in store_data

  out <- list(
    model = m_out,
    anova = stats::anova(m_out, baseline),
    parameter_names = param_names,
    k_n = KN
  )

  if(!is.null(store_var)){
    attr(out, "store_var") <- mlm_storevar(m_out, .envir = .envir, store_var)
  }

  class(out) <- "metalm"
  out

}

mlm_storevar <- function(model, .envir, instructions){

  model_data <- eval(model$call$data, envir = .envir)
  model_data <- data.table(model_data[row.names(model$data),])

  # Process within function by cluster

  within_results <-
    model_data[, list(. = eval(parse(text = as.character(
      instructions$within
    )[2]))), by = eval(as.character(model$call$cluster))]

  # Process between function across clusters
  out <- within_results[, eval(parse(text = as.character(instructions$between)[2]))]
  if (length(out) > 1) {
    cli::cli_alert_danger(glue::glue("The length of store_var is {length(out)}"))

   cli::cli_abort(
        "The result of store_var must be length 1"
    )
  }
  out
}

#' print.metalm
#'
#' Print function for metalm
#' @param x metalm object

print.metalm = function(x){
  print(x$model)
}

model_slopes <- function(x, rename = TRUE){
  if (methods::is(x, "meta")) {
    return(stats::coef(summary(x)))
  }

  out <- summary(x$model)
  if(rename){
    rownames(out$coefficients)[grepl("Slope", rownames(out$coefficients))] = x$parameter_names
  }
  stats::coef(out)
}

#' summary.metalm
#'
#' Summary function for metalm
#' @param object metalm object
#' @param ... additional arguments passed to model_slopes
#' @export

summary.metalm = function(object, ...) {
  out <- summary(object$model)

  out$coefficients <- model_slopes(object, ...)

  out <- list(summary = out,
              anova = object$anova)

  class(out) <- "mlm_summary"
  out

}

#' print.mlm_summary
#'
#' Print function for mlm_summary
#' @param x mlm object
#' @param ... additional arguments passed to print
#' @export

print.mlm_summary = function(x, ...){

  print(x$summary,...)

  cat("\n")
  cat(crayon::underline("ANOVA results"))
  cat("\n")
  print(x$anova, ...)
}

`%~%` = function(lhs, rhs){

  lhs <- eval(substitute(lhs))
  rhs <- deparse(substitute(rhs))
  rhs = gsub("\\~","", rhs)

  mlm(m = lhs, formula = rhs)

}

#' model.matrix2
#'
#' model.matrix but does not assign a prefix to dummy variables
#' @param object a formula
#' @param data data.frame object
#' @param contrasts.args argument passed to model.matrix
#' @param xlev argument passed to model.matrix
#' @param ... other arguments passed to model.matrix

model.matrix2 <- function(object, data = environment, contrasts.args = NULL, xlev = NULL, ...){

  ans <- stats::model.matrix(object, data, contrasts.args, xlev,...)
  factor_names <- names(attr(ans, "contrasts"))
  object_names <- trimws(unlist(strsplit(as.character(object), "\\+")))[-1]

  for(cn in seq_along(colnames(ans))){

    if(!colnames(ans)[cn] %in% object_names){
      for(fn in factor_names){
        colnames(ans)[cn] <- gsub(fn, "", colnames(ans)[cn], fixed = TRUE)
      }
    }
  }
  return(ans)
}

#' dummy_matrix
#'
#' creates a predictor matrix when cells can contain multiple tags
#' @param x a character vector
#' @param levels a vector. Provides the level (and order) of colnames.
#' @param pattern pattern provided to str_split
#' @importFrom stringr str_split
#' @importFrom dplyr %>%
#' @export

dummy_matrix = function(x, levels = NULL, pattern = ",") {
  if (!(is.factor(x))) {
    x = factor(x)
  }
  if(all(is.na(x))){
    x = as.matrix(rep(1, length(x)))
    colnames(x) = "NA"
    return(x)
  }  # return singular result

  split = x %>%
    stringr::str_split(pattern) %>% #split based on pattern
    lapply(., trimws) #remove whitespace

  contents = split %>%
    unlist() %>%
    unique() %>%
    stats::na.omit() #don't record NAs

  out = lapply(seq_along(contents), function(c) {
    lapply(seq_along(split), function(s) {
      tag = contents[c]
      current = split[s][[1]]
      out = ifelse(tag %in% current, 1, 0)

      if (all(is.na(current))) {
        out = NA
      }
      out
    }) %>% unlist
  }) %>% do.call(cbind, .)

  colnames(out) = contents

  matrix_levels = levels(droplevels(x))
  levels_length_same <- length(matrix_levels) == ncol(out)

  if (levels_length_same & all(matrix_levels %in% colnames(out))) {
    out = out[, matrix_levels, drop = FALSE] #reorder matrix if possible
  }

  if (!is.null(levels)) {
    if (!all(colnames(out) %in% levels)) {
      current_colnames =  paste0("(", paste(colnames(out), collapse = ","), ")")
      current_levels = paste0("(", paste(levels, collapse = ","), ")")
      stop(
        paste0(
          "colnames ",
          current_colnames,
          " do not match supplied levels: ",
          current_levels
        )
      )
    } else {
      levels_missing = levels[!levels %in% colnames(out)]

      if (length(levels_missing) > 0) {
        find_data_rows = lapply(seq_along(out[, 1]), function(r) {
          all(!is.na(out[r, ]))
        }) %>% unlist
        warning(
          paste0(
            "The following levels were not found in the vector: ",
            paste(levels_missing, collapse = ", "),
            ". They have been added."
          )
        )

        for (l in levels_missing) {
          missing = matrix(rep(NA, nrow(out)))
          colnames(missing) = l


          missing[, 1][find_data_rows] = 0
          out =  cbind(out, missing)
        }
      }
    }
    out = out[, as.character(levels)]
  }

  return(out)
}

#' meta3_OK
#'
#' Check whether meta3 status is OK
#' @param m meta3 model
#' @details This function looks up the OpenMx status code. If 0 or 1 returns TRUE, else returns FALSE
#' @export

meta3_OK <- function(m) {
  methods::is(m, "meta3")
  m$mx.fit$output$status$code[[1]] %in% c(0, 1)
}
