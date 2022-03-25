#' moderation_instructions
#'
#' Takes instructions, can then be provided to moderate
#'
#' @param ... moderation instructions
#' @export moderation_instructions

moderation_instructions = function(...){
  sapply(rlang::enexprs(...), deparse)
}

#' moderate
#'
#' Create a list of meta3 objects moderated by predictors
#' @param m a meta3 object
#' @param ... a named list of moderators
#' @param moderators You can feed a named character vector for consistent moderators across models
#' @param debug If true, ... contents are returned
#' @export

moderate = function(m,...,moderators = NULL, debug = FALSE){
  call = match.call()
  m$call$model.name = "Baseline"
  attr(m, "Baseline") <- TRUE
  elip = sapply(rlang::enexprs(...), deparse)
  if(length(elip) > 0){
  if(!is.null(moderators)) elip <- c(elip, moderators)
  }else{
    elip <- moderators
  }

  if(debug){
    return(elip)
  }
  names(elip)[names(elip) == ""] = gsub("\\.{1,}","_",make.names(elip[names(elip) == ""]))

  parent_envir <- parent.frame()

  models <- lapply(seq_along(elip), function(x) {
    moderated_model <-
      mlm(m, formula = elip[x], model.name = names(elip)[x], .envir = parent_envir)
    attr(moderated_model, "Baseline") <- FALSE
    moderated_model
  })

  models <- append(list(m), models)
  names(models) = c("Baseline",names(elip))

  models <- models[!is.na(models)] # Remove models without valid predictor matrices


  out <- list(
    call = call,
    models = models
  )
  class(out) <- "meta_list"

  model_info <- do.call(rbind,lapply(out$models, mlm_overview))

  #Flag to the users if models have a mx status of greater than one (which indicates issues)

  if(any(!model_info$mx_status %in% c(0,1))){

    mx_results <- model_info$mx_status
    names(mx_results) <- model_info$moderation

    warning("The following models had Mx status other than zero or one which could indicate issues with convergence: ", paste(names(mx_results)[!mx_results %in% c(0,1)], collapse = ", "))
  }

  if(any(!is.na(model_info$fixed_tau))){
    warning("At least one tau were < 0.001 and constrained to zero to assist with estimating standard errors for predictors",call. = FALSE)
  }

  out

}

#' print.meta_list
#' @param x object to print
#' @param ... additional arguments. Not currently used.
#' @export

print.meta_list = function(x, ...) {

  if(length(x$models) == 1) return(cat(crayon::red("<meta_list> with no moderators")))
  baseline_summary <- summary(x$models$Baseline)
  I2_2 = digits(baseline_summary$I2.values[1,"Estimate"]*100,1)
  I2_3 = digits(baseline_summary$I2.values[2,"Estimate"]*100,1)


  tab <- lapply(2:length(x$models), function(i){
    mlm_overview(x$models[[i]])
  })

  tab <- do.call(rbind, tab)
  tab$LRT <- NULL
  tab$R2_2 <- digits(tab$R2_2, 2)
  tab$R2_3 <- digits(tab$R2_3, 2)
  tab$p.value <- round_p(tab$p.value,3, stars = 0.05)

  tab$p.value <- sapply(tab$p.value, function(i){
    if(!grepl("\\*",i)){
      i = paste0(i," ") #if no star, add a space to keep things nicely lined up
    }
    i
  })

problem_models <- tab$moderation[!tab$mx_status %in% c(0,1)]
fixed_tau <- tab$moderation[!is.na(tab$fixed_tau)]
if(length(fixed_tau) > 0){
  fixed_tau <- glue::glue("`{fixed_tau}` [{tab$fixed_tau[!is.na(tab$fixed_tau)]}]")
  fixed_tau <- paste(fixed_tau, collapse = ", ")
}

tab$fixed_tau <- NULL

 out = utils::capture.output(tab[, !names(tab) %in% "mx_status"])

  header = out[1]
  width = max(nchar(out))
  bar = paste(rep(crayon::silver("-"), width),collapse = "")
  text = gsub("\\*", crayon::yellow("*"),out[-1])

  tab = paste(c(bar,
                header,
                bar,
                text,
                bar), collapse = "\n")

  cat(crayon::underline("Moderation results:\n"))
  cat("\n")

  cat(paste0("I2(2): ", I2_2, "%"))
  cat("\n")
  cat(paste0("I2(3): ", I2_3, "%"))
  cat("\n")

  cat(tab)
  cat("\n")

  if(length(problem_models) > 0){
    mx_message = paste0("Did not converge: " ,paste(problem_models, collapse = ", "),".") %>%
      crayon::red()
  } else {
    mx_message = crayon::cyan("All models converged.")
  }

  cat(mx_message)
  cat("\n\n")
  if(length(fixed_tau) > 0){
    tau_message <- glue::glue("Tau values for the following models were < 0.001 and needed to be fixed to zero to estimate standard errors for at least one predictor:\n")
    cli::cli_alert_info(crayon::yellow(tau_message))
    cat(fixed_tau)
  }

  # removed_moderators = names(x$removed_moderators)[x$removed_moderators]
  #
  # cat("\n\n")
  # if(length(removed_moderators) > 0){
  #   removed_moderator_message = paste0(length(removed_moderators) %>% papyr::as_word(T),
  #                                      " moderator(s) were removed due to no variance:\n",
  #                                      paste(removed_moderators, collapse = ", "),".")
  #   cli::cli_alert_warning(crayon::red(removed_moderator_message))
  # }

}

c_sentence <- function(x, and = "and"){
  last = x[length(x)]
  first = paste0(x[-length(x)], collapse = ", ")
  return(paste0(first, ", ",and," ", last))
}


digits <- function(x, n = 2) {
  x = round(x, n)
  x[] = sapply(x, function(i) {
    ifelse(!is.na(i), trimws(format(round(as.numeric(as.character(i)), n), nsmall = n)),NA)
  })
  return(x)
}

#' round_p
#'
#' Rounds p using APA formatting
#' @param p numeric
#' @param n number of digits
#' @param stars numeric vector for thresholds after which stars will be added. e.g. 0.05
#' @param leading.zero If FALSE the leading zero is removed
#' @param apa_threshold numeric. The number after which results are summarised with less than symbols
#' @param simplify numeric. Defines the threshold above which a digit is dropped
#' @export

round_p <- function(p, n = 3, stars = c(), leading.zero = T, apa_threshold = 0.001, simplify = .1){
  rounded = digits(p,n)
  out <- lapply(seq_along(rounded), function(x){

    if(!is.na(rounded[x])){
      #message(x)
      original = p[x]
      r_original = rounded[x]
      r = rounded[x]

      if(as.numeric(r) == 0){
        r = strsplit(r,split="")[[1]]
        r[length(r)] = 1
        r = paste(r,collapse = "")
      }

      #  add stars --------------
      stars_to_add = c()
      if(!is.null(stars)){
        stars_to_add = lapply(stars,function(s){
          if(as.numeric(original) < s){
            return("*")
          }else{
            return(NA)
          }
        })

        stars_to_add <- paste(stats::na.omit(unlist(stars_to_add)), collapse = "")

      }

      if(!leading.zero){
        r = sub("^(-)?0[.]",
                "\\1.", r)
      }

      if(r_original < as.numeric(r)){
        r = paste0("< ",r)
      }

      if(original < apa_threshold){
        r = paste0("< ", apa_threshold)
      }

      if(original >= simplify){
        r = digits(original, 2)
      }

      r = paste0(r,stars_to_add)

      return(r)

    }else{
      NA
    }
  })
  unlist(out)

}

mlm_overview = function(x, include_tau2 = FALSE) {
  if (methods::is(x, "metalm")) {
    model_summary = summary(x)$summary
    Pval = x$anova$p[2]
    LRT = as.character(glue::glue("$\\chi^2$({x$anova$diffdf[2]}) = {digits(x$anova$diffLL[2],2)}, $p$ = {round_p(Pval)}"))
    LRT <- gsub("= <", "<", LRT)
    r2_val <- model_summary$R2.values
    R2_2 <- r2_val["R2", "Level 2"]
    R2_3 <- r2_val["R2", "Level 3"]
    model_name = x$model$call$model.name
    fixed_tau <- attr(x$model, "fixed_tau")
    fixed_tau.original <- attr(x$model, "fixed_tau.originalvalue")

  }else{
    if(methods::is(x, "meta_list")) x <- x$models[[1]]

  model_summary = summary(x)
  Pval = NA
  R2_2 <- NA
  R2_3 <- NA
  LRT <- NA_character_
  model_name <- x$call$model.name
  fixed_tau <- attr(x, "fixed_tau")
  fixed_tau.original <- attr(x, "fixed_tau.originalvalue")
  }

  Tau2_2 <- model_summary$coefficients["Tau2_2", "Estimate"]
  Tau2_3 <- model_summary$coefficients["Tau2_3", "Estimate"]

  if(!is.null(fixed_tau)){
    if(fixed_tau == "Tau2_2"){
      Tau2_2 <- fixed_tau.original
    }
    if(fixed_tau == "Tau2_3"){
      Tau2_3 <- fixed_tau.original
    }

  }else{
    fixed_tau = NA
  }

  k <- model_summary$no.studies
  n <- model_summary$obsStat

  if(is.null(model_name)) model_name <- "Unnamed model"

  out <- data.frame(
    moderation = model_name,
    k = k,
    n = n,
    Tau2_2 = Tau2_2,
    Tau2_3 = Tau2_3,
    R2_2 = R2_2,
    R2_3 = R2_3,
    p.value = Pval,
    mx_status = model_summary$Mx.status1,
    fixed_tau = fixed_tau,
    LRT = LRT
  )

  if(!include_tau2){
    out$Tau2_2 <- NULL
    out$Tau2_3 <- NULL
  }

  out

}

#' summary.meta_list
#' @param object meta_list object
#' @param ... additional arguments passed to format_nicely
#' @export

summary.meta_list = function(object, ...){

  tab <- format_nicely(object, slope_p = TRUE, ...)

  tab$Moderation[attr(tab, "indent")] <- paste0("--",tab$Moderation[attr(tab, "indent")])
  names(tab) <-gsub("[$^{}_]","", names(tab))
  tab$LRT_p <- gsub(".*p\\$\\s=?\\s?", "", tab$`Likelihood Ratio Test`)
  tab$`Likelihood Ratio Test` <- NULL

  capture_simply = function(cols, right) utils::capture.output(print.data.frame(tab[,cols, drop = FALSE, with = FALSE], right = right, row.names = FALSE))
  col1 <- trimws(capture_simply(1, FALSE), which = "left")
  others <- capture_simply(-1, TRUE)
  out <- paste0(col1,others)

  out <- gsub("--","  ",out)
  width <- nchar(out[1])
  bar <- paste(rep(crayon::silver("-"), width),collapse = "")
  out <- gsub("\\-(?![0-9])",crayon::silver("-"),out, perl = TRUE)
  vred <- crayon::make_style(grDevices::rgb(1,.2,.2))
  out <- gsub("\\-(?=[0-9])",vred("-"),out, perl = TRUE)

  cat(crayon::blue("Moderation results"))
  cat("\n")
  cat(bar)
  cat("\n")
  cat(out[1])
  cat("\n")
  cat(bar)
  cat("\n")
  cat(paste(out[-1], collapse = "\n"))
  cat("\n")
  cat(bar)
  cat("\n")
}

#' coef.meta_list
#'
#' @param object meta_list object
#' @param ... unused
#' @export

coef.meta_list <- function(object, ...){

  crit_val <- stats::qnorm(.975)

 effects <-
    data.table::data.table(object$models$Baseline$data)[, .(
      cluster,
      est = y,
      SE = sqrt(v),
      lower = y - crit_val * sqrt(v),
      upper = y + crit_val * sqrt(v),
      setting = "Effect sizes",
      type = "Effect sizes"
    )]

  moderators <- lapply(object$models[-1], function(m) {
    info <- moderator_info(m)
    info$mod <- attr(info, "model.name")
    out <- info[moderator_level == TRUE, .(
      cluster = moderation,
      moderation = mod,
      est = Estimate,
      SE = SE,
      lower = lbound,
      upper = ubound,
      type = "moderator level",
      setting = "Pooled"
    )]
    out[, model_p := info$p.value[[1]]]
  })

  moderators <- do.call(rbind, moderators)

  baseline <-
    data.table::data.table(cbind(
      mlm_overview(object$models$Baseline),
      just_estimates(object$models$Baseline)
    ))[, .(
      cluster = "Baseline",
      moderation = "Baseline",
      est = Estimate,
      SE,
      lower = lbound,
      upper = ubound,
      type = "Baseline",
      setting = "Pooled"
    ),]

data.table::rbindlist(list(baseline, moderators, effects), fill = TRUE)

}

#' print.KIN_summary
#' @param x meta_list object
#' @param ... additional arguments passed to format_nicely
#' @export

print.KIN_summary = function(x, ...){
  tab <- x
  tab$Predictor[attr(tab, "indent")] <- paste0("--",tab$Predictor[attr(tab, "indent")])
  names(tab) <-gsub("[$^{}_]","", names(tab))

  capture_simply = function(cols, right) utils::capture.output(print.data.frame(tab[,cols, drop = FALSE], right = right, row.names = FALSE))

  col1 <- trimws(capture_simply(1, FALSE), which = "left")
  others <- capture_simply(-1, TRUE)
  out <- paste0(col1,others)

  out <- gsub("--","  ",out)
  width <- nchar(out[1])
  bar <- paste(rep(crayon::silver("-"), width),collapse = "")
  out <- gsub("\\-(?![0-9])",crayon::silver("-"),out, perl = TRUE)
  vred <- crayon::make_style(grDevices::rgb(1,.2,.2))
  out <- gsub("\\-(?=[0-9])",vred("-"),out, perl = TRUE)

  cat(crayon::blue(attr(x,"title")))
  cat("\n")
  cat(bar)
  cat("\n")
  cat(out[1])
  cat("\n")
  cat(bar)
  cat("\n")
  cat(paste(out[-1], collapse = "\n"))
  cat("\n")
  cat(bar)
  cat("\n")
}

