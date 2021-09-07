#' just_estimates
#'
#' Obtain estimates for a model
#' @param x the meta object

just_estimates <- function(x){
  mod_lvl <- !attr(x, "Baseline")
  data.table::data.table(model_slopes(x),
                         keep.rownames = "moderation")[
                           !moderation %in% c("Tau2_2", "Tau2_3"),
                           .(moderation, Estimate, SE = Std.Error, lbound, ubound,
                             moderator_level = mod_lvl)]
}

#' moderator_info
#'
#' moderator_info
#' @param x meta_list model

moderator_info <- function(x){
  overview <- data.table::data.table(mlm_overview(x))[,moderator_level := FALSE]
  slopes <- just_estimates(x)
  slopes$moderator_level = TRUE
  slopes <- merge(slopes, x$k_n, by = "moderation", all.x = TRUE, sort = FALSE)
  out <- data.table::rbindlist(list(overview, slopes), fill = TRUE)
  out$model <- overview$moderation
  attr(out, "model.name") <- x$model$call$model.name
  out

}

#' format_nicely
#'
#' Formats a list of moderated meta3 objects into a data.frame which can be used as a table in a LaTex document
#' @param meta_list a meta_list object
#' @param round a scalar.
#' @param transf a function. If provided will transform effects and confidence intervals.
#' @param effect_name a string. If provided, will rename Estimate column with string provided.
#' @param transf_name a character string. If provided, will name the transformed column.
#' @param hide_insig a bool.
#' @param escape_pc a bool. If TRUE, \% symbols will be escaped in header, captions and notes.
#' @param p_digits a scalar. The number of digits to round p to.
#' @param leading_zero a bool. If TRUE, p-values will have leading zeros
#' @param ci_sep separator for confidence intervals
#' @param include_i2 A bool, should i2 be included next to baseline?
#' @param replace a vector with names included. gsub will be applied to the moderation column such that the vector's names are replaced with the vector's contents
#' @import data.table
#' @export

format_nicely = function(meta_list,
                         round = 2,
                         transf = NULL,
                         effect_name = "Estimate",
                         transf_name = NULL,
                         hide_insig = TRUE,
                         escape_pc = FALSE,
                         p_digits = 3,
                         leading_zero = FALSE,
                         ci_sep = ", ",
                         include_i2 = FALSE,
                         replace = c("_" = " ")) {
  call <- match.call()

  if(is.null(transf)){
    transf <- function(x) x
    is_transf = FALSE
  } else{
    is_transf = TRUE
  }

  baseline <-
    cbind(mlm_overview(meta_list$models[[1]]),
          just_estimates(meta_list$models[[1]]))

  moderators <-
    lapply(meta_list$models[2:length(meta_list$models)], moderator_info)

  tab <-
    data.table::rbindlist(append(list(baseline), moderators), fill = TRUE)

  if (hide_insig) {
    sig_mods <- get_sig_moderators(meta_list)
    tab <-
      tab[(tab$model %in% c("Baseline", sig_mods) |
             moderator_level == FALSE),]
  }
  tab$Transformed_estimate = digits(transf(tab$Estimate), round)
  tab$Estimate = digits(tab$Estimate, round)
  tab$SE = digits(tab$SE, round)
  tab$lbound = digits(transf(tab$lbound), round)
  tab$ubound = digits(transf(tab$ubound), round)

  tab$Estimate_formatted = glue::glue_data(tab, "{Transformed_estimate} [{lbound}{ci_sep}{ubound}]")
  nullreplace = glue::glue_data(tab, "{NA} [{NA}{ci_sep}{NA}]")
  tab$Estimate_formatted[tab$Estimate_formatted == nullreplace] <-
    NA

  indent <- tab$moderator_level

  tab <- tab[, .(
    Moderation = moderation,
    k = digits(k, 0),
    n = digits(n, 0),
    Estimate_formatted,
    Estimate,
    SE,
    "$R^2_{(2)}$" = digits(R2_2, round),
    "$R^2_{(3)}$" = digits(R2_3, round),
    `$p$` = round_p(p.value, p_digits)
  )]

  if (!is_transf) {
    tab$Estimate <- NULL
    if (is.null(transf_name)) {
      transf_name <- paste0(effect_name, " [95% CI]")
    }
  }

  if (is.null(transf_name)) {
    transf_name <- as.character(call$transf)
    transf_name <- paste0(transf_name, " [95% CI]")

  }

  if (length(replace) > 0 & !identical(replace, FALSE)) {
    for (i in seq_along(replace)) {
      tab$Moderation <-
        gsub(names(replace)[i], replace[i], tab$Moderation)
    }
  }

  data.table::setnames(tab, "Estimate", effect_name, skip_absent = TRUE)
  data.table::setnames(tab, "Estimate_formatted", transf_name)
  if (all(is.na(tab$"$R^2_{(3)}$")))
    tab$"$R^2_{(3)}$" <- NULL
  tab[is.na(tab)] <- "-"
  attr(tab, "indent") <- indent
  tab
}
