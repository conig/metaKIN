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
                             moderator_level = mod_lvl,
                             p_value = `Pr(>|z|)`)]
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
#' @param p_digits a scalar. The number of digits to round p to.
#' @param ci_sep separator for confidence intervals
#' @param include_i2 A bool, should i2 be included next to baseline?
#' @param stars should significance stars be included for factor levels?
#' @param slope_p if TRUE slope p-values are included
#' @param replace a vector with names included. gsub will be applied to the moderation column such that the vector's names are replaced with the vector's contents
#' @import data.table
#' @export

format_nicely = function(meta_list,
                         effect_name = "Estimate",
                         transf = NULL,
                         transf_name = NULL,
                         slope_p = TRUE,
                         hide_insig = FALSE,
                         round = 2,
                         p_digits = 3,
                         ci_sep = ", ",
                         include_i2 = FALSE,
                         stars = FALSE,
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

  # replace SE with "-" when NA and moderator level
  tab$SE[is.na(tab$SE) & tab$moderator_level] <- "-"

  tab$lbound = digits(transf(tab$lbound), round)
  tab$ubound = digits(transf(tab$ubound), round)

  tab$Estimate_formatted = glue::glue_data(tab, "{Transformed_estimate} [{lbound}{ci_sep}{ubound}]")
  nullreplace = glue::glue_data(tab, "{NA} [{NA}{ci_sep}{NA}]")
  tab$Estimate_formatted[tab$Estimate_formatted == nullreplace] <-
    NA
  tab$Estimate_formatted <- gsub(glue::glue(" \\[NA{ci_sep}NA\\]"), "", tab$Estimate_formatted)

  if(stars){
    tab$Estimate_formatted[which(tab$p_value < 0.05)] <- paste0(tab$Estimate_formatted[which(tab$p_value < 0.05)], "*")
  }

  indent <- tab$moderator_level

  tab <- tab[, .(
    Moderation = moderation,
    k = digits(k, 0),
    n = digits(n, 0),
    Estimate_formatted,
    Estimate,
    SE,
    `$p$` = round_p(p_value, p_digits),
    "$R^2_{(2)}$" = digits(R2_2 * 100, round),
    "$R^2_{(3)}$" = digits(R2_3 * 100, round),
    `Likelihood Ratio Test` = LRT
  )]

  if (!is_transf) {
    tab$Estimate <- NULL
    if (is.null(transf_name)) {
      transf_name <- paste0(effect_name, " [95\\% CI]")
    }
  }

  if (is.null(transf_name)) {
    transf_name <- as.character(call$transf)
    transf_name <- paste0(transf_name, " [95\\% CI]")

  }

  if (length(replace) > 0 & !identical(replace, FALSE)) {
    for (i in seq_along(replace)) {
      tab$Moderation <-
        gsub(names(replace)[i], replace[i], tab$Moderation)
    }
  }

  if(!slope_p) tab$`$p$` <- NULL
  data.table::setnames(tab, "Estimate", effect_name, skip_absent = TRUE)
  data.table::setnames(tab, "Estimate_formatted", transf_name)
  if (all(is.na(tab$"$R^2_{(3)}$")))
    tab$"$R^2_{(3)}$" <- NULL
  tab[is.na(tab)] <- ""
  attr(tab, "indent") <- indent
  tab
}
