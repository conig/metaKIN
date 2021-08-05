#' PET_PEESE
#'
#' Calculates a PET-PEESEE adjusted estimate
#' @param m a metaSEM or meta_list model
#' @param transf function with which to transform results
#' @param round number of digits for rounding
#' @param alpha the alpha threshold used in PET-PEESE. Defaults to 0.1 \(10\% as per Stanley 2017\)
#' @param type Should PET-PEESE, PET or PEESE be used?
#' @details PET-PEESE uses meta-regression in order to adjust estimates for publication bias. When there is publication bias the sampling variance becomes correlated with effect size. By using the standard error (PET) or the sampling variance (PEESE) as moderators in a meta-regression, estimates (the intercept) can be made which partial out the correlation between effect size and variance. PET-PEESE first tests whether the variance component is a significant (p < 0.05) predictor of effect size. If it is, PEESE is used. Otherwise, PET is used.

PET_PEESE = function(m, transf = function(x) x, round = 2, alpha = .1, type = c("PET-PEESE","PET","PEESE")){
  if(methods::is(m, "meta_list")){
  m <- m$models$Baseline
  }


 if(!methods::is(m, "meta3")){
   stop("This function will only work with meta3 or meta_ninja objects as inputs")
 }

  type = type[1]
  types = c("PET-PEESE","PET","PEESE")

  if(!type %in% c("PET-PEESE","PET","PEESE")){
    stop("type must be one of: '", paste(types, collapse = "', '"), "'.")
  }
  call = m$call
  pet_call = call
  pet_call$x = str2lang(glue::glue("sqrt({call$v})"))

  peese_call = call
  peese_call$x = str2lang(glue::glue("{call$v}"))

  if(type == "PET"){
    adj_m <- eval(pet_call)
  }
  if(type == "PEESE"){
    adj_m <- eval(peese_call)
  }

  if(type == "PET-PEESE") {
    adj_m <- eval(pet_call)
    pet_p <- summary(adj_m)$coefficients["Slope_1","Pr(>|z|)"]
    if (pet_p < alpha) {
      adj_m <- eval(peese_call)
    }
  }

  if(is.null(transf)){
    transf = function(x) x
  }

  unname(get_val(adj_m, "estimate95", transf = transf))

}


#' threePSM
#'
#' Perform a weight-function model to test whether the pvalue of estimates affects population estimates
#' @param model the meta3 model
#' @param ... arguments passed to weight::weightfunct
#' @param transf function to transform adjusted estimate

threePSM = function(model, ..., transf = NULL){

  if(methods::is(model, "meta_list")){
    model = model$models[[1]]
  }

  result <- weightr::weightfunct(
    effect = model$data$y,
    v = model$data$v,
    ...
  )

  df <- length(result[[2]]$par) - length(result[[1]]$par)
  lrchisq <- 2 * (abs(result[[1]]$value - result[[2]]$value))
  pvalue <- 1 - stats::pchisq(lrchisq, df)
  adj_t = threePSM_table(result)

  if(is.null(transf)) transf <- function(x) x
  est = digits(transf(adj_t$estimate[1]),2)
  lower = digits(transf(adj_t$ci.lb[1]), 2)
  upper = digits(transf(adj_t$ci.ub[1]), 2)

  adjusted_result = glue::glue("{est} [95% CI {lower}, {upper}]")

  apa = glue::glue("$\\chi^2$({df}) = {digits(lrchisq,2)}, $p$ = {round_p(pvalue)}")

  list(chisq = lrchisq,
       df = df,
       pvalue = pvalue,
       adjusted_table = adj_t,
       adjusted_result = adjusted_result,
       apa = apa,
       raw = result)
}

# This function is the modified print function from weightr.
# It is used internally in threePSM

threePSM_table = function(x){
  if (x$fe == FALSE) {
    if (is.null(x$weights)) {
      adj_int_est <- cbind(x$adj_est[2:((x$nsteps - 1) +
                                          (x$npred + 2))])
      adj_int_se <- cbind(x$adj_se[2:((x$nsteps - 1) +
                                        (x$npred + 2))])
    }
    else {
      adj_int_est <- cbind(c(round(x$adj_est[2:((x$npred +
                                                   2))], digits = 4), sprintf("%.4f", x$weights[2:length(x$weights)])))
      adj_int_se <- cbind(rep("---", length(x[[2]]$par[2:length(x[[2]]$par)])))
    }
  }
  if (x$fe == TRUE) {
    if (is.null(x$weights)) {
      adj_int_est <- cbind(x$adj_est[1:((x$nsteps - 1) +
                                          (x$npred + 1))])
      adj_int_se <- cbind(x$adj_se[1:((x$nsteps - 1) +
                                        (x$npred + 1))])
    }
    else {
      adj_int_est <- cbind(c(round(x$adj_est[1:((x$npred +
                                                   1))], digits = 4), sprintf("%.4f", x$weights[2:length(x$weights)])))
      adj_int_se <- cbind(rep("---", length(x[[2]]$par[1:length(x[[2]]$par)])))
    }
  }
  if (is.null(x$weights)) {
    z_stat_int <- adj_int_est/adj_int_se
    p_val_int <- (2 * stats::pnorm(-abs(z_stat_int)))
    ci.lb_int <- adj_int_est - stats::qnorm(0.975) * adj_int_se
    ci.ub_int <- adj_int_est + stats::qnorm(0.975) * adj_int_se
  }
  else {
    if (x$fe == FALSE) {
      length_a <- length(x[[2]]$par[2:length(x[[2]]$par)])
      z_stat_int <- rep("---", length_a)
      p_val_int <- rep("---", length_a)
      ci.lb_int <- rep("---", length_a)
      ci.ub_int <- rep("---", length_a)
    }
    else {
      length_aF <- length(x[[2]]$par[1:length(x[[2]]$par)])
      z_stat_int <- rep("---", length_aF)
      p_val_int <- rep("---", length_aF)
      ci.lb_int <- rep("---", length_aF)
      ci.ub_int <- rep("---", length_aF)
    }
  }
  res.table <- data.frame(matrix(c(adj_int_est, adj_int_se,
                                   z_stat_int, p_val_int, ci.lb_int, ci.ub_int), nrow = (x$npred +
                                                                                           1 + (x$nsteps - 1)), byrow = F), stringsAsFactors = FALSE)
  rowlabels1 <- rep(0, (x$npred + 1))
  rowlabels1[1] <- "Intercept"
  if (x$npred > 0) {
    for (i in 2:length(rowlabels1)) {
      rowlabels1[i] <- paste(c(colnames(x$XX)[i]))
    }
  }
  rowlabels2 <- rep(0, (x$nsteps - 1))
  for (i in 1:(length(rowlabels2))) {
    rowlabels2[i] <- paste(c(x$steps[i], "< p <", x$steps[i +
                                                            1]), collapse = " ")
  }
  row.names(res.table) <- c(rowlabels1, rowlabels2)
  colnames(res.table) <- c("estimate", "std.error",
                           "z-stat", "p-val", "ci.lb", "ci.ub")
  if (is.null(x$weights)) {
    res.table[, "p-val"] <- format.pval(res.table[,
                                                  "p-val"])
  }
  res.table
}

