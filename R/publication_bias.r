#' FPP
#'
#' Calculates a FAT-PET-PEESEE adjusted estimate
#' @param m a metaSEM or meta_list model
#' @param transf function with which to transform results
#' @param round number of digits for rounding
#' @param alpha the alpha threshold used in PET-PEESE. Defaults to 0.1 \(10\% as per Stanley 2017\)
#' @details PET-PEESE uses meta-regression in order to adjust estimates for publication bias. When there is publication bias the sampling variance becomes correlated with effect size. By using the standard error (PET) or the sampling variance (PEESE) as moderators in a meta-regression, estimates (the intercept) can be made which partial out the correlation between effect size and variance. PET-PEESE first tests whether the intercept is a significant when controlling for the standard error of effects (p < 0.05) predictor of effect size. If it is, PEESE is used. Otherwise, PET is used.

FPP = function(m, transf = function(x) x, round = 2, alpha = .05){
  if(methods::is(m, "meta_list")){
  m <- m$models$Baseline
  }

 if(!methods::is(m, "meta3")){
   stop("This function will only work with meta3 or metaKIN objects as inputs")
 }

  dat <- m$data

  call = m$call
  call$y <- as.name("y")
  call$v <- as.name("v")
  call$cluster = as.name("cluster")
  call$data = as.name("dat")


  pet_call = call
  pet_call$model.name = "PET"
  pet_call$x = str2lang(glue::glue("sqrt({call$v})"))

  peese_call = call
  peese_call$x = str2lang(glue::glue("{call$v}"))

  # PET
  pet_m <- eval(pet_call)

  pet_p <- summary(pet_m)$coefficients["Intercept","Pr(>|z|)"]

    if (pet_p < alpha) {
    peese_m <- eval(peese_call)
      adj_m <- peese_m
      meth <- "PEESE"
    }else{
      adj_m <- pet_m
      meth <- "PET"
    }

  # Publication bias test (FAT)

  pub_bias <- stats::anova(pet_m, m)

  if(is.null(transf)){
    transf = function(x) x
  }

  out <-
    list(
      anova = pub_bias,
      diffLL = glue::glue(
        "$\\chi^2$({pub_bias$diffdf[2]}) = {digits(pub_bias$diffLL[2])}, $p$ = {round_p(pub_bias$p[2])}"
      ),
      adj_est = unname(get_val(adj_m, "estimate95", transf = transf)),
      adj_method = meth,
      adj_model = adj_m
    )
  class(out) = "fpp"
  out

}

#' aggregate_to_cluster
#'
#' Aggregate effect sizes to cluster
#' @param model a meta3 model
#' @details multiple effect sizes per cluster are aggregated using a fixed effects meta-analysis

aggregate_to_cluster <- function(model){
if(methods::is(model, "meta_list")) model <- model$models[[1]]

  model_data <- data.table(model$data)

  model_data[, {

    if(nrow(.SD) == 1){
      out <- list(y = .SD$y, v = .SD$v)
    }else{
    temp_mod <- metafor::rma(y, v, data = .SD, method = "FE")
    out <- list(y = as.numeric(temp_mod$b), v = as.numeric(temp_mod$se^2))
    }
    out
  }
       , by = cluster]

}

#' trim_and_fill
#'
#' Perform trim and fill for meta3 objects.
#' @param model a model of class meta3
#' @param aggregate should data be aggregated first?
#' @param ... additional arguments passed to metafor::rma.uni
#' @details uses data aggregated by cluster.
#' @return rma.uni.trimfill

trim_and_fill <- function(model, aggregate = TRUE, ...){
  if(methods::is(model, "meta_list")) model <- model$models[[1]]

  if(aggregate){
  data <- aggregate_to_cluster(model)
  }else{
    data <- model$data
    warning("If there is clustering within the data the assumptions of trim and fill are violated.")
  }

  rma_model <- metafor::rma.uni(y, v, data = data, ...)

  metafor::trimfill(rma_model)
}

