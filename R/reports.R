#reports

report_q = function(x,
                    rmarkdown = FALSE,
                    round = 2) {
  call = match.call()
  envir = sys.parent()
  if (!methods::is(x, "name"))
    x <- call$x
  stat_q = rmarkdown_wrap(
    glue::glue('get_val({x}, "q", round = {round})'),
    rmarkdown = rmarkdown,
    envir = envir
  )
  stat_df = rmarkdown_wrap(
    glue::glue('get_val({x}, "q_df")'),
    rmarkdown = rmarkdown,
    envir = envir
  )
  stat_p = rmarkdown_wrap(
    glue::glue('round_p(get_val({x}, "q_p"),3)'),
    rmarkdown = rmarkdown,
    envir = envir
  )

  if (stat_p < 0.05) {
    mess = glue::glue(
      "Inspecting the Q statistic revealed significant heterogeneity $Q$({stat_df}) = {stat_q}, $p$ {stat_p}."
    )
  } else{
    mess = glue::glue("Evidence for heterogeneity was not found $Q$({stat_df}) = {stat_q}, $p$ = {stat_p}.")
  }

  return(mess)

}

#' report_n
#'
#' Get number of studies and effect sizes of a baseline model
#' @param x meta3L model
#' @param rmarkdown if TRUE, rmarkdown code is produced which can be edited

report_n = function(x, rmarkdown = FALSE){
  if(methods::is(x, "meta_list")) x <- x$models[[1]]

  call = match.call()
  envir = sys.parent()

  if(!methods::is(x,"name")) x <- call$x


  stat_k = rmarkdown_wrap(glue::glue('metaKIN:::as_word(get_val({x},"k"), TRUE)'), rmarkdown = rmarkdown, envir = envir)
  stat_n = rmarkdown_wrap(glue::glue('metaKIN:::as_word(get_val({x},"n"), FALSE)'), rmarkdown = rmarkdown, envir = envir)

  mess = glue::glue("{stat_k} studies (including {stat_n} effect sizes) reported data which could be pooled.")

  mess

}

report_baseline = function(x, rmarkdown = FALSE, round = 2, transf = function(x) x){
  call = match.call()
  envir = sys.parent()


  if(!rmarkdown){
    stat_est <- get_val(x, "estimate95", round = round, transf = transf)
  }else{
    if(is.null(transf)) transf <- as.name("function(x) x")
    if (!methods::is(x, "name"))
      x <- call$x
    if (!methods::is(transf, "name"))
      transf <- call$transf
    stat_est = rmarkdown_wrap(
      glue::glue('get_val({x}, "estimate95", round = {round}, transf = {transf})'),
      rmarkdown = rmarkdown,
      envir = envir)
  }


  mess = glue::glue("The pooled effect size was {stat_est}.")
  mess
}

report_i2 = function(x, rmarkdown = FALSE){
  call = match.call()
  envir = sys.parent()
  if(!methods::is(x, "name")) x <- call$x
  stat_i2_2 = rmarkdown_wrap(glue::glue('get_val({x}, "I2_2%")'), rmarkdown = rmarkdown, envir = envir)
  stat_i2_3 = rmarkdown_wrap(glue::glue('get_val({x}, "I2_3%")'), rmarkdown = rmarkdown, envir = envir)
  mess = glue::glue("The heterogeneity at level 2 was {stat_i2_2}. The heterogeneity at level 3 was {stat_i2_3}.")
  return(mess)
}


report_moderators = function(x, rmarkdown = FALSE, digits = 2){
  call = match.call()
  #return(call)
  mods = get_sig_moderators(x)
  if(!methods::is(x, "name")) x <- call$x
  envir = sys.parent()

  if(length(mods) == 0){
    return("No covariate was found to be a significant moderator of the baseline model.")
    }


  r2_2.md = rmarkdown_wrap(glue::glue("get_val({x}, 'R2_2%', '{mods}')"),rmarkdown, envir=envir)
  r2_3.md = rmarkdown_wrap(glue::glue("get_val({x}, 'R2_3%', '{mods}')"),rmarkdown, envir=envir)

  if (length(mods) == 1) {
    sent_mods = glue::glue(
      "'<<mods>>' ($R^2_{(2)}$ = <<r2_2.md>>; $R^2_{(3)}$ = <<r2_3.md>>)",
      .open = "<<",
      .close = ">>"
    )
  } else{
    sent_mods = c_sentence(
      glue::glue(
        "'<<mods>>' ($R^2_{(2)}$ = <<r2_2.md>>; $R^2_{(3)}$ = <<r2_3.md>>)",
        .open = "<<",
        .close = ">>"
      )
    )
  }

  if(length(mods) ==1){
    mess = glue::glue("The covariate which significantly moderated the baseline model was {sent_mods}.")
  }


  if(length(mods) > 1){
  mess = glue::glue("The covariates which significantly moderated the baseline model were {sent_mods}.")
  }
return(mess)

}

#' report
#'
#' Constructs written reports about the contents of models
#' @param meta_list the meta_list object
#' @param ... things to report. One o
#' @param rmarkdown return results in rmarkdown?
#' @param digits the number of digits to return
#' @param transf you can supply a function to transform baseline pooled estimates
#' @export report

report = function(meta_list,..., rmarkdown = FALSE, digits = 2, transf = function(x) x){
call = match.call()



elip = unname(sapply(rlang::enexprs(...), deparse))

options <- c("n","q","baseline","i2","moderators")

if(length(elip) == 0) elip <- options

options <- options[options %in% elip]

mess = list(
  n = report_n(meta_list, rmarkdown = rmarkdown),
  q = report_q(meta_list, rmarkdown = rmarkdown, round = 2),
  baseline = report_baseline(
    meta_list,
    rmarkdown = rmarkdown,
    round = 2,
    transf = transf
  ),
  i2 = report_i2(meta_list, rmarkdown = rmarkdown),
  moderators = report_moderators(meta_list, rmarkdown = rmarkdown, digits = 2)
)

mess = paste(mess[options], collapse = " ")
mess = gsub("call\\$transf","NULL", mess)

if(rmarkdown){
  mess <- gsub("meta_list", call$meta_list, mess)
  mess <-
    gsub("transf = transf",
         glue::glue("transf = {deparse(call$transf)}"),
         mess)
  clipr::write_clip(mess)
  return(cat(crayon::blue("<script sent to clipboard>")))
}else{
  return(mess)
}

}



rmarkdown_wrap = function(code, rmarkdown = FALSE, envir = sys.parent()){
call = match.call()

  if(!rmarkdown){

    if(length(code) >1){

      return(unlist(lapply(seq_along(code), function(i) eval(parse(text = code[[i]]), envir = envir))))
    }else{
    return(eval(parse(text = code), envir = envir))
    }

  }else{
    return(glue::glue('`r {code}`'))

  }

}


#' as_word
#'
#' Convert numbers to words
#' @param x integer to convert to a word
#' @param sentence if TRUE capitalise the first letter
#' @param hyphenate should hyphens be included?

as_word <- function(x = NULL,
                    sentence = F,
                    hyphenate = T) {
  x = as.character(english::english(as.numeric
                                    (x)))
  if (sentence == T) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  }

  if (hyphenate) {
    compounds = paste("y",english::english(1:9))
    compounds_hyphen = paste("y",english::english(1:9), sep = "-")

    for (n in seq_along(compounds)) {
      x <- gsub(compounds[n], compounds_hyphen[n], x)
    }
  }

  x
}

default_note = function(){
  "k = number of studies; n = numbers of effect sizes; Estimate = population average; SE = standard error. $R^2^~(2)~$ = the proportion of within-cluster heterogeneity explained by the covariate; $R^2^~(3)~$ = the proportion of between-cluster heterogeneity explained by the covariate; p-value = ANOVA p-value; * indicates p < 0.05"
}


#describing models

#' describe_baseline
#' @param obj a meta_ninja
#' @importFrom dplyr %>%

describe_q = function(obj){
  meta_ninja = get(obj)
  q_info = meta_ninja$models$Baseline %>%
    summary %>%
    .$Q.stat
  if(q_info$pval < 0.05){
    starting_message = "Inspecting the Q statistic revealed significant heterogeneity"
  } else{
    starting_message = "Inspecting the Q statistic did not reveal significant heterogeneity"
  }
  q = paste0("`r summary(",obj,"$models$Baseline)$Q.stat$Q %>% papyr::digits(2)`")
  df = paste0("`r summary(",obj,"$models$Baseline)$Q.stat$Q.df`")
  p = paste0("`r summary(",obj,"$models$Baseline)$Q.stat$pval %>% papyr::round_p(2)`")
  stats_text = paste0(" (Q(df = ",df,") = ",q, ", *p* = ", p,").")
  paste0(starting_message, stats_text)
}

get_sig_moderators = function(x, p = 0.05){
  if(methods::is(x,"name")) x <- eval(x)
  if(!methods::is(x, "meta_list")) stop("x must be a meta_list")
  tab <- do.call(rbind, lapply(x$models[2:length(x$models)], moderator_info))
  unname(unlist(tab[moderator_level == FALSE & p.value < p, "moderation"]))
}

#' get_val
#'
#' Return values
#' @param model the meta3L or meta_list model to extract
#' @param val character vector of values to extract
#' @param moderator character string of moderator to extract from. Defaults to NULL
#' @param round 2 number of digits
#' @param transf function by which to transform relevant estimates (only applies to estimate and estimate95)
#' @export

get_val <- function(model, val, moderator = NULL, round = 2, transf = function(x) x){
  if(is.null(moderator)){
    if(methods::is(model,"meta_list")){
      model <- model$models[[1]]
    }
    return(get_val_baseline(model, val,round, transf))
  }

  model <- model$models[[moderator]]$model
  if(is.null(model)) stop("Moderator not found, check spelling")
  summ <- summary(model)

  retrieve_tau2_2 <- function(x){
    if(!is.null(attr(x, "fixed_tau"))){
      if(attr(x, "fixed_tau") == "Tau2_2"){
        return(attr(x, "fixed_tau.originalvalue"))
      }
    }
    summ$coefficients["Tau2_2", "Estimate"]
  }

   retrieve_tau2_3 <- function(x){
    if(!is.null(attr(x, "fixed_tau"))){
      if(attr(x, "fixed_tau") == "Tau2_3"){
        return(attr(x, "fixed_tau.originalvalue"))
      }
    }
    summ$coefficients["Tau2_3", "Estimate"]
  }

  inst <- list(
    "k" = "summ$no.studies",
    "n" = "summ$obsStat",
    "q" = "digits(summ$Q.stat$Q,round)",
    "q_p" = "summ$Q.stat$pval",
    "q_df" = "summ$Q.stat$Q.df",
    "R2_2" = 'digits(summ$R2.values["R2","Level 2"],round)',
    "R2_2%" = 'paste0(digits(summ$R2.values["R2","Level 2"]*100,round),"%")',
    "R2_3" = 'digits(summ$R2.values["R2","Level 3"], round)',
    "R2_3%" = 'paste0(digits(summ$R2.values["R2","Level 3"]*100, round),"%")',
    "Tau2_2" = 'digits(retrieve_tau2_2(model), round)',
    "Tau2_3" = 'digits(retrieve_tau2_3(model), round)'
  )

  res <- lapply(inst[names(inst) %in% val], function(x) eval(parse(text = x)))
  if(length(res) == 0) stop("value not found, check spelling")
  if(length(res) == 1) res <- unlist(res)
  res

}

get_val_baseline <- function(model, val, round, transf){

  summ <- summary(model)

  estimate95 <- "with(
    summ,
    glue::glue(
      \"{digits(transf(coefficients['Intercept', 'Estimate']), round)} [95% CI {digits(transf(coefficients['Intercept', 'lbound']), round)}, {digits(transf(coefficients['Intercept', 'ubound']), round)}]\"
    )
  )"

  inst <- list(
    "estimate" = 'digits(transf(summ$coefficients["Intercept", "Estimate"]),round)',
    "estimate95" = estimate95,
    "k" = "summ$no.studies",
    "n" = "summ$obsStat",
    "q" = "digits(summ$Q.stat$Q,round)",
    "q_p" = "summ$Q.stat$pval",
    "q_df" = "summ$Q.stat$Q.df",
    "I2_2" = 'digits(summ$I2.values["I2_2 (Typical v: Q statistic)", "Estimate"],round)',
    "I2_3" = 'digits(summ$I2.values["I2_3 (Typical v: Q statistic)", "Estimate"],round)',
    "I2_2%" = 'paste0(digits(summ$I2.values["I2_2 (Typical v: Q statistic)", "Estimate"]*100,round),"%")',
    "I2_3%" = 'paste0(digits(summ$I2.values["I2_3 (Typical v: Q statistic)", "Estimate"]*100,round),"%")',
    "Tau2_2" = 'digits(summ$coefficients["Tau2_2", "Estimate"], round)',
    "Tau2_3" = 'digits(summ$coefficients["Tau2_3", "Estimate"], round)'
  )

  res <- lapply(inst[names(inst) %in% val], function(x) eval(parse(text = x)))
  if(length(res) == 0) stop("value not found, check spelling")
  if(length(res) == 1) res <- unlist(res)
  res


}
