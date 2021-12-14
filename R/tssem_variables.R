#' cormat_list
#'
#' Generate a list of correlation matrices for tssem
#' @param yi correlation effect
#' @param vi correlation variance
#' @param var1 the name of the column containing the first variable
#' @param var2 the name of the column containing the second variable
#' @param cluster the name of the clustering variable
#' @param data data.frame
#' @param transf a function to transform individual correlations
#' @param ni the name of the column containing the sample sizes
#' @param ... optional moderators e.g. moderator = max(variable, na.rm = TRUE). Variables will be created by data.table by cluster.
#' @import data.table
#' @details N is calculated as the average sample size per correlation

cormat_list <- function(yi, vi, ni, var1, var2, cluster, data, transf = NULL, ...){

  data = data.table::data.table(data)
  inspect_vars = data[,c(yi,vi,var1,var2,cluster), with = FALSE]
  missing = unlist(lapply(1:nrow(inspect_vars), function(x) any(is.na(inspect_vars[x,]))))
  data = data[!missing, ]
  vars <- unique(unlist(c(data[,..var1], data[,..var2])))

  make_matrix = function(vars) {
    m = matrix(nrow = length(vars), ncol = length(vars))
    rownames(m) = vars
    colnames(m) = vars

    diag(m) = 1
    m
  }

  get_matrix = function(id, vars, data, transform = NULL){

    m = make_matrix(vars)
    dat = data[which(data[,cluster, with = FALSE] == id),]

    for(col in seq_along(colnames(m))){
      for(row in seq_along(rownames(m))){
        #message("row = ",row, "col = ", col)
        if(col == row){
          m[row, col] = 1
          next
        }

        current_vars = c(colnames(m)[col], rownames(m)[row])
        valid_rows = unlist(dat[,var1, with = FALSE]) %in% current_vars &
          unlist(dat[,var2, with = FALSE]) %in% current_vars &
          !is.na(unlist(dat[,yi, with = FALSE])) &
          !is.na(unlist(dat[,vi, with = FALSE]))

        if(sum(valid_rows) == 0){
          m[row,col] = NA
          next
        }

        y <- dat[valid_rows, yi, with = FALSE]
        v <- dat[valid_rows, vi, with = FALSE]

        res <- stats::weighted.mean(y, 1/v)

        if(!is.null(transf)){
          res <- transf(res)
        }
        m[row,col] = res
      }

    }
    rownames(m) = rownames(m)
    colnames(m) = colnames(m)
    m

  }

  unique_ids = unique(unlist(data[,cluster, with = FALSE]))

  dat_list = list()
  dat_list$data = lapply(unique_ids, function(x){
    get_matrix(x, vars = vars, data, transform)})

  names(dat_list$data) <- unique_ids

  # Apply moderators

  moderators <- as.list(substitute(list(...)))[-1L]

  n_dat <- data[,
                list(n = mean(eval(parse(text = ni), envir = .SD), na.rm = TRUE))
                , by = cluster]

  n_order <- order(match(n_dat[[1]], names(dat_list$data)))

  dat_list <- append(dat_list, n_dat[n_order, -1])

  if(length(moderators) > 0){

  moderators <- data[,
       {
       lapply(moderators, function(x) eval(x, envir = .SD))
         }
       , by = cluster]

  order_mods <- order(match(moderators[[1]], names(dat_list$data)))
  dat_list <- append(dat_list, moderators[order_mods,-1])

  }

  dat_list

}


star_matrix <- function(m, stars) {
  get_stars = function(p, stars) {
    if (is.na(p))
      p <- 1
    n_stars = sum(p < stars)
    paste(rep("*", n_stars), collapse = "")
  }

  s_matrix = m
  s_matrix[] =  sapply(m, function(p)
    get_stars(p, stars = stars))
  return(s_matrix)
}

#' osmasem_table

report_tssem2 <- function(x, pattern = NULL) {
  if(!methods::is(x, "wls")) {
    stop("Can only be used with wls objects")
  }
  if(is.null(pattern)) pattern <- "The model had {assessment} fit $\\chi^2$({df}) = {chi}; $p$ = {p}; RMSEA = {RMSEA}; SRMR = {SRMR}; TLI = {TLI}"
  x <- summary(x)
  l = as.list(x$stat)
  names(l) = rownames(x$stat)

  if (l$`DF of target model` == 0) {
    mess = "The model had zero degrees of freedom and so had perfect fit."
    return(mess)
  }

  RMSEA <- l$RMSEA
  RMSEA_cat <- dplyr::case_when(RMSEA <= 0.01  ~ 1, # excellent (MacCallum, Browne, Sugawara)
                                RMSEA <= 0.05  ~ 2, # good
                                RMSEA <= 0.08  ~ 3, # mediocre
                                RMSEA > 0.08   ~ 4) # poor
  TLI <- l$TLI
  TLI_cat <- dplyr::case_when(TLI > 0.95  ~ 1, # excellent
                              TLI >= .95  ~ 2, # good
                              TLI >= 0.9  ~ 3, # mediocre
                              TLI <  0.9  ~ 4) # poor
  CFI <- l$CFI
  CFI_cat <- dplyr::case_when(CFI > 0.95  ~ 1, # excellent
                              CFI >= .95  ~ 2, # good
                              CFI >= 0.9  ~ 3, # mediocre
                              CFI <  0.9  ~ 4) # poor
  SRMR <- l$SRMR
  SRMR_cat <- dplyr::case_when(SRMR <= 0.04  ~ 1, # excellent (MacCallum, Browne, Sugawara)
                               SRMR <= 0.08  ~ 2, # good
                               SRMR <= 0.1   ~ 3, # mediocre
                               SRMR > .1     ~ 4) # poor



  assessment <- round(mean(c(TLI_cat, RMSEA_cat, SRMR_cat)),0)
  assessment <- c("excellent", "good","mediocre", "poor")[assessment]

  chi = digits(l$`Chi-square of target model`,2)
  df = digits(l$`DF of target model`,0)
  p = round_p(l$`p value of target model`)

  RMSEA = digits(RMSEA, 3)
  TLI = digits(TLI, 3)
  CFI = digits(CFI, 3)
  SRMR = digits(SRMR,3)

  mess <- glue::glue(pattern)
  mess
}

to_rowhead <- function(data, x, italics = FALSE) {
  row_head = as.character(unlist(data[, x]))
  new_data = data[, !names(data) %in% x, drop = FALSE]
  new_data$indent_ = T

  new_head = lapply(seq_along(row_head), function(i) {
    new_head = T
    i = unlist(i)
    if (i > 1) {
      if (row_head[[i]] == row_head[[i - 1]]) {
        new_head = F
      }
    }
    return(new_head)
  })

  new_head = unlist(new_head)

  table_out = lapply(seq_along(new_head), function(i) {
    if (new_head[i]) {
      new_row = new_data[i, , drop = FALSE]
      new_row[, 1] = unlist(data[i, x])

      if(italics){
        new_row[,1] = paste0("*",new_row[,1],"*")
      }

      new_row[, 2:ncol(new_row)] = ""
      new_row$indent_ = F

      return(rbind(new_row, new_data[i,]))

    } else{
      return(new_data[i, , drop = FALSE])
    }
  })

  out <- do.call(rbind, table_out)
  attr(out, "indent") <- out$indent_
  out$indent_ <- NULL
  class(out) <- c("KIN_summary","data.frame")
  out
}
