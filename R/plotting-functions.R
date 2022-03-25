#' diamond_df
#'
#' This function creates diamonds for moderators and summary statistics.
#' @param plot a ggplot2
#' @param data a data.frame
#' @param fill ggplot2 fill value
#' @param colour character. A hex code. If provided, gives diamond coloured border
#' @importFrom ggplot2 aes
add_diamond = function(plot, data, fill = "grey20", colour = NA) {
  diamond_shape = data.frame(
    x = c(data$lower, data$est, data$upper, data$est),
    y = c(data$position,data$position + .4,data$position,data$position - .4),
    names = c("xmin", "ymax", "xmax", "ymax"),
    setting = "Pooled"
  )
  plot = plot + ggplot2::geom_polygon(
    data = diamond_shape,
    ggplot2::aes(x = x, y = y),
    fill = fill,
    colour = colour,
    inherit.aes = F
  )
  plot
}

#' forest_plot
#'
#' Produce a forest plot for meta_list objects
#'
#' @param model a meta_list object
#' @param ... moderators to plot character strings only
#' @param xlab a string.
#' @param effect_label a string. If provided relabels "effect size" in facetting.
#' @param transf a function. If supplied effect sizes are transformed by this function
#' @param baseline_name a string. The label for the baseline model
#' @param factor.levels  a character vector. If supplied, only the factor.levels specified will be plotted.
#' @param facet_by a colname. Facets effect sizes by a supplied variable.
#' @param vline a scalar. Dictates the x-intercept (the dashed line).
#' @param author the name of the author column
#' @param year the name of the year column
#' @param moderator.shape a scalar ggplot2 geom_point shape value
#' @param moderator.size a scalar. ggplot2 geom_point size value
#' @param summary.shape a scalar ggplot2 geom_point shape value
#' @param summary.size a scalar ggplot2 geom_point size value
#' @param moderator_diamond a bool. If true, diamond is created for moderators
#' @param font A string. The name of a font family. Defaults to serif
#' @param envir the calling environment to retrieve model data by name
#' @import ggplot2
#' @export

forest_plot <- function(model,
                        ...,
                        xlab = "Effect size",
                        effect_label = NULL,
                        transf = NULL,
                        baseline_name = "Pooled estimate",
                        factor.levels = NULL,
                        facet_by = NULL,
                        vline = 0,
                        author = NULL,
                        year = NULL,
                        moderator.shape = 23,
                        moderator.size = 3,
                        summary.shape = 23,
                        summary.size = 4,
                        moderator_diamond = FALSE,
                        font = "serif",
                        envir = parent.frame()) {

  if (is.null(transf)) {
    transf <- function(x)
      x
  }
  if (methods::is(model, "meta3"))
    model <- moderate(model)

  .internalDat <- stats::coef(model)
  .internalDat$est <- transf(.internalDat$est)
  .internalDat$lower <- transf(.internalDat$lower)
  .internalDat$upper <- transf(.internalDat$upper)
  vline   <- transf(vline)
  .internalDat$order <- seq_along(.internalDat$est)

  if (!methods::is(model, "meta_list"))
    stop("model is a ", class(model), ". model class must be a 'meta_list'")

  # Filter to requested moderators
  mod <- unlist(list(...))
  .internalDat <-
    .internalDat[is.na(moderation) |
                   moderation %in% c(mod, "Baseline"),]

  source_data <- eval(model$models$Baseline$call$data, envir = envir)

  if (is.null(author))
    author <- find_author(source_data)
  if (is.null(year))
    year <- find_year(source_data)
  cluster <- as.character(model$models$Baseline$call$cluster)

  key <-
    data.table::data.table(source_data)[, c(cluster, year, author), with = FALSE]
  names(key) = c("cluster", "year", "author")
  key <- unique(key[, cluster := as.character(cluster)])
  if (max(table(key$cluster)) > 1)
    warning("Inconsistent author and year for at least one cluster. Check data.")

  plot_dat <- merge(.internalDat, key, all.x = TRUE, by = "cluster")

  plot_dat$year[is.na(plot_dat$year)] <-
    plot_dat$order[is.na(plot_dat$year)]


  plot_dat[type == "Effect sizes"]$cluster <-
    with(plot_dat[type == "Effect sizes"], paste(author, year))

  # Order plot data
  plot_dat <- plot_dat[order(plot_dat$year),]

  plot_dat$cluster <- factor(plot_dat$cluster, levels = unique(plot_dat$cluster))

  mod_levels <-
    levels(droplevels(plot_dat[type == "moderator level"]$cluster))
  effect_levels <-
    levels(droplevels(plot_dat[type == "Effect sizes"]$cluster))
  plot_dat$cluster <-
    factor(plot_dat$cluster,
           levels = c("Baseline", mod_levels, effect_levels))
  levels(plot_dat$cluster)[1] <- baseline_name


  plot_dat$position <- as.numeric(plot_dat$cluster)


  p <- ggplot2::ggplot(plot_dat,
                       ggplot2::aes(
                         y = cluster,
                         x = est,
                         xmin = lower,
                         xmax = upper
                       )) +
    ggplot2::labs(x = xlab, y = "") +
    ggplot2::facet_grid(rows = ggplot2::vars(setting),
                        scales = 'free',
                        space = 'free_y') +
    ggplot2::theme_classic() +
    ggplot2::geom_point(data = plot_dat) +
    ggplot2::geom_errorbar(data = plot_dat[type == "Effect sizes",], width = .1) +
    ggplot2::geom_vline(xintercept = vline,
                        #add horizontal line
                        color = 'black',
                        linetype = 'dashed')

  p <- add_diamond(p, plot_dat[type == "Baseline"])

  if (moderator_diamond) {
    types = plot_dat$cluster[plot_dat$type == "moderator level"]
    for (i in types) {
      p <-
        add_diamond(p, plot_dat[plot_dat$cluster == i,], fill = "white", colour = "grey20")
    }
  } else {
    p = p + ggplot2::geom_errorbar(data = plot_dat[type == "moderator level",] , width = .25) + ggplot2::geom_point(
      #add summary points
      data = plot_dat[type == "moderator level", ],
      color = 'black',
      shape = moderator.shape,
      size = moderator.size,
      fill = "white"
    )
  }

      if (!is.null(font))
    p <-
    p + ggplot2::theme(text = ggplot2::element_text(family = font))

  if (nrow(plot_dat[plot_dat$setting == "Pooled",]) < 2) {
    p <-
      p + ggplot2::theme(
        strip.text.y = ggplot2::element_text(angle = 0),
        strip.background.y = ggplot2::element_blank(),
        text = element_text(family = font)
      )
  }


  p

}

#' find_year
#'
#' detects the string name of a column with year information
#' @param data data object

find_year = function(data) {
  data = data.frame(data)
  count_chars = lapply(seq_along(names(data)), function(i) {
    suppressWarnings(var <-
                       as.numeric(as.character(data[, names(data)[i]]))) #strip factors, make numeric
    year = strsplit(as.character(Sys.Date()), split = "-")[[1]][1]
    var = ifelse(var > as.numeric(year) + 1 , NA, var)
    var = ifelse(var < 1800, NA, var)
    n = nchar(as.character(var)) == 4
    sum(n, na.rm = T) / length(n)
  })
  string = names(data)[which.max(count_chars)]
  message(paste0("year was not specified, using: '", string, "'."))
  return(string)
}


#' find_author
#'
#' detects the string name of a column with year information
#' @param data data object

find_author = function(data) {
  data = data.frame(data)
  vars = names(data)
  has_author_title =  as.numeric(agrepl("author", tolower(vars)))
  et_al = lapply(vars, function(i) {
    num = grepl("et al", tolower(data[, i]))
    journal = grepl("journal", tolower(data[, i])) * 2
    sum(num, na.rm = T) - sum(journal, na.rm = T)
  })
  score = has_author_title + unlist(et_al)
  string = vars[which.max(score)]
  message(paste0(
    "author was not specified, using: '",
    string,
    "'."
  ))
  string
}


#' funnel_plot
#'
#' This function is used for plotting funnel plots
#'
#' @param model an object belonging to the class 'meta_list'. These objects are created by the function 'meta3_moderation'.
#' @param xlab a character string. Label for the x-axis.
#' @param ylab a character string. Label for the y-axis.
#' @param font a character string. Set's font family. Defaults to times new roman ('serif')
#' @param trimfill should trimfill be performed?
#' @param aggregate should data be aggregated to cluster level?
#' @param alpha a scalar. Sets alpha transparency.
#' @param size a scalar. The ggplot2 size value for the points
#' @param CI_linetype ggplot2 linetype for confidence intervals
#' @param CI_size ggplot2 size value for confidence intervals
#' @param pool_linetype ggplot2 linetype for pooled estimate
#' @param pool_size ggplot2 size value for pooled estimate
#' @param density a bool. If True, plots densit of points.
#' @param ... extra arguments passed to rma.uni for the purpose of trimfill.
#' @details produces a funnel plot using ggplot2. If aggregate is set to TRUE, data is aggregated using fixed-effects meta-analysis to each cluster level (the pooled effect is not updated). If trimfill is set to TRUE, extra observations will be created as needed to achieve symmetry using metafor::trimfill. If trimfill results in new studies, the pooled effect will be updated based on the results of that meta-analysis—which is not multi-level. Note that funnel plots are often difficult to interpret, especially for multi-level meta-analysis. Note that the assumptions of trim and fill are violated if there are dependencies within the data (e.g., if data is clustered and you used multi-level meta-analysis to address this clustering).
#'
#' @importFrom ggplot2 aes theme xlab geom_point coord_flip scale_x_reverse geom_line geom_segment geom_errorbar element_text stat_density_2d stat scale_fill_gradient scale_shape_discrete theme_bw
#' @export funnel_plot

#test arguments:
funnel_plot <- function(model,
                       xlab = "Estimate",
                       ylab = "Standard error",
                       font = "serif",
                       aggregate = FALSE,
                       trimfill = FALSE,
                       alpha = 1,
                       size = 1.5,
                       CI_linetype = "dashed",
                       CI_size = .7,
                       pool_linetype = "dotted",
                       pool_size = .5,
                       density = FALSE,
                       ...
                       ) {
   t_model <- NULL
  if (methods::is(model, "meta_list")) {
    t_model <- model$models[[1]]
  }
  if (methods::is(model, "meta")) {
    t_model = model
  }

  estimate <- summary(t_model)$coefficients["Intercept","Estimate"]
  se <- summary(t_model)$coefficients["Intercept","Std.Error"]

  requireNamespace("metafor", quietly = TRUE)

  funnel_data <- t_model$data

    if (!is.null(funnel_data$x1)) {
    stop(
      "Moderated models cannot be used to create funnel plots as they have no single estimate. Please use a baseline model.",
      call. = F
    )
  }

  funnel_data$type <- "effect_size"

  # Obtain aggregated and trim fill data
  if(aggregate){
  funnel_data <- aggregate_to_cluster(t_model)
  funnel_data$type <- "aggregate"
  }
  funnel_data$trimfill <- "Studies"
  n_obs <- nrow(funnel_data)

  if(trimfill) {
    trimfill_m <- trim_and_fill(t_model, aggregate = aggregate, ...)

    if (length(trimfill_m$yi) > n_obs) {
      extra_rows <- (n_obs + 1):length(trimfill_m$yi)
      extra_data <- data.frame(cluster = "trimfill",
                               y = trimfill_m$yi[extra_rows],
                               v = trimfill_m$vi[extra_rows])
      extra_data$type <-
        ifelse(aggregate, "aggregate", "effect_size")
      extra_data$trimfill <- "Filled studies"
      funnel_data <- rbindlist(list(funnel_data, extra_data), fill = TRUE)

      estimate <- trimfill_m$b[[1]]
      se <- trimfill_m$se
    }
  }

  crit_val <- stats::qnorm(0.975)

  funnel_data$se <- sqrt(funnel_data$v)


  se.seq <- seq(0, max(funnel_data$se), 0.001)

  if(trimfill){

  }

  ll95 <- estimate - (crit_val * se.seq)
  ul95 <- estimate + (crit_val * se.seq)
  meanll95 <- estimate - crit_val * se
  meanul95 <- estimate + crit_val * se

  funnel_data$lower = funnel_data$y  - crit_val * funnel_data$se
  funnel_data$upper = funnel_data$y + crit_val * funnel_data$se
  dfCI = data.frame(ll95, ul95, se.seq, estimate, meanll95, meanul95)

  # Start ggplot2

  funnel_data$trimfill <- factor(funnel_data$trimfill, levels = c("Studies", "Filled studies"))

  fp <-
    ggplot2::ggplot(funnel_data,
                    ggplot2::aes(
                      x = se,
                      y = y,
                      ymin = lower,
                      ymax = upper,
                      shape = trimfill
                    )) + scale_shape_manual(values = c(16, 1)) +
    theme_bw() + coord_flip() +
    labs(x = ylab, y = xlab, shape = "") +
    ggplot2::scale_x_reverse()

  if (density) {
    fp <- fp +
      ggplot2::stat_density_2d(
        ggplot2::aes(fill = stat(level)),
        bins = 7,
        colour = "white",
        geom = "polygon",
        alpha = 0.2,
        show.legend = F
      ) +
      ggplot2::scale_fill_gradient(low = "grey60", high = "grey30")
  }

  # Apply studies and lines

  fp <- fp + geom_point(size = size, alpha = alpha) + ggplot2::geom_line(
    ggplot2::aes(x = se.seq, y = ll95),
    linetype = CI_linetype,
    data = dfCI,
    size = CI_size,
    inherit.aes = FALSE
  ) +
    ggplot2::geom_line(
      ggplot2::aes(x = se.seq, y = ul95),
      linetype = CI_linetype,
      data = dfCI,
      size = CI_size,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min(se.seq),
        y = estimate,
        xend = max(se.seq),
        yend = estimate
      ),
      linetype = pool_linetype,
      data = dfCI,
      size = pool_size,
      inherit.aes = FALSE
    )


  if (!is.null(font)) {
    fp <-
      fp + ggplot2::theme(text = ggplot2::element_text(family = font))
  }


  if (sum(funnel_data$trimfill == "Filled studies") == 0) {
    fp <- fp + theme(legend.position = "none")

  }

  fp
}

#' forest_height
#'
#' Use this function to resize forest plots in rmarkdown.
#'
#' @param meta3_plot a moderated meta3 plot
#' @param slope the numeric value to multiple number of rows by
#' @param intercept the numeric constant
#' @export forest_height

forest_height = function(meta3_plot, slope = .12, intercept = .52){
  length(unique(meta3_plot$data$cluster)) * slope + intercept
}

#' moderation_matrix
#'
#' Plots a moderation matrix
#' @param ... Named meta_list models. All common moderators should have the same name.
#' @param effect_size a string
#' @param moderators a list of moderators to include, order retained.
#' @param null_value a scalar indicating non-significance (where the dashed line will be drawn).
#' @param transf function to transform values
#' @param leading_zero when true, leading zeros are allowed on the x-axis
#' @param black list of outcome names containing their moderators
#' @param replace a vector with names included. gsub will be applied to the moderation column such that the vector's names are replaced with the vector's contents
#' @import ggplot2 data.table
#' @export

moderation_matrix <- function(..., effect_size = "Effect size", moderators = NULL,
                              null_value = NULL, transf = NULL, leading_zero = TRUE,
                              black = NULL, replace = c("_" = " ")){

  models <- list(...)
  #return(models)
  if(is.null(names(models))) stop("All models must be named within the moderation_matrix function call")

  dat_list <- lapply(models, function(x) stats::coef(x))

  DL <- lapply(seq_along(dat_list), function(x){
    dat_list[[x]]$outcome = names(dat_list)[x]
    dat_list[[x]]
  })

  DL <- data.table::rbindlist(DL)

  if(is.null(transf)){
    transf = function(x) x
  }

  if(is.null(null_value)) null_value <- transf(0)

  DL <- DL[type != "Effect sizes", .(
    y = transf(est),
    lower = transf(lower),
    upper = transf(upper),
    cluster = factor(cluster, levels = unique(cluster)),
    moderation,
    outcome = factor(outcome, levels = names(models)),
    type,
    model_p
  )]

  graph_dat <- DL

  if(is.null(moderators)){
    mod_levels <- unique(graph_dat$moderation)
    mod_levels <- mod_levels[mod_levels != "Baseline"]
  }else{
    mod_levels = moderators
  }
  graph_dat$moderation[graph_dat$moderation == "Baseline"] = ""
  mod_levels = c(mod_levels, "")

  graph_dat$moderation = factor(graph_dat$moderation, levels = mod_levels)
  graph_dat = graph_dat[!is.na(graph_dat$moderation),]
  graph_dat$outcome = factor(graph_dat$outcome)
  # ------- fix cluster order

  cluster_levels <- levels(as.factor(graph_dat$cluster))
  cluster_levels <- cluster_levels[cluster_levels != "Baseline"]
  cluster_levels <- c(cluster_levels, "Baseline")

  graph_dat$cluster = factor(graph_dat$cluster, levels = rev(cluster_levels))

  # prepare significance -------

  final_dat <- graph_dat # give baseline a p value
  final_dat$model_p <- tidyr::replace_na(final_dat$model_p, 0)  # Baseline p is NA. I want it
  final_dat <- stats::na.omit(final_dat)

  setkey(final_dat, outcome, moderation)
  sig_dat  <- final_dat[CJ(outcome,moderation, unique = TRUE), .(p = mean(model_p, na.rm = T), y = 0, cluster = NA) , by = .EACHI]
  sig_dat$p[is.na(sig_dat$p)] = 1

  sig_dat$cluster = factor(sig_dat$cluster, levels = levels(final_dat$cluster))

  if(!is.null(black)){

    black <- unlist(sapply(seq_along(black), function(i){
      paste(names(black)[i], black[[i]])
    }))

  sig_dat$black <- paste(sig_dat$outcome, sig_dat$moderation) %in% black

  }else{
    sig_dat$black <- FALSE
  }

    if (length(replace) > 0 & !identical(replace, FALSE)) {
      for (i in seq_along(replace)) {
        final_dat$cluster <-
          gsub(names(replace)[i], replace[i], final_dat$cluster)
      }
      final_dat$cluster <- factor(final_dat$cluster, levels = unique(final_dat$cluster))
    }

  p <- ggplot(final_dat, aes(
    x = y,
    y = cluster,
    xmin = lower,
    xmax = upper
  )) +
    ggplot2::geom_rect(data = sig_dat[sig_dat$p >= 0.05 | is.na(sig_dat$p),], # add in grey rectangles if not sig.
                       fill = "black",
                       xmin = -Inf, xmax = Inf,
                       ymin = -Inf, ymax = Inf,
                       alpha = 0.15, inherit.aes = F) +    # inherit.aes caused factors to lose order

    # Add in black rectangles

    ggplot2::geom_rect(data = sig_dat[sig_dat$black,],
                       fill = "black",
                       xmin = -Inf, xmax = Inf,
                       ymin = -Inf, ymax = Inf,
                       alpha = 1, inherit.aes = F) +



    geom_vline(xintercept = null_value, linetype = 2) + # add in vertical line at 0
    ggplot2::geom_point() + geom_errorbarh(height = .1)+
    ggplot2::facet_grid( # grid by moderation and outcome
      rows = vars(moderation),
      cols = vars(outcome),
      scales = "free_y",
      space = "free_y"
    ) +
    labs(y = " ", x = effect_size) + theme(text = element_text(family = "serif")) +
    ggplot2::scale_y_discrete(labels = c("Baseline" = expression(bold(Baseline)), parse = T)) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                   strip.background.y = ggplot2::element_blank(),
                   text = ggplot2::element_text(family = "serif"))
  if(!leading_zero){
    p <- p + scale_x_continuous(labels = function(x) gsub("^0\\.",".",x))
  }
  p
}

utils::globalVariables(
  c(
    "x",
    "y",
    ".",
    "cluster",
    "y",
    "v",
    "moderator_level",
    "moderation",
    "mode",
    "Estimate",
    "SE",
    "lbound",
    "ubound",
    "model_p",
    "type",
    "est",
    "lower",
    "upper",
    "setting",
    "k",
    "n",
    "Estimate_formatted",
    "R2_2",
    "R2_3",
    "p.value",
    "level",
    "na.pass",
    "mod",
    "Std.Error",
    "(Intercept)",
    "outcome"
  )
)
