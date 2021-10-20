## Interface to rbacon

# Bacon wrapper ------


#' A Wrapper for rbacon::Bacon
#'
#' @param depth A vector of depths
#' @param obs_age A vector of observed ages
#' @param obs_err A vector of errors on the observed ages
#' @description Wraps the Bacon function from rbacon so that it can be used in a
#'   more typical "R" way. Returns age-depth models in a format to match hamstr
#'   output. Most Bacon functionality is accessible, for example hiatuses and
#'   boundaries, different calibration curves, postbomb curves.
#'   "accept.suggestions" defaults to TRUE so that the Bacon's accumulation rate
#'   prior suggestion is used by default, however, the section thickness
#'   parameter "thick" is not changed.
#' @inheritParams rbacon::Bacon
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' hb1 <- hamstr_bacon(id = "sdf",
#'                     depth = MSB2K$depth,
#'                     obs_age = MSB2K$age,
#'                     obs_err = MSB2K$error
#'                     )
#'
#' plot(hb1)
#' predict(hb1)
#' summary(hb1)
#' }
hamstr_bacon <- function(
  id = "default",
  depth,
  obs_age, obs_err,
  thick = 5,
  d.min = NA, d.max = NA,
  d.by = NULL,
  seed = NA,
  acc.shape = 1.5,  acc.mean = 20,
  mem.strength = 10, mem.mean = 0.5,
  boundary = NA,
  hiatus.depths = NA, hiatus.max = 10000,
  add = c(),
  cc = 1, cc1 = "IntCal20", cc2 = "Marine20", cc3 = "SHCal20",
  cc4 = "ConstCal", ccdir = "", postbomb = 0, delta.R = 0,
  delta.STD = 0, t.a = 3, t.b = 4, normal = FALSE,
  suggest = TRUE, accept.suggestions = TRUE,
  reswarn = c(10, 200),
  ask = FALSE,
  slump = c(),
  remove = FALSE, ssize = 2000, th0 = c(),
  burnin = min(500, ssize), MinAge = c(), MaxAge = c(),
  plot.pdf = FALSE,
  close.connections = FALSE,
  verbose = FALSE, suppress.plots = TRUE
){
  if(packageVersion("rbacon") < "2.5.2")
    stop("hamstr_bacon requires rbacon version 2.5.2 or higher")

  # check inputs
  if (any(length(depth) == 0, length(obs_age) == 0, length(obs_err)== 0))
    stop("depth, age, or obs_err are missing")

  #if (length(cc) > 1) stop("Argument cc can only be length 1")


  if (is.na(d.min)) d.min <- min(depth)
  if (is.na(d.max)) d.max <- max(depth)

  # use temp directory to store Bacon input and output
  tmpdir <- tempdir()

  dirbase <- basename(tmpdir)
  dirnm <- dirname(tmpdir)

  # create datafile in format required for Bacon
  datfl <- tempfile(tmpdir = tmpdir)

  bacon_dat <- dplyr::tibble(
      id = id,
      age = obs_age,
      error = obs_err,
      depth = depth,
      cc = cc,
      delta.R = delta.R,
      delta.STD = delta.STD
    )


  #fl <- suppressWarnings(normalizePath(paste0(tmpdir, "//", dirbase, ".csv")))

  fl <- paste0(tmpdir, "//", dirbase, ".csv")

  utils::write.csv(bacon_dat,
                   file = fl,
                   row.names = FALSE, quote = FALSE)



  # capture the passed arguments
  pars <- c(as.list(environment()))

  # hard code remember
  pars$remember <- FALSE

  # set name and location for output
  pars$core = dirbase
  pars$coredir = normalizePath(paste0(dirnm, "/"))

  # subset of args for Bacon
  to.keep <- names(pars) %in% methods::formalArgs(Bacon2)
  pars <- pars[to.keep]


  # get defaults values of arguments to Bacon
  default.args <- formals(Bacon2)
  default.arg.nms <- names(default.args)

  # Overwrite the defaults with non-null passed arguments
  #dpars <- pars
  non.null.pars <- pars[lapply(pars, is.null) == FALSE]

  merged.pars <- default.args
  merged.pars[names(non.null.pars)] <- non.null.pars



  # set d.by to thick unless specified. We will do interpolation separately.
  if (is.null(d.by)){
    pars$d.by <- merged.pars$thick
  }


  # Call Bacon


  do.call(Bacon2, pars)

  par_list <- info[names(info) %in% names(pars)]

  #posterior <- info$output

  # create output, add class attributes and return
  out <- list(pars = par_list, data = bacon_dat, info = info)
  class(out) <- append("hamstr_bacon_fit", class(out))

  rbacon::Bacon.cleanup()


  return(out)
}


#' Stack iterations
#' @param x extracted bacon age models
#'
#' @return
#' @keywords internal
.StackIterations <- function(x){
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  n.row <- nrow(x)

  x <- utils::stack(x, stringsAsFactors = FALSE)
  x$depth.index <- (1:n.row) -1
  x
}


#' Extract and reconstruct bacon age models
#'
#' @param hamstr_bacon_fit
#'
#' @return
#'
#' @keywords internal
return_bacon_age_mods <- function(hamstr_bacon_fit){

  posterior <- hamstr_bacon_fit$info$output
  pars <- hamstr_bacon_fit$pars

  n.col <- ncol(posterior)
  log.lik <- posterior[, n.col]

  posterior[, 2:(n.col - 2)] <- posterior[, 2:(n.col - 2)] * pars$thick

  age_mods <- apply(posterior[, 1:(n.col - 2)], 1, cumsum)

  age_mods <- .StackIterations(age_mods)

  age_mods$depth <- pars$d.min + age_mods$depth.index * pars$thick

  age_mods <- age_mods[, c("ind", "depth", "values")]
  names(age_mods) <- c("iter", "depth", "age")

  age_mods$iter <- readr::parse_number(as.character(age_mods$iter))

  return(dplyr::as_tibble(age_mods))

}




#' Title
#'
#' @param hamstr_bacon_fit
#' @param depth
#'
#' @return
#'
#' @keywords internal
interpolate_bacon_age_models <- function(hamstr_bacon_fit, depth){

  if (is.numeric(depth) == FALSE) {

    depth <- match.arg(depth, choices = c("modelled", "data"))

    if (depth == "modelled") {

      # get posterior age models
      out <- return_bacon_age_mods(hamstr_bacon_fit)

      class(out) <- append("hamstr_bacon_interpolated_ages", class(out))

      return(out)

    } else if (depth == "data") {

      depth <- hamstr_bacon_fit$data$depth

    }

    }

  # get posterior age models
  pst_age <- return_bacon_age_mods(hamstr_bacon_fit)

  # use base list, split methods, much faster than dplyr::do
  pst_age_lst <- split(pst_age, pst_age$iter)

  new_pst_age <- lapply(pst_age_lst, function(x) {
    stats::approx(x$depth, x$age, depth)$y
  })

  out <- expand.grid(depth = depth,
                     iter = 1:length(pst_age_lst)
  )

  out$age <- unlist(new_pst_age)

  out <- dplyr::as_tibble(out)
  out <- out[,c(2,1,3)]

  class(out) <- append("hamstr_bacon_interpolated_ages", class(out))

  return(out)
}



#' Summarise Bacon age models
#'
#' @param hamstr_bacon_fit
#'
#' @return
#'
#' @keywords internal
summarise_bacon_age_models <- function(hamstr_bacon_fit){

  age_mods <- return_bacon_age_mods(hamstr_bacon_fit)

  age_summary <- age_mods %>%
    dplyr::group_by(depth) %>%
    dplyr::summarise(mean = mean(age),
              sd = stats::sd(age),
              `2.5%` = stats::quantile(age, probs = 0.025, na.rm = T),
              `25%` = stats::quantile(age, probs = 0.25, na.rm = T),
              `50%` = stats::quantile(age, probs = 0.5, na.rm = T),
              `75%` = stats::quantile(age, probs = 0.75, na.rm = T),
              `97.5%` = stats::quantile(age, probs = 0.975, na.rm = T),
              )

  return(age_summary)
}


#' Summarise interpolated Bacon age models
#'
#' @param interpolated_age_mods
#'
#' @return
#' @keywords internal
summarise_interpolated_bacon_age_models <- function(interpolated_age_mods){

  age_summary <- interpolated_age_mods %>%
    dplyr::group_by(depth) %>%
    dplyr::summarise(mean = mean(age),
                     sd = stats::sd(age),
                     `2.5%` = stats::quantile(age, probs = 0.025, na.rm = T),
                     `25%` = stats::quantile(age, probs = 0.25, na.rm = T),
                     `50%` = stats::quantile(age, probs = 0.5, na.rm = T),
                     `75%` = stats::quantile(age, probs = 0.75, na.rm = T),
                     `97.5%` = stats::quantile(age, probs = 0.975, na.rm = T),
    )

  return(age_summary)
}



#' Title
#'
#' @param hamstr_bacon_fit
#' @param summarise
#' @param n.iter
#'
#' @return
#' @keywords internal
plot_hamstr_bacon_fit <- function(hamstr_bacon_fit, summarise = TRUE, n.iter = 1000) {

  if (summarise == TRUE){
    p.fit <- plot_summary_bacon_age_models(hamstr_bacon_fit)
  } else if (summarise == FALSE){
    p.fit <- plot_bacon_age_models(hamstr_bacon_fit, n.iter = n.iter)
  }
  return(p.fit)
}



#' Title
#'
#' @param hamstr_bacon_fit
#'
#' @return
#' @keywords internal
get_bacon_obs_ages <- function(hamstr_bacon_fit){

  pdf.sums <- lapply(hamstr_bacon_fit$info$calib$probs,
                     function(x) {
                       # suppress warnings about modes as mode not used anyway
                       suppressWarnings(
                         SummariseEmpiricalPDF(x[,1], x[,2])
                         )
                       })

  obs_ages <- dplyr::bind_rows(
    lapply(pdf.sums, function(x) c(age = x[["median"]], err = x[["sd"]]))
  )    %>%
    dplyr::mutate(age_upr = age + 2*err,
                  age_lwr = age - 2*err) %>%
    dplyr::mutate(depth = hamstr_bacon_fit$info$calib$d)

  return(obs_ages)

}



#' Plot individual Bacon age models
#'
#' @param hamstr_bacon_fit
#' @param n.iter
#'
#' @return
#' @keywords internal
plot_bacon_age_models <- function(hamstr_bacon_fit, n.iter = 1000){


  pdf.sums <- lapply(hamstr_bacon_fit$info$calib$probs,
                     function(x) {
                       # suppress warnings about modes as mode not used anyway
                       suppressWarnings(
                       SummariseEmpiricalPDF(x[,1], x[,2])
                       )
                       })


  obs_ages <- get_bacon_obs_ages(hamstr_bacon_fit)


  posterior_ages <- return_bacon_age_mods(hamstr_bacon_fit)

  p.fit <- posterior_ages %>%
    dplyr::filter(iter %in% sample(unique(iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(ggplot2::aes(x = depth, y = age, group = iter))


  p.fit <- p.fit +
    ggplot2::geom_line(alpha = 0.5 / sqrt(n.iter))

  p.fit <- p.fit +
    ggplot2::geom_linerange(data = obs_ages,
                            ggplot2::aes(x = depth,
                                         ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                            colour = "Blue", size = 1.25) +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(x = depth, y = age), inherit.aes = FALSE,
                        colour = "Blue") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "Depth", y = "Age")

  p.fit
}

#' Plot summarised Bacon age models
#'
#' @param hamstr_bacon_fit
#'
#' @return
#' @keywords internal
plot_summary_bacon_age_models <- function(hamstr_bacon_fit){


  age_summary <- summarise_bacon_age_models(hamstr_bacon_fit)

  obs_ages <- get_bacon_obs_ages(hamstr_bacon_fit)

  p.age.sum <- age_summary %>%
    plot_downcore_summary(.) +
    ggplot2::labs(x = "Depth", y = "Age")

  p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(data = obs_ages,
                            ggplot2::aes(x = depth,
                                         ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                            colour = "Blue", size = 1.25) +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age),
                        colour = "Blue")

  p.age.sum
}


# Methods -------

#' Interpolate Age Models at Given Depths
#' @description Method for generic function predict. Returns the posterior age
#' models interpolated to new depths given in depth.
#' @param object
#' @param depth
#' @inheritParams interpolate_bacon_age_models
#' @return
#'
#' @examples
#' @export
#' @method predict hamstr_bacon_fit
predict.hamstr_bacon_fit <- function(object, depth = "modelled"){

  interpolate_bacon_age_models(object, depth)

}

#' Title
#'
#' @param object a hamstr_bacon_fit object
#' @return A ggplot object
#'
#' @examples
#' @export
#' @method plot hamstr_bacon_fit
plot.hamstr_bacon_fit <- function(hamstr_bacon_fit,
                                  type = "default",
                                  summarise = TRUE,
                                  ...){

  plot_hamstr_bacon_fit(hamstr_bacon_fit,
                        summarise = summarise, ...)

}

#' Summarise hamstr_bacon Age Models
#'
#' @param object
#' @param type
#'
#' @return
#'
#' @examples
#' @export
#' @method summary hamstr_bacon_fit
summary.hamstr_bacon_fit <- function(object){

  summarise_bacon_age_models(object)

}


#' Summarise hamstr_bacon Age Models
#'
#' @param object
#' @param type
#'
#' @return
#'
#' @examples
#' @export
#' @method summary hamstr_bacon_interpolated_ages
summary.hamstr_bacon_interpolated_ages <- function(object){

  summarise_interpolated_bacon_age_models(object)

}


