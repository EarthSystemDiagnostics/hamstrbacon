# Methods -------

#' Interpolate Age Models at Given Depths
#' @description Method for generic function predict. Returns the posterior age
#' models interpolated to new depths given in new_depth.
#' @param object 
#' @param new_depth
#' @inheritParams interpolate_bacon_age_models
#' @return
#'
#' @examples
#' @export
#' @method predict hamstr_bacon_fit
predict.hamstr_bacon_fit <- function(object, new_depth = NULL){
  
  interpolate_bacon_age_models(object, new_depth)
  
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




# Bacon wrapper ------


#' A Wrapper for rbacon::Bacon
#' 
#' @param depth 
#' @param obs_age 
#' @param obs_err 
#' @param ... 
#' @description Wraps the Bacon function from rbacon so that it can be used in a
#' more typical "R" way. Returns age-depth models in a format to match hamstr 
#' output. Not all Bacon functionality is accessible, for example there is no 
#' hiatuses and boundaries.
#' @inheritParams rbacon::Bacon
#' @return
#' @export
#'
#' @examples
hamstr_bacon <- function(depth, obs_age, obs_err,
                         thick = NULL,
                         d.min = NULL, d.max = NULL,
                         d.by = NULL,
                         acc.shape = 1.5, acc.mean = 20,
                         mem.strength = 10, mem.mean = 0.5,
                         plot.pdf = FALSE,
                         ask = FALSE,
                         suggest = TRUE, accept.suggestions = TRUE,
                         cc = 1,
                         verbose = FALSE){
  
  if(packageVersion("rbacon") < "2.5.2") stop("hamstr_bacon requires rbacon version 2.5.2 or higher")
  
  # check inputs
  if (any(length(depth) == 0, length(obs_age) == 0, length(obs_err)== 0)) 
    stop("depth, age, or obs_err are missing")
  
  
  if (is.null(d.min)) d.min <- min(depth)
  
  if (is.null(d.max)) d.max <- max(depth)
 
  # use temp directory to store Bacon input and output
  tmpdir <- tempdir()
  
  dirbase <- basename(tmpdir)
  dirnm <- dirname(tmpdir)
  
  # create datafile in format required for Bacon
  datfl <- tempfile(tmpdir = tmpdir)
  
  bacon_dat <- dplyr::tibble(
      id = "test",
      age = obs_age,
      error = obs_err,
      depth = depth
    )
  
  utils::write.csv(bacon_dat, file = paste0(tmpdir, "\\", dirbase, ".csv"),
              row.names = FALSE, quote = FALSE)
  
  
  
  # capture the passed arguments
  pars <- c(as.list(environment()))
  
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
  
  #do.call(rbacon::Bacon, pars)
  do.call(Bacon2, pars)
  
  par_list <- info[names(info) %in% names(pars)]
  
  # read the produced settings file 
  #settings.file <- utils::read.csv(file = paste0(tmpdir, "\\", dirbase, "_settings.txt"), header = FALSE)[,1]

  # create list of used parameters
  # settings.file <- settings.file %>%
  #   dplyr::as_tibble() %>%
  #   tidyr::separate(value, into = c("value", "par"), sep = "#") %>%
  #   dplyr::mutate(value = readr::parse_number(value)) %>%
  #   dplyr::select(par, value)

  #par_list <- as.list(settings.file$value)
  #names(par_list) <- settings.file$par

  #par_list$thick <- info$thick

  #par_list$d.by <- info$d.by
  #par_list$d.min <- info$d.min
  #par_list$d.max <- info$d.max

  # construct the file name
  #K = length(seq(floor(d.min), ceiling(d.max), by = merged.pars$thick))
  #K <- info$K
  
  #outflnm <- paste0(tmpdir, "\\", dirbase, "_", K, ".out")

  # read the posterior
  #posterior <- utils::read.table(outflnm, header = FALSE)

  posterior <- info$output
  
  # create output, add class attributes and return
  out <- list(pars = par_list, data = bacon_dat, posterior = posterior, info = info)
  class(out) <- append("hamstr_bacon_fit", class(out))

  return(out)
}



#' @param x 
#'
#' @return
#' @keywords internal
#'
#' @examples
.StackIterations <- function(x){
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  n.row <- nrow(x)
  
  x <- utils::stack(x, stringsAsFactors = FALSE)
  x$depth.index <- (1:n.row) -1
  x
}


#' Title
#'
#' @param hamstr_bacon_fit 
#'
#' @return
#' @export
#'
#' @examples
return_bacon_age_mods <- function(hamstr_bacon_fit){
  
  posterior <- hamstr_bacon_fit$posterior
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
#' @param new_depth 
#'
#' @return
#' @export
#'
#' @examples
interpolate_bacon_age_models <- function(hamstr_bacon_fit, new_depth){
  
  if (is.null(new_depth)) {
    new_depth <- hamstr_bacon_fit$data$depth
  }
  
  
  # get posterior age models
  pst_age <- hamstr::return_bacon_age_mods(hamstr_bacon_fit)
  
  # use base list, split methods, much faster than dplyr::do
  pst_age_lst <- split(pst_age, pst_age$iter)
  
  new_pst_age <- lapply(pst_age_lst, function(x) {
    stats::approx(x$depth, x$age, new_depth)$y
  })
  
  out <- expand.grid(depth = new_depth,
                     iter = 1:length(pst_age_lst)
  )
  
  out$age <- unlist(new_pst_age)
  
  out <- dplyr::as_tibble(out) 
  out <- out[,c(2,1,3)]
  
  class(out) <- append("hamstr_bacon_interpolated_ages", class(out))
  
  return(out)
}



summarise_bacon_age_models <- function(hamstr_bacon_fit){
  
  age_mods <- hamstr::return_bacon_age_mods(hamstr_bacon_fit)
  
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


summarise_interpolated_bacon_age_models <- function(age_mods){
  
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




plot_hamstr_bacon_fit <- function(hamstr_bacon_fit, summarise = TRUE, n.iter = 1000) {
  
  if (summarise == TRUE){
    p.fit <- plot_summary_bacon_age_models(hamstr_bacon_fit)
  } else if (summarise == FALSE){
    p.fit <- plot_bacon_age_models(hamstr_bacon_fit, n.iter = n.iter)
  }
  return(p.fit)
}


plot_bacon_age_models <- function(hamstr_bacon_fit, n.iter = 1000){
  
  
  # obs_ages <- hamstr_bacon_fit$info$dets 
  # 
  # if (info$cc != 0){
  #  obs_ages <- hamstr::calibrate_14C_age(obs_ages, age.14C = "age", age.14C.se = "error")
  #  obs_ages <- obs_ages %>%
  #   dplyr::select(-age, -error) %>% 
  #   dplyr::rename(.,
  #                 age = age.14C.cal,
  #                 err  = age.14C.cal.se) %>% 
  #   dplyr::mutate(age_upr = age + 2*err,
  #                 age_lwr = age - 2*err)
  # }  else{
  #   obs_ages <- obs_ages %>%
  #     dplyr::mutate(age_upr = age + 2*error,
  #                   age_lwr = age - 2*error)
  # }

  
  pdf.sums <- lapply(hamstr_bacon_fit$info$calib$probs,
                     function(x) SummariseEmpiricalPDF(x[,1], x[,2]))
 
  obs_ages <- dplyr::bind_rows(
    lapply(pdf.sums, function(x) c(age = x[["median"]], err = x[["sd"]]))
    )    %>%
    dplyr::mutate(age_upr = age + 2*err,
                  age_lwr = age - 2*err) %>% 
    dplyr::mutate(depth = hamstr_bacon_fit$info$calib$d)
  

  
  
  posterior_ages <- hamstr::return_bacon_age_mods(hamstr_bacon_fit)
  
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

plot_summary_bacon_age_models <- function(hamstr_bacon_fit){
  
 
  age_summary <- summarise_bacon_age_models(hamstr_bacon_fit)
  
 
  pdf.sums <- lapply(hamstr_bacon_fit$info$calib$probs,
                     function(x) SummariseEmpiricalPDF(x[,1], x[,2]))
  
  obs_ages <- dplyr::bind_rows(
    lapply(pdf.sums, function(x) c(age = x[["median"]], err = x[["sd"]]))
  )    %>%
    dplyr::mutate(age_upr = age + 2*err,
                  age_lwr = age - 2*err) %>% 
    dplyr::mutate(depth = hamstr_bacon_fit$info$calib$d)
  
  
  
  
  
  p.age.sum <- age_summary %>% 
    ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`), fill = "Lightgrey") +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`), fill = "Darkgrey") +
    ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = `50%`), colour = "Green") +
    ggplot2::labs(x = "Depth", y = "Age") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(data = obs_ages,
                            ggplot2::aes(x = depth, 
                                         ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                            colour = "Blue", size = 1.25) +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age),
                        colour = "Blue")
  
  p.age.sum
}

