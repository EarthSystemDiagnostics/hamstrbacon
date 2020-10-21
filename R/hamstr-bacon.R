#' Run Wrapped Bacon
#'
#' @param depth 
#' @param obs_age 
#' @param obs_err 
#' @param ... 
#'
#' @inheritParams rbacon::Bacon
#' @return
#' @export
#'
#' @examples
hamstr_bacon <- function(depth, obs_age, obs_err,
                         thick = NULL,
                         d.min = NA, d.max = NA, d.by = NA,
                         plot.pdf = FALSE){
  
  
  pars <- c(as.list(environment()))
  
  tmpdir <- tempdir()
  
  #unlink(tmpdir, recursive = TRUE, force = TRUE)
  
  dirbase <- basename(tmpdir)
  dirnm <- dirname(tmpdir)
  
  datfl <- tempfile(tmpdir = tmpdir)
  
   d <- dplyr::tibble(
      id = "test",
      age = obs_age,
      error = obs_err,
      depth = depth
    )
  
  utils::write.csv(d, file = paste0(tmpdir, "\\", dirbase, ".csv"),
              row.names = FALSE, quote = FALSE)
  
  
  # call Bacon
  
  pars$core = dirbase
  pars$coredir = normalizePath(paste0(dirnm, "/"))
  pars$ask = FALSE

  # subset of args for Bacon
  to.keep <- names(pars) %in% methods::formalArgs(rbacon::Bacon)
  pars <- pars[to.keep]
  
  # get defaults
  default.args <- formals(rbacon::Bacon)
  default.arg.nms <- names(default.args)
  
  # Overwrite the defaults with non-null passed arguments 
  dpars <- pars
  dpars <- dpars[lapply(dpars, is.null) == FALSE]
  
  default.args[names(dpars)] <- dpars
  
  dpars <- default.args
  
  # set d.by to thick unless specified. We will do interpolation separately
  if (is.na(d.by)){
    pars$d.by <- dpars$thick
  }
  
  do.call(rbacon::Bacon, pars)
  #do.call(envibacon::Bacon2, pars)
  
  settings.file <- utils::read.csv(file = paste0(tmpdir, "\\", dirbase, "_settings.txt"), header = FALSE)[,1]
  
  settings.file <- settings.file %>% 
    dplyr::as_tibble() %>% 
    tidyr::separate(value, into = c("value", "par"), sep = "#") %>% 
    dplyr::mutate(value = readr::parse_number(value)) %>% 
    dplyr::select(par, value)
  
  par_list <- as.list(settings.file$value)
  names(par_list) <- settings.file$par
  
  par_list$thick <- pars$thick
  
  d.by <- par_list[["d.by"]]
  d.min <- par_list[["d.min"]]
  d.max <- par_list[["d.max"]]

  
  # contruct the file name
  K = length(seq(floor(d.min), ceiling(d.max), by = dpars$thick))
 
  outflnm <- paste0(tmpdir, "\\", dirbase, "_", K, ".out")

  posterior <- utils::read.table(outflnm, header = FALSE)
  
  #age_mods <- return_age_mods(outflnm, thick = dpars$thick, d.min = d.min)

  out <- list(pars = par_list, posterior = posterior, data = d)
  
  class(out) <- append("hamstr_bacon_fit", class(out))
  
  return(out)
}


.StackIterations <- function(x){
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  n.row <- nrow(x)
  
  x <- utils::stack(x, stringsAsFactors = FALSE)
  x$depth.index <- (1:n.row) -1
  x
}


return_age_mods <- function(hamstr_bacon_fit){
  
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
  
  return(dplyr::as_tibble(age_mods))
  
}


interpolate_bacon_age_models <- function(hamstr_bacon_fit, new_depth){
  
  # get posterior age models
  pst_age <- hamstr::return_age_mods(hamstr_bacon_fit)
  
  new_age <- pst_age %>% 
    dplyr::group_by(iter) %>% 
    dplyr::do({
      tibble::tibble(
        iter = .$iter[1],
        depth = new_depth,
        age = stats::approx(.$depth, .$age, new_depth)$y
      )
    }) %>% 
    dplyr::ungroup()
  
  class(new_age) <- append(class(new_age), "hamstr_bacon_interpolated_ages")
  
  return(new_age)
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


summarise_bacon_age_models <- function(hamstr_bacon_fit){
  
  age_mods <- hamstr::return_age_mods(hamstr_bacon_fit)
  
  age_summary <- age_mods %>% 
    dplyr::group_by(depth) %>% 
    dplyr::summarise(mean = mean(age),
              sd = stats::sd(age),
              `2.5%` = stats::quantile(age, probs = 0.025),
              `25%` = stats::quantile(age, probs = 0.25),
              `50%` = stats::quantile(age, probs = 0.5),
              `75%` = stats::quantile(age, probs = 0.75),
              `97.5%` = stats::quantile(age, probs = 0.975),
              )
  
  return(age_summary)
}




plot_summary_bacon_age_models <- function(hamstr_bacon_fit){
  
 
  age_summary <- hamstr::summarise_bacon_age_models(hamstr_bacon_fit)
  
  obs_ages <- hamstr::calibrate_14C_age(hamstr_bacon_fit$data, age.14C = "age", age.14C.se = "error")
  
  obs_ages <- obs_ages %>%
    dplyr::select(-age, -error) %>% 
    dplyr::rename(.,
                  age = age.14C.cal,
                  err  = age.14C.cal.se) %>% 
    dplyr::mutate(age_upr = age + 2*err,
               age_lwr = age - 2*err)
  
  
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


