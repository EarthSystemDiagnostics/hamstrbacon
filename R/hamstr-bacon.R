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
  
  unlink(tmpdir, recursive = TRUE, force = TRUE)
  
  dirbase <- basename(tmpdir)
  dirnm <- dirname(tmpdir)
  
  datfl <- tempfile(tmpdir = tmpdir)
  
   d <- tibble(
      id = "test",
      age = obs_age,
      error = obs_err,
      depth = depth
    )
  
  write.csv(d, file = paste0(tmpdir, "\\", dirbase, ".csv"),
              row.names = FALSE, quote = FALSE)
  
  
  # call Bacon
  
  pars$core = dirbase
  pars$coredir = normalizePath(paste0(dirnm, "/"))
  pars$ask = FALSE

  # subset of args for Bacon
  to.keep <- names(pars) %in% formalArgs(rbacon::Bacon)
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
  
  settings.file <- read.csv(file = paste0(tmpdir, "\\", dirbase, "_settings.txt"), header = FALSE)[,1]
  
  settings.file <- settings.file %>% 
    as_tibble() %>% 
    separate(value, into = c("value", "par"), sep = "#") %>% 
    mutate(value = readr::parse_number(value)) %>% 
    select(par, value)
  
  par_list <- as.list(settings.file$value)
  names(par_list) <- settings.file$par
  
  par_list$thick <- pars$thick
  
  d.by <- par_list[["d.by"]]
  d.min <- par_list[["d.min"]]
  d.max <- par_list[["d.max"]]

  
  # contruct the file name
  K = length(seq(floor(d.min), ceiling(d.max), by = dpars$thick))
 
  outflnm <- paste0(tmpdir, "\\", dirbase, "_", K, ".out")

  posterior <- read.table(outflnm, header = FALSE)
  
  #age_mods <- return_age_mods(outflnm, thick = dpars$thick, d.min = d.min)

  out <- list(pars = par_list, posterior = posterior, data = d)
  
  class(out) <- append(class(out), "hamstr_bacon_fit")
  
  return(out)
}


.StackIterations <- function(x){
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  n.row <- nrow(x)
  
  x <- stack(x, stringsAsFactors = FALSE)
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
  
  return(as_tibble(age_mods))
  
}


interpolate_bacon_age_models <- function(hamstr_bacon_fit, new_depth){
  
  # get posterior age models
  pst_age <- return_age_mods(hamstr_bacon_fit)
  
  new_age <- pst_age %>% 
    dplyr::group_by(iter) %>% 
    dplyr::do({
      tibble::tibble(
        iter = .$iter[1],
        depth = new_depth,
        age = stats::approx(.$depth, .$age, new_depth)$y
      )
    }) %>% 
    ungroup()
  
  class(new_age) <- append(class(new_age), "hamstr_bacon_interpolated_ages")
  
  return(new_age)
}

summarise_bacon_age_models <- function(hamstr_bacon_fit){
  
  age_mods <- return_age_mods(hamstr_bacon_fit)
  
  age_summary <- age_mods %>% 
    group_by(depth) %>% 
    summarise(mean = mean(age),
              sd = sd(age),
              `2.5%` = quantile(age, probs = 0.025),
              `25%` = quantile(age, probs = 0.25),
              `50%` = quantile(age, probs = 0.5),
              `75%` = quantile(age, probs = 0.75),
              `97.5%` = quantile(age, probs = 0.975),
              )
  
  return(age_summary)
}




plot_summary_bacon_age_models <- function(hamstr_bacon_fit){
  
 
  age_summary <- summarise_bacon_age_models(hamstr_bacon_fit)
  
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

# MSB2K_cal <- hamstr::calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")
# ham1 <- hamstr(depth = MSB2K_cal$depth, obs_age = MSB2K_cal$age.14C.cal, obs_err = MSB2K_cal$age.14C.cal.se)
# 
# 
# hambac1 <- hamstr_bacon(depth = MSB2K$depth, obs_age = MSB2K$age, obs_err = MSB2K$error, thick = 3, d.min = 0)
# 
# return_age_mods(hambac1)
# summarise_bacon_age_models(hambac1)
# 
# plot_summary_bacon_age_models(hambac1)
# plot_summary_age_models(ham1)
# 
# hamstr::get_posterior_ages(ham1)
# 
# #rbacon:::.Bacon.settings
# 
# 
# #rbacon::Bacon(core = "RtmpU1bx5j", coredir = "C:\\Users\\Andrew\\AppData\\Local\\Temp\\")
# 
# ham_sum <- summarise_age_models(ham1)
# bac_sum <- summarise_bacon_age_models(hambac1)
# 
# hb_sum <- bind_rows(ham = ham_sum, bac = bac_sum, .id = "model")
# 
# p.age.sum <- hb_sum %>% 
#   ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
#   ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`), fill = "Lightgrey") +
#   ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`), fill = "Darkgrey") +
#   ggplot2::geom_line() +
#   ggplot2::geom_line(ggplot2::aes(y = `50%`), colour = "Green") +
#   ggplot2::labs(x = "Depth", y = "Age") +
#   ggplot2::theme_bw() +
#   ggplot2::theme(panel.grid = ggplot2::element_blank()) +
#   ggplot2::facet_wrap(~model)
# 
# p.age.sum <- p.age.sum +
#   ggplot2::geom_linerange(data = MSB2K_cal,
#                           ggplot2::aes(x = depth, 
#                                        ymax = age.14C.cal + age.14C.cal.se,
#                                        ymin = age.14C.cal - age.14C.cal.se),
#                           inherit.aes = FALSE,
#                           colour = "Blue", size = 1.25) +
#   ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age.14C.cal),
#                       colour = "Blue")
# 
# 
# ### compare coverage
# 
# join_mod_dat <- function(hamstr_fit, dat){
#   i1 <- hamstr::interpolate_age_models(hamstr_fit, dat$depth)
#   s1 <- hamstr::summarise_age_models(i1)
#   c1 <- left_join(dat, s1)
#   return(c1)
# }
# 
# c0 <- join_mod_dat(h0, acc)
# c1 <- join_mod_dat(h1, acc)
# c2 <- join_mod_dat(h2, acc)
# c2b <- join_mod_dat(h2b, acc)
# c3 <- join_mod_dat(h3, acc)
# 
# c0123 <- bind_rows(c0=c0, c1=c1, c2=c2, c2b=c2b, c3=c3, .id = "exp")
# 
# c0123 %>% 
#   ggplot(aes(x = depth, y = mean - cal_age, colour = exp)) +
#   geom_line() +
#   expand_limits(y = c(-100, 100))
# 
# c0123 %>% 
#   ggplot(aes(x = depth, y = mean - age.14C.cal, colour = exp)) +
#   geom_line()+
#   expand_limits(y = c(-100, 100))
# 
# AgeModCoverage <- function(hamstr_fit, dat){
#   comb <- join_mod_dat(hamstr_fit, dat)
#   
#   comb %>% 
#     filter(complete.cases(`75%`)) %>% 
#     mutate(c95 = (cal_age >= `2.5%`) & (cal_age <= `97.5%`),
#            c50 = (cal_age >= `25%`) & (cal_age <= `75%`)) %>% 
#     #select(c95)
#     summarise(c95 = sum(c95) / n(),
#               c50 = sum(c50) / n(),
#               RMSE = sqrt(mean((cal_age - mean)^2)))
#   
# }
# 
# AgeModCoverage(h0, acc)
# AgeModCoverage(h1, acc)
# AgeModCoverage(h2, acc)
# AgeModCoverage(h2b, acc)
# AgeModCoverage(h3, acc)
# 

