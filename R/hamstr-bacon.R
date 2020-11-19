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
                         d.min = NA, d.max = NA, d.by = NA,
                         acc.shape = 1.5, acc.mean = 20,
                         mem.strength = 4, mem.mean = 0.7,
                         plot.pdf = FALSE,
                         ask = FALSE, suggest = FALSE){
  
 
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
  if (is.na(d.by)){
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


#' Title
#'
#' @param object a hamstr_bacon_fit object
#' @return A ggplot object
#'
#' @examples
#' @export
#' @method plot hamstr_bacon_fit
plot.hamstr_bacon_fit <- function(object,
                            type = c("default"),
                            summarise = TRUE,
                            ...){
  
  type <- match.arg(type)
  
  switch(type,
         default = plot_summary_bacon_age_models(object, ...)
         )
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


#' Bacon2, A Wrapper for rbacon::Bacon so that it plays nicely when run in
#' parallel. Suppresses calls to graphics devices.
#'
#' @param suppress.plots If TRUE it prevents calls to on screen graphics devices.
#' @inheritParams rbacon::Bacon
#' @param ...
#'
#' @return
#' @import rbacon
#' @keywords internal
#'
#' @examples
Bacon2 <- function (suppress.plots = TRUE,
                    core = "MSB2K", thick = 5, coredir = "", prob = 0.95,
                    d.min = NA, d.max = NA, d.by = 1, seed = NA, depths.file = FALSE,
                    depths = c(), depth.unit = "cm", age.unit = "yr", unit = depth.unit,
                    acc.shape = 1.5, acc.mean = 20, mem.strength = 4, mem.mean = 0.7,
                    boundary = NA, hiatus.depths = NA, hiatus.max = 10000, add = c(),
                    after = 1e-04/thick, cc = 1, cc1 = "IntCal20", cc2 = "Marine20",
                    cc3 = "SHCal20", cc4 = "ConstCal", ccdir = "", postbomb = 0,
                    delta.R = 0, delta.STD = 0, t.a = 3, t.b = 4, normal = FALSE,
                    suggest = c("accept", FALSE, TRUE), reswarn = c(10, 200), remember = TRUE, ask = FALSE,
                    run = TRUE, defaults = "defaultBacon_settings.txt", sep = ",",
                    dec = ".", runname = "", slump = c(), BCAD = FALSE, ssize = 2000,
                    th0 = c(), burnin = min(500, ssize), MinAge = c(), MaxAge = c(),
                    MinYr = MinAge, MaxYr = MaxAge, cutoff = 0.001, plot.pdf = TRUE,
                    dark = 1, date.res = 100, age.res = 200, yr.res = age.res,
                    close.connections = TRUE, verbose = TRUE, ...)
{
  
 
  suggest <- match.arg(as.character(suggest), 
                       choices = c("accept", "FALSE", "TRUE"))
  
  # Temporarily attach all internal functions from package rbacon
  attach(loadNamespace("rbacon"), name = "rbacon_all")
  
  coredir <- assign_coredir(coredir, core, ask)
  if (core == "MSB2K" || core == "RLGH3") {
    dir.create(paste(coredir, core, "/", sep = ""), showWarnings = FALSE,
               recursive = TRUE)
    fileCopy <- system.file(paste("extdata/Cores/", core,
                                  sep = ""), package = "rbacon")
    file.copy(fileCopy, coredir, recursive = TRUE, overwrite = FALSE)
  }
  if (ccdir == "") 
    ccdir <- system.file("extdata", package = "IntCal")
  ccdir <- .validateDirectoryName(ccdir)
  defaults <- system.file("extdata", defaults, package = "rbacon")
  dets <- .read.dets(core, coredir, sep = sep, dec = dec,
                     cc = cc)
  if (ncol(dets) > 4 && length(cc) > 0) {
    cc.csv <- unique(dets[, 5])
    if (verbose) {
      if (length(cc.csv) == 1) {
        if (cc.csv != cc)
          message(" Using calibration curve specified within the .csv file,",
                  cc[cc.csv], "\n")
      }
      else if (min(cc.csv) == 0)
        message(" Using a mix of cal BP and calibrated C-14 dates\n")
      else message(" Using several C-14 calibration curves\n")
    }
  }
  if (suggest == "TRUE") {
    sugg <- sapply(c(1, 2, 5), function(x) x * 10^(-1:2))
    ballpacc <- lm(dets[, 2] * 1.1 ~ dets[, 4])$coefficients[2]
    ballpacc <- abs(sugg - ballpacc)
    ballpacc <- ballpacc[ballpacc > 0]
    sugg <- sugg[order(ballpacc)[1]]
    if (!sugg %in% acc.mean) {
      ans <- readline(message(" Ballpark estimates suggest changing the prior for acc.mean to ",
                              sugg, " ", age.unit, "/", depth.unit, ". OK? (y/N) "))
      if (tolower(substr(ans, 1, 1)) == "y")
        acc.mean <- sugg
      else message(" No problem, using the provided prior")
    }
  } else if (suggest == "accept"){
    sugg <- sapply(c(1, 2, 5), function(x) x * 10^(-1:2))
    ballpacc <- lm(dets[, 2] * 1.1 ~ dets[, 4])$coefficients[2]
    ballpacc <- abs(sugg - ballpacc)
    ballpacc <- ballpacc[ballpacc > 0]
    sugg <- sugg[order(ballpacc)[1]]
    
    message("Setting acc.mean to Ballpark estimate of ", sugg)
    acc.mean <- sugg
    
    }
  if (!is.na(boundary[1]))
    boundary <- sort(unique(boundary))
  if (!is.na(hiatus.depths[1])) {
    hiatus.depths <- sort(unique(hiatus.depths))
    if (length(acc.mean) == 1)
      acc.mean <- rep(acc.mean, length(hiatus.depths) +
                        1)
  }
  info <- .Bacon.settings(core = core, coredir = coredir,
                          dets = dets, thick = thick, remember = remember, d.min = d.min,
                          d.max = d.max, d.by = d.by, depths.file = depths.file,
                          slump = slump, acc.mean = acc.mean, acc.shape = acc.shape,
                          mem.mean = mem.mean, mem.strength = mem.strength, boundary = boundary,
                          hiatus.depths = hiatus.depths, hiatus.max = hiatus.max,
                          BCAD = BCAD, cc = cc, postbomb = postbomb, cc1 = cc1,
                          cc2 = cc2, cc3 = cc3, cc4 = cc4, depth.unit = depth.unit,
                          normal = normal, t.a = t.a, t.b = t.b, delta.R = delta.R,
                          delta.STD = delta.STD, prob = prob, defaults = defaults,
                          runname = runname, ssize = ssize, dark = dark, MinAge = MinAge,
                          MaxAge = MaxAge, cutoff = cutoff, age.res = age.res,
                          after = after, age.unit = age.unit)
  .assign_to_global("info", info)
  info$coredir <- coredir
  info$seed <- seed
  info$isplum <- FALSE
  if (any(info$acc.shape == info$acc.mean))
    stop("acc.shape cannot be equal to acc.mean", call. = FALSE)
  if (info$t.b - info$t.a != 1)
    stop("t.b - t.a should always be 1, check the manual",
         call. = FALSE)
  if (min(acc.shape) < 1)
    cat("\nWarning, using values <1 for acc.shape might cause unexpected results\n")
  if (info$cc > 0)
    if (info$postbomb == 0 && ((ncol(info$dets) == 4 &&
                                min(info$dets[, 2]) < 0) || ncol(info$dets) > 4 &&
                               max(info$dets[, 5]) > 0 && min(info$dets[info$dets[,
                                                                                  5] > 0, 2]) < 0))
      stop("you have negative C14 ages so should select a postbomb curve",
           call. = FALSE)
  info$calib <- .bacon.calib(dets, info, date.res, ccdir = ccdir)
  info$rng <- c()
  for (i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[which(tmp[, 2] > cutoff),
                                    1])
  }
  if (length(th0) == 0)
    info$th0 <- round(rnorm(2, max(MinAge, dets[1, 2]),
                            dets[1, 3]))
  info$th0[info$th0 < info$MinAge] <- info$MinAge
  if (length(depths) == 0)
    depths <- seq(info$d.min, info$d.max, by = d.by)
  if (depths.file) {
    dfile <- paste0(info$coredir, info$core, "/", info$core,
                    "_depths.txt")
    if (!file.exists(dfile))
      stop("I cannot find the file ", paste0(info$coredir,
                                             info$core, "/", info$core, "_depths.txt"), call. = FALSE)
    depths <- read.table(dfile, header = FALSE)[, 1]
    if (!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers",
           call. = FALSE)
  }
  info$depths <- depths
  if (min(depths) < info$d.min)
    info$d.min <- min(depths)
  if (max(depths) > info$d.max)
    info$d.max <- max(depths)
  info$elbows <- seq(floor(info$d.min), ceiling(info$d.max),
                     by = thick)
  info$K <- length(info$elbows)
  info$cK <- info$d.min + (info$thick * info$K)
  if (!is.na(info$hiatus.depths[1]) || !is.na(info$boundary[1])) {
    ifelse(is.na(info$boundary[1]), hd <- info$hiatus.depths,
           hd <- info$boundary)
    if (min(hd) < info$d.min)
      stop("cannot have hiatus above the core's top depth. Adapt hiatus.depths or d.min.",
           call. = FALSE)
    if (max(hd) + info$thick > info$d.max)
      stop("the age-depth model should have at least one section below the one containing the deepest hiatus. Adapt thick or d.max?",
           call. = FALSE)
    if (length(hd) > 1) {
      above <- c()
      for (i in hd) above <- c(above, max(which(info$elbows <=
                                                  i)))
      if (any(diff(above) < 2))
        stop("we need at least 2 section elbows between hiatuses. Choose fewer hiatuses, different depths, more sections (decrease thick) or a different d.min.\n ",
             call. = FALSE)
    }
  }
  
  ans <- "n"
  if (suggest == "TRUE"){
    if (length(reswarn) == 2){
      if (info$K < min(reswarn)) {
        sugg <- pretty(thick * (info$K/min(reswarn)),
                       10)
        sugg <- min(sugg[sugg > 0])
        ans <- readline(message(" Warning, the current value for thick, ",
                                thick, ", will result in very few age-model sections (",
                                info$K, ", not very flexible). Suggested maximum value for thick: ",
                                sugg, " OK? (y/n) "))
      }
  else if (info$K > max(reswarn)) {
    sugg <- max(pretty(thick * (info$K/max(reswarn))))
    ans <- readline(message(" Warning, the current value for thick, ",
                            thick, ", will result in very many age-model sections (",
                            info$K, ", possibly hard to run). Suggested minimum value for thick: ",
                            sugg, " OK? (y/n) "))
  }}}
  if (suggest == "accept"){
      if (length(reswarn) == 2)
        if (info$K < min(reswarn)) {
          sugg <- pretty(thick * (info$K/min(reswarn)),
                         10)
          sugg <- min(sugg[sugg > 0])
          ans <- "y"
        } else if (info$K > max(reswarn)) {
      sugg <- max(pretty(thick * (info$K/max(reswarn))))
      ans <- "y"
    }
   
  }
  
  if (tolower(substr(ans, 1, 1)) == "y") {
    message(" OK, setting thick to ", sugg, "\n")
    thick <- sugg
    info$thick = thick
    info$d.by = thick
    info$elbows <- seq(floor(info$d.min), ceiling(info$d.max),
                       by = thick)
    if (length(info$slump) > 0)
      info$elbows <- seq(floor(info$d.min), toslump(ceiling(info$d.max),
                                                    info$slump), by = thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min + (info$thick * info$K)
  }

 
  
  
  if (length(slump) > 0) {
    if (length(slump)%%2 == 1)
      stop("slumps need both upper and lower depths. Please check the manual",
           call. = FALSE)
    slump <- matrix(sort(slump), ncol = 2, byrow = TRUE)
    info$slump <- slump
    slumpdmax <- toslump(ceiling(info$d.max), slump)
    info$elbows <- seq(floor(info$d.min), slumpdmax, by = thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min + (info$thick * info$K)
    info$slumpfree <- toslump(depths, slump)
    info$slumphiatus <- toslump(info$hiatus.depths, slump)
    if (!is.na(info$boundary[1])) {
      info$slumpboundary <- toslump(info$boundary, slump)
      info$slumphiatus <- info$slumpboundary
    }
    slumpdets <- info$dets
    slumpdets[, 4] <- toslump(slumpdets[, 4], slump, remove = FALSE)
    info$slumpdets <- slumpdets[!is.na(slumpdets[, 4]),
    ]
  }
  info$prefix <- paste(coredir, core, "/", core, runname,
                       "_", info$K, sep = "")
  info$coredir <- coredir
  info$bacon.file <- paste(info$prefix, ".bacon", sep = "")
  if (!file.exists(outfile <- paste(info$prefix, ".out", sep = "")))
    file.create(outfile)
  if (BCAD)
    info$BCAD <- TRUE
  if (!is.na(boundary[1])) {
    if (length(slump) > 0)
      boundary <- info$slumpboundary
    info$hiatus.depths <- boundary
    if (length(add) == 0)
      add <- info$acc.mean
    info$hiatus.max <- add
  }
  .assign_to_global("info", info)
  prepare <- function() {
    if (suppress.plots == FALSE){
      pn <- c(1, 2, 3, 3)
      if (!is.na(info$hiatus.depths[1]))
        if (is.na(info$boundary[1]))
          pn <- c(1, 2, 3, 4, 4, 4)
        layout(matrix(pn, nrow = 2, byrow = TRUE), heights = c(0.3,
                                                               0.7))
        oldpar <- par(mar = c(3, 3, 1, 1), mgp = c(1.5, 0.7,
                                                   0), bty = "l")
        on.exit(par(oldpar))
        .PlotAccPrior(info$acc.shape, info$acc.mean, depth.unit = depth.unit,
                      age.unit = age.unit)
        .PlotMemPrior(info$mem.strength, info$mem.mean, thick)
        if (!is.na(info$hiatus.depths)[1])
          if (is.na(info$boundary)[1])
            .PlotHiatusPrior(info$hiatus.max, info$hiatus.depths)
        calib.plot(info, BCAD = BCAD)
        legend("topleft", core, bty = "n", cex = 1.5)
    }
  }
  cook <- function() {
    txt <- paste(info$prefix, ".bacon", sep = "")
    bacon(txt, outfile, ssize, ccdir)
    scissors(burnin, info)
    if (suppress.plots == FALSE){
      agedepth(info, BCAD = BCAD, depths.file = depths.file,
               depths = depths, verbose = TRUE, age.unit = age.unit,
               depth.unit = depth.unit, ...)
    }
    if (plot.pdf){
      pdf(file = paste0(info$prefix, ".pdf"))
      agedepth(info, BCAD = BCAD, depths.file = depths.file,
               depths = depths, verbose = FALSE, age.unit = age.unit,
               depth.unit = depth.unit, ...)
      dev.off()
    }
  }
  .write.Bacon.file(info)
  if (!run)
    prepare()
  else if (!ask)
    cook()
  else {
    prepare()
    ans <- readline(message(" Run ", core, " with ", info$K,
                            " sections? (Y/n) "))
    ans <- tolower(substr(ans, 1, 1))[1]
    if (ans == "y" || ans == "")
      cook()
    else message("  OK. Please adapt settings")
  }
  if (close.connections)
    closeAllConnections()
  
  # detach internal rbacon functions
  detach("rbacon_all")
  
  # cleanup global
  rm(info)
}


