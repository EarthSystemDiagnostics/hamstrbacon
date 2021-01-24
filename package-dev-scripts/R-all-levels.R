#' Make the data object required by the Stan program
#'
#' @inheritParams hamstr
#'
#' @return a list of data and parameters to be passed as data to the Stan sampler
#' @export
#'
#' @examples
#' make_stan_dat_hamstr(depth = MSB2K$depth,
#'               obs_age = MSB2K$age,
#'               obs_err = MSB2K$error)
make_stan_dat_hamstr_dev <- function(depth=NULL, obs_age=NULL, obs_err=NULL,
                                 top_depth=NULL, bottom_depth=NULL,
                                 pad_top_bottom=NULL,
                                 K=NULL, nu=NULL,
                                 acc_mean_prior=NULL,
                                 shape=NULL,
                                 mem_mean=NULL, mem_strength=NULL,
                                 #scale_R=NULL,
                                 inflate_errors=NULL,
                                 infl_sigma_sd=NULL, 
                                 infl_shape_shape=NULL, infl_shape_mean=NULL,
                                 ...) {
  
  l <- c(as.list(environment()))
  
  # get defaults
  default.args <- formals(hamstr)
  default.arg.nms <- names(default.args)
  
  
  # Overwrite the defaults with non-null passed arguments 
  
  l <- l[lapply(l, is.null) == FALSE]
  
  default.args[names(l)] <- l
  
  l <- default.args
  
  if (is.null(l$acc_mean_prior)){
    
    d <- data.frame(depth = depth, obs_age = obs_age)
    acc_mean <- stats::coef(MASS::rlm(obs_age~depth, data = d))[2]
    
    # if negative replace with 20
    if (acc_mean <= 0) {
      warning("Estimated mean accumulation rate is negative - using value = 20")
      acc_mean <- 20
    }
    l$acc_mean_prior <- acc_mean
  }
  
  
  ord <- order(depth)
  
  l$depth <- depth[ord]
  l$obs_age <- obs_age[ord]
  l$obs_err <- obs_err[ord]
  
  
  if (is.null(infl_sigma_sd)){
    l$infl_sigma_sd <- 10 * mean(obs_err)
  }
  
  
  if (l$pad_top_bottom == TRUE){
    # Set start depth to 5% less than first depth observation, and DO allow negative depths
    depth_range <- diff(range(l$depth))
    buff <- 0.05 * depth_range
  } else {
    buff <- 0
  }
  
  if (is.null(top_depth)) l$top_depth <- l$depth[1] - buff
  
  if (is.null(bottom_depth)) l$bottom_depth <- utils::tail(l$depth, 1) + buff
  
  depth_range <- l$bottom_depth - l$top_depth
  
  if(l$top_depth > min(l$depth)) stop("top_depth must be above or equal to the shallowest data point")
  if(l$bottom_depth < max(l$depth)) stop("bottom_depth must be deeper or equal to the deepest data point")
  
  
  # Transformed arguments
  l$N <- length(l$depth)
  
  stopifnot(l$N == length(obs_err), l$N == length(obs_age))
  
  alpha_idx <- alpha_indices_dev(l$K)
  
  l$K_tot <- sum(alpha_idx$nK)
  l$K_fine <- utils::tail(alpha_idx$nK, 1)
  l$c <- 1:l$K_fine
  
  l$mem_alpha = mem_strength * mem_mean
  l$mem_beta = mem_strength * (1-mem_mean)
  
  l$mem_mean = mem_mean
  l$mem_strength = mem_strength
  
  l$delta_c = depth_range / alpha_idx$nK
  l$c_depth_bottom = tail(l$delta_c, 1) * l$c + l$top_depth
  l$c_depth_top = c(l$top_depth, l$c_depth_bottom[1:(l$K_fine-1)])
  
  l$modelled_depths <- c(l$c_depth_top[1], l$c_depth_bottom)
  
  # Index for which sections the target depth is in
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
  
  l <- append(l, alpha_idx)
  
  return(l)
}

# Make index creating functions for K levels
#' Create alpha level indices
#'
#' @inheritParams hamstr
#'
#' @return a list
#' @keywords internal
alpha_indices_dev <- function(K){
  
  K <- eval(K)
  
  # prepend 1 for the single overall mean alpha 
  K <- c(1, K)
  
  # number of sections at each level
  nK <- cumprod(K)
  
  K_lvls <- length(K)
  
  K <- c(0, K)
  nK <- c(0, nK)
  
  alpha_idx <- 1:sum(nK)
  
  lvl_SE <- cbind(head(cumsum(nK), -1)+1, tail(cumsum(nK), -1))
  
  finest_idx <- lvl_SE[K_lvls, 1]:lvl_SE[K_lvls, 2]
  
  nLevels <- length(K)-1
  
  # which level is each parameter
  lvl <- unlist(lapply(seq_along(nK[-1]), function(i) rep(i, times = nK[i+1])))
  
  # index the parent of each parameter
  parent <- c(rep(0, K[2]), unlist(lapply(alpha_idx[1:sum(nK[1:nLevels])], function(i) rep(i, K[lvl[i]+2]))))
  
  list(alpha_idx=alpha_idx, lvl=lvl, lvl_SE=lvl_SE, finest_idx=finest_idx, parent=parent, nK = nK[-1], K_lvls=K_lvls)
}



sd1 <- make_stan_dat_hamstr(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
                         K = c(10, 10), mem_mean = 0.7, mem_strength = 4,
                         cores = 3)


sd2 <- make_stan_dat_hamstr_dev(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
                         K = c(10, 10, 5), mem_mean = 0.7, mem_strength = 4,
                         cores = 3)




### attempt to compile hierarchical model with AR1 at all levels

hd0 <- hamstr(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
                         K = c(10, 10), mem_mean = 0.7, mem_strength = 4,
                         cores = 3)

hd0b <- hamstr_bacon(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
                     mem.strength = 4, mem.mean = 0.7, acc.mean = 50,
                     thick = 4.8)

hd0b2 <- hamstr_bacon(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
                     mem.strength = 4, mem.mean = 0.7, acc.mean = 50,
                     thick = 48)

hd0b3 <- hamstr_bacon(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
                      mem.strength = 4, mem.mean = 0.7, acc.mean = 50,
                      thick = 1)


hd2 <- stan(file = "package-dev-scripts/stan-models/hamstr_R_all.stan", data = sd2, cores = 3, chains = 3)

hd2 <- list(fit=hd2, data=sd2)

class(hd2) <- append("hamstr_fit", class(hd2))


s1 <- as_tibble(summary(hd2$fit)$summary, rownames = "parameter")
hist(s1$Rhat)


plot(hd0b)
plot(hd0b2)
plot(hd0b3)

plot(hd0, type = "age_models")
plot(hd2, type = "age_models")

plot(hd2)


s1 %>% 
  filter(parameter %in% c("a", "R", "acc_mean"))

s1 %>% 
  filter(grepl("w", parameter))








