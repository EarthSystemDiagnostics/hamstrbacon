#devtools::install(quick = FALSE, dependencies = FALSE)
# install.packages("../hamstr", repos = NULL, type = "source")
library(hamstr)
library(tidyverse)

MSB2K_cal <- hamstr::calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")


###
# simulate correlated gamma?

SimulateAgeDepth <- function(top, bottom, d_depth, gamma_shape,
                             gamma_mean, gamma_breaks = NULL,
                             ar1){

  pars <- c(as.list(environment()))


  depth <- seq(top, bottom, by = d_depth)

  if (is.null(gamma_breaks) == FALSE){
    if(length(gamma_mean) - length(gamma_breaks) != 1)
      stop("gamma_mean should be 1 longer than gamma_breaks")

    foo <- stepfun(gamma_breaks, gamma_mean)

    gamma_mean <- foo(depth[-1])
  }


  n <- length(depth)

  gamma_scale = gamma_mean / gamma_shape

  y <- as.numeric(arima.sim(model = list(ar = ar1^d_depth), n-1))
  y <- y / sd(y)
  y <- y - mean(y)
  sd_y <- sd(y)

  py <- (pnorm(y, mean = 0, sd = sd_y))

  # translate to gamma
  acc.rates_yr_cm <- qgamma(py, shape = gamma_shape, scale = gamma_scale)

  min.age <- gamma_mean[1] * top

  age <- cumsum(c(min.age, acc.rates_yr_cm) * d_depth)


  core <- dplyr::tibble(depth = depth, age = age, acc.rates_yr_cm = c(acc.rates_yr_cm, NA))

  return(list(pars = pars, core = core))

  }



tmp <- SimulateAgeDepth(top = 100, bottom = 700, d_depth = 1.1,
                        gamma_shape = 1.5,
                        gamma_mean = c(50, 10, 50), gamma_breaks = c(200, 600),
                        ar1 = 0.7)

tmp$core %>%
  ggplot(aes(depth, age)) +
  geom_line()



SimulateCarbonDates <- function(age_depth_obj, sampling_interval, cal_curve,
                                sample_gap = NULL){

  pars <- c(as.list(environment()))

  obs <- age_depth_obj$core %>%
    mutate(sampled = (rank(depth) %% sampling_interval == 0)) %>%
    filter(sampled) %>%
    mutate(age_14C = Bchron::unCalibrate(age, type = "ages")) %>%
    mutate(age_14C_sigma = 20 + age * 0.02) %>%
    mutate(age_14C_hat = rnorm(n(), age_14C, age_14C_sigma))

  if (is.null(sample_gap) == FALSE){
    obs <- obs %>%
      filter(depth < sample_gap[1] | depth > sample_gap[2])
  }

  obs <- obs %>%
    select(depth, age_14C_hat, age_14C_sigma) %>%
    calibrate_14C_age(.,
                      "age_14C_hat", "age_14C_sigma",
                      cal_curve = cal_curve)


  return(list(obs = obs,
              core = age_depth_obj$core,
              pars = pars))

}

## should be comparing coverage of simulated "true" ages not the data points
## as data might have extra error added
join_mod_dat <- function(hamstr_fit, dat){

  #browser()
  i1 <- predict(hamstr_fit, dat$depth)
  s1 <- summary(i1)

  # try a bind rather than join, should be much faster.
  c1 <- dplyr::left_join(dat, s1)

  #c1 <- bind_cols(dat, s1)


  return(c1)
}


AgeModCoverage <- function(hamstr_fit, dat){
  comb <- join_mod_dat(hamstr_fit, dat)

  comb %>%
    filter(complete.cases(`75%`)) %>%
    mutate(c95 = (age >= `2.5%`) & (age <= `97.5%`),
           c50 = (age >= `25%`) & (age <= `75%`)) %>%
    #select(c95)
    summarise(c95 = sum(c95) / n(),
              c50 = sum(c50) / n(),
              RMSE = sqrt(mean((age - `50%`)^2)))

}


### Compare hamstr and Bacon

CompareHamBac <- function(age_depth_obj,
                          ar_coefs, acc_mean,
                          mem_mean, mem_strength, acc_shape,
                          scale_shape = TRUE,
                          K_hamstr, K_bacon,
                          inflate_errors = FALSE){



  ham1 <-  hamstr(depth = age_depth_obj$obs$depth, obs_age = age_depth_obj$obs$age.14C.cal,
                  obs_err = age_depth_obj$obs$age.14C.cal.se,
                  top_depth = min(age_depth_obj$core$depth),
                  bottom_depth = max(age_depth_obj$core$depth),
                  K = K_hamstr,
                  acc_shape = acc_shape,
                  mem_mean = mem_mean, mem_strength = mem_strength,
                  inflate_errors = inflate_errors,
                  scale_shape = scale_shape,
                  chains = 3, cores = 3)


  bac1 <- hamstr_bacon(depth = age_depth_obj$obs$depth,
                       obs_age = age_depth_obj$obs$age_14C_hat,
                       obs_err = age_depth_obj$obs$age_14C_sigma,
                       acc.mean = acc_mean, acc.shape = acc_shape,
                       mem.mean = mem_mean, mem.strength = mem_strength,
                       d.min = min(age_depth_obj$core$depth),
                       d.max = max(age_depth_obj$core$depth),
                       thick = diff(range(age_depth_obj$core$depth)) / K_bacon,
                       suggest = "FALSE")

  return(list(age_depth_obj = age_depth_obj,
              hamstr = ham1, bacon = bac1
  ))
}



GetCoverages <- function(sim){

  hamstr = AgeModCoverage(sim$hamstr, sim$age_depth_obj$core)
  bacon = AgeModCoverage(sim$bacon, sim$age_depth_obj$core)

  bind_rows(hamstr=hamstr, bacon=bacon, .id = "method")

}

PlotHamstrBacon <- function(sim){

  obs_ages <- hamstr::calibrate_14C_age(sim$bacon$data, age.14C = "age", age.14C.se = "error")

  obs_ages <- obs_ages %>%
    dplyr::select(-age, -error) %>%
    dplyr::rename(.,
                  age = age.14C.cal,
                  err  = age.14C.cal.se) %>%
    dplyr::mutate(age_upr = age + 2*err,
                  age_lwr = age - 2*err)

  bacon_age_summary <- hamstr:::summarise_bacon_age_models(sim$bacon)

  hamstr_age_summary <- hamstr:::summarise_age_models(sim$hamstr)

  alpha.lvl <- 0.5
  #pal <- c("Hamstr" = "Blue", "Bacon" = "Red", "Simulated" = "Yellow")
  pal <- c("Hamstr" = "#7570b3", "Bacon" = "#e7298a", "Simulated" = "green")

  p.age.sum <- bacon_age_summary %>%
    ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`, fill = "Bacon"), alpha = alpha.lvl/1.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`, fill = "Bacon"), alpha = alpha.lvl) +
    #ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = `50%`, colour = "Bacon"), lwd = 1) +
    ggplot2::labs(x = "Depth", y = "Age") #+
  #ggplot2::theme(legend.key = element_rect(fill  = "grey40", colour = "grey40"))#+
  #ggplot2::theme_bw()# +
  #ggplot2::theme(panel.grid = ggplot2::element_blank()) +


  p.age.sum <- p.age.sum +
    ggplot2::geom_ribbon(data = hamstr_age_summary,
                         ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`, fill = "Hamstr"), alpha = alpha.lvl/1.5) +
    ggplot2::geom_ribbon(data = hamstr_age_summary,
                         ggplot2::aes(ymax = `75%`, ymin = `25%`, fill = "Hamstr"), alpha = alpha.lvl) +
    #ggplot2::geom_line(data = hamstr_age_summary ) +
    ggplot2::geom_line(data = hamstr_age_summary, ggplot2::aes(y = `50%`, colour = "Hamstr"), lwd = 1)

  p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(data = obs_ages,
                            ggplot2::aes(x = depth,
                                         ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                            colour = "green") +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age),
                        colour = "green") +
    geom_line(data = sim$age_depth_obj$core, aes(x = depth, y = age, colour = "Simulated"), lwd = 1)+

    scale_fill_manual("", values = pal) +
    scale_color_manual("", values = pal)


  p.age.sum

}

CompareMedAgeMod <- function(sim){

  hamstr = join_mod_dat(sim$hamstr, sim$sim_core)
  bacon = join_mod_dat(sim$bacon, sim$sim_core)

  hb <- bind_rows(hamstr=hamstr, bacon=bacon, .id = "method") %>%
    mutate(error = `50%` - age)

  p <- hb %>%
    ggplot(aes(x = depth, y = `50%` - age, colour = method)) +
    geom_line() +
    geom_hline(yintercept = 0)

  return(list(data = hb, plot = p))

}




ad.obj <- SimulateAgeDepth(top = 100, bottom = 700, d_depth = 1.1,
                           gamma_shape = 1.5,
                           gamma_mean = c(50, 10, 50), gamma_breaks = c(200, 600),
                           ar1 = 0.5)

tmp <- SimulateCarbonDates(ad.obj, sampling_interval = 24, sample_gap = c(201, 349),
                    cal_curve = "intcal20")

tst1 <- CompareHamBac(tmp,
                       K_hamstr = c(10,10),
                       K_bacon = 100,
                      scale_shape = FALSE,
                       acc_shape = 1.5,
                      acc_mean = mean(tmp$pars$age_depth_obj$pars$gamma_mean),
                       mem_mean = 0.5, mem_strength = 10)

tst2 <- CompareHamBac(tmp,
                      K_hamstr = c(10,10),
                      K_bacon = 100,
                      scale_shape = TRUE,
                      acc_shape = 1.5,
                      acc_mean = mean(tmp$pars$age_depth_obj$pars$gamma_mean),
                      mem_mean = 0.5, mem_strength = 10)


PlotHamstrBacon(tst1)
GetCoverages(tst1)

PlotHamstrBacon(tst2)
GetCoverages(tst2)




## Simulations ------

# hbc1 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
#               sampling_interval = 12)
#
# hbc2 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
#                       sampling_interval = 1)

hbc3a <- CompareHamBac(100, 400, 1, sim_gamma_shape = 1.5,
                      sim_gamma_mean = 80,
                      ar_coefs = 0.5,
                      sampling_interval = 12,
                      sample_gap = c(201, 349),
                      #K_hamstr = hamstr:::optimal_K(100, 10),
                      K_hamstr = 100,
                      K_bacon = 100,
                      acc_shape = 1.5,
                      mem_mean = 0.5, mem_strength = 10)

hbc3b <- CompareHamBac(100, 400, 1,
                       sim_gamma_shape = 1.5,
                       sim_gamma_mean = 40,
                      ar_coefs = 0.5,
                      sampling_interval = 12,
                      sample_gap = c(201, 349),
                      K_hamstr = hamstr:::optimal_K(100, 10),
                      K_bacon = 100,
                      acc_shape = 1.5,
                      mem_mean = 0.8, mem_strength = 10)

hbc3c <- CompareHamBac(100, 400, 1, gamma_shape = 1.5,
                      gamma_mean = 40,
                      ar_coefs = 0.5,
                      sampling_interval = 12,
                      sample_gap = c(201, 349),
                      K_hamstr = hamstr:::optimal_K(100, 10),
                      K_bacon = 100,
                      acc_shape = 1.5,
                      mem_mean = 0.2, mem_strength = 10)


GetCoverages(hbc3a)
GetCoverages(hbc3b)
GetCoverages(hbc3c)

PlotHamstrBacon(hbc3a)
PlotHamstrBacon(hbc3b)
PlotHamstrBacon(hbc3c)

hbc4 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5,
                      gamma_mean =  c(20, 90), gamma_breaks = 250,
                      ar_coefs = 0.7,
                      sampling_interval = 24,
                      sample_gap = c(225, 300),
                      K_hamstr = hamstr:::optimal_K(100, 10),
                      K_bacon = 100)

hbc5 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5,
                      gamma_mean =  c(90, 20), gamma_breaks = 250,
                      ar_coefs = 0.7,
                      sampling_interval = 24,
                      sample_gap = c(225, 300),
                      K_hamstr = hamstr:::optimal_K(100, 10),
                      K_bacon = 100)

hbc6 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5,
                      gamma_mean =  c(90, 20), gamma_breaks = 250,
                      ar_coefs = 0.7,
                      sampling_interval = 24,
                      sample_gap = c(225, 300),
                      K_hamstr = hamstr:::optimal_K(100, 10),
                      K_bacon = 100,
                      inflate_errors = TRUE)



GetCoverages(hbc3)
GetCoverages(hbc4)
GetCoverages(hbc5)
GetCoverages(hbc6)


PlotHamstrBacon(hbc3)
PlotHamstrBacon(hbc4)
PlotHamstrBacon(hbc5)
PlotHamstrBacon(hbc6)

CompareMedAgeMod(hbc4)

plot(hbc6$hamstr)

comps <- list(hbc1, hbc2, hbc3, hbc4, hbc5, hbc6)

covr <- plyr::ldply(comps, GetCoverages, .progress = TRUE)

covr %>%
  gather(percentile, coverage, -RMSE, -method) %>%
  mutate(percentile_num = readr::parse_number(percentile)) %>%
  ggplot(aes(x = percentile, y = coverage, colour = method)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_hline(aes(yintercept = percentile_num/100)) +
  facet_wrap(~percentile, scales = "free_y")

if (parallel::detectCores() >= 4) options(mc.cores = 4)
hbc_list <- replicate(64, CompareHamBac(100, 700, 1,
                                       gamma_shape = 1.5,
                                       gamma_mean = c(rep(10, 200), rep(50, 200), rep(10, 200)),
                                       ar_coefs = 0.7,
                      sampling_interval = 48, K_hamstr = c(10, 10), K_bacon = 600/5),
                      simplify = FALSE)

PlotHamstrBacon(hbc_list[[1]])
GetCoverages(hbc_list[[1]])

names(hbc_list) <- 1:length(hbc_list)

#saveRDS(hbc_list, file = gsub(" ", "_", (paste0("package-dev-scripts/cached-runs/ar0.5_gs3_gm10_50_10_si48_", Sys.time(),".RDS"))))
# Big

lapply(hbc_list[[1]]$bacon, function(x) format(object.size(x), units = "Mb"))

covr <- plyr::ldply(hbc_list, GetCoverages, .progress = "text")

covr %>%
  gather(percentile, coverage, -RMSE, -method, -.id) %>%
  mutate(percentile_num = readr::parse_number(percentile)) %>%
  ggplot(aes(x = percentile, y = coverage, colour = method)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_hline(aes(yintercept = percentile_num/100)) #+
  #facet_wrap(~percentile, scales = "free_y")


covr %>%
  gather(percentile, coverage, -RMSE, -method, -.id) %>%
  mutate(percentile_num = readr::parse_number(percentile),
         method_num = as.numeric(as.factor(method))) %>%
  tbl_df() %>%

  ggplot(aes(x = method_num, y = coverage, group = .id)) +
  geom_point() +
  geom_line(alpha = 0.5) +
  geom_hline(aes(yintercept = percentile_num/100)) +
  facet_wrap(~percentile, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 2), labels = c("bacon", "hamstr"))


covr %>%
  gather(percentile, coverage, -RMSE, -method, -.id) %>%
  mutate(percentile_num = readr::parse_number(percentile)) %>%
  ggplot(aes(x = percentile, y = coverage, fill = method)) +
  geom_violin() +
  #geom_point(position = position_jitter(width = 0.2)) +
  geom_hline(aes(yintercept = percentile_num/100)) +
  facet_wrap(~percentile, scales = "free")


covr %>%
  ggplot(aes(x = method, y = RMSE)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))


covr %>%
  mutate(method_num = as.numeric(as.factor(method))) %>%
  ggplot(aes(x = method_num, y = RMSE, group = .id)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = c(1, 2), labels = c("bacon", "hamstr"))




covr %>%
  group_by(method) %>%
  summarise(across(everything(), mean))



get_bacon_memory <- function(hamstr_bacon_fit){

  b1 <- hamstr_bacon_fit
  b1$posterior <- b1$info$output
  w <- b1$posterior[,ncol(b1$posterior) -1]
  thick <- b1$info$thick
  R <- w^(1/thick)

  dplyr::tibble(iter = 1:nrow(b1$posterior),
                R = R, w = w)

}

get_bacon_memory(hbc3$bacon) %>%
  pivot_longer(cols = c("R", "w")) %>%
  ggplot(aes(x = value, fill = name)) +
    geom_density() +
  scale_x_continuous(limits = c(0, 1))


plot(hbc4$hamstr)


CompareMedAgeMod <- function(sim){

  hamstr = join_mod_dat(sim$hamstr, sim$sim_core)
  bacon = join_mod_dat(sim$bacon, sim$sim_core)

  hb <- bind_rows(hamstr=hamstr, bacon=bacon, .id = "method") %>%
    mutate(error = `50%` - age)

  p <- hb %>%
    ggplot(aes(x = depth, y = `50%` - age, colour = method)) +
    geom_line() +
    geom_hline(yintercept = 0)

  return(list(data = hb, plot = p))

}

err <- CompareMedAgeMod(hbc4)

err$plot

err$data %>%
  ggplot(aes(x = depth, y = abs(`50%` - age), colour = method)) +
  geom_line() +
  geom_hline(yintercept = 0)

err$data %>%
  ggplot(aes(x = (`50%` - age), fill = method)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~method, ncol = 1)


err$data %>%
  filter(complete.cases(mean)) %>%
  group_by(method) %>%
  summarise(n = n())





PlotHamList <- function(mod.list){

  obs_ages <- mod.list[[1]]$data[c("depth", "obs_age", "obs_err")]

  obs_ages <- dplyr::as_tibble(obs_ages)

  obs_ages <- obs_ages %>%
    dplyr::rename(.,
                  age = obs_age,
                  err  = obs_err) %>%
    dplyr::mutate(age_upr = age + 2*err,
                  age_lwr = age - 2*err)


  s <- lapply(mod.list, function(x) {
    summary(x)
  })

  s <- bind_rows(s, .id = "model")

  alpha.lvl <- 0.5


  p.age.sum <- s %>%
    ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`, fill = model), alpha = alpha.lvl/1.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`, fill = model), alpha = alpha.lvl) +
    ggplot2::geom_line(ggplot2::aes(y = `50%`, colour = model), lwd = 1) +
    ggplot2::labs(x = "Depth", y = "Age") +
    theme_bw()


  p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(data = obs_ages,
                            ggplot2::aes(x = depth,
                                         ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                            colour = "black") +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age),
                        colour = "black")

  p.age.sum

}

PlotHamList(tail(hbc4, -1)) +
  geom_line(data = hbc4$sim_core,
            aes(x = depth, y = age), colour = "black")

PlotHamList(list("hamstr_10_10" = hamstr_10_10, hamstr_50 = hamstr_50, "hambac_5" = hambac_5)) #+
# facet_wrap(~model)
