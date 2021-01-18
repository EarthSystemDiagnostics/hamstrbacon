#devtools::install(quick = FALSE, dependencies = FALSE)
# install.packages("../hamstr", repos = NULL, type = "source")
library(hamstr)

MSB2K_cal <- hamstr::calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")
ham1 <- hamstr(depth = MSB2K_cal$depth, obs_age = MSB2K_cal$age.14C.cal, obs_err = MSB2K_cal$age.14C.cal.se)

hambac1 <- hamstr_bacon(depth = MSB2K$depth, obs_age = MSB2K$age, obs_err = MSB2K$error, thick = 5, d.min = 0,
                        acc.mean = 30, ask = FALSE, suggest = "accept")




###
# simulate correlated gamma?

SimulateAgeDepth <- function(top, bottom, d_depth, gamma_shape, gamma_mean,
                             ar1){
  
  depth <- seq(top, bottom, by = d_depth)
  
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
  
  data.frame(depth = depth, age = age, acc.rates_yr_cm = c(acc.rates_yr_cm, NA))
}



library(tidyverse)

tmp <- SimulateAgeDepth(100, 700, 1, 1.5, c(rep(10, 200), rep(50, 200), rep(10, 200)), 0.7) %>% 
  tbl_df()

tmp %>% 
  ggplot(aes(depth, age)) +
  geom_line()


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
              RMSE = sqrt(mean((age - mean)^2)))
  
}


### Compare hamstr and Bacon

CompareHamBac <- function(top, bottom, d_depth, gamma_shape, gamma_mean, ar_coefs,
                          sampling_interval, K_hamstr, K_bacon){
  
  ad1 <- SimulateAgeDepth(top = top, bottom = bottom, d_depth = d_depth,
                          gamma_shape = gamma_shape, gamma_mean = gamma_mean,
                          ar1 = ar_coefs) %>% 
    tbl_df()
  
  ad1 <- ad1 %>% 
    #rowwise() %>% 
    mutate(rad_age = Bchron::unCalibrate(age, type = "ages")) %>% 
    mutate(sampled = (rank(depth) %% sampling_interval == 0), 
           rad_age_sigma = 20 + age * 0.02) %>% 
    mutate(rad_age_hat = rnorm(n(), rad_age, rad_age_sigma))
  
  ad1_samp <- ad1 %>% 
    filter(sampled)
  
  ad1_samp <- calibrate_14C_age(ad1_samp, "rad_age_hat", "rad_age_sigma")
  
  
  ham1 <-  hamstr(depth = ad1_samp$depth, obs_age = ad1_samp$age.14C.cal,
                     obs_err = ad1_samp$age.14C.cal.se,
                     top_depth = min(ad1$depth), bottom_depth = max(ad1$depth),
                  K = K_hamstr)
  
 
  bac1 <- hamstr_bacon(depth = ad1_samp$depth, obs_age = ad1_samp$rad_age_hat,
                       obs_err = ad1_samp$rad_age_sigma,
                       acc.mean = 30, acc.shape = 1.5,
                       d.min = min(ad1$depth), d.max = max(ad1$depth),
                       thick = diff(range(ad1$depth)) / K_bacon,
                       suggest = "FALSE")
  
  return(list(
    sim_core = ad1, hamstr = ham1, bacon = bac1
  ))
}

hbc1 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
              sampling_interval = 12)

hbc2 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
                      sampling_interval = 1)

hbc3 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
                      sampling_interval = 12, K_hamstr = optimal_K(300, 10), K_bacon = 100)


hbc4 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
                      sampling_interval = 12, K_hamstr = optimal_K(300, 10), K_bacon = 300)

hbc5 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
                      sampling_interval = 12, K_hamstr = optimal_K(100, 10), K_bacon = 100)


hbc6 <- CompareHamBac(100, 400, 1, gamma_shape = 1.5, gamma_mean = 20, ar_coefs = 0.7,
                      sampling_interval = 24, K_hamstr = optimal_K(100, 10), K_bacon = 100)


GetCoverages <- function(sim){
  
  hamstr = AgeModCoverage(sim$hamstr, sim$sim_core)
  bacon = AgeModCoverage(sim$bacon, sim$sim_core)
  
  bind_rows(hamstr=hamstr, bacon=bacon, .id = "method")
  
  }

PlotHamstrBacon <- function(sim){
  
  p1 <- plot(sim$hamstr, type = "age_models")+
  geom_line(data = sim$sim_core, aes(y = age), colour = "red")
  
  p2 <- plot(sim$bacon) +
  geom_line(data = sim$sim_core, aes(y = age), colour = "red")
  
  egg::ggarrange(p1, p2, ncol = 2)

}

GetCoverages(hbc1)
GetCoverages(hbc2)
GetCoverages(hbc3)
GetCoverages(hbc4)
GetCoverages(hbc5)
GetCoverages(hbc6)

dim(hbc6$bacon$posterior)


PlotHamstrBacon(hbc3)
PlotHamstrBacon(hbc4)
PlotHamstrBacon(hbc5)
PlotHamstrBacon(hbc6)

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


######


# Simulated with power law spec

set.seed(1)
ar.mod <- arima0(PaleoSpec::SimPowerlaw(1.5, 1e05), order = c(10,0,0))
ar.cfs <- head(ar.mod$coef, -1)
ar.cfs


ad1 <- SimulateAgeDepth(top = 0, bottom = 10000, d_depth = 1, gamma_shape = 3, gamma_mean = 50, ar1 = ar.cfs)

#plot(acc.rates_yr_cm~depth, type = "l", data = ad1, ylim = c(0, max(ad1$acc.rates_yr_cm, na.rm = T)))

plot(age~depth, type = "l", data = ad1, ylim = c(0, max(ad1$age, na.rm = T)))
abline(0, 50, col = "red")

rbacon::Bacon

library(tidyverse)
exp <- tibble(r = 1:10) %>% 
  group_by(r) %>% 
  do({
    SimulateAgeDepth(top = 0, bottom = 500, d_depth = 1,
                     gamma_shape = 3, gamma_mean = 50,
                     ar1 = 0.7)
  })

exp %>% 
  ggplot(aes(x = depth, y = age, group = r)) +
  geom_line(alpha = 0.25) +
  geom_abline(intercept = 0, slope = 50, colour = "red", lwd = 2) +
  theme_bw()





## simulated with big gap
library(hamstr)

sim.dat <- SimulateAgeDepth(top = 0, bottom = 500, d_depth = 1,
                            gamma_shape = 1.5, gamma_mean = 50,
                            ar1 = 0.7) %>% 
  as_tibble() %>% 
  mutate(sigma.age = 30 + age * 0.025,
         rad.age = Bchron::unCalibrate(age, type = "ages"))

obs.dat <- sim.dat %>% 
  filter(depth %in% c(seq(20, 200, by = 20), seq(350, 500, by = 50))) 

sim.dat %>% 
  ggplot(aes(x = depth, y = age)) +
  geom_line(aes(colour = "cal.age")) +
  geom_line(aes(y = rad.age, colour = "rad.age"))


h1 <- hamstr(depth = obs.dat$depth, obs_age = obs.dat$age, obs_err = obs.dat$sigma.age,
             K = c(10, 10), mem_mean = 0.9,
             cores = 3)

b1 <- hamstr_bacon(depth = obs.dat$depth, obs_age = obs.dat$rad.age, obs_err = obs.dat$sigma.age,
                   thick = 1, acc.mean = 50)

plot(h1, type = "age_models") +
  geom_line(data = sim.dat, aes(x = depth, y = age, colour = "Simulated"))


plot(b1) +
  geom_line(data = sim.dat, aes(x = depth, y = age, colour = "Simulated"))



h1$data

