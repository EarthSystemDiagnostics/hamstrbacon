#devtools::install(quick = FALSE, dependencies = FALSE)

library(hamstr)

MSB2K_cal <- hamstr::calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")
ham1 <- hamstr(depth = MSB2K_cal$depth, obs_age = MSB2K_cal$age.14C.cal, obs_err = MSB2K_cal$age.14C.cal.se)


hambac1 <- hamstr_bacon(depth = MSB2K$depth, obs_age = MSB2K$age, obs_err = MSB2K$error, thick = 5, d.min = 0,
                        acc.mean = 50, ask = FALSE)


extract_hamstr_fit(ham1)
phb <- predict(hambac1)
summarise_bacon_age_models(hambac1)

plot_summary_bacon_age_models(hambac1)

plot(ham1)

summarise_bacon_age_models(phb)

extract_hamstr_fit(ham1)

predict(ham1)

ham_sum <- summary(ham1)
bac_sum <- summary(hambac1)

hb_sum <- bind_rows(ham = ham_sum, bac = bac_sum, .id = "model")

p.age.sum <- hb_sum %>%
  ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`), fill = "Lightgrey") +
  ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`), fill = "Darkgrey") +
  ggplot2::geom_line() +
  ggplot2::geom_line(ggplot2::aes(y = `50%`), colour = "Green") +
  ggplot2::labs(x = "Depth", y = "Age") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = ggplot2::element_blank()) +
  ggplot2::facet_wrap(~model)

p.age.sum <- p.age.sum +
  ggplot2::geom_linerange(data = MSB2K_cal,
                          ggplot2::aes(x = depth,
                                       ymax = age.14C.cal + age.14C.cal.se,
                                       ymin = age.14C.cal - age.14C.cal.se),
                          inherit.aes = FALSE,
                          colour = "Blue", size = 1.25) +
  ggplot2::geom_point(data = MSB2K_cal, ggplot2::aes(y = age.14C.cal),
                      colour = "Blue")

p.age.sum

### compare coverage

join_mod_dat <- function(hamstr_fit, dat){
  i1 <- predict(hamstr_fit, dat$depth)
  s1 <- summary(i1)
  c1 <- dplyr::left_join(dat, s1)
  return(c1)
}

plot(ham1)
plot_summary_bacon_age_models(hambac1)

cham1 <- join_mod_dat(ham1, MSB2K_cal)

cbac1 <- join_mod_dat(hambac1, MSB2K_cal)

c0123 <- bind_rows(cham1=cham1, cbac1=cbac1, .id = "exp")

c0123 %>%
  ggplot(aes(x = depth, y = mean - age.14C.cal, colour = exp)) +
  geom_line()+
  expand_limits(y = c(-100, 100))


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
  
  min.age <- gamma_mean * top
  
  age <- cumsum(c(min.age, acc.rates_yr_cm) * d_depth)
  
  data.frame(depth = depth, age = age, acc.rates_yr_cm = c(acc.rates_yr_cm, NA))
}

set.seed(1)
#ar.mod <- arima0(PaleoSpec::SimPowerlaw(1.5, 1e05), order = c(10,0,0))
#ar.cfs <- head(ar.mod$coef, -1)

ar.cfs <- c(0.75)

ad1 <- SimulateAgeDepth(top = 100, bottom = 500, d_depth = 1, gamma_shape = 0.5, gamma_mean = 30, ar1 = ar.cfs) %>% 
  tbl_df()

ad1 <- ad1 %>% 
  #rowwise() %>% 
  mutate(rad_age = Bchron::unCalibrate(age, type = "ages")) %>% 
  mutate(sampled = (depth %% 48 == 0), 
         rad_age_sigma = 20 + age * 0.02) %>% 
  mutate(rad_age_hat = rnorm(n(), rad_age, rad_age_sigma))


#plot(acc.rates_yr_cm~depth, type = "l", data = ad1, ylim = c(0, max(ad1$acc.rates_yr_cm, na.rm = T)))

plot(age~depth, type = "l", data = ad1)


ad1_samp <- ad1 %>% 
  filter(sampled)

ad1_samp <- calibrate_14C_age(ad1_samp, "rad_age_hat", "rad_age_sigma")

lm(age ~ depth, data = ad1_samp)

if (parallel::detectCores() >= 3) options(mc.cores = 3)



ham1 <- hamstr(depth = ad1_samp$depth, obs_age = ad1_samp$age.14C.cal, obs_err = ad1_samp$age.14C.cal.se,
               top_depth = min(ad1$depth), bottom_depth = max(ad1$depth),
               K = optimal_K(400, 5))


ham2 <- hamstr(depth = ad1_samp$depth, obs_age = ad1_samp$age.14C.cal, obs_err = ad1_samp$age.14C.cal.se,
               top_depth = min(ad1$depth), bottom_depth = max(ad1$depth))

ham3 <- hamstr(depth = ad1_samp$depth, obs_age = ad1_samp$age.14C.cal, obs_err = ad1_samp$age.14C.cal.se,
               top_depth = min(ad1$depth), bottom_depth = max(ad1$depth),
               #mem_mean = 0.5, mem_strength = 2,
               K = 400)

#ham3 <- ham2

system.time(
hb2 <- hamstr_bacon(depth = ad1_samp$depth, obs_age = ad1_samp$rad_age_hat, obs_err = ad1_samp$rad_age_sigma,
                    thick = 5, acc.mean = 30, acc.shape = 1.5,
                    d.min = min(ad1$depth), d.max = max(ad1$depth))
)


plot(ham1, type = "age_models") +
  geom_line(data = ad1, aes(y = age), colour = "red")

plot(ham2, type = "age_models") +
  geom_line(data = ad1, aes(y = age), colour = "red")

plot(ham3, type = "age_models") +
  geom_line(data = ad1, aes(y = age), colour = "red")

plot(ham3, summarise  =FALSE, n.iter = 100)

plot_summary_bacon_age_models(hb2)+
  geom_line(data = ad1, aes(x = depth, y = age), colour = "red")

## should be comparing coverage of simulated "true" ages not the data points
## as data might have extra error added
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

plot(ham1, type = "mem")
plot(ham2, type = "mem")
plot(ham3, type = "mem")

AgeModCoverage(ham1, ad1)
AgeModCoverage(ham2, ad1)
AgeModCoverage(ham3, ad1)

AgeModCoverage(hb2, ad1)

