pacman::p_load(tidyr, rstanarm, simsurv, foreach, doParallel, tibble, purrr, dplyr)

generate_data <- function(iteration, max_ss){
  tibble(id = 1:max_ss,
         treatment = rbinom(n = max_ss, size = 1, p= 0.5),
         x1 = rbinom(n = max_ss, size = 1, p= 0.5),
         x2 = rbinom(n = max_ss, size = 1, p= 0.5),
         x3 = rnorm(n = max_ss),
         x4 = x3**2,
         x5 = rnorm(n = max_ss),
         ### x6, x7, x8 are just noise
         x6 = rbinom(n = max_ss, size = 1, p= 0.5),
         x7 = rnorm(n = max_ss),
         x8 = rnorm(n = max_ss))
  
}

generate_outcomes <- function(data, effect_treatment, beta_1, beta_2, beta_3, beta_4, beta_5){

  time_fixed_covariate_effects <- tibble(treatment = effect_treatment, #lh_c,
                                         x1 = beta_1,
                                         x2 = beta_2,
                                         x3 = beta_3,
                                         x4 = beta_4,
                                         x5 = beta_5) %>%
    uncount(nrow(data)) # repeats the rows 
  
  time_data <- simsurv(dist = "exponential",
                       lambdas = 0.05, # mean of 1/0.05 =20 time units
                       betas = time_fixed_covariate_effects,
                       x = data) %>%
    select(-status)
  
  data %>%
    left_join(time_data, by = "id")
}

## Raw estimate of marginal hazard ratio at time t=75

get_marginal_hr <- function(iteration, lh_c){
  data_y <- generate_outcomes(generate_data(1, max_ss = 5000),
                              effect_treatment = lh_c,
                              beta_1 = 1,
                              beta_2 = -0.5,
                              beta_3 = 1,
                              beta_4 = -0.1, 
                              beta_5 = 0.5)
  
  res <- data_y %>%
    mutate(survived = if_else(eventtime > 75, 1, 0)) %>%
    group_by(treatment) %>%
    summarize(p_surv = mean(survived, na.rm = TRUE))
  
  p_surv_ctr <- res %>% filter(treatment == 0) %>% pull(p_surv)
  p_surv_trt <- res %>% filter(treatment == 1) %>% pull(p_surv)
  
  exp(log(-log(p_surv_trt)) - log(-log(p_surv_ctr)))
  
}


## The values of the conditional log hazard ratio (lh_c) are selected as those where
## the unadjuste model achieves approximately 30% and 40% power (N=100), or 50% and 80%
## power (all other N)

cores <- detectCores()
# N=100
registerDoParallel(cores = cores)
hr_100_30p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.6250000)
}
hr_100_40p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c =  -0.84)
}


# N= 200

hr_200_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.53)
}
hr_200_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.8333333)
}


# N = 500

hr_500_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.3500000)
}
hr_500_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.5293333)
}


# N =1000

hr_1000_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.24)
}
hr_1000_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_hr(i, lh_c = -0.3611111)
}
stopImplicitCluster()

res <- tibble(max_ss = c(100, 100, 200, 200, 500, 500, 1000, 1000),
       lh_c = c(-0.6250000, -0.84, -0.53, -0.8333333, -0.3500000, -0.5293333, -0.24, -0.3611111), 
       approx_power = c(30, 40, 50, 80, 50, 80, 50, 80),
       hr = c(mean(hr_100_30p), mean(hr_100_40p), mean(hr_200_50p), mean(hr_200_80p),
              mean(hr_500_50p), mean(hr_500_80p), mean(hr_1000_50p), mean(hr_1000_80p)))

saveRDS(res, "PATH/FILENAME.RDS")