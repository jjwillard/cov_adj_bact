### COVID Trial - CER = 0.07

pacman::p_load(tidyr, rstanarm, foreach, doParallel, tibble, purrr, dplyr)

r_age <- function(max_ss){
  age <- tibble(age = ceiling(rnorm(250*max_ss, mean = 62, sd = 40))) %>%
    filter(age >= 18, age <= 90) 
  int_1 <- age %>% filter(age < 39) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  int_2 <- age %>% filter(age >= 39, age < 55) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  int_3 <- age %>% filter(age >= 55, age < 70) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  int_4 <- age %>% filter(age >= 70) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  collection <- c(int_1, int_2, int_3, int_4)
  sample(collection, length(collection))
}

r_resp_rate <- function(max_ss){
  rr <- tibble(rr = ceiling(rnorm(250*max_ss, mean = 30, sd = 6))) %>%
    filter(rr >= 12, rr <= 40) 
  int_1 <- rr %>% filter(rr < 18) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  int_2 <- rr %>% filter(rr >= 18, rr < 20) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  int_3 <- rr %>% filter(rr >= 20, rr < 22) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  int_4 <- rr %>% filter(rr >= 22) %>% slice_sample(n = max_ss/4) %>% flatten_dbl()
  collection <- c(int_1, int_2, int_3, int_4)
  sample(collection, length(collection))
}

generate_data <- function(iteration, max_ss){
  tibble(treatment = rbinom(max_ss, 1, 0.5),
         x1 = r_age(max_ss), #age
         x2 = r_resp_rate(max_ss), #rr
         x3 = rbinom(n = max_ss, size = 1, prob = 0.478), #female
         x4 = rbinom(n = max_ss, size = 1, prob = 0.216), #chest pain, 
         x5 = rbinom(n = max_ss, size = 1, prob = 0.403)) # arrival police/ambulance
}

logit <- function(p){
  log(p/(1-p))
}

logit_inverse <- function(eta){
  exp(eta)/(1+exp(eta))
}

generate_outcomes <- function(data, effect_treatment, b_0_c, beta_1, beta_2, beta_3, beta_4, beta_5, control_risk, max_ss){

  data %>%
    mutate(eta = b_0_c + effect_treatment*treatment + beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4 + beta_5*x5,
           p = logit_inverse(eta),
           y = rbinom(n = max_ss, size = 1, p))
  
}

## Raw estimates of marginal relative risk
get_marginal_rr <- function(iteration, b_0_c, lo_c){
  data_y <- generate_outcomes(generate_data(1, 5000),
                              effect_treatment = lo_c,
                              b_0_c = b_0_c,
                              beta_1 = 0.09193548, 
                              beta_2 = 0.09666667, 
                              beta_3 = -0.61, 
                              beta_4 = -0.8, 
                              beta_5 = 0.63,
                              max_ss = 5000)
  
  res <- data_y %>%
    group_by(treatment) %>%
    summarize(mean_y = mean(y, na.rm = TRUE))
  
  risk_ctr <- res %>% filter(treatment == 0) %>% pull(mean_y)
  risk_trt <- res %>% filter(treatment == 1) %>% pull(mean_y)
  
  risk_trt / risk_ctr
  exp(log(risk_trt) - log(risk_ctr))
}

cores <- detectCores()
# N=1000, CER = 0.07
registerDoParallel(cores = cores)
rr_20ev_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -10.76457, lo_c = -0.80086207)
}
rr_20ev_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -1.18965517)
}

rr_30ev_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -0.790)
}
rr_30ev_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -1.180)
}

rr_40ev_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -0.79586207)
}
rr_40ev_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -1.20)
}

rr_50ev_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -0.82758621)
}
rr_50ev_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c =  -10.76457, lo_c = -1.204655)
}
stopImplicitCluster()

res <- tibble(max_ss = 1000,
              events = c(20, 20, 30, 30, 40, 40, 50, 50),
              lo_c = c(-0.80086207, -1.18965517, -0.790, -1.180, -0.79586207, -1.20,
                       -0.82758621, -1.204655),
              approx_power = c(50, 80, 50, 80, 50, 80, 50, 80),
              rr = c(mean(rr_20ev_50p), mean(rr_20ev_80p), mean(rr_30ev_50p), mean(rr_30ev_80p),
                     mean(rr_40ev_50p), mean(rr_40ev_80p), mean(rr_50ev_50p), mean(rr_50ev_80p)))

saveRDS(res, "PATH/FILENAME.RDS")
