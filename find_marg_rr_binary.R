### Finds marginal relative risk given values of b_0_c and conditional log odds (lo_c)
### Control event risk of 0.3

pacman::p_load(tidyr, rstanarm, foreach, doParallel, tibble, purrr, dplyr)

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

logit <- function(p){
  log(p/(1-p))
}

logit_inverse <- function(eta){
  exp(eta)/(1+exp(eta))
}

generate_outcomes <- function(data, effect_treatment, b_0_c, beta_1, beta_2, beta_3, beta_4, beta_5, control_risk){
  
  max_ss <- nrow(data)
  
  data %>%
    mutate(eta = b_0_c + effect_treatment*treatment + beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4 + beta_5*x5,
           p = logit_inverse(eta),
           y = rbinom(n = max_ss, size = 1, p))
  
}

get_marginal_rr <- function(iteration, b_0_c, lo_c){
  data_y <- generate_outcomes(generate_data(1, 5000),
                              effect_treatment = lo_c,
                              b_0_c = b_0_c,
                              beta_1 = 1,
                              beta_2 = -0.5,
                              beta_3 = 1,
                              beta_4 = -0.1, 
                              beta_5 = 0.5)
  
  res <- data_y %>%
    group_by(treatment) %>%
    summarize(mean_y = mean(y, na.rm = TRUE))
  
  risk_ctr <- res %>% filter(treatment == 0) %>% pull(mean_y)
  risk_trt <- res %>% filter(treatment == 1) %>% pull(mean_y)
  
  risk_trt / risk_ctr
}

# b_0_c obtained from conditional_b0_binary.R
# lo_c obtained as those values where the unadjusted model achieves 50% and 80% power
# except for N=100, where it's 30% and 40% power

cores <- detectCores()
## CER = 0.3
registerDoParallel(cores = cores)
rr_100_30p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.265619, lo_c = -1.03448276)
}
rr_100_40p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.265619, lo_c = -1.24137931)
}

rr_200_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.264603, lo_c = -0.87931034)
}
rr_200_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.264603, lo_c = -1.36)
}

rr_500_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.267382, lo_c = -0.56)
}
rr_500_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.267382, lo_c = -0.82758621)
}

rr_1000_50p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.265569, lo_c = -0.3925)
}
rr_1000_80p <- foreach(i = 1:5000, .combine = 'c', .errorhandling = 'remove') %dopar% {
  get_marginal_rr(i, b_0_c = -1.265569, lo_c = -0.55172414)
}
stopImplicitCluster()

res <- tibble(max_ss = c(100, 100, 200, 200, 500, 500, 1000, 1000),
              lo_c = c(-1.03448276, -1.24137931, -0.87931034, -1.36, -0.56, -0.82758621, -0.3925,
                       -0.55172414),
              approx_power = c(30, 40, 50, 80, 50, 80, 50, 80),
              rr = c(mean(rr_100_30p), mean(rr_100_40p), mean(rr_200_50p), mean(rr_200_80p),
                     mean(rr_500_50p), mean(rr_500_80p), mean(rr_1000_50p), mean(rr_1000_80p)))

saveRDS(res, "PATH/FILENAME.RDS")
