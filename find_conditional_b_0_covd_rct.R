pacman::p_load(tidyr, foreach, doParallel, tibble, purrr, dplyr)

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

generate_lin_pred <- function(data, beta_1, beta_2, beta_3, beta_4, beta_5, type = NULL){
  
  if(is.null(type)){
    data %>%
      mutate(lin_pred = beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4 + beta_5*x5)
  } else if(type == "control"){
    data %>%
      filter(treatment == 0) %>%
      mutate(lin_pred_ctr = beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4 + beta_5*x5)
  } else if(type == "treatment"){
    data %>%
      filter(treatment == 1) %>%
      mutate(lin_pred_trt = beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4 + beta_5*x5)
  }
}

find_b_0_c <- function(b_0_c, marginal_control_risk, lin_pred_control){
  l <- length(lin_pred_control)
  sum(logit_inverse(b_0_c + lin_pred_control)) - l*marginal_control_risk
}

get_conditional_pars <- function(iteration, data, beta_1, beta_2, beta_3, beta_4, beta_5, 
                                 marginal_control_risk, max_ss){
  ### find b_0_c only
  
  lin_pred_ctr <- generate_lin_pred(data, beta_1, beta_2, beta_3, beta_4, beta_5, type = "control")$lin_pred_ctr
  b_0_c <- uniroot(find_b_0_c, interval = c(-500, 500), 
                   marginal_control_risk = marginal_control_risk, 
                   lin_pred_control = lin_pred_ctr)$root

  tibble(iteration = iteration,
         max_ss = max_ss,
         b_0_c = b_0_c,
         beta_1 = beta_1, 
         beta_2 = beta_2, 
         beta_3 = beta_3, 
         beta_4 = beta_4, 
         beta_5 = beta_5, 
         marginal_control_risk = marginal_control_risk)
}

### Only interested in max-ss of 1000 given small CER = 0.07
cores <- parallel::detectCores()
registerDoParallel(cores = cores)
max_ss <- 1000
all_res <- foreach(k = 1:length(max_ss), 
                   .combine = "rbind", 
                   .errorhandling = "remove") %do% {
  
  set.seed(123, kind = "L'Ecuyer-CMRG")
  data_structure <- foreach(i = 1:5000,
                            .errorhandling = "remove") %dopar% {
                              generate_data(i, 5000) # generate huge dataset under the settings
                            }
  
  marginal_control_risk <- 0.07
  
  res <- foreach(i = 1:length(data_structure), .combine = "rbind", .errorhandling = "remove") %dopar% {
    get_conditional_pars(i, data_structure[[i]], beta_1 = 0.09193548, beta_2 = 0.09666667, beta_3 = -0.61, 
                         beta_4 =  -0.8, beta_5 = 0.63, marginal_control_risk = marginal_control_risk, 
                         max_ss = max_ss[k])
  }
  
  res %>%
    summarise(max_ss = mean(max_ss),
              marginal_control_risk = mean(marginal_control_risk),
              b_0_c = mean(b_0_c))
  
} 

stopImplicitCluster()
saveRDS(all_res, "PATH/FILENAME.RDS")
