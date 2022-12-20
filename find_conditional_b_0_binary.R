pacman::p_load(tidyr, foreach, doParallel, tibble, purrr, dplyr)

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

logit <- function(p){
  log(p/(1-p))
}

### also called expit()
logit_inverse <- function(eta){
  exp(eta)/(1+exp(eta))
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


cores <- parallel::detectCores()
registerDoParallel(cores = cores)
max_ss <- c(100, 200, 500, 1000)
all_res <- foreach(k = 1:length(max_ss), 
                   .combine = "rbind", 
                   .errorhandling = "remove") %do% {
  
  set.seed(123, kind = "L'Ecuyer-CMRG")
  data_structure <- foreach(i = 1:5000,
                            .errorhandling = "remove") %dopar% {
                              generate_data(i, max_ss[k])
                            }
  
  marginal_control_risk <- 0.3
  
  res <- foreach(i = 1:length(data_structure), .combine = "rbind", .errorhandling = "remove") %dopar% {
    get_conditional_pars(i, data_structure[[i]], beta_1 = 1, beta_2 = -0.5, beta_3 = 1, 
                         beta_4 = -0.1, beta_5 = 0.5, marginal_control_risk = marginal_control_risk, 
                         max_ss = max_ss[k])
  }
  
  res %>%
    summarise(max_ss = mean(max_ss),
              marginal_control_risk = mean(marginal_control_risk),
              b_0_c = mean(b_0_c))
  
} 

stopImplicitCluster()
saveRDS(all_res, "PATH/FILENAME.RDS")

