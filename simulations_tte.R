pacman::p_load(simsurv, tidyr, ggplot2, rstanarm, foreach, doParallel, tibble, purrr, dplyr)

### This is to generate the data

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


### This is for time-FIXED outcomes only
generate_outcomes <- function(data, effect_treatment, beta_1, beta_2, beta_3, beta_4, beta_5, max_ss){
  # must have the same names as the df used above for the data to be used in the
  # the simsurv function below

  conditional_params <- tibble(max_ss = c(100, 100, 100, 200, 200, 200, 500, 500, 500, 1000, 1000, 1000),
                               lh_m = log(c(0.7317716, 0.6517841, 1, 0.7690409, 0.6536824, 1, 0.8437124,
                                          0.7697936, 1, 0.8907880, 0.8389165, 1)),
                               lh_c = c(-0.625, -0.84, 0, -0.53, -0.8333333, 0, -0.35, -0.5293333, 0,
                                        -0.24, -0.3611111, 0))
  lh_c <- conditional_params %>% 
    filter(max_ss == max_ss, lh_m == effect_treatment) %>%
    pull(lh_c) 
  
  
  time_fixed_covariate_effects <- tibble(treatment = lh_c,
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

calculate_events <- function(data, type = "total"){
  events <-  data %>%
    group_by(treatment) %>%
    summarise(n_events = sum(event)) 
  n_events_control <- events %>% filter(treatment == 0) %>% pull(n_events)
  n_events_treatment <- events %>% filter(treatment == 1) %>% pull(n_events)
  n_events_total <- n_events_control + n_events_treatment
  n_events_min <- min(n_events_control, n_events_treatment)
  
  if(type == "total"){
    n_events_total
  } else if (type == "control"){
    n_events_control
  } else if (type == "treatment"){
    n_events_treatment
  } else if (type == "min"){
    n_events_min
  }
}

## marginalizes conditional posterior samples
get_marginal_log_hr <- function(fitted_model, time, draws = 3000){
  new_data <- model.frame(fitted_model)
  new_data$treatment <- 1
  surv_ps_trt <- posterior_survfit(fitted_model, 
                                   newdata = new_data,
                                   type = "surv", 
                                   standardise = TRUE, #averages over covariates in sample, essentially a rowMeans() call on the matrix
                                   draws = draws,
                                   times = time,
                                   extrapolate = FALSE,
                                   return_matrix = TRUE)
  
  new_data$treatment <- 0
  surv_ps_ctr <- posterior_survfit(fitted_model, 
                                   newdata = new_data,
                                   type = "surv", 
                                   standardise = TRUE, #averages over covariates in sample
                                   draws = draws,
                                   times = time,
                                   extrapolate = FALSE,
                                   return_matrix = TRUE)
  
  as.numeric(log(-log(surv_ps_trt[[1]])) - log(-log(surv_ps_ctr[[1]])))
}

### Models

## correct
model_correct <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x1 + x2 + x3 + x4 + x5,
                                 data = sim_data,
                                 basehaz = "ms", #"exp"
                                 chains = 3))

## no quad
model_incorrect <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x1 + x2 + x3 + x5,
                                   data = sim_data,
                                   basehaz = "ms", #"exp"
                                   chains = 3))

## correct noise
model_noise <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
                               data = sim_data,
                               basehaz = "ms", #"exp"
                               chains = 3))


## no strong prog
model_no_strong_prog <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x2 + x5,
                                        data = sim_data,
                                        basehaz = "ms", #"exp"
                                        chains = 3))

## no strong prog noise
model_no_strong_prog_noise <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x2 + x5 + x6 + x7 + x8,
                                              data = sim_data,
                                              basehaz = "ms", #"exp"
                                              chains = 3))

## unadjusted
model_unadjusted <- quote(stan_surv(formula = Surv(t, event) ~ treatment,
                                    data = sim_data,
                                    basehaz = "ms", #"exp"
                                    chains = 3))

## correct prior
model_correct_prior <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x1 + x2 + x3 + x4 + x5,
                                       data = sim_data,
                                       basehaz = "ms", #"exp"
                                       prior = prior_correct,
                                       chains = 3))

## correct strong prior
model_correct_prior_strong <- quote(stan_surv(formula = Surv(t, event) ~ treatment + x1 + x2 + x3 + x4 + x5,
                                              data = sim_data,
                                              basehaz = "ms", #"exp"
                                              prior = prior_correct_strong,
                                              chains = 3))

### results output

results_output <- quote(
  tibble(iteration = iteration,
         p_sup_thresh = p_sup_thresh,
         batch_size = batch_size,
         max_ss = max_ss,
         initial_ss = initial_ss,
         max_time = max_time,
         initial_events_req_ia = initial_events_req_ia,
         new_events_req_ia = new_events_req_ia,
         model = !!mod_name,
         effect_treatment = effect_treatment,
         beta_1 = beta_1,
         beta_2 = beta_2,
         beta_3 = beta_3,
         beta_4 = beta_4,
         beta_5 = beta_5,
         runtime_sec = toc - tic,
         n_total = n_total, 
         n_events = n_events,
         n_analyses_total = n_analyses_total,
         time_total = time,
         p_sup = p_sup, 
         superiority = if_else(p_sup > p_sup_thresh, 1, 0),
         reach_max_ss = if_else(n_total >= max_ss, 1, 0), #ge as failsafe in case bug introduce to increment past intended size
         reach_max_time = if_else(time >= max_time, 1, 0),
         trt_est_mean = mean(results, na.rm = TRUE),
         trt_est_mean_hazard = mean(exp(results), na.rm = TRUE),
         trt_est_median = median(results, na.rm = TRUE),
         trt_est_median_hazard = median(exp(results), na.rm = TRUE),
         bias = mean(results - effect_treatment, na.rm = TRUE), # updated 2022-03-09
         bias_hazard = mean(exp(results) - exp(effect_treatment), na.rm = TRUE), # updated 2022-03-09
         relative_bias = if_else(effect_treatment != 0, 
                                 mean((results - effect_treatment)/effect_treatment, na.rm = TRUE), # updated 2022-03-09
                                 NA_real_),
         relative_bias_hazard = mean((exp(results) - exp(effect_treatment))/exp(effect_treatment), na.rm = TRUE),# updated 2022-03-09
         rmse = sqrt(mean((results - effect_treatment)**2)),
         rmse_hazard = sqrt(mean((exp(results) - exp(effect_treatment))**2)),
         mae = mean(abs(results - effect_treatment)),
         mae_hazard = mean(abs(exp(results) - exp(effect_treatment))),
         post_var = var(results),
         post_var_hazard = var(exp(results)),
         message = "Adequate number of events for interim analyses.",
         trt_posterior = list(results %>% as_tibble_col(column_name = "treatment")),
         stan_summary = list(as_tibble(fit_mod$stan_summary, rownames = "parameter"))
))

### results output in event of not enough events within max time
results_output_na <- quote(
  tibble(iteration = iteration, 
         p_sup_thresh = p_sup_thresh,
         batch_size = batch_size,
         max_ss = max_ss,
         initial_ss = initial_ss,
         max_time = max_time,
         initial_events_req_ia = initial_events_req_ia,
         new_events_req_ia = new_events_req_ia,
         model = !!mod_name,
         effect_treatment = effect_treatment,
         beta_1 = beta_1,
         beta_2 = beta_2,
         beta_3 = beta_3,
         beta_4 = beta_4,
         beta_5 = beta_5,
         runtime_sec = toc - tic,
         n_total = n_total, 
         n_events = n_events, 
         n_analyses_total = n_analyses_total,
         time_total = time,
         p_sup = NA_real_, 
         superiority = NA_real_,
         reach_max_ss = NA_real_, #ge as failsafe in case bug introduce to increment past intended size
         reach_max_time = NA_real_,
         trt_est_mean = NA_real_,
         trt_est_mean_hazard = NA_real_,
         trt_est_median = NA_real_,
         trt_est_median_hazard = NA_real_,
         bias = NA_real_,
         bias_hazard = NA_real_,
         relative_bias = NA_real_,
         relative_bias_hazard = NA_real_,
         rmse = NA_real_,
         rmse_hazard = NA_real_,
         mae = NA_real_,
         mae_hazard = NA_real_,
         post_var = NA_real_,
         post_var_hazard = NA_real_,
         message = "Insufficient number of events for interim analyses.",
         trt_posterior = list(tibble(treatment = NA_real_)),
         stan_summary = list(tibble(stan_summary = NA_real_))
))

### single trial simulation
run_single_sim_tte <- function(iteration, 
                               batch_size = 25, # average recruitment rate between analyses
                               max_ss = 100,
                               initial_ss = 25,
                               max_time = 50,
                               enrollment_time_cutoff = 25,
                               initial_events_req_ia = 10, # initial number of events for first interim analysis
                               new_events_req_ia = 10, # number of additional events per additional interim analysis
                               p_sup_thresh = 0.975,
                               effect_treatment, 
                               beta_1 = NULL, beta_2 = NULL, beta_3 = NULL, beta_4 = NULL, beta_5 = NULL,
                               full_data,
                               model){
  
  ### Tracks name of the models used
  mod_name <- rlang::as_label(enquo(model))
  
  ### Pull out one full dataset from list
  iteration_data <- full_data[[iteration]]
  
  ### create enrollment_time variable, set baseline enrollment at 0
  sim_data <- iteration_data[1:initial_ss, ] %>%
    mutate(enrollment_time = 0)
  
  ### initial events required for interim analysis to be performed
  events_req_ia <- initial_events_req_ia
  
  ### time is taken to be at the end of the time interval 
  ### (i,e. 1 is at the end of the first week)
  time <- 0
  n_total <- initial_ss
  p_sup <- 0 
  n_events <- 0
  n_analyses_total <- 0
  
  tic <- tictoc::tic()
  
  
  ## This design will wait for pre-specified number of events before performing interim analysis
  ## And will only enroll them up to a specific time cut-off point
  
  ## set up two while loops
  while(time < max_time){ #n_total < max_ss & 
    
    ### Check for specified number of events before performing interim analysis
    ### increments time but not number of participants
    while(n_events < events_req_ia & time < max_time){
      ## increment time by one and recalculate
      time <- time + 1
      sim_data <- sim_data %>%
        # recalculate event and t under staggered entry
        mutate(event = if_else((enrollment_time + eventtime) <= time, 1, 0),
               t = if_else(event == 1, enrollment_time + eventtime, time))
      n_events <- calculate_events(sim_data, type = "total")
    } 
    
    ### sets a break if not enough events occur for a single interim analysis 
    ### within the required time frame
    
    if(n_events < initial_events_req_ia & time >= max_time){
      results_output <- results_output_na
      break
    }

    sd_treatment <- sd(sim_data$treatment)
    sd_x1 <- sd(sim_data$x1)
    sd_x2 <- sd(sim_data$x2)
    sd_x3 <- sd(sim_data$x3)
    sd_x4 <- sd(sim_data$x4)
    sd_x5 <- sd(sim_data$x5)
    
    # centered at DGM except for trt
    prior_correct <- normal(location = c(0, 1, -0.5, 1, -0.1, 0.5),
                            scale = 2.5,
                            autoscale = TRUE)
    # centered and scaled
    prior_correct_strong <- normal(location = c(0, 1, -0.5, 1, -0.1, 0.5),
                                   scale = c(2.5/sd_treatment, 1/sd_x1, 1/sd_x2, 1/sd_x3, 1/sd_x4, 1/sd_x5),
                                   autoscale = FALSE)
    
    n_analyses_total <- n_analyses_total + 1
    fit_mod <- eval(model)
    results <- get_marginal_log_hr(fitted_model = fit_mod, time = time) #log hazard scale
    
    p_sup <- mean(exp(results) < 1)
    n_total <- as.numeric(length(sim_data$treatment))
    
    ## Determine whether or not to stop for superiority based on stopping rule
    if(p_sup > p_sup_thresh){
      break
    }
    
    ## Stop if maximum time
    if(time >= max_time){ #n_total >= max_ss | 
      break
    }
    
    ### If not stopping, then enroll more participants only if time is less than enrollment cutoff
    ### must append to old sim_data and calculate new enrollment time
    ### batch_size is average recruitment between analyses
    if(n_total < max_ss & time < enrollment_time_cutoff){
      new_data <- iteration_data[(n_total + 1):(n_total + batch_size),] %>%
        mutate(enrollment_time = time)
      sim_data <- bind_rows(sim_data, new_data)
    }
    
    ## Increase number of events over current number to be required at next interim analysis
    events_req_ia <- n_events + new_events_req_ia
    
  }
  
  toc <- tictoc::tic()
  
  ### produce final results
  eval(results_output)
  
}

### Runs .n_iterations trial simulations

complete_sim_tte <- function(.n_iterations,
                             .max_ss, 
                             .initial_ss,
                             .max_time,
                             .enrollment_time_cutoff,
                             .initial_events_req_ia,
                             .new_events_req_ia,
                             .batch_size,
                             .effect_treatment,
                             .beta_1, .beta_2, .beta_3, .beta_4, .beta_5,
                             .p_sup_thresh = 0.99,
                             .n_cores = parallel::detectCores(),
                             .seed = 123){
  
  set.seed(.seed, kind = "L'Ecuyer-CMRG")
  .data_structure <- foreach(i = 1:.n_iterations,
                             .errorhandling = "remove") %do% {
                               generate_data(i, .max_ss)
                             }
  
  set.seed(.seed, kind = "L'Ecuyer-CMRG")
  .data <- foreach(j = 1:.n_iterations) %do% {
    generate_outcomes(.data_structure[[j]], .effect_treatment, .beta_1, .beta_2, 
                      .beta_3, .beta_4, .beta_5, .max_ss)
  }
  
  registerDoParallel(cores = .n_cores)
  
  ## correct
  .model_correct <- foreach(i=1:.n_iterations,
                            .combine='bind_rows',
                            .inorder = FALSE,
                            .errorhandling = "remove") %dopar% {
                              run_single_sim_tte(i,
                                                 batch_size = .batch_size,
                                                 max_ss = .max_ss,
                                                 initial_ss = .initial_ss,
                                                 max_time = .max_time,
                                                 enrollment_time_cutoff = .enrollment_time_cutoff,
                                                 initial_events_req_ia = .initial_events_req_ia,
                                                 new_events_req_ia = .new_events_req_ia,
                                                 effect_treatment = .effect_treatment,
                                                 beta_1 = .beta_1,
                                                 beta_2 = .beta_2,
                                                 beta_3 = .beta_3,
                                                 beta_4 = .beta_4,
                                                 beta_5 = .beta_5,
                                                 p_sup_thresh = .p_sup_thresh,
                                                 full_data = .data,
                                                 model = model_correct)
                            }
  
  ## no quad
  .model_incorrect <- foreach(i=1:.n_iterations,
                              .combine='bind_rows',
                              .inorder = FALSE,
                              .errorhandling = "remove") %dopar% {
                                run_single_sim_tte(i,
                                                   batch_size = .batch_size,
                                                   max_ss = .max_ss,
                                                   initial_ss = .initial_ss,
                                                   max_time = .max_time,
                                                   enrollment_time_cutoff = .enrollment_time_cutoff,
                                                   initial_events_req_ia = .initial_events_req_ia,
                                                   new_events_req_ia = .new_events_req_ia,
                                                   effect_treatment = .effect_treatment,
                                                   beta_1 = .beta_1,
                                                   beta_2 = .beta_2,
                                                   beta_3 = .beta_3,
                                                   beta_4 = .beta_4,
                                                   beta_5 = .beta_5,
                                                   p_sup_thresh = .p_sup_thresh,
                                                   full_data = .data,
                                                   model = model_incorrect)
                              }
  
  ## correct noise
  .model_noise <- foreach(i=1:.n_iterations,
                          .combine='bind_rows',
                          .inorder = FALSE,
                          .errorhandling = "remove") %dopar% {
                            run_single_sim_tte(i,
                                               batch_size = .batch_size,
                                               max_ss = .max_ss,
                                               initial_ss = .initial_ss,
                                               max_time = .max_time,
                                               enrollment_time_cutoff = .enrollment_time_cutoff,
                                               initial_events_req_ia = .initial_events_req_ia,
                                               new_events_req_ia = .new_events_req_ia,
                                               effect_treatment = .effect_treatment,
                                               beta_1 = .beta_1,
                                               beta_2 = .beta_2,
                                               beta_3 = .beta_3,
                                               beta_4 = .beta_4,
                                               beta_5 = .beta_5,
                                               p_sup_thresh = .p_sup_thresh,
                                               full_data = .data,
                                               model = model_noise)
                          }
  
  ## unadjusted
  .model_unadjusted <- foreach(i=1:.n_iterations,
                               .combine='bind_rows',
                               .inorder = FALSE,
                               .errorhandling = "remove") %dopar% {
                                 run_single_sim_tte(i,
                                                    batch_size = .batch_size,
                                                    max_ss = .max_ss,
                                                    initial_ss = .initial_ss,
                                                    max_time = .max_time,
                                                    enrollment_time_cutoff = .enrollment_time_cutoff,
                                                    initial_events_req_ia = .initial_events_req_ia,
                                                    new_events_req_ia = .new_events_req_ia,
                                                    effect_treatment = .effect_treatment,
                                                    beta_1 = .beta_1,
                                                    beta_2 = .beta_2,
                                                    beta_3 = .beta_3,
                                                    beta_4 = .beta_4,
                                                    beta_5 = .beta_5,
                                                    p_sup_thresh = .p_sup_thresh,
                                                    full_data = .data,
                                                    model = model_unadjusted)
                               }
  
  ## correct prior
  .model_correct_prior <- foreach(i=1:.n_iterations,
                            .combine='bind_rows',
                            .inorder = FALSE,
                            .errorhandling = "remove") %dopar% {
                              run_single_sim_tte(i,
                                                 batch_size = .batch_size,
                                                 max_ss = .max_ss,
                                                 initial_ss = .initial_ss,
                                                 max_time = .max_time,
                                                 enrollment_time_cutoff = .enrollment_time_cutoff,
                                                 initial_events_req_ia = .initial_events_req_ia,
                                                 new_events_req_ia = .new_events_req_ia,
                                                 effect_treatment = .effect_treatment,
                                                 beta_1 = .beta_1,
                                                 beta_2 = .beta_2,
                                                 beta_3 = .beta_3,
                                                 beta_4 = .beta_4,
                                                 beta_5 = .beta_5,
                                                 p_sup_thresh = .p_sup_thresh,
                                                 full_data = .data,
                                                 model = model_correct_prior)
                            }
  
  ## correct strong prior
  .model_correct_prior_strong <- foreach(i=1:.n_iterations,
                                  .combine='bind_rows',
                                  .inorder = FALSE,
                                  .errorhandling = "remove") %dopar% {
                                    run_single_sim_tte(i,
                                                       batch_size = .batch_size,
                                                       max_ss = .max_ss,
                                                       initial_ss = .initial_ss,
                                                       max_time = .max_time,
                                                       enrollment_time_cutoff = .enrollment_time_cutoff,
                                                       initial_events_req_ia = .initial_events_req_ia,
                                                       new_events_req_ia = .new_events_req_ia,
                                                       effect_treatment = .effect_treatment,
                                                       beta_1 = .beta_1,
                                                       beta_2 = .beta_2,
                                                       beta_3 = .beta_3,
                                                       beta_4 = .beta_4,
                                                       beta_5 = .beta_5,
                                                       p_sup_thresh = .p_sup_thresh,
                                                       full_data = .data,
                                                       model = model_correct_prior_strong)
                                  }
  
  stopImplicitCluster()
  
  # return tibble which includes data the model was run on and its results
  .model_correct %>%
    bind_rows(.model_incorrect) %>%
    bind_rows(.model_noise) %>%
    bind_rows(.model_unadjusted) %>%
    bind_rows(.model_correct_prior) %>%
    bind_rows(.model_correct_prior_strong)
}

# on LOG-hazard scale
# make separate scripts fo N=200-1000 due to computation/time constraints
trt_effect_100 <- log(c(0.7317716, 0.6517841, 1))


model_pars <- tibble(max_ss = c(rep(100, length(trt_effect_100))),
                                # rep(200, length(trt_effect_200)),
                                # rep(500, length(trt_effect_500)),
                                # rep(1000, length(trt_effect_1000))),
                     batch_size = c(rep(20, length(trt_effect_100))),
                                    # rep(40, length(trt_effect_200)),
                                    # rep(100, length(trt_effect_500)),
                                    # rep(200, length(trt_effect_1000))),
                     initial_ss = c(rep(20, length(trt_effect_100))),
                                    # rep(40, length(trt_effect_200)),
                                    # rep(100, length(trt_effect_500)),
                                    # rep(200, length(trt_effect_1000))),
                     initial_events_req_ia = c(rep(10, length(trt_effect_100))),
                                               # rep(20, length(trt_effect_200)),
                                               # rep(50, length(trt_effect_500)),
                                               # rep(100, length(trt_effect_1000))),
                     new_events_req_ia = c(rep(10, length(trt_effect_100))),
                                           # rep(20, length(trt_effect_200)),
                                           # rep(50, length(trt_effect_500)),
                                           # rep(100, length(trt_effect_1000))),
                     beta_1 = 1,
                     beta_2 = -0.5,
                     beta_3 = 1,
                     beta_4 = -0.1,
                     beta_5 = 0.5)

model_pars <- model_pars %>%
  bind_cols(effect_treatment  = c(trt_effect_100)) #, trt_effect_200, trt_effect_500, trt_effect_1000))

### Run the full simulation

sim_res <- foreach(j = 1:nrow(model_pars),
                   .errorhandling = "remove",
                   .combine = 'rbind') %do% {
                     sim_res <- complete_sim_tte(.n_iterations = 1000,
                                      .max_ss = pull(model_pars[j, "max_ss"]),
                                      .initial_ss =pull(model_pars[j, "initial_ss"]),
                                      .max_time = 75,
                                      .enrollment_time_cutoff = 50,
                                      .initial_events_req_ia = pull(model_pars[j, "initial_events_req_ia"]),
                                      .new_events_req_ia = pull(model_pars[j, "new_events_req_ia"]),
                                      .batch_size = pull(model_pars[j, "batch_size"]),
                                      .effect_treatment = pull(model_pars[j, "effect_treatment"]),
                                      .beta_1 = pull(model_pars[j, "beta_1"]),
                                      .beta_2 = pull(model_pars[j, "beta_2"]),
                                      .beta_3 = pull(model_pars[j, "beta_3"]),
                                      .beta_4 = pull(model_pars[j, "beta_4"]),
                                      .beta_5 = pull(model_pars[j, "beta_5"]),
                                      .p_sup_thresh = 0.99)
        }

saveRDS(sim_res, "PATH/FILENAME.RDS")

