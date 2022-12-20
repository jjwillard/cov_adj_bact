pacman::p_load(tidyr, rstanarm, foreach, doParallel, tibble, purrr, dplyr)

generate_data <- function(iteration, max_ss){
  tibble(treatment = rbinom(n = max_ss, size = 1, p= 0.5),
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

generate_outcomes <- function(data, effect_treatment, b_0_c, beta_1, beta_2, beta_3, beta_4, beta_5, control_risk, max_ss){

  # corresponds to values found through optimization
  conditional_params <- tibble(max_ss = c(100, 100, 100, 200, 200, 200, 500, 500, 500, 1000, 1000, 1000),
                               rr_m = c(0.5180120, 0.4447808, 1, 0.5773397, 0.4075299, 1,
                                        0.7136076, 0.5976614, 1, 0.7947651, 0.7196724, 1),
                               lo_c = c(-1.0344828, -1.2413793, 0, -0.8793103, -1.3600000, 0,
                                        -0.5600000, -0.8275862, 0, -0.3925000, -0.5517241, 0))

  lo_c <- conditional_params %>%
    filter(max_ss == max_ss, rr_m == effect_treatment) %>%
    pull(lo_c)
  
  data %>%
    mutate(eta = b_0_c + lo_c*treatment + beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4 + beta_5*x5,
           p = logit_inverse(eta),
           y = rbinom(n = max_ss, size = 1, p))
  
}

calculate_events <- function(data, type = "total"){
  events <-  data %>%
    group_by(treatment) %>%
    summarise(n_events = sum(y)) 
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

### Models

## correct
model_correct <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5,
                                family = binomial(link = "logit"),
                                data = sim_data,
                                chains = 3))

## no quad
model_incorrect <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x5,
                                  family = binomial(link = "logit"),
                                  data = sim_data,
                                  chains = 3))

## correct noise
model_noise <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
                              family = binomial(link = "logit"),
                              data = sim_data,
                              chains = 3))

## no strong prog
model_no_strong_prog <- quote(stan_glm(y ~ treatment + x2 + x5,
                                       family = binomial(link = "logit"),
                                       data = sim_data,
                                       chains = 3))

## no strong prog noise
model_no_strong_prog_noise <- quote(stan_glm(y ~ treatment + x2 + x5 + x6 + x7 + x8,
                                             family = binomial(link = "logit"),
                                             data = sim_data,
                                             chains = 3))

## unadjusted
model_unadjusted <- quote(stan_glm(y ~ treatment,
                                   family = binomial(link = "logit"),
                                   data = sim_data,
                                   chains = 3))

## correct prior
model_correct_prior <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5,
                                family = binomial(link = "logit"),
                                data = sim_data,
                                prior = prior_correct,
                                chains = 3))

## correct strong prior
model_correct_prior_strong <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5,
                                      family = binomial(link = "logit"),
                                      data = sim_data,
                                      prior = prior_correct_strong,
                                      chains = 3))

## marginalize conditional posterior samples
get_relative_risk <- function(fitted_model){
  
  mod_df <- model.frame(fitted_model)
  mod_df$treatment <- 0
  probs_control <- rowMeans(posterior_epred(fitted_model, newdata = mod_df))
  mod_df$treatment <- 1
  probs_treatment <- rowMeans(posterior_epred(fitted_model, newdata = mod_df))
  
  # posterior of relative risk
  probs_treatment / probs_control
}

results_output <- quote(
  tibble(iteration = iteration,
         p_sup_thresh = p_sup_thresh,
         batch_size = batch_size,
         max_ss = max_ss,
         initial_ss = initial_ss,
         initial_events_req_ia = initial_events_req_ia,
         new_events_req_ia = new_events_req_ia,
         model = !!mod_name,
         effect_treatment = effect_treatment,
         b_0_c = b_0_c,
         beta_1 = beta_1,
         beta_2 = beta_2,
         beta_3 = beta_3,
         beta_4 = beta_4,
         beta_5 = beta_5,
         control_risk = control_risk,
         runtime_sec = toc - tic,
         n_total = n_total, 
         n_analyses_total = n_analyses_total,
         n_events_total = n_events,
         p_sup = p_sup, 
         superiority = if_else(p_sup > p_sup_thresh, 1, 0),
         reach_max_ss = if_else(n_total >= max_ss, 1, 0), #ge as failsafe in case bug introduce to increment past intended size
         stopped_early = if_else(reach_max_ss == 0, 1, 0),
         trt_est_mean = mean(results, na.rm = TRUE),
         trt_est_mean_log = mean(log(results), na.rm = TRUE),
         trt_est_median = median(results, na.rm = TRUE),
         trt_est_median_log = median(log(results), na.rm = TRUE),
         rmse = sqrt(mean((results - effect_treatment)**2)),
         rmse_log = sqrt(mean((log(results) - log(effect_treatment))**2)),
         mae = mean(abs(results - effect_treatment)),
         post_var = var(results),
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
         initial_events_req_ia = initial_events_req_ia,
         new_events_req_ia = new_events_req_ia,
         model = !!mod_name,
         effect_treatment = effect_treatment,
         b_0_c = b_0_c,
         beta_1 = beta_1,
         beta_2 = beta_2,
         beta_3 = beta_3,
         beta_4 = beta_4,
         beta_5 = beta_5,
         control_risk = control_risk,
         runtime_sec = toc - tic,
         n_total = n_total, 
         n_analyses_total = n_analyses_total,
         n_events_total = n_events,
         p_sup = NA_real_, 
         superiority = NA_real_,
         reach_max_ss = NA_real_, #ge as failsafe in case bug introduce to increment past intended size
         stopped_early = NA_real_,
         trt_est_mean = NA_real_,
         trt_est_mean_log = NA_real_,
         trt_est_median = NA_real_,
         trt_est_median_log = NA_real_,
         rmse = NA_real_,
         rmse_log = NA_real_,
         mae = NA_real_,
         post_var = NA_real_,
         message = "Insufficient number of events for interim analyses.",
         trt_posterior = list(tibble(treatment = NA_real_)),
         stan_summary = list(tibble(stan_summary = NA_real_))
))

## single trial simulation
run_single_sim_binary <- function(iteration, 
                                  batch_size = 20, # this is recruitment rate per unit time
                                  max_ss = 100,
                                  initial_ss = 20,
                                  initial_events_req_ia = 10, # initial number of events for first interim analysis
                                  new_events_req_ia = 10, # number of additional events per additional interim analysis
                                  p_sup_thresh = 0.99,
                                  effect_treatment, 
                                  b_0_c,
                                  control_risk,
                                  beta_1 = NULL, beta_2 = NULL, beta_3 = NULL, beta_4 = NULL, beta_5 = NULL,
                                  full_data,
                                  model){
  
  ### Tracks name of the models used
  mod_name <- rlang::as_label(enquo(model))
  
  ### Pull out one full dataset from list
  iteration_data <- full_data[[iteration]]
  
  ### Enroll initial participants
  sim_data <- iteration_data[1:initial_ss, ] 
  
  ### initial events required for interim analysis to be performed
  events_req_ia <- initial_events_req_ia
  
  ### assuming that outcomes/events are observed immediately
  n_total <- initial_ss
  p_sup <- 0 
  n_events <- calculate_events(sim_data, type = "total")
  n_analyses_total <- 0
  
  tic <- tictoc::tic()
  
  ## set up two while loops to keep trial within proper size and only perform interim analysis after
  ## pre-specified number of events
  while(n_total < max_ss){
    
    ### Check for specified number of events before performing interim analysis
    while(n_events < events_req_ia & n_total < max_ss){
      ## enroll more participants, since assuming outcomes observed immediately
      ## batch_size = recruitment rate per unit time
      new_data <- iteration_data[(n_total + 1):(n_total + batch_size),]
      sim_data <- bind_rows(sim_data, new_data)
      n_total <- as.numeric(length(sim_data$treatment))
      n_events <- calculate_events(sim_data, type = "total")
    } 
    
    ### sets a break if 0 events occur by end
    
    if(n_events == 0 & n_total >= max_ss){
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
    results <- get_relative_risk(fit_mod)
    
    ## on relative risk scale
    p_sup <- mean(results < 1)

    ## Determine whether or not to stop based on stopping rule
    if(p_sup > p_sup_thresh){
      break
    }
    
    ## Stop if reach maximium sample size 
    if(n_total >= max_ss){
      break
    }
    
    ### If not stopping, then increase number of events required and repeat
    ## Increase number of events over current number to be required at next interim analysis
    events_req_ia <- n_events + new_events_req_ia
    
  }
  toc <- tictoc::tic()
  
  ### produce final results
  eval(results_output)
  
}

### Wrapper function to run .n_iterations trial simulations

complete_sim_binary <- function(.n_iterations,
                                .max_ss, 
                                .initial_ss,
                                .initial_events_req_ia,
                                .new_events_req_ia,
                                .batch_size,
                                .effect_treatment,
                                .b_0_c,
                                .beta_1, .beta_2, .beta_3, .beta_4, .beta_5, .control_risk,
                                .p_sup_thresh = 0.99,
                                .seed = 123,
                                .n_cores = parallel::detectCores()){
  
  set.seed(.seed, kind = "L'Ecuyer-CMRG")
  .data_structure <- foreach(i = 1:.n_iterations,
                             .errorhandling = "remove") %do% {
                               generate_data(i, .max_ss)
                             }
  
  set.seed(.seed, kind = "L'Ecuyer-CMRG")
  .data <- foreach(j = 1:.n_iterations) %do% {
    generate_outcomes(.data_structure[[j]], .effect_treatment, .b_0_c, .beta_1, .beta_2, 
                      .beta_3, .beta_4, .beta_5, .control_risk, .max_ss)
  }
  
  registerDoParallel(cores = .n_cores)
  
  ## correct
  .model_correct <- foreach(i=1:.n_iterations,
                            .combine='bind_rows',
                            .inorder = FALSE,
                            .errorhandling = "remove") %dopar% {
                              run_single_sim_binary(iteration = i,
                                                    batch_size = .batch_size,
                                                    max_ss = .max_ss,
                                                    initial_ss = .initial_ss,
                                                    initial_events_req_ia = .initial_events_req_ia,
                                                    new_events_req_ia = .new_events_req_ia,
                                                    effect_treatment = .effect_treatment,
                                                    b_0_c = .b_0_c,
                                                    beta_1 = .beta_1,
                                                    beta_2 = .beta_2,
                                                    beta_3 = .beta_3,
                                                    beta_4 = .beta_4,
                                                    beta_5 = .beta_5,
                                                    control_risk = .control_risk,
                                                    p_sup_thresh = .p_sup_thresh,
                                                    full_data = .data,
                                                    model = model_correct)
                            }

  ## no quad
  .model_incorrect <- foreach(i=1:.n_iterations,
                              .combine='bind_rows',
                              .inorder = FALSE,
                              .errorhandling = "remove") %dopar% {
                                run_single_sim_binary(iteration = i,
                                                      batch_size = .batch_size,
                                                      max_ss = .max_ss,
                                                      initial_ss = .initial_ss,
                                                      initial_events_req_ia = .initial_events_req_ia,
                                                      new_events_req_ia = .new_events_req_ia,
                                                      effect_treatment = .effect_treatment,
                                                      b_0_c = .b_0_c,
                                                      beta_1 = .beta_1,
                                                      beta_2 = .beta_2,
                                                      beta_3 = .beta_3,
                                                      beta_4 = .beta_4,
                                                      beta_5 = .beta_5,
                                                      control_risk = .control_risk,
                                                      p_sup_thresh = .p_sup_thresh,
                                                      full_data = .data,
                                                      model = model_incorrect)
                              }
  
  ## correct noise
  .model_noise <- foreach(i=1:.n_iterations,
                          .combine='bind_rows',
                          .inorder = FALSE,
                          .errorhandling = "remove") %dopar% {
                            run_single_sim_binary(iteration = i,
                                                  batch_size = .batch_size,
                                                  max_ss = .max_ss,
                                                  initial_ss = .initial_ss,
                                                  initial_events_req_ia = .initial_events_req_ia,
                                                  new_events_req_ia = .new_events_req_ia,
                                                  effect_treatment = .effect_treatment,
                                                  b_0_c = .b_0_c,
                                                  beta_1 = .beta_1,
                                                  beta_2 = .beta_2,
                                                  beta_3 = .beta_3,
                                                  beta_4 = .beta_4,
                                                  beta_5 = .beta_5,
                                                  control_risk = .control_risk,
                                                  p_sup_thresh = .p_sup_thresh,
                                                  full_data = .data,
                                                  model = model_noise)
                          }

  ## no strong prog
  .model_no_strong_prog <- foreach(i=1:.n_iterations,
                                   .combine='bind_rows',
                                   .inorder = FALSE,
                                   .errorhandling = "remove") %dopar% {
                                     run_single_sim_binary(iteration = i,
                                                           batch_size = .batch_size,
                                                           max_ss = .max_ss,
                                                           initial_ss = .initial_ss,
                                                           initial_events_req_ia = .initial_events_req_ia,
                                                           new_events_req_ia = .new_events_req_ia,
                                                           effect_treatment = .effect_treatment,
                                                           b_0_c = .b_0_c,
                                                           beta_1 = .beta_1,
                                                           beta_2 = .beta_2,
                                                           beta_3 = .beta_3,
                                                           beta_4 = .beta_4,
                                                           beta_5 = .beta_5,
                                                           control_risk = .control_risk,
                                                           p_sup_thresh = .p_sup_thresh,
                                                           full_data = .data,
                                                           model = model_no_strong_prog)
                                   }

  ## no strong prog noise
  .model_no_strong_prog_noise <- foreach(i=1:.n_iterations,
                                         .combine='bind_rows',
                                         .inorder = FALSE,
                                         .errorhandling = "remove") %dopar% {
                                           run_single_sim_binary(iteration = i,
                                                                 batch_size = .batch_size,
                                                                 max_ss = .max_ss,
                                                                 initial_ss = .initial_ss,
                                                                 initial_events_req_ia = .initial_events_req_ia,
                                                                 new_events_req_ia = .new_events_req_ia,
                                                                 effect_treatment = .effect_treatment,
                                                                 b_0_c = .b_0_c,
                                                                 beta_1 = .beta_1,
                                                                 beta_2 = .beta_2,
                                                                 beta_3 = .beta_3,
                                                                 beta_4 = .beta_4,
                                                                 beta_5 = .beta_5,
                                                                 control_risk = .control_risk,
                                                                 p_sup_thresh = .p_sup_thresh,
                                                                 full_data = .data,
                                                                 model = model_no_strong_prog_noise)
                                         }
  
  ## unadjusted
  .model_unadjusted <- foreach(i=1:.n_iterations,
                               .combine='bind_rows',
                               .inorder = FALSE,
                               .errorhandling = "remove") %dopar% {
                                 run_single_sim_binary(iteration = i,
                                                       batch_size = .batch_size,
                                                       max_ss = .max_ss,
                                                       initial_ss = .initial_ss,
                                                       initial_events_req_ia = .initial_events_req_ia,
                                                       new_events_req_ia = .new_events_req_ia,
                                                       effect_treatment = .effect_treatment,
                                                       b_0_c = .b_0_c,
                                                       beta_1 = .beta_1,
                                                       beta_2 = .beta_2,
                                                       beta_3 = .beta_3,
                                                       beta_4 = .beta_4,
                                                       beta_5 = .beta_5,
                                                       control_risk = .control_risk,
                                                       p_sup_thresh = .p_sup_thresh,
                                                       full_data = .data,
                                                       model = model_unadjusted)
                               }
  
  ## correct prior
  .model_correct_prior <- foreach(i=1:.n_iterations,
                                    .combine='bind_rows',
                                    .inorder = FALSE,
                                    .errorhandling = "remove") %dopar% {
                                      run_single_sim_binary(iteration = i,
                                                            batch_size = .batch_size,
                                                            max_ss = .max_ss,
                                                            initial_ss = .initial_ss,
                                                            initial_events_req_ia = .initial_events_req_ia,
                                                            new_events_req_ia = .new_events_req_ia,
                                                            effect_treatment = .effect_treatment,
                                                            b_0_c = .b_0_c,
                                                            beta_1 = .beta_1,
                                                            beta_2 = .beta_2,
                                                            beta_3 = .beta_3,
                                                            beta_4 = .beta_4,
                                                            beta_5 = .beta_5,
                                                            control_risk = .control_risk,
                                                            p_sup_thresh = .p_sup_thresh,
                                                            full_data = .data,
                                                            model = model_correct_prior)
                                    }

  ## correct strong prior
  .model_correct_prior_strong <- foreach(i=1:.n_iterations,
                                       .combine='bind_rows',
                                       .inorder = FALSE,
                                       .errorhandling = "remove") %dopar% {
                                         run_single_sim_binary(iteration = i,
                                                               batch_size = .batch_size,
                                                               max_ss = .max_ss,
                                                               initial_ss = .initial_ss,
                                                               initial_events_req_ia = .initial_events_req_ia,
                                                               new_events_req_ia = .new_events_req_ia,
                                                               effect_treatment = .effect_treatment,
                                                               b_0_c = .b_0_c,
                                                               beta_1 = .beta_1,
                                                               beta_2 = .beta_2,
                                                               beta_3 = .beta_3,
                                                               beta_4 = .beta_4,
                                                               beta_5 = .beta_5,
                                                               control_risk = .control_risk,
                                                               p_sup_thresh = .p_sup_thresh,
                                                               full_data = .data,
                                                               model = model_correct_prior_strong)
                                       }

  stopImplicitCluster()
  
  # return tibble which includes data the model was run on and its results
  .model_correct %>%
    bind_rows(.model_incorrect) %>%
    bind_rows(.model_noise) %>%
    bind_rows(.model_no_strong_prog) %>%
    bind_rows(.model_no_strong_prog_noise) %>%
    bind_rows(.model_unadjusted) %>%
    bind_rows(.model_correct_prior) %>%
    bind_rows(.model_correct_prior_strong)
}

# on relative risk scale
trt_effect_100 <- c(1, 0.5180120, 0.4447808) #null, 30% and 40% power
trt_effect_200 <- c(1, 0.5773397, 0.4075299) #null, 50 and 80% power for rest
trt_effect_500 <- c(1, 0.7136076, 0.5976614)
trt_effect_1000 <- c(1, 0.7947651, 0.7196724)

model_pars <- tibble(max_ss = c(rep(100, length(trt_effect_100)),
                                rep(200, length(trt_effect_200)),
                                rep(500, length(trt_effect_500)),
                                rep(1000, length(trt_effect_1000))),
                     # for CER=0.3
                     b_0_c = c(rep(-1.265619, length(trt_effect_100)),
                               rep(-1.264603, length(trt_effect_200)),
                               rep(-1.267382, length(trt_effect_500)),
                               rep(-1.265569, length(trt_effect_1000))),
                     batch_size = c(rep(20, length(trt_effect_100)),
                                    rep(40, length(trt_effect_200)),
                                    rep(100, length(trt_effect_500)),
                                    rep(200, length(trt_effect_1000))),
                     initial_ss = c(rep(20, length(trt_effect_100)),
                                    rep(40, length(trt_effect_200)),
                                    rep(100, length(trt_effect_500)),
                                    rep(200, length(trt_effect_1000))),
                     initial_events_req_ia = c(rep(10, length(trt_effect_100)),
                                               rep(20, length(trt_effect_200)),
                                               rep(50, length(trt_effect_500)),
                                               rep(100, length(trt_effect_1000))),
                     new_events_req_ia = c(rep(10, length(trt_effect_100)),
                                           rep(20, length(trt_effect_200)),
                                           rep(50, length(trt_effect_500)),
                                           rep(100, length(trt_effect_1000))),
                     beta_1 = 1,
                     beta_2 = -0.5,
                     beta_3 = 1,
                     beta_4 = -0.1, 
                     beta_5 = 0.5,
                     control_risk = 0.3)

model_pars <- model_pars %>%
  bind_cols(effect_treatment  = c(trt_effect_100, trt_effect_200, trt_effect_500, trt_effect_1000))


### Run the full simulation
sim_res <- foreach(j = 1:nrow(model_pars),
                   .errorhandling = "remove",
                   .combine = 'rbind') %do% {
                     complete_sim_binary(.n_iterations = 1000,
                                         .max_ss = pull(model_pars[j, "max_ss"]),  
                                         .initial_ss = pull(model_pars[j, "initial_ss"]), 
                                         .initial_events_req_ia = pull(model_pars[j, "initial_events_req_ia"]), 
                                         .new_events_req_ia = pull(model_pars[j, "new_events_req_ia"]), 
                                         .batch_size = pull(model_pars[j, "batch_size"]), 
                                         .effect_treatment = pull(model_pars[j, "effect_treatment"]), 
                                         .b_0_c = pull(model_pars[j, "b_0_c"]), 
                                         .beta_1 = pull(model_pars[j, "beta_1"]),
                                         .beta_2 = pull(model_pars[j, "beta_2"]),
                                         .beta_3 = pull(model_pars[j, "beta_3"]),
                                         .beta_4 = pull(model_pars[j, "beta_4"]),
                                         .beta_5 = pull(model_pars[j, "beta_5"]),
                                         .control_risk = pull(model_pars[j, "control_risk"]),
                                         .p_sup_thresh = 0.99)
                   }

saveRDS(sim_res, "PATH/FILENAME.RDS")

