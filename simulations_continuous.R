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

generate_outcomes <- function(data, effect_treatment, beta_1, beta_2, beta_3, beta_4, beta_5, max_ss){
  
  data %>%
    mutate(mu = effect_treatment*treatment + beta_1*x1 + beta_2*x2 + beta_3*x3 +
             beta_4*x4 + beta_5*x5,
           y = rnorm(n = max_ss, mean = mu)) #sd = 1
} 



### Models

## correct
model_correct <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5,
                                family = gaussian,
                                data =  sim_data,
                                chains = 3))

## correct prior 
model_correct_prior <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5,
                                             family = gaussian,
                                             data =  sim_data,
                                             prior = prior_correct,
                                             chains = 3))
## correct strong prior
model_correct_prior_strong <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5,
                                      family = gaussian,
                                      data =  sim_data,
                                      prior = prior_correct_strong,
                                      chains = 3))

### no quad 
model_incorrect <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x5,
                                  family = gaussian,
                                  data = sim_data,
                                  chains = 3))

## correct noise
model_noise <- quote(stan_glm(y ~ treatment + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
                              family = gaussian,
                              data =  sim_data,
                              chains = 3))

## drops strong prognostic variables (not in manuscript)
model_no_strong_prog <- quote(stan_glm(y ~ treatment + x2 + x5,
                                       family = gaussian,
                                       data =  sim_data,
                                       chains = 3))

## drops strong prognostic variables and adds noise (not in manuscript)
model_no_strong_prog_noise <- quote(stan_glm(y ~ treatment + x2 + x5 + x6 + x7 + x8,
                                             family = gaussian,
                                             data =  sim_data,
                                             chains = 3))

## unadjusted
model_unadjusted <- quote(stan_glm(y ~ treatment,
                                   family = gaussian,
                                   data = sim_data,
                                   chains = 3))

##marginalizes conditional samples 
get_marginal_effect <- function(fitted_model){
  
  mod_df <- model.frame(fitted_model)
  mod_df$treatment <- 0
  exp_y_control <- rowMeans(posterior_epred(fitted_model, newdata = mod_df))
  mod_df$treatment <- 1
  exp_y_treatment <- rowMeans(posterior_epred(fitted_model, newdata = mod_df))
  # posterior of relative risk
  exp_y_treatment - exp_y_control
}

results_output <- quote(
  tibble(iteration = iteration, 
         p_sup_thresh = p_sup_thresh,
         batch_size = batch_size,
         max_ss = max_ss,
         model = !!mod_name,
         effect_treatment = effect_treatment,
         beta_1 = beta_1,
         beta_2 = beta_2,
         beta_3 = beta_3,
         beta_4 = beta_4,
         beta_5 = beta_5,
         runtime_sec = toc - tic,
         n_total = n_total, 
         n_analyses_total = n_analyses_total,
         p_sup = p_sup, 
         superiority = if_else(p_sup > p_sup_thresh, 1, 0),
         reach_max_ss = if_else(n_total == max_ss, 1, 0),
         stopped_early = if_else(reach_max_ss == 0, 1, 0), # same info as above just flipped
         trt_est_mean = mean(results, na.rm = TRUE), # using posterior MEAN
         trt_est_median = median(results, na.rm = TRUE),
         rmse = sqrt(mean((results - effect_treatment)**2)),
         post_var = var(results),
         mae = mean(abs(results - effect_treatment)),
         trt_posterior = list(results %>% as_tibble_col(column_name = "treatment")),
         stan_summary = list(as_tibble(fit_mod$stan_summary, rownames = "parameter"))
))


## single trial simulation
run_single_sim <- function(iteration, 
                           batch_size,
                           max_ss = 100,
                           p_sup_thresh = 0.99,
                           effect_treatment,
                           beta_1, beta_2, beta_3, beta_4, beta_5,
                           full_data,
                           model){
  
  ### Tracks name of the models used
  mod_name <- rlang::as_label(enquo(model))
  
  ### Pull out one full dataset from list
  iteration_data <- full_data[[iteration]]
  
  ### Enroll initial participants
  sim_data <- iteration_data[1:batch_size, ]
  
  n_total <- batch_size
  p_sup <- 0 
  n_analyses_total <- 0
  
  tic <- tictoc::tic()
  while(n_total < max_ss){
    
    # fits the model 
    sd_y <- sd(sim_data$y)
    sd_treatment <- sd(sim_data$treatment)
    sd_x1 <- sd(sim_data$x1)
    sd_x2 <- sd(sim_data$x2)
    sd_x3 <- sd(sim_data$x3)
    sd_x4 <- sd(sim_data$x4)
    sd_x5 <- sd(sim_data$x5)

    # centered at DGM except for trt
    prior_correct <- normal(location = c(0, 0.5, -0.25, 0.5, -0.05, 0.25),
                            scale = 2.5,
                            autoscale = TRUE)
    # centered and scaled
    prior_correct_strong <- normal(location = c(0, 0.5, -0.25, 0.5, -0.05, 0.25),
                                   scale = c(2.5/sd_treatment, sd_y/sd_x1, sd_y/sd_x2, sd_y/sd_x3,
                                             sd_y /sd_x4, sd_y/sd_x5),
                                   autoscale = FALSE)
    
    n_total <- length(sim_data$y)
    n_analyses_total <- n_analyses_total + 1
    fit_mod <- eval(model)
    
    # posterior samples
    results <- get_marginal_effect(fit_mod)
    # probability of superiority
    p_sup <- mean(results < 0)
    
    # p_sup_thresh = u in manuscript
    if(p_sup > p_sup_thresh){
      break
    }
    
    if(n_total >= max_ss){
      break
    }
    
    # appends new data
    sim_data <- iteration_data[1:(n_total + batch_size), ]
  }
  
  toc <- tictoc::tic()
  
  # produce final results
  
  eval(results_output)
  
}

## Wrapper to run .n_iterations trial simulations

complete_sim_continuous <- function(.n_iterations,
                                      .max_ss, 
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
                              run_single_sim(iteration = i,
                                             batch_size = .batch_size,
                                             max_ss = .max_ss,
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
  
  ## correct prior
  .model_correct_prior <- foreach(i=1:.n_iterations,
                            .combine='bind_rows',
                            .inorder = FALSE,
                            .errorhandling = "remove") %dopar% {
                              run_single_sim(iteration = i,
                                             batch_size = .batch_size,
                                             max_ss = .max_ss,
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
                                    run_single_sim(iteration = i, 
                                                   batch_size = .batch_size, 
                                                   max_ss = .max_ss,
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

  
  ## no quad
  .model_incorrect <- foreach(i=1:.n_iterations,
                              .combine='bind_rows',
                              .inorder = FALSE,
                              .errorhandling = "remove") %dopar% {
                                run_single_sim(iteration = i,
                                               batch_size = .batch_size,
                                               max_ss = .max_ss,
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
                            run_single_sim(iteration = i,
                                           batch_size = .batch_size,
                                           max_ss = .max_ss,
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

  ## no strong prog
  .model_no_strong_prog <- foreach(i=1:.n_iterations,
                                   .combine='bind_rows',
                                   .inorder = FALSE,
                                   .errorhandling = "remove") %dopar% {
                                     run_single_sim(iteration = i,
                                                    batch_size = .batch_size,
                                                    max_ss = .max_ss,
                                                    effect_treatment = .effect_treatment,
                                                    beta_1 = .beta_1,
                                                    beta_2 = .beta_2,
                                                    beta_3 = .beta_3,
                                                    beta_4 = .beta_4,
                                                    beta_5 = .beta_5,
                                                    p_sup_thresh = .p_sup_thresh,
                                                    full_data = .data,
                                                    model = model_no_strong_prog)
                                   }

  ## no strong prog noise
  .model_no_strong_prog_noise <- foreach(i=1:.n_iterations,
                                         .combine='bind_rows',
                                         .inorder = FALSE,
                                         .errorhandling = "remove") %dopar% {
                                           run_single_sim(iteration = i,
                                                          batch_size = .batch_size,
                                                          max_ss = .max_ss,
                                                          effect_treatment = .effect_treatment,
                                                          beta_1 = .beta_1,
                                                          beta_2 = .beta_2,
                                                          beta_3 = .beta_3,
                                                          beta_4 = .beta_4,
                                                          beta_5 = .beta_5,
                                                          p_sup_thresh = .p_sup_thresh,
                                                          full_data = .data,
                                                          model = model_no_strong_prog_noise)
                                         }
  
  ## unadjusted
  .model_unadjusted <- foreach(i=1:.n_iterations,
                               .combine='bind_rows',
                               .inorder = FALSE,
                               .errorhandling = "remove") %dopar% {
                                 run_single_sim(iteration = i,
                                                batch_size = .batch_size,
                                                max_ss = .max_ss,
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

## trt effects, null, approx 50% power and 80% power for unadjusted model
trt_effect_100 <- c(0, -0.73, -0.53)
trt_effect_200 <- c(0, -0.52, -0.35)
trt_effect_500 <- c(0, -0.33, -0.22)
trt_effect_1000 <- c(0, -0.23, -0.16)

model_pars <- tibble(max_ss = c(rep(100, length(trt_effect_100)),
                                rep(200, length(trt_effect_200)),
                                rep(500, length(trt_effect_500)),
                                rep(1000, length(trt_effect_1000))),
                     batch_size = c(rep(20, length(trt_effect_100)),
                                    rep(40, length(trt_effect_200)),
                                    rep(100, length(trt_effect_500)),
                                    rep(200, length(trt_effect_1000))),
                     beta_1 = 0.5,
                     beta_2 = -0.25,
                     beta_3 = 0.5,
                     beta_4 = -0.05, 
                     beta_5 = 0.25)

model_pars <- model_pars %>%
  bind_cols(trt_effect = c(trt_effect_100, trt_effect_200, trt_effect_500, trt_effect_1000))


### Run the full simulation

sim_res <- foreach(j = 1:nrow(model_pars),
                   .errorhandling = "remove",
                   .combine = 'rbind') %do% {
                     complete_sim_continuous(
                       .n_iterations = 1000,
                       .max_ss = pull(model_pars[j, "max_ss"]), 
                       .batch_size = pull(model_pars[j, "batch_size"]), 
                       .effect_treatment = pull(model_pars[j, "trt_effect"]), 
                       .beta_1 = pull(model_pars[j, "beta_1"]),
                       .beta_2 = pull(model_pars[j, "beta_2"]),
                       .beta_3 = pull(model_pars[j, "beta_3"]),
                       .beta_4 = pull(model_pars[j, "beta_4"]),
                       .beta_5 = pull(model_pars[j, "beta_5"]),
                       .p_sup_thresh = 0.99)
                   }

saveRDS(sim_res, "PATH")


