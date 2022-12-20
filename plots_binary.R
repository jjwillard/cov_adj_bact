pacman::p_load(dplyr, foreach, ggplot2, cowplot, tidyr, purrr, rlang, kableExtra)

## Read in data
full_sim_res <- readRDS("PATH for sim_res from simulations_binary.R") %>%
  mutate(model = factor(model, levels = c("model_correct", "model_incorrect", "model_noise",
                                          "model_no_strong_prog", "model_no_strong_prog_noise", 
                                          "model_correct_prior", "model_correct_prior_strong", "model_unadjusted"),
                        labels = c("correct", "no quad", "correct noise", "no strong prog",
                                   "no strong prog noise", "correct prior", "correct strong prior",  "unadjusted")))

summarize_results <- function(complete_sim_results){
  complete_sim_results %>%
    group_by(model, max_ss, effect_treatment, beta_1, beta_2, beta_3, beta_4, beta_5, control_risk) %>%
    summarise(power = mean(superiority, na.rm = TRUE),
              bias_post_median = mean(trt_est_median - effect_treatment, na.rm = TRUE),
              bias_post_mean = mean(trt_est_mean - effect_treatment, na.rm = TRUE),
              bias_post_median_lrr = mean(trt_est_median_log - log(effect_treatment), na.rm = TRUE),
              bias_post_mean_lrr = mean(trt_est_mean_log - log(effect_treatment), na.rm = TRUE),
              p_superiority = mean(superiority, na.rm = TRUE),
              exp_sample_size = mean(n_total, na.rm = TRUE),
              p_stop_early = mean(stopped_early, na.rm = TRUE),
              .groups = "keep") %>%
    ungroup()
}

full_res_summary <- summarize_results(full_sim_res)


plot_summary_metric <- function(.data_list, 
                                metric = metric,
                                ylim = c(NA, NA), 
                                round_digits = 2,
                                xlab = 'relative risk',
                                xscale = "response",
                                ylab = NULL,
                                scales = 'free_x',
                                dodge_width = 0.5,
                                legend_position = "top",
                                n_rows_legend = 2){ 
  if(is.null(ylab)){
    metric_name <- rlang::as_name(metric)
  } else {
    metric_name <- ylab
  }
  
  if(xscale == "link"){
    .data_list <- .data_list %>%
      mutate(effect_treatment = log(effect_treatment))
  }
  
  metric <- sym(metric)
  
  max_ss_names <- c(`100` = "max ss = 100",
                    `200` = "max ss = 200",
                    `500` = "max ss = 500",
                    `1000` = "max ss = 1000")
  
  ggplot(.data_list, aes(x = factor(round(effect_treatment, round_digits)), 
                         y = !!metric, 
                         group = model, 
                         color = model, 
                         shape = model)) +
    coord_cartesian(ylim = ylim) +
    facet_wrap(~ max_ss, scales = scales,  nrow = 2,
               labeller = as_labeller(max_ss_names)) +
    geom_point(size = 1.5, position = position_dodge(width = dodge_width)) +
    scale_color_brewer(palette = 'Set2') +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),#,face="bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text.y = element_text(size = 8),
          strip.text.x = element_text(size = 8),
          legend.position = legend_position,
          legend.text = element_text(size = 10),
          legend.title = element_blank()) +
    guides(colour = guide_legend(nrow = n_rows_legend))
}

### Manuscript plot with power and pse ###-------------------------------------------

plot_legend <- get_legend(plot_summary_metric(full_res_summary %>%
                                                filter(model %in% c("correct", "no quad", "correct noise", 
                                                             "correct prior", "correct strong prior", "unadjusted")),
                                              metric = "power",
                                              xlab = 'difference in means',
                                              ylab = 'power',
                                              dodge_width = 0.5,
                                              legend_position = "top",
                                              n_rows_legend = 1))
 
plot(plot_legend)
prow <- plot_grid(plot_summary_metric(full_res_summary %>%
                                        filter(effect_treatment != 1,
                                               model %in% c("correct", "no quad", "correct noise", 
                                                            "correct prior", "correct strong prior", "unadjusted")),
                                      metric = "power",
                                      xlab = 'relative risk',
                                      ylab = 'power',
                                      scales = "free",
                                      dodge_width = 0.5,
                                      round_digits = 2,
                                      legend_position = "none"),
                  plot_summary_metric(full_res_summary %>%
                                        filter(effect_treatment != 1, 
                                               model %in% c("correct", "no quad", "correct noise", 
                                                            "correct prior", "correct strong prior", "unadjusted")),
                                      metric = "p_stop_early",
                                      xlab = 'relative risk',
                                      ylab = 'probability of stopping early',
                                      scales = "free",
                                      dodge_width = 0.5,
                                      round_digits = 2,
                                      legend_position = "none"),
                  labels = LETTERS[1:2],
                  align = "vh",
                  hjust = -1,
                  nrow = 1
)


power_pse <- plot_grid(prow, plot_legend, ncol = 1, rel_heights = c(1, .1))

ggsave2(filename = "PATH/FILENAME.eps",
        plot = power_pse,
        width = 8,
        height = 3,
        units = "in",
        dpi = 800)

#----------------------------------------------------------------------------------#

### False positive rate table
fpr <- full_res_summary %>%
  filter(effect_treatment == 1) %>%
  arrange(max_ss) %>%
  select(model, max_ss, power) %>%
  pivot_wider(names_from = "max_ss",
              values_from = "power")

## Bias under null
bias_null <- full_res_summary %>%
  filter(effect_treatment == 1) %>%
  arrange(max_ss) %>%
  select(model, max_ss, bias_post_median) %>%
  pivot_wider(names_from = "max_ss",
              values_from = "bias_post_median") %>%
  mutate(across(where(is.double), round, 3))

### Expected sample size
ess_100 <- full_res_summary %>%
  filter(max_ss == 100) %>%
  select(model, max_ss, effect_treatment, exp_sample_size) %>%
  pivot_wider(names_from = "effect_treatment",
              values_from = "exp_sample_size")

ess_200 <- full_res_summary %>%
  filter(max_ss == 200) %>%
  select(model, max_ss, effect_treatment, exp_sample_size) %>%
  pivot_wider(names_from = "effect_treatment",
              values_from = "exp_sample_size")

ess_500 <- full_res_summary %>%
  filter(max_ss == 500) %>%
  select(model, max_ss, effect_treatment, exp_sample_size) %>%
  pivot_wider(names_from = "effect_treatment",
              values_from = "exp_sample_size")

ess_1000 <- full_res_summary %>%
  filter(max_ss == 1000) %>%
  select(model, max_ss, effect_treatment, exp_sample_size) %>%
  pivot_wider(names_from = "effect_treatment",
              values_from = "exp_sample_size")

### Build table
### null is null, one is lowest power, two is highest power 
n_100 <- fpr %>% select(model, frp = `100`) %>%
  left_join(bias_null %>% select(model, bias_null = `100`), by = 'model') %>%
  left_join(ess_100 %>% select(-max_ss, null = `1`, one = `0.518012`, two = `0.4447808`) %>% 
              relocate(model, null, one, two), by = "model")

n_200 <- fpr %>% select(model, frp = `200`) %>%
  left_join(bias_null %>% select(model, bias_null = `200`), by = 'model') %>%
  left_join(ess_200 %>% select(-max_ss, null = `1`, one = `0.5773397`, two = `0.4075299`) %>% 
              relocate(model, null, one, two), by = "model")

n_500 <- fpr %>% select(model, frp = `500`) %>%
  left_join(bias_null %>% select(model, bias_null = `500`), by = 'model') %>%
  left_join(ess_500 %>% select(-max_ss, null = `1`, one = `0.7136076`, two = `0.5976614`) %>% 
              relocate(model, null, one, two), by = "model")

n_1000 <- fpr %>% select(model, frp = `1000`) %>%
  left_join(bias_null %>% select(model, bias_null = `1000`), by = 'model') %>%
  left_join(ess_1000 %>% select(-max_ss, null = `1`, one = `0.7947651`, two = `0.7196724`)  %>% 
              relocate(model, null, one, two), by = "model")

tbl_summary <- bind_rows(n_100  %>%
                           left_join(n_200, by = "model"),
                         n_500 %>%
                           left_join(n_1000, by = "model"))

tbl_summary <- tbl_summary %>%
  filter(model %in% c("correct", "no quad", "correct noise", "correct prior",
                      "correct strong prior", "unadjusted")) 

kbl(tbl_summary, booktabs = T, format = "latex") %>%
  add_header_above(c(" "=3, "Expected sample size"=3,
                     " "=2, "Expected sample size"=3)) %>%
  add_header_above(c(" "=1, "Maximum sample size = 100" = 5,
                     "Maximum sample size = 200" = 5)) %>%
  kable_styling(latex_options = c("scale_down"))


### Visualize rmse as boxplots

# https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

plot_trt_estimates_boxplot <- function(.data_list, 
                                       metric, 
                                       xlab = "treatment effect",
                                       ylab = NULL,
                                       scales = "free",
                                       round_digits = 2,
                                       ylim = c(NA, NA),
                                       legend_position = "top",
                                       n_rows_legend = 2){
  
  n_models <- length(unique(.data_list$model))
  
  if(is.null(ylab)){
    metric_name <- rlang::as_name(metric)
  } else {
    metric_name <- ylab
  }
  
  metric <- sym(metric)
  
  max_ss_names <- c(`100` = "max ss = 100",
                    `200` = "max ss = 200",
                    `500` = "max ss = 500",
                    `1000` = "max ss = 1000")
  
  
  ggplot(.data_list, aes(x = factor(round(effect_treatment, round_digits)), 
                         y = !!metric,
                         fill = model,
                         color = model)) + 
    stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", position = "dodge",
                 lwd = 0.25) + 
    facet_wrap(~ max_ss, scales = scales,  nrow = 2,
               labeller = as_labeller(max_ss_names)) +
    xlab(xlab) +
    ylab(metric_name) +
    scale_color_manual(values = rep("#003B46", n_models)) +
    scale_fill_brewer(palette = 'Set2') +
    theme_bw()+
    theme(axis.text.x= element_text(size = 10),#, angle = -90, hjust = 0),
          axis.title.x= element_text(size = 10),#, face = 'bold'),
          axis.title=element_text(size=10),#,face="bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text.y = element_text(size = 8),
          strip.text.x = element_text(size = 8),
          legend.position = legend_position,
          legend.text = element_text(size = 10),
          legend.title = element_blank()) +
    guides(colour = guide_legend(nrow = n_rows_legend))
}

### Manuscript appendix plot for bias and rmse

plot_legend <- get_legend(plot_summary_metric(full_res_summary %>%
                                                filter(effect_treatment != 1,
                                                       model %in% c("correct", "no quad", "correct noise", "correct prior",
                                                                    "correct strong prior", "unadjusted")),
                                              metric = "bias_post_median",
                                              xlab = 'difference in means',
                                              ylab = 'power',
                                              dodge_width = 0.5,
                                              legend_position = "top",
                                              n_rows_legend = 1))

plot(plot_legend)
prow <- plot_grid(plot_summary_metric(full_res_summary %>%
                                        filter(effect_treatment != 1,
                                               model %in% c("correct", "no quad", "correct noise", "correct prior",
                                                            "correct strong prior", "unadjusted")),
                                      metric = "bias_post_median",
                                      xlab = 'relative risk',
                                      ylab = 'posterior median bias',
                                      scales = "free",
                                      dodge_width = 0.5,
                                      round_digits = 2,
                                      legend_position = "none"),
                  plot_trt_estimates_boxplot(full_sim_res %>%
                                               filter(model %in% c("correct", "no quad", "correct noise", "correct prior",
                                                                   "correct strong prior", "unadjusted")), 
                                             metric = "rmse", 
                                             ylab = "root mean squared error",
                                             xlab = "relative risk",
                                             scales = "free",
                                             legend_position = "none"),
                  labels = LETTERS[1:2],
                  align = "vh",
                  hjust = -1,
                  nrow = 1
)
prow

bias_rmse <- plot_grid(prow, plot_legend, ncol = 1, rel_heights = c(1, .1))

ggsave2(filename = "PATH/FILENAME.eps",
        plot = bias_rmse,
        width = 8,
        height = 3,
        units = "in",
        dpi = 800)






















plot_trt_estimates_boxplot(full_sim_res %>%
                             filter(model %in% c("correct", "no quad", "correct noise", "no strong prog",
                                                           "no strong prog noise", "unadjusted")), metric = "rmse", 
                           ylab = "root mean squared error", 
                           ylim = c(0,NA),
                           xlab = "relative risk")

plot_grid(plot_trt_estimates_boxplot(full_sim_res, metric = "bias", xlab = "treatment effect (relative risk)"),
          plot_trt_estimates_boxplot(full_sim_res, metric = "rmse", 
                                     ylab = "root mean squared error", 
                                     ylim = c(0,NA),
                                     xlab = "treatment effect (relative risk)"),
          ncol = 2, nrow = 1, byrow = TRUE
)

### Manuscript plot 
### RR Scale

plot_grid(plot_power_overall(full_res_summary,
                             ylim = c(0,NA),
                             xlab = "relative risk"),
          plot_doc_barplot(full_res_summary, 
                           metric = "p_stop_early",
                           ylab = "probability of stopping early",
                           xlab = "relative risk"),
          # plot_doc_barplot(full_res_summary,
          #                  metric = "bias_post_mean",
          #                  ylab = "posterior mean bias",
          #                  xlab = "relative risk"),
          plot_doc_barplot(full_res_summary,
                           metric = "bias_post_median",
                           ylab = "posterior median bias",
                           xlab = "relative risk"),
          plot_trt_estimates_boxplot(full_sim_res, metric = "rmse", 
                                     ylab = "root mean squared error", 
                                     ylim = c(0,1.2),
                                     xlab = "relative risk"),
          labels = LETTERS[1:4],
          ncol = 2, nrow = 2, byrow = TRUE
)


### Manuscript plot 
### Log RR Scale

### Need to get rmse on lrr scale

#sqrt(mean((results - effect_treatment)**2))

library(foreach)
library(doParallel)
registerDoParallel(cores = 4)
rmse_lrr <- foreach(i=1:nrow(full_sim_res), .combine = 'rbind',
                    .errorhandling = 'remove') %dopar% {
                      full_sim_res %>%
                        filter(row_number() == i) %>%
                        select(effect_treatment, trt_posterior) %>%
                        unnest(trt_posterior) %>%
                        rename(results = treatment) %>%
                        summarize(rmse_lrr = sqrt(mean((log(results) - log(effect_treatment))**2)))
                    }
stopImplicitCluster()

full_sim_res2 <- full_sim_res %>%
  bind_cols(rmse_lrr)

names(full_sim_res)

plot_grid(plot_power_overall(full_res_summary,
                             ylim = c(0,NA),
                             xscale = "link",
                             axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                             xlab = "log relative risk"),
          plot_doc_barplot(full_res_summary, 
                           metric = "p_stop_early",
                           ylab = "probability of stopping early",
                           xscale = "link",
                           xlab = "log relative risk"),
          # plot_doc_barplot(full_res_summary,
          #                  metric = "bias_post_mean_lrr",
          #                  ylab = "posterior mean bias",
          #                  xlab = "treatment effect (log relative risk)"),
          plot_doc_barplot(full_res_summary,
                           metric = "bias_post_median_lrr",
                           xscale = "link",
                           ylab = "posterior median bias",
                           xlab = "log relative risk"),
          plot_trt_estimates_boxplot(full_sim_res2, metric = "rmse_lrr", 
                                     ylab = "root mean squared error", 
                                     xscale = "link",
                                     ylim = c(0, NA),
                                     xlab = "log relative risk"),
          labels = LETTERS[1:4],
          ncol = 2, nrow = 2, byrow = TRUE
)
names(full_sim_res)





full_res_summary %>%
  filter(model == "unadjusted") %>%
  select(effect_treatment, bias_post_median)
### Need plot showing spuriour correlation being induced between outcome 
### and noise variables X6-X8 due to the small sample size.  Since these are now 
### correlated, adjustment for them will push the conditonal estimate further from
### the null due to the non-collapsibility of the OR in the logistic regression model


get_cor_with_y <- function(data, i){
  cor_mat <- data %>%
    filter(row_number() == i) %>%
    pull(data) %>%
    flatten_df() %>%
    select(y, x6, x7, x8) %>%
    cov() %>%
    cov2cor() 
  
  ### correlations of these variables with the outcome
  bind_rows(cor_mat[1,2:4]) %>%
    rename(cor_y_x6 = x6,
           cor_y_x7 = x7,
           cor_y_x8 = x8)
}

get_cor_with_y(data = full_sim_res, i = 34)

### Now apply this for each of 500 datsets under different parameter settings

all_datasets <- full_sim_res %>%
  filter(model == "model_correct") 
all_corrs <- foreach(i=1:nrow(all_datasets), .combine = 'rbind') %do% get_cor_with_y(data = all_datasets, i = i)

dataset_with_corrs <- all_datasets %>%
  bind_cols(all_corrs) %>%
  group_by(effect_treatment) %>%
  pivot_longer(cols = c(cor_y_x6:cor_y_x8),
               names_to = "correlation_variable",
               values_to = "correlations_with_y") %>%
  select(correlation_variable, correlations_with_y)

 
ggplot(data = dataset_with_corrs, aes(x = correlations_with_y)) +
  geom_histogram() +
  facet_grid(effect_treatment ~ correlation_variable)



### Distribution of interim analysis stops, sample size

sim_res_ia_stop_distn <- readRDS("ia_cov_adjustment/binary/2022-03-11/Data/sim_res_ia_stop_distn.RDS") %>%
  mutate(model = factor(model, levels = c("model_correct", "model_correct_t_prior", "model_incorrect", "model_noise",
                                          "model_no_strong_prog", "model_no_strong_prog_noise", "model_unadjusted",
                                          "model_unadjusted_t_prior"),
                        labels = c("COR", "CORT", "INC", "NOISE", "NSP", "NSPN", "UNADJ", "UNADJT"))) %>%
  filter(!(model %in% c("CORT", "UNADJT")))

  
        


plot_ia_stop_distn <- function(.data_list){
  
  #factor(n_total, levels = sort(unique(n_total), decreasing = TRUE))
  ggplot(.data_list, aes(x = model, 
                         fill = factor(n_interim_analyses, levels = sort(unique(n_interim_analyses), decreasing = TRUE)))) +
    facet_grid(.~round(effect_treatment,3)) +
    geom_bar(aes(y = n), position = "fill", stat = "identity") +
    scale_fill_brewer(palette = 'Set3') +
    theme_bw() +
    #"total ss"
    labs(fill = "# analyses") +
    ylab("proportion") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=10, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text.y = element_text(size = 8, face = "bold"),
          strip.text.x = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"))
}

plot_ia_stop_distn(sim_res_ia_stop_distn)

### Tables

fn_sum_metric <- function(metric, data = full_res_summary, type = "response"){
  metric <- sym(metric)
  
  if(type == "response"){
    data <- data %>%
      mutate(effect_treatment = exp(effect_treatment))
  }
  
  data %>%
    select(model, effect_treatment, !!metric) %>%
    mutate(across(where(is.double), round, 3)) %>%
    pivot_wider(names_from = effect_treatment,
                values_from = !!metric)
}


table_power <- fn_sum_metric("power", type = "link") %>%
  arrange(model)
saveRDS(table_power, "ia_cov_adjustment/binary/2022-03-11/Data/table_power.RDS")

### Bias
fn_summary_5 <- function(metric){
  tibble(min = min(metric, na.rm = TRUE),
         q1 = quantile(metric, probs = 0.25, na.rm = TRUE),
         mean = mean(metric, na.rm = TRUE),
         median = quantile(metric, probs = 0.5, na.rm = TRUE),
         q3 = quantile(metric, probs = 0.75, na.rm = TRUE),
         max = max(metric, na.rm = TRUE))
}

fn_sum_metric_boxplot <- function(metric, data = full_sim_res, type = "response"){
  metric <- sym(metric)
  
  if(type == "response"){
    data <- data %>%
      mutate(effect_treatment = exp(effect_treatment))
  }
  
  data %>%
    mutate(across(where(is.double), round, 3)) %>%
    group_by(model, effect_treatment) %>%
    select(!!metric) %>%
    summarise(fn_summary_5(!!metric), .groups = "drop") %>%
    arrange(effect_treatment) %>%
    mutate(across(min:max, ~ if_else(is.infinite(.x), NA_real_, .x)),
           across(min:max, ~ if_else(is.nan(.x), NA_real_, .x)))
  
}

### Bias
table_bias <- fn_sum_metric_boxplot("bias", type = "link") %>%
  arrange(effect_treatment, model)
saveRDS(table_bias, "ia_cov_adjustment/binary/2022-03-11/Data/table_bias.RDS")

### RMSE
table_rmse <- fn_sum_metric_boxplot("rmse", type = "link") %>%
  arrange(effect_treatment, model)
saveRDS(table_rmse, "ia_cov_adjustment/binary/2022-03-11/Data/table_rmse.RDS")

### Expected sample size
table_exp_ss <- fn_sum_metric("exp_sample_size", type = "link") %>%
  arrange(model)
saveRDS(table_exp_ss, "ia_cov_adjustment/binary/2022-03-11/Data/table_exp_ss.RDS")

### Probability of stopping early
table_p_stop_early <- fn_sum_metric("p_stop_early", type = "link") %>%
  arrange(model)
saveRDS(table_p_stop_early, "ia_cov_adjustment/binary/2022-03-11/Data/table_p_stop_early.RDS")

### Distn of IA by model
table_sim_res_ia_stop_distn <- sim_res_ia_stop_distn %>%
  select(effect_treatment, model, n_interim_analyses, n) %>%
  pivot_wider(names_from = n_interim_analyses,
              values_from = n) %>%
  arrange(effect_treatment, model)
saveRDS(table_sim_res_ia_stop_distn, "ia_cov_adjustment/binary/2022-03-11/Data/table_sim_res_ia_stop_distn.RDS")


#-----------------------------------------------------------------------------#
#
# Adjustment models comparing t distn priors
#
#-----------------------------------------------------------------------------#

full_res_summary_t <- readRDS("ia_cov_adjustment/binary/2022-03-11/Data/res_summary.RDS") %>%
  mutate(model = factor(model, levels = c("model_correct", "model_correct_t_prior", "model_incorrect", "model_noise",
                                          "model_no_strong_prog", "model_no_strong_prog_noise", "model_unadjusted",
                                          "model_unadjusted_t_prior"))) %>%
  filter(model %in% c("model_correct",
                      "model_correct_t_prior",
                      "model_unadjusted",
                      "model_unadjusted_t_prior"))

plot_power_overall(full_res_summary_t,
                   ylim = c(0,NA),
                   xlab = "treatment effect (relative risk)")

plot_grid(plot_doc_barplot(full_res_summary_t, 
                           metric = "exp_sample_size",
                           ylab = "expected sample size",
                           xlab = "treatment effect (relative risk)"),
          plot_doc_barplot(full_res_summary_t, 
                           metric = "p_stop_early",
                           ylab = "probability of stopping early",
                           xlab = "treatment effect (relative risk)"),
          ncol = 2, nrow = 1, byrow = TRUE)


full_sim_res_t <- readRDS("ia_cov_adjustment/binary/2022-03-11/Data/sim_res.RDS") %>%
  mutate(model = factor(model, levels = c("model_correct", "model_correct_t_prior", "model_incorrect", "model_noise",
                                          "model_no_strong_prog", "model_no_strong_prog_noise", "model_unadjusted",
                                          "model_unadjusted_t_prior")))  %>%
  filter(model %in% c("model_correct",
                      "model_correct_t_prior",
                      "model_unadjusted",
                      "model_unadjusted_t_prior"))

plot_grid(plot_trt_estimates_boxplot(full_sim_res_t, metric = "bias", xlab = "treatment effect (relative risk)"),
          plot_trt_estimates_boxplot(full_sim_res_t, metric = "rmse", 
                                     ylab = "root mean squared error", 
                                     ylim = c(0,NA),
                                     xlab = "treatment effect (relative risk)"),
          ncol = 2, nrow = 1, byrow = TRUE
)


sim_res_ia_stop_distn <- readRDS("ia_cov_adjustment/binary/2022-03-11/Data/sim_res_ia_stop_distn.RDS") %>%
  mutate(model = factor(model, levels = c("model_correct", "model_correct_t_prior", "model_incorrect", "model_noise",
                                          "model_no_strong_prog", "model_no_strong_prog_noise", "model_unadjusted",
                                          "model_unadjusted_t_prior"),
                        labels = c("COR", "CORT", "INC", "NOISE", "NSP", "NSPN", "UNADJ", "UNADJT"))) %>%
  filter(model %in% c("COR", "CORT", "UNADJ", "UNADJT"))


plot_ia_stop_distn(sim_res_ia_stop_distn)
