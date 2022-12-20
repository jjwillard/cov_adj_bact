pacman::p_load(dplyr, foreach, ggplot2, cowplot, tidyr, purrr, rlang, kableExtra)

### Read in data
full_res <- readRDS("PATH for sim_res N = 100 from simulations_tte.R") %>%
  bind_rows(readRDS("PATH for sim_res N = 200 from simulations_tte.R")) %>%
  bind_rows(readRDS("PATH for sim_res N = 500 from simulations_tte.R")) %>%
  bind_rows(readRDS("PATH for sim_res N = 1000 from simulations_tte.R")) %>%
  mutate(model = factor(model, levels = c("model_correct", "model_incorrect", "model_noise",
                                          "model_correct_prior", "model_correct_prior_strong", "model_unadjusted"),
                        labels = c("correct", "no quad", "correct noise", "correct prior",
                        "correct strong prior", "unadjusted")))

## effect_treatment on log-hazard scale
summarize_results <- function(complete_sim_results){
  complete_sim_results %>%
    group_by(model, max_ss, effect_treatment, beta_1, beta_2, beta_3, beta_4, beta_5) %>%
    summarise(power = mean(superiority, na.rm = TRUE),
              exp_sample_size = mean(n_total, na.rm = TRUE),
              p_stop_early = 1 - mean(reach_max_time, na.rm = TRUE),
              bias_post_median =  mean(trt_est_median_hazard - exp(effect_treatment), na.rm = TRUE), 
              bias_post_mean =  mean(trt_est_mean_hazard - exp(effect_treatment), na.rm = TRUE), 
              bias_post_median_lh =  mean(trt_est_median - effect_treatment, na.rm = TRUE), 
              bias_post_mean_lh =  mean(trt_est_mean - effect_treatment, na.rm = TRUE), 
              .groups = "keep") %>%
    ungroup()
}

full_res_summary <- summarize_results(full_res)

### Power

plot_summary_metric <- function(.data_list, 
                                metric = metric,
                                ylim = c(NA, NA), 
                                round_digits = 2,
                                xlab = 'hazard ratio',
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
  
  if(xscale == "response"){
    .data_list <- .data_list %>%
      mutate(effect_treatment = exp(effect_treatment))
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
                         #linetype = model,
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

### Manuscript plot with just power and pse ###-------------------------------------------

plot_legend <- get_legend(plot_summary_metric(full_res_summary %>%
                                                filter(effect_treatment != 0),
                                              metric = "power",
                                              xlab = 'hazard ratio',
                                              ylab = 'power',
                                              scales = "free",
                                              dodge_width = 0.5,
                                              legend_position = "top",
                                              n_rows_legend = 1))
  


plot(plot_legend)
prow <- plot_grid(plot_summary_metric(full_res_summary %>%
                                        filter(effect_treatment != 0),
                                      metric = "power",
                                      xlab = 'hazard ratio',
                                      ylab = 'power',
                                      scales = "free",
                                      dodge_width = 0.5,
                                      round_digits = 2,
                                      legend_position = "none"),
                  plot_summary_metric(full_res_summary %>%
                                        filter(effect_treatment != 0),
                                      metric = "p_stop_early",
                                      xlab = 'hazard ratio',
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
prow

power_pse <- plot_grid(prow, plot_legend, ncol = 1, rel_heights = c(1, .1))

ggsave2(filename = "PATH/FILENAME.eps",
        plot = power_pse,
        width = 8,
        height = 3,
        units = "in",
        dpi = 800)


####-------------------------------------------------------------------####
### Plot bias and RMSE

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
                                       xscale = "response",
                                       scales = "free",
                                       round_digits = 2,
                                       ylim = c(NA, NA),
                                       legend_position = "top",
                                       n_rows_legend = 2){
  
  n_models <- length(unique(.data_list$model))
  
  if(xscale == "response"){
    .data_list <- .data_list %>%
      mutate(effect_treatment = exp(effect_treatment))
  }
  
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
                                                filter(effect_treatment != 0),
                                              metric = "bias_post_median",
                                              xlab = 'hazard ratio',
                                              ylab = 'posterior median bias',
                                              scales = "free",
                                              dodge_width = 0.5,
                                              round_digits = 2,
                                              legend_position = "top",
                                              n_rows_legend = 1))
                          

plot(plot_legend)
prow <- plot_grid(plot_summary_metric(full_res_summary %>%
                                        filter(effect_treatment != 0),
                                      metric = "bias_post_median",
                                      xlab = 'hazard ratio',
                                      ylab = 'posterior median bias',
                                      scales = "free",
                                      dodge_width = 0.5,
                                      round_digits = 2,
                                      legend_position = "none"),
                  plot_trt_estimates_boxplot(full_res,
                                             metric = "rmse", 
                                             ylab = "root mean squared error",
                                             xlab = "hazard ratio",
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

#----------------------------------------------------------------------------------#

### False positive rate table

## change scale of effect treatment to HAZARD
full_res_summary <- full_res_summary %>%
  mutate(effect_treatment = exp(effect_treatment))

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
  left_join(ess_100 %>% select(-max_ss, null = `1`, one = `0.7317716`, two = `0.6517841`) %>% 
              relocate(model, null, one, two), by = "model")

n_200 <- fpr %>% select(model, frp = `200`) %>%
  left_join(bias_null %>% select(model, bias_null = `200`), by = 'model') %>%
  left_join(ess_200 %>% select(-max_ss, null = `1`, one = `0.7690409`, two = `0.6536824`) %>% 
              relocate(model, null, one, two), by = "model")

n_500 <- fpr %>% select(model, frp = `500`) %>%
  left_join(bias_null %>% select(model, bias_null = `500`), by = 'model') %>%
  left_join(ess_500 %>% select(-max_ss, null = `1`, one = `0.8437124`, two = `0.7697936`) %>% 
              relocate(model, null, one, two), by = "model")

n_1000 <- fpr %>% select(model, frp = `1000`) %>%
  left_join(bias_null %>% select(model, bias_null = `1000`), by = 'model') %>%
  left_join(ess_1000 %>% select(-max_ss, null = `1`, one = `0.890788`, two = `0.8389165`)  %>% 
              relocate(model, null, one, two), by = "model")

tbl_summary <- bind_rows(n_100  %>%
                           left_join(n_200, by = "model"),
                         n_500 %>%
                           left_join(n_1000, by = "model"))

kbl(tbl_summary, booktabs = T, format = "latex") %>%
  add_header_above(c(" "=3, "Expected sample size"=3,
                     " "=2, "Expected sample size"=3)) %>%
  add_header_above(c(" "=1, "Maximum sample size = 100" = 5,
                     "Maximum sample size = 200" = 5)) %>%
  kable_styling(latex_options = c("scale_down"))



