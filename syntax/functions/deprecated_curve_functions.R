plot_age_curve <- function(covar, outcome, b_size = 12, min_at_risk = 10, w5_sample = FALSE, weights = FALSE, trim_weights = NULL, intervals = TRUE, font_fam = "sans", ylims = NULL){
  require(survival)
  require(icenReg)
  require(patchwork)
  outcome_text <- switch(outcome,
                         concealed_carry = "Concealed carry", 
                         used_gun = "Used gun",
                         shot_shotat = "Shot or shot at someone",
                         gun_use = "Any gun use")
  pals <- unname(list("race" = rev(khroma::color("high contrast")(3)),
                      "sex" = khroma::color("vibrant")(2),
                      "cohort" = rev(khroma::color("bright")(4)))[[covar]])
  input_form <- formula(paste0("Surv(age_left, age_right, type = 'interval2') ~ ", covar))
  
  plot_data <- survival_data[[outcome]] |>
    mutate(sex    = fct_relevel(fct_recode(sex, Male = "m", Female = "f"), "Female", "Male"),
           race   = fct_relevel(race, "Black", "Hispanic", "White"),
           cohort = fct_recode(cohort, 
                               `1996` = "0", 
                               `1981` = "15", 
                               `1984` = "12", 
                               `1987` = "9")) 
  if(w5_sample){
    plot_data <- plot_data |> filter(w5_sample)
  }
  survfit_model <- survfit(input_form, 
                           data = plot_data)
  plot_weight <- NULL
  if(weights){
    plot_weight <- plot_data$weight
  }
  if(!is.null(trim_weights)){
    plot_weight <- scale(case_when(
      plot_weight >= quantile(plot_weight, trim_weights) ~ quantile(plot_weight, trim_weights),
      plot_weight <= quantile(plot_weight, 1-trim_weights) ~ quantile(plot_weight, 1-trim_weights), 
      TRUE ~ plot_weight), center = TRUE, scale = FALSE) + 1
  }
  risk_index <-  survfit_model |>
    summary(times = seq(22, 45, by = 0.01)) |>
    with(tibble(
      max_time = time,
      at_risk = n.risk,
      strata = str_remove(as.character(strata), "^.*="))) |>
    group_by(strata) |>
    filter(at_risk >= min_at_risk) |>
    slice_max(max_time) |>
    select(-at_risk)
  
  plot_curves <- ic_np(formula = formula(input_form), data = plot_data, weights = plot_weight) |>
    extract_curve_df() |>
    filter(is.finite(time)) |>
    rename(age = time) |>
    mutate(strata = str_remove(strata, ".*=")) |>
    left_join(risk_index, by = "strata") |>
    filter(age <= max_time) |>
    select(-max_time) |>
    bind_rows(risk_index |>
                transmute(age = max_time, cdf = NA, strata = strata)) |>
    group_by(strata) |>
    arrange(age) |>
    mutate(cdf = ifelse(is.na(cdf), lag(cdf), cdf)) |>
    ungroup()
  
  int_data <- plot_data |>
    filter(event == 1) |>
    rename(strata = {{covar}}) |>
    mutate(id = row_number(),
           y = -0.005 + as.numeric(as.factor(strata)) * -0.02,
           y = y + runif(n(), -0.008,0.008))
  text_size <- 4.5
  title_text <- str_c(outcome_text, " by ", covar)
  if(covar == "sex"){
    plot_curves <- plot_curves |> mutate(strata = fct_relevel(strata, "Male"))
  }
  curve_plot <- plot_curves |>
    mutate(cdf = cdf*100) |>
    ggplot(aes(x = age, y = cdf, color = strata), linewidth = 0.25) +
    geom_line(linewidth = 1) + 
    labs(x=NULL, y = NULL, title = title_text, color = NULL) +
    scale_x_continuous(limits = c(0, 42.5)) +
    scale_y_continuous(labels = ~ paste0(round(., 2), "%"), limits = ylims) +
    scale_color_manual(values = pals) +
    theme_classic(base_size = b_size) +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.key.height = unit(7, "pt"),
          legend.margin=margin(0,4,4,4),
          panel.grid.major.y = element_line(linetype = "solid", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          text = element_text(family = font_fam),
          axis.text = element_text(color = "black"),
          axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5, "pt"),
          title = element_text(size = rel(0.8)),
          legend.box.background = element_rect(colour = "black"))
  if(intervals){
    int_plot <- ggplot(int_data) +
      geom_errorbarh(aes(xmin = age_left, xmax = age_right, y = y, color = strata), height = 0.005, alpha = 0.35,  inherit.aes=FALSE) +
      labs(x=NULL, y = NULL, title = NULL, color = NULL) +
      scale_x_continuous(limits = c(0, 42.5)) +
      scale_color_manual(values = pals) +
      theme_void() +
      theme(legend.position = "none")
    plot_out <- curve_plot / int_plot &
      plot_layout(heights = c(5, 1)) &
      theme(plot.background = element_rect(fill = "white", color = "white"),
            plot.margin = margin(0, 5.5, 5.5, 5.5, "pt"))
  } else {
    plot_out <- curve_plot + 
      theme(plot.background = element_rect(fill = "white", color = "white"),
            plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"))
  }
  return(plot_out)
}

plot_year_curve <- function(covar, outcome, b_size = 12, min_at_risk = 10, ret = FALSE, w5_sample = FALSE, weights = FALSE, trim_weights = NULL, intervals = TRUE){
  outcome_text <- switch(outcome,
                         concealed_carry = "Concealed carry", 
                         used_gun = "Used gun",
                         shot_shotat = "Shot or shot at someone",
                         gun_use = "Any gun use")
  pals <- unname(list("race" = rev(khroma::color("high contrast")(3)),
                      "sex" = khroma::color("vibrant")(2),
                      "cohort" = rev(khroma::color("bright")(4)))[[covar]])
  input_form <- formula(paste0("Surv(age_left, age_right, type = 'interval2') ~ ", covar))
  
  plot_data <- survival_data[[outcome]] |>
    mutate(sex    = fct_relevel(fct_recode(sex, Male = "m", Female = "f"), "Male", "Female"),
           race   = fct_relevel(race, "Black", "Hispanic", "White"),
           cohort = fct_recode(cohort, 
                               `1996` = "0", 
                               `1981` = "15", 
                               `1984` = "12", 
                               `1987` = "9"))  |>
    left_join(survey_index |> distinct(subid, dob), by = "subid") |>
    mutate(year_left = decimal_date(dob + round(365.25*age_left)),
           year_right = ifelse(!is.finite(age_right), Inf, decimal_date(dob + round(365.25*age_right))))|>
    mutate(age_left = year_left, age_right = year_right)
  
  if(w5_sample){
    plot_data <- plot_data |> filter(w5_sample == 1)
  }
  plot_weight <- NULL
  if(weights){
    plot_weight <- plot_data$weight
  }
  if(!is.null(trim_weights)){
    plot_weight <- scale(case_when(
      plot_weight >= quantile(plot_weight, trim_weights) ~ quantile(plot_weight, trim_weights),
      plot_weight <= quantile(plot_weight, 1-trim_weights) ~ quantile(plot_weight, 1-trim_weights), 
      TRUE ~ plot_weight), center = TRUE, scale = FALSE) + 1
  }
  risk_index <-  survfit(input_form, 
                         data = plot_data) |>
    summary(times = seq(2020, 2023, by = 0.001)) |>
    with(tibble(
      max_time = time,
      at_risk = n.risk,
      strata = str_remove(as.character(strata), "^.*="))) |>
    group_by(strata) |>
    filter(at_risk >= min_at_risk) |>
    slice_max(max_time) |>
    select(-at_risk)
  
  plot_curves <- ic_np(formula = formula(input_form), data = plot_data, weights = plot_weight) |>
    extract_curve_df() |>
    rename(age = time) |>
    left_join(risk_index, by = "strata") |>
    filter(age <= max_time) |>
    select(-max_time) |>
    bind_rows(risk_index |>
                transmute(age = max_time, cdf = NA, strata = strata)) |>
    group_by(strata) |>
    arrange(age) |>
    mutate(cdf = ifelse(is.na(cdf), lag(cdf), cdf)) |>
    ungroup() |>
    mutate(cdf = cdf*100)
  
  int_data <- plot_data |>
    filter(event == 1) |>
    rename(strata = {{covar}}) |>
    mutate(id = row_number(),
           y = -0.005 + as.numeric(as.factor(strata)) * -0.02,
           y = y + runif(n(), -0.008,0.008),
           age_left = ifelse(age_left < 1992.75, 1992.75, age_left)) |>
    filter(age_right > 1992.75)
  text_size <- 4.5
  title_text <- str_c(outcome_text, " by ", covar)
  if(covar == "sex"){
    plot_curves <- plot_curves |> mutate(strata = fct_relevel(strata, "Male"))
  }
  if(covar == "cohort" & outcome == "concealed_carry"){
    age_21 <- plot_curves |>
      arrange(age) |>
      filter( (age-as.numeric(strata)) >= 21) |>
      group_by(strata) |>
      slice(1) |>
      mutate(age = as.numeric(strata) + 21)
  }
  
  
  curve_plot <- plot_curves |>
    ggplot(aes(x = age, y = cdf, color = strata), linewidth = 0.25) +
    geom_line(linewidth = 0.75) + 
    labs(x=NULL, y = NULL, title = title_text, color = NULL) +
    scale_x_continuous(breaks = seq(1995, 2020, by = 5)) +
    scale_y_continuous(labels = ~ paste0(round(., 2), "%")) +
    scale_color_manual(values = pals) +
    coord_cartesian(xlim = c(1994, 2022), expand = TRUE ) +
    theme_classic(base_size = b_size) +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.key.height = unit(7, "pt"),
          legend.margin=margin(0,4,4,4),
          panel.grid.major.y = element_line(linetype = "solid", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          text = element_text(family = font_fam),
          axis.text = element_text(color = "black"),
          axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5, "pt"),
          title = element_text(size = rel(0.8)),
          legend.box.background = element_rect(colour = "black"))
  if(outcome == "concealed_carry"){
    curve_plot <- curve_plot +
      geom_vline(xintercept = 2014, linetype = "dashed", alpha = 0.5) +
      annotate("label", x = 2014, y = max(plot_curves$cdf), label = "Legalized", size = 4, hjust = 0.5, fill = "white", label.size = 0)
  }
  if(covar == "cohort" & outcome == "concealed_carry"){
    curve_plot <- curve_plot +
      annotate("text", x = 2007, y = 14.5, label = "Age 21", size = 4) +
      geom_segment(data = age_21 |> filter(strata == 1984), aes(x = 2006, xend = age, y = 13, yend = cdf), inherit.aes = FALSE) +
      geom_point(data = age_21)
  }
  if(intervals){
    int_plot <- ggplot(int_data) +
      geom_errorbarh(aes(xmin = age_left, xmax = age_right, y = y, color = strata), height = 0.005, alpha = 0.35,  inherit.aes=FALSE) +
      labs(x=NULL, y = NULL, title = NULL, color = NULL) +
      scale_color_manual(values = pals) +
      coord_cartesian(xlim = c(1994, 2022), expand = TRUE, clip = "off") +
      theme_void() +
      theme(legend.position = "none")
    plot_out <- curve_plot / int_plot &
      plot_layout(heights = c(5, 1)) &
      theme(plot.background = element_rect(fill = "white", color = "white"),
            plot.margin = margin(0, 5.5, 5.5, 5.5, "pt"))
  } else {
    plot_out <- curve_plot + 
      theme(plot.background = element_rect(fill = "white", color = "white"),
            plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"))
  }
  if(ret){
    return(plot_out)
  } else {
    print(plot_out)
  }
}