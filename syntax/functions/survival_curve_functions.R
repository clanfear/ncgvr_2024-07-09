# These functions generate multi-panel survival curve plots by
# age or year using the Turnbull NPMLE

plot_curve <- 
  function(covar, 
           outcome, 
           xaxis = "age", 
           title = "strata_xaxis", 
           b_size = 12, 
           min_at_risk = 10, 
           w5_sample = FALSE, 
           weights = FALSE, 
           trim_level = NULL, 
           font_fam = "sans", 
           ylims = NULL, 
           year_mark = FALSE,
           title_letters = FALSE,
           risk_table = TRUE){
  require(survival)
  require(icenReg)
  require(patchwork)
    title_text <- switch(title,
                       strata = stringr::str_to_title(covar),
                       strata_xaxis = str_c(stringr::str_to_title(covar), " by ", xaxis),
                       none = NULL)
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
  if(xaxis == "year"){
    plot_data <- plot_data |>
      left_join(survey_index |> distinct(subid, dob), by = "subid") |>
      mutate(year_left = decimal_date(dob + round(365.25*age_left)),
             year_right = ifelse(!is.finite(age_right), Inf, decimal_date(dob + round(365.25*age_right))))|>
      mutate(age_left = year_left, age_right = year_right)
  }
  if(w5_sample){
    plot_data <- plot_data |> filter(w5_sample)
  }
  survfit_model <- survfit(input_form, 
                           data = plot_data)
  plot_weight <- NULL
  if(weights){
    plot_weight <- plot_data$weight
  }
  if(!is.null(trim_level)){
    plot_weight <- trim_weights(plot_weight, 1-trim_level, trim_level)
  }
  risk_times <- switch(xaxis,
                       age = seq(22, 45, by = 0.005),
                       year = seq(2020, 2023, by = 0.001)
  )
  risk_index <-  survfit_model |>
    summary(times = risk_times) |>
    with(tibble(
      max_time = time,
      at_risk = n.risk,
      strata = str_remove(as.character(strata), "^.*="))) |>
    group_by(strata) |>
    filter(at_risk >= min_at_risk) |>
    slice_max(max_time) |>
    select(-at_risk)
  xbreaks <- switch(xaxis,
                    age = seq(0, 40, by = 10),
                    year = seq(1996, 2020, by = 6))
  xlims <- switch(xaxis,
                  age = c(0, 42.5),
                  year = c(1996, 2022))
  if(risk_table){
    covar_levels <- levels(plot_data[[covar]])
    at_risk_label_x <- switch(xaxis, age = -12.5, year = 1988.5)
    at_risk_title_x <- switch(xaxis, age = -15, year = 1987)
    surv_counts <- survfit_model |>
      summary(xbreaks) %>%
      with(.,
                        tibble(
                          time = time,
                          at_risk = n.risk,
                          covar = strata)) |>
      mutate(covar = str_remove(covar, "^.*="),
             at_risk = as.character(round(at_risk, 0))) |>
      pivot_wider(names_from = time, values_from = at_risk) |>
      mutate(across(-covar, ~ifelse(is.na(.), 0, .)))
      
    text_df <- surv_counts |>
      mutate(covar = factor(covar, levels = covar_levels),
             y = -as.numeric(covar)) |>
      pivot_longer(-y, names_to = "x", values_transform = as.character) |>
      mutate(x = as.integer(ifelse(x != "covar", x, at_risk_label_x))) |>
      bind_rows(data.frame(y = 0, x = at_risk_title_x, value = "No. at risk")) 
    
    risk_plot <- ggplot(text_df, aes(x = x, y = y, label = value)) + 
      coord_cartesian(ylim = c(-4, 1), xlim = xlims, clip = "off") +
      theme_void() + 
      geom_text(data = ~. |> filter(x >= xlims[1]), size = 2.75, hjust = 0.5) +
      geom_text(data = ~. |> filter(x < xlims[1]), size = 2.75, hjust = 0) +
      theme(plot.margin = margin(-5.5, 5.5, 0, 5.5*4, "pt"))
  }
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
    ungroup() |>
    mutate(cdf = cdf*100)
  
  text_size <- 4.5
  
  if(covar == "sex"){
    plot_curves <- plot_curves |> mutate(strata = fct_relevel(strata, "Male"))
  }
  if(xaxis == "year" & covar == "cohort"){
    age_21 <- plot_curves |>
      arrange(age) |>
      filter( (age-as.numeric(strata)) >= 21) |>
      group_by(strata) |>
      slice(1) |>
      mutate(age = as.numeric(strata) + 21)
  }
  
  curve_plot <- plot_curves |>
    ggplot(aes(x = age, y = cdf, color = strata)) +
    coord_cartesian(clip = "on", xlim = xlims, ylim = ylims)
  if(xaxis == "year" & year_mark){
    curve_plot <- curve_plot +
      annotate("text", x = c(2016, 2020), y = c(0.75*ylims[2], 0.82*ylims[2]), label = c(2016, 2020), size = 4, vjust = 0, color = "white") +
      annotate("segment", x = c(2016, 2020), xend = c(2016, 2020), y = c(0.75*ylims[2], 0.82*ylims[2]) - (ylims[2]/100), yend = c(-Inf,-Inf), linetype = "dashed", alpha = 0.75, linewidth = 0.5, color = "white")
  }
  if(xaxis == "age" & year_mark){
    curve_plot <- curve_plot +
      annotate("text", x = 21, y = 1.03*(0.7*ylims[2]), label = "Age 21", size = 4, vjust = 0) +
      annotate("segment", x = 21, xend = 21, y = 0.7*ylims[2], yend = c(-Inf,-Inf), linetype = "dashed", alpha = 0.75, linewidth = 0.5)
  }
  curve_plot <- curve_plot +
    geom_line(linewidth = 0.75) +
    labs(x=NULL, y = NULL, 
         title = NULL, 
         color = NULL)
  if(!title_letters){
    curve_plot <- curve_plot +
      labs(title = title_text)
  }
  curve_plot <- curve_plot +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(labels = ~ paste0(round(., 2), "%")) +
    scale_color_manual(values = pals) +
    theme_classic() +
    theme(legend.position = "inside",
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = "white"),
          legend.position.inside = c(0, 1),
          legend.justification = c(0, 1),
          legend.key.height = unit(7, "pt"),
          legend.margin=margin(4,4,4,4),
          legend.box.background = element_rect(colour = "black"),
          panel.grid.major.y = element_line(linetype = "solid", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          text = element_text(family = font_fam),
          axis.text = element_text(color = "black"),
          axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
          title = element_text(size = rel(0.8)))
  if(xaxis == "year" & covar == "cohort"){
    curve_plot <- curve_plot +
      annotate("text", x = 2000, y = 0.55*ylims[2], label = "Age 21", size = 4, color = "white") +
      geom_segment(data = age_21 |> filter(strata == 1981), aes(x = 2000, xend = age, y = 0.52*ylims[2], yend = cdf), inherit.aes = FALSE, color = "white") +
      geom_point(data = age_21, shape = 21, fill = "white")
  }
  if(title_letters){
    title_letter_mat <- matrix(c("A", "B", "C", "D", "E", "F"), ncol = 3, byrow = TRUE)
    rownames(title_letter_mat) <- c("age", "year")
    colnames(title_letter_mat) <- c("race", "sex", "cohort")
    title_letter <- title_letter_mat[xaxis, covar]
    panel_label <- ggplot() +
      annotate("label",  x = 0 , y = 0.5, label = title_letter, size = 3.5, label.r = unit(0, "pt")) +
      annotate("text",  x =  .75, hjust = 0, y = 0.5, label = title_text, size = 3.5) +
      theme_void() +
      coord_cartesian(expand=FALSE, xlim = c(0,10), ylim = c(0,1), clip = "off") +
      theme(text = element_text(family = "Open Sans"),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = "white"),
            plot.margin = margin(5.5, 5.5,3,0))
    curve_plot <- curve_plot + 
      theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"))
    if(risk_table){
      curve_plot <- curve_plot + 
        labs(y = "Cumulative onset, %", x = str_to_title(xaxis)) +
        theme(plot.margin = margin(5.5, 5.5, 0, 5.5*4, "pt"))
      plot_out <- free(panel_label) +  
        curve_plot + 
        risk_plot + 
        plot_layout(ncol = 1, nrow = 3, widths = c(20, 20, 20), heights = c(1,22,6))
    } else {
      plot_out <- (free(panel_label) / curve_plot ) + 
        plot_layout(widths = c(20,20), heights = c(1,20))
    }
  } else {
    plot_out <- curve_plot 
  }
  return(plot_out)
}

render_curve_plot <- function(filename, 
                              plot_vars, 
                              xaxis = "age",
                              title = "full",
                              year_mark = FALSE,
                              ret = FALSE, 
                              min_at_risk = 10, 
                              weights = FALSE, 
                              w5_sample = FALSE, 
                              trim_level = NULL, 
                              horizontal = FALSE, 
                              fig_scale = 0.75,
                              fig_width  = 6,
                              fig_height = 9,
                              ylims = NULL,
                              font_fam = "sans",
                              title_letters = TRUE,
                              risk_table = TRUE){
  plots_out <- map2(plot_vars$x, plot_vars$y, \(x, y) plot_curve(x, y, xaxis = xaxis, title = title, min_at_risk = min_at_risk, weights = weights, w5_sample = w5_sample, trim_level = trim_level, ylims = ylims, year_mark = year_mark, title_letters = title_letters, risk_table = risk_table))
  if(horizontal){
    plot_list <- wrap_plots(plots_out, nrow = length(unique(plot_vars$y)), ncol = length(unique(plot_vars$x)), byrow = FALSE)
    fw <- fig_width
    fh <- fig_height
    fig_width  <- fh
    fig_height <- fw
  } else {
    plot_list <- wrap_plots(plots_out, ncol = length(unique(plot_vars$y)), nrow = length(unique(plot_vars$x)))
  }
  if(ret){
    return(plot_list)
  } else {
    filename <- glue::glue("./docs/draft/fig/{filename}.pdf")
    ggsave(plot = plot_list, filename = filename, 
           device = pdf, 
           width = unit(fig_width/fig_scale, units = "in"), height = unit(fig_height/fig_scale, units = "in"))
    plot_pdf <- pdftools::pdf_render_page(filename, dpi = 320)
    png::writePNG(plot_pdf, str_replace(filename, "\\.pdf", ".png" ), dpi = 320)
    message(glue::glue("Plot saved to {filename}"))
  }
}
