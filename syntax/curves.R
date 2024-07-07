library(tidyverse)
library(survival)
library(lubridate)
library(interval)
library(khroma)
library(showtext)
library(pdftools)
library(ggtext)
library(patchwork)
library(icenReg)
library(survey)

source("./syntax/project_functions.R")
source("./syntax/functions/survival_curve_functions.R")
font_fam <- "Open Sans"
font_add_google(name = font_fam)
showtext_auto()

load("./data/survival_data.RData")
load("./data/survey_index.RData") # Needed to get dobs

# EXPOSURE
plot_vars <- expand_grid(x = c("race", "sex"),
                         y = c("seen_shot")) 

seen_shot_plot <- render_curve_plot("age-seen-shot-curves-weighted", plot_vars, 
                              xaxis = "age",
                              title = "strata_xaxis",
                              ret = TRUE,
                              year_mark = TRUE,
                              weights = TRUE, trim_level = 0.95,
                              fig_width  = 3, horizontal = TRUE, 
                              font_fam = font_fam,
                              ylims = c(0, 80),
                              title_letters = TRUE,
                              risk_table = FALSE) 
seen_shot_plot

plot_vars <- expand_grid(x = c("race", "sex"),
                         y = c("been_shot")) 

been_shot_plot <- render_curve_plot("age-been-shot-curves-weighted", plot_vars, 
                               xaxis = "age",
                               title = "strata_xaxis",
                               year_mark = TRUE,
                               ret = TRUE,
                               weights = TRUE, trim_level = 0.95,
                               fig_width  = 3, horizontal = TRUE, 
                               font_fam = font_fam,
                               ylims = c(0, 20),
                               title_letters = TRUE,
                               risk_table = FALSE) 

plot_list <- seen_shot_plot  / plot_spacer() / been_shot_plot + plot_layout(nrow = 3, heights = c(20,1,20))
plot_list
filename <- glue::glue("./img/exposure-curves-weighted.pdf")
fig_width <- 9
fig_height <- 6
fig_scale <- 0.8
ggsave(plot = plot_list, filename = filename, 
       device = pdf, 
       width = unit(fig_width/fig_scale, units = "in"), height = unit(fig_height/fig_scale, units = "in"))
plot_pdf <- pdftools::pdf_render_page(filename, dpi = 320)
png::writePNG(plot_pdf, str_replace(filename, "\\.pdf", ".png" ), dpi = 320)
message(glue::glue("Plot saved to {filename}"))

# CARRY
plot_vars <- expand_grid(x = c("race", "sex", "cohort"),
                         y = c("concealed_carry")) 

age_plot <- render_curve_plot("age-carry-curves-weighted", plot_vars, 
                              xaxis = "age",
                              title = "strata_xaxis",
                              ret = TRUE,
                              year_mark = TRUE,
                              weights = TRUE, trim_level = 0.95,
                              fig_width  = 3, horizontal = TRUE, 
                              font_fam = font_fam,
                              ylims = c(0, 50),
                              title_letters = TRUE,
                              risk_table = FALSE) 

year_plot <- render_curve_plot("year-carry-curves-weighted", plot_vars, 
                               xaxis = "year",
                               title = "strata_xaxis",
                               year_mark = TRUE,
                               ret = TRUE,
                               weights = TRUE, trim_level = 0.95,
                               fig_width  = 3, horizontal = TRUE, 
                               font_fam = font_fam,
                               ylims = c(0, 50),
                               title_letters = TRUE,
                               risk_table = FALSE) 

plot_list <- age_plot  / plot_spacer() / year_plot + plot_layout(nrow = 3, heights = c(20,1,20))
plot_list
filename <- glue::glue("./img/carry-curves-weighted.pdf")
fig_width <- 9
fig_height <- 6
fig_scale <- 0.8
ggsave(plot = plot_list, filename = filename, 
       device = pdf, 
       width = unit(fig_width/fig_scale, units = "in"), height = unit(fig_height/fig_scale, units = "in"))
plot_pdf <- pdftools::pdf_render_page(filename, dpi = 320)
png::writePNG(plot_pdf, str_replace(filename, "\\.pdf", ".png" ), dpi = 320)
message(glue::glue("Plot saved to {filename}"))