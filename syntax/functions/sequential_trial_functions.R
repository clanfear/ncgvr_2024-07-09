# These functions generate and analyze data for sequential trial models

# The basic idea of the trial models is to pick a window time to analyze and 
# then define exposure as being any exposure occurring before the window starts
# and the outcome as being the target outcome occurring within the window. This
# is trivial with normal data but weird with interval censored data because the
# censoring intervals of the exposure and outcome could overlap. So, I impute
# exact times then reconstruct the data for each imputation. This means you can
# get different sample sizes every time, because you have to drop people if they
# experienced the outcome prior to the window!

# First, one needs to generate data to estimate the trial models on. This takes:
# data: A single dataframe of combined PHDCN data
# exposure_max_age: Max age at which exposure is considered (window start)
# outcome_max_age: Max at at which the outcome is considered (window end) 
#   If NULL, there is no max age and a Cox model will eventually get run
# outcome: Which outcome ("carry" or "gun_use")
# exposure: Which exposure ("exposed" or "carry")

create_trial_data <- function(data, exposure_max_age, outcome_max_age = NULL, outcome = "carry", exposure = "exposed"){
  # If carry is exposure, just replace the vars
  if(exposure == "carry"){
    data <- data |>
      mutate(exposed = concealed_carry,
             exp_age = cc_age)
  }
  if(outcome == "carry"){
    carry_max_age <- outcome_max_age
    # Only include people who responded to concealed carry at least once
    data <- data |> filter(!is.na(concealed_carry)) 
    # NULL max age means window goes to infinity; right censored data
    if(is.null(carry_max_age)){
      
      trial_data <- data |>
        mutate(
          # If exposed before start of window, you're exposed, otherwise not exposed.
          exposed = case_when(
            exposed == 0 ~ "Not exposed",
            exp_age >= exposure_max_age ~ "Not exposed",
            exp_age < exposure_max_age ~ "Exposed",
            TRUE ~ "ERROR"),
          # If carried before start of window, you're ineligible, and if you started carrying after window opens, you carried, otherwise not
          carry = case_when(
            concealed_carry == 0 ~ "No carry",
            cc_age >= exposure_max_age  ~ "Carry",
            cc_age < exposure_max_age ~ "Prior carry",
            TRUE ~ "ERROR"
          )) |> 
        # Drop anyone who carried prior to start of window or was not observed in window
        filter(carry != "Prior carry" & cc_age >= exposure_max_age) |>
        mutate(carry = as.numeric(fct_relevel(carry, "No carry")),
               exposed = as.numeric(exposed == "Exposed")) |>
        mutate(trial_age = exposure_max_age)
    } else if (!is.null(carry_max_age)) {
      trial_data <- data |>
        mutate(
          exposed = case_when(
            exposed == 0 ~ "Not exposed",
            exp_age >= exposure_max_age ~ "Not exposed",
            exp_age < exposure_max_age ~ "Exposed",
            TRUE ~ "ERROR"),
          # Same deal as above but now there's a max value for the window; people who started carrying after window didn't carry
          carry = case_when(
            concealed_carry == 0 ~ "No carry",
            cc_age >= carry_max_age ~ "No carry",
            cc_age >= exposure_max_age & cc_age < carry_max_age ~ "Carry",
            cc_age < exposure_max_age ~ "Prior carry",
            TRUE ~ "ERROR"
          )) |> 
        # Drop anyone who carried prior to start of window or was not observed in window
        filter(carry != "Prior carry" & cc_age >= exposure_max_age) |>
        mutate(carry = as.numeric(carry == "Carry"),
               exposed = as.numeric(exposed == "Exposed")) |>
        mutate(trial_age = exposure_max_age)
    }
    # If gun use is outcome, basically the same process as exposure
  } else if (outcome == "gun_use"){
    gun_use_max_age <- outcome_max_age
    data <- data |> filter(!is.na(gun_use)) 
    if(is.null(gun_use_max_age)){
      trial_data <- data |>
        mutate(
          exposed = case_when(
            exposed == 0 ~ "Not exposed",
            exp_age >= exposure_max_age ~ "Not exposed",
            exp_age < exposure_max_age ~ "Exposed",
            TRUE ~ "ERROR"),
          gun_use = case_when(
            gun_use == 0 ~ "No use",
            gun_use_age >= exposure_max_age  ~ "Used gun",
            gun_use_age < exposure_max_age ~ "Prior use",
            TRUE ~ "ERROR"
          )) |> 
        filter(gun_use != "Prior use" & gun_use_age >= exposure_max_age) |>
        mutate(gun_use = as.numeric(gun_use == "Used gun"),
               exposed = as.numeric(exposed == "Exposed")) |>
        mutate(trial_age = exposure_max_age)
    } else if (!is.null(gun_use_max_age)) {
      trial_data <- data |>
        mutate(
          exposed = case_when(
            exposed == 0 ~ "Not exposed",
            exp_age >= exposure_max_age ~ "Not exposed",
            exp_age < exposure_max_age ~ "Exposed",
            TRUE ~ "ERROR"),
          gun_use = case_when(
            gun_use == 0 ~ "No use",
            gun_use_age >= gun_use_max_age ~ "No use",
            gun_use_age >= exposure_max_age & gun_use_age < gun_use_max_age ~ "Used gun",
            gun_use_age < exposure_max_age ~ "Prior use",
            TRUE ~ "ERROR"
          )) |> 
        filter(gun_use != "Prior use" & gun_use_age >= exposure_max_age) |>
        mutate(gun_use = as.numeric(gun_use == "Used gun"),
               exposed = as.numeric(exposed == "Exposed")) |>
        mutate(trial_age = exposure_max_age)
    }
  }
  # Recode reference categories, standardize neighb vars, rename exposure var
  trial_data <- trial_data |>
    mutate(fam_crim = ifelse(fam_crim == "none", 0, 1),
           race = fct_relevel(race, "White"),
           cohort = fct_relevel(cohort, "15")) |>
    mutate(across(matches("12|18"), ~standardize(.))) |>
    rename({{exposure}} := exposed) 
  return(trial_data)
}

# Run sequential trial model
# Arguments:
# exposure_max_age: Same as above
# outcome_max_age: Same as above
# outcome: Same as above
# covars: Age at which control covariates are measured (12 or 18)
# weighted: Use survey weights? (TRUE or FALSE)
# weight_trim: Trim weights? (center percentile to include, e.g., 95)
# exposure: : Same as above
# predictions: Return predicted probabilities instead of estimates (TRUE or FALSE)
# pred_var: For predicted probs, strata to generate probs across
# pred_vals: Values of strata var (e.g., c(0,9,12,15) )
# interaction_term: Extra text to bolt on to formula to run interactions or whatever

run_trial_model <- function(exposure_max_age,
                            outcome_max_age = NULL,
                            outcome = "carry",
                            covars = 12,
                            weighted = FALSE,
                            weight_trim = NULL,
                            exposure = "exposed",
                            predictions = FALSE,
                            pred_var = "cohort",
                            pred_vals = c(0,9,12,15),
                            interaction_term = NULL,
                            data = concealed_carry_exposure_imputed,
                            contrast = FALSE,
                            cox_max_age = 40){
  # First, generate trial data. This will spit out a list of trial datasets that you'll need to pool results from
  data_list  <- map(data, ~ create_trial_data(.x, exposure_max_age, outcome_max_age, outcome = outcome, exposure = exposure))
  # If weighted and trimming weights... trim the weights; 0.95 will give you middle 90%
  if(weighted & !is.null(weight_trim)){
    data_list <- data_list |> 
      map(\(x) x |> 
            mutate(weight = trim_weights(weight, 1-weight_trim, weight_trim)))
  }
  # If there's no max outcome age, you're running a survival model.
  if(is.null(outcome_max_age)){
    # Choose outcome
    outcome_age <- c("carry" = "cc_age", "gun_use" = "gun_use_age")[outcome]
    # Generate formula
    trial_form <- glue::glue("Surv({outcome_age}, {outcome}) ~ {exposure} + sex + race + cohort + ses + pc_fborn + fam_crim + percblack_{covars} + poverty_{covars} + pop_density_{covars} + homicide_{covars}")
    # If no exposure, delete exposure from formula
    if(exposure == "none"){
      trial_form <- str_remove(trial_form, glue("{exposure} \\+ "))
    }
    # Tack on interactions
    if(!is.null(interaction_term)){
      trial_form <- paste0(trial_form, " + ", interaction_term)
    }
    # If weighting and not doing predicted probs, we can use survey::svycoxph()
    if(weighted & !predictions){
      model_list <- map(data_list, \(x) svycoxph(formula(trial_form), design = svydesign(id = ~1, weights = ~weight, data = x)))
    # {marginaleffects} doesn't support svycoxph, so use survival with weights here. CIs might be off but aren't needed
    } else if(weighted & predictions){
      model_list <- map(data_list, \(x) survival::coxph(formula(trial_form), weights = x$weight, data = x, robust = TRUE))
    # If not weighting, use normal Cox model
    } else {
      model_list <- map(data_list, \(x) survival::coxph(formula(trial_form), data = x))
    }
  # If there IS a max age, we're running a logit model
  } else {
    # Generate formula
    trial_form <- glue::glue("{outcome} ~  {exposure} + sex + race + cohort + ses + pc_fborn + fam_crim + percblack_{covars} + poverty_{covars} + pop_density_{covars} + homicide_{covars}")
    # If no exposure, delete exposure from formula
    if(exposure == "none"){
      trial_form <- str_remove(trial_form, glue("{exposure} \\+ "))
    }
    # Tack on interactions
    if(!is.null(interaction_term)){
      trial_form <- paste0(trial_form, " + ", interaction_term)
    }
    # If weighting, we can use svyglm() with bias-reducing penalized likelihood from brglm
    if(weighted){
      model_list <- map(data_list, \(x) svyglm(formula(trial_form), design = svydesign(id = ~1, weights = ~weight, data = x), family = binomial(logit), method = "brglmFit"))
    # Otherwise, use glm with brglm
    } else {
      model_list <- map(data_list, \(x) glm(formula(trial_form), family = binomial(logit), method = "brglmFit", data = x))
    }
  }
  # If not predicting probs, we can just pool to get correct standard errors
  if (!predictions){
    pooled     <- mice::pool(model_list) |> summary()
  # If predicting, we're just after means so not using usual pooling; CIs won't be valid here!
  } else if (predictions){
    # If not using a strata var, get average predictions over everything
    if(is.null(pred_var)){
      # If no max age, Cox model; get the predicted survival
      if(!contrast){
        if(is.null(outcome_max_age)){
        # marginaleffects seems to have broken this way of doing this:
        # pooled <- lapply(seq_along(model_list), \(i) 
        #                  avg_predictions(model_list[[i]], newdata = data_list[[i]], 
        #                                  variables = setNames(list(c(0,1), 40), c(exposure, outcome_age)),
        #                                  type = "survival"))
        
          pooled <- lapply(seq_along(model_list), \(i) 
                           predictions(model_list[[i]], newdata = datagrid(model = model_list[[i]], exposed = 0:1, cc_age = cox_max_age, grid_type = "counterfactual"), 
                                       type = "survival") |> 
                             tidy() |> 
                             group_by(exposed) |> 
                             summarize(estimate = mean(estimate), conf.low = NA, conf.high = NA))         
        
      # If there is a max age, it is a log model; just leave default for probs
        } else if(!is.null(outcome_max_age)) {
          pooled <- lapply(seq_along(model_list), \(i) 
                           avg_predictions(model_list[[i]], newdata = data_list[[i]], variables = setNames(list(c(0,1)), exposure)))
        }
        # For either model, just take the means of the estimates across all datasets to get PE; no uncertainty needed here
        pooled <- pooled |> 
          list_rbind(names_to = "m") |>
          group_by(pick({{exposure}})) |>
          summarize(estimate = mean(estimate), ub = mean(conf.high), lb = mean(conf.low),  .groups = "drop")
        # If Cox model, flip survival to exposure; again, CIs here prob wrong, so upper/lower bounds are likely wrong
        if(is.null(outcome_max_age)){
          pooled <-  pooled |> 
            mutate(new_estimate = 1-estimate, new_ub = 1-lb, new_lb = 1-ub) |>
            mutate(estimate = new_estimate, ub = new_ub, lb = new_lb) |>
            select(-starts_with("new_"))
        }
      } else if (contrast) {
        # Contrasts need proper uncertainty calculation, but we can do this with avg_comparisons()
          if(is.null(outcome_max_age)){
            pooled <- lapply(seq_along(model_list), \(i) 
                             avg_comparisons(model_list[[i]], newdata = datagrid(model = model_list[[i]], cc_age = cox_max_age, grid_type = "counterfactual"), 
                                           type = "survival", variable = "exposed")) |> 
              mice::pool() |> 
              tidy(conf.int = TRUE) |> 
              # If Cox model, flip survival to exposure
              mutate(estimate = -1*estimate, ub = -1*conf.low, lb = -1*conf.high, statistic = -1*statistic) |>
              select(term, estimate, ub, lb, statistic, std.error, p.value)
          } else {
            pooled <- lapply(seq_along(model_list), \(i) 
                             avg_comparisons(model_list[[i]], newdata = data_list[[i]], variable = "exposed", type = "response")) |> 
              mice::pool() |> 
              tidy(conf.int = TRUE) |>
              select(term, estimate, ub = conf.high, lb = conf.low, statistic, std.error, p.value)
          }
        }
    # IF there is a strata var, do predictions across the strata; note this is only implemented for the logit models and will fail for survival
    } else {
      pooled <- lapply(seq_along(model_list), \(i) 
                       avg_predictions(model_list[[i]], newdata = data_list[[i]], variables = setNames(list(c(0,1), pred_vals), c(exposure, pred_var))))
      pooled <- pooled |> 
        list_rbind(names_to = "m") |>
        group_by(pick({{pred_var}}), pick({{exposure}})) |>
        summarize(estimate = mean(estimate), ub = mean(conf.high), lb = mean(conf.low), .groups = "drop")
    }
  }
  return(pooled)
}
