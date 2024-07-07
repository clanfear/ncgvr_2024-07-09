standardize <- function(input, scale=1){
  x <- as.numeric(input)
  return((x- mean(x, na.rm=T))/(scale*sd(x, na.rm=T)))
}

valcheck <- function(df){
  numerics <- df[,sapply(df, is.numeric)]
  sapply(numerics, function(x){
    x_no_na <- x[!is.na(x)]
    c("min"     = min(x_no_na),
      "mean"    = mean(x_no_na),
      "median"  = median(x_no_na),
      "max"     = max(x_no_na),
      "sd"      = sd(x_no_na),
      "n_missing" = sum(is.na(x)),
      "pr_missing" = sum(is.na(x))/length(x))
  })
}

`%!in%` <- Negate(`%in%`)

x_notin_y <- function(x, y){
  setNames(list(
    unique(x)[unique(x) %!in% unique(y)],
    unique(y)[unique(y) %!in% unique(x)]
  ), 
  c(paste0(deparse1(substitute(x)), " not in ", deparse1(substitute(y))),
    paste0(deparse1(substitute(y)), " not in ", deparse1(substitute(x)))
  )
  )
}

x_in_y <- function(x, y){
  setNames(list(
    unique(x)[unique(x) %in% unique(y)],
    unique(y)[unique(y) %in% unique(x)]
  ), 
  c(paste0(deparse1(substitute(x)), " in ", deparse1(substitute(y))),
    paste0(deparse1(substitute(y)), " in ", deparse1(substitute(x)))
  )
  )
}


list_missing <- function(x){
  missings <- lapply(x, function(x) sum(is.na(x)))
  data.frame(missing = unlist(missings)) |> 
    arrange(desc(missing)) |>
    mutate(pr_missing = round(missing / nrow(x), 3))
}

invert <- function (x) # This is from searchable package which is deprecated
{
  if (is.null(names(x))) 
    stop("vector does not have names.")
  v <- names(x)
  names(v) <- as.character(x)
  return(v)
}

lme_reliability_2lvl <- function(x){
  var_components <- insight::get_variance(x)
  t00 <- var_components$var.intercept[[1]]
  s2  <- var_components$var.residual
  n_count <- table(insight::get_random(x))
  J <- length(n_count) # number of neighbs
  sum(t00 / (t00 + s2 / n_count)) / J
}

lme_reliability_3lvl <- function(x){
  # This is RRK 1991 calculation
  var_components <- insight::get_variance(x)
  tb <- var_components$var.intercept[[2]]
  tp <- var_components$var.intercept[[1]]
  s2  <- var_components$var.residual
  rel <- x@frame %>% 
    count(year_ward, id) %>%
    group_by(year_ward) %>%
    summarize(jk = n(),
              njk = sum(n)) %>%
    mutate(reliability = tb / (tb + tp/jk + s2/(njk))) %>% 
    pull(reliability)
  return(mean(rel))
}

# read_cohort() takes a wave and battery name and grabs the relevant data from
# the folders, strips out duplicate info from the master, and returns a single
# dataframe with a row per respondent ready to bind to others. All data are
# returned as character data. Anything that was labelled will have the format
# value_number:value_label so it can be split up.
read_cohort <- function(wave, battery, master_vars = NULL){
  battery_pattern <- paste0("W", wave, "_", battery)
  battery_dir <- list.files(paste0(cohort_study_path, "W", wave, "/"), pattern = battery_pattern, full.names = TRUE)
  battery_files <- list.files(battery_dir, pattern = "\\.por", recursive=TRUE, full.names = TRUE)
  if(battery == "Master"){
    drop_vars <-   ""
  } else {
    drop_vars <- c(master_vars, "BASISID", "WAVE", "MONTH", "DAY", "YEAR")
  }
  return(map_dfr(battery_files, ~ read_spss(.) |> 
                   mutate(across(where(is.labelled), 
                                 ~as.character(as_factor(.)))) |>
                   mutate(across(everything(), ~as.character(.))) |>
                   select(-any_of(drop_vars))))
}

# Function to estimate SPT models with arbitrary right hand side, accepts either
# a formula or text of formula and can filter data arbitrarily

estimate_ic_sp <- function(form, outcome, dlist, .filter = TRUE, bs_samples = 500, weights = FALSE){
  if(!exists("cl")){
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
  }
  if(is_formula(form)){
    rhs <- deparse1(form)
  } else if (is.character(form)) {
    rhs <- form
  }
  lhs <- glue::glue("Surv({outcome}_left, {outcome}_right, type = 'interval2')")
  input_form <- formula(paste0(lhs, rhs))
  na_check_var <- glue::glue("{outcome}_left")
  input_df <- dlist[[outcome]] |>
    filter(if_all(matches("_left|_right"), ~!is.na(.))) |>
    filter({{.filter}}) |>
    droplevels()
  input_weights <- NULL
  if(weights){
    input_df <- input_df |> 
      filter(!is.na(weight))
    input_weights <- input_df |> 
      pull(weight) %>% 
      {(./sum(.))*length(.)}
  }
  output_estimates <- ic_sp(input_form, 
                            bs_samples = bs_samples, useMCores = TRUE,
                            data = input_df, model = "ph", weights = input_weights)
  return(output_estimates)
}

estimate_ic_sp_imp <- function(form, df, .filter = TRUE, bs_samples = 500, weights = FALSE, np = TRUE, dist = "lnorm"){
  input_df <- df %>% 
    filter({{.filter}}) %>%
    droplevels()
  input_weights <- NULL
  if(weights){
    input_df <- input_df %>% filter(!is.na(weight))
    input_weights <- input_df %>% 
      pull(weight) %>% 
      {(./sum(.))*length(.)}
  }
  if(np){
    if(!exists("cl")){
      require(doParallel)
      cl <- makeCluster(detectCores())
      registerDoParallel(cl)
    }
    output_estimates <- ic_sp(form, 
                              bs_samples = bs_samples, useMCores = TRUE,
                              data = input_df, model = "ph", weights = input_weights)
  } else if (!np){
    output_estimates <- ic_par(form, 
                               dist = dist,
                               data = input_df, model = "aft", weights = input_weights)
  }
  return(output_estimates)
}
# A crude tidy() for icenReg models
get_spt_estimates <- function(model_list, outcome, zval = 1.96){
  estimates <- summary(model_list[[outcome]])$summaryParameters %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    janitor::clean_names() %>%
    transmute(outcome = outcome,
              term,
              estimate, 
              std.error = std_error, 
              conf.low  = estimate - zval*std.error, 
              conf.high = estimate + zval*std.error) %>%
    mutate(across(c(estimate, conf.low, conf.high), ~exp(.)))
  return(estimates)
}

# Interval censoring plots
interval_plot <- function(df, title = NULL){
  mutate(df, midpoint = (age_left + age_right) / 2 ,
         exact_age = ifelse(near(age_right - age_left,1 ), "Exact Age", "Interval")) %>% 
    arrange(midpoint, cohort) %>% 
    mutate(row = row_number()) %>%
    ggplot(aes(y = row, color = cohort)) + 
    geom_errorbarh(aes(xmin = age_left, xmax = age_right), height = 0) +
    geom_line(aes(x = midpoint, y = row), inherit.aes = FALSE, alpha = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_viridis_d(option = "D") +
    theme_minimal() +
    labs(x = "Age", y = NULL, color = "Cohort", title = title) +
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          axis.text.y = element_blank(),
          text = element_text(family = "serif"))
}

# Nice proportion tabulations
make_prop_table <- function(tab, denominator = "row"){
  results_table <- the_counts <- tab %>%
    adorn_totals(where = c("row", "col"))
  the_percents <- the_counts %>%
    adorn_percentages(denominator)
  results_table[, -1] <- sapply(2:ncol(the_counts), function(x){
    paste0(the_counts[[x]], "\n(", round(the_percents[[x]], 3)*100, "%)")
  })
  return(results_table)
}

# Extracting curves from icfit and ic_np objects for use in ggplot
# These do not return equivalent dfs, but they will plot identically.

extract_lines <-  function(x, range=NULL){
  intmap <- x$intmap
  pf <- x$pf
  if(!is.null(range)){
    intmap <- intmap[,range]
    pf <- pf[range]
  }
  time <- c(0, as.vector(intmap))
  S <- c(1, 1 - cumsum(pf))
  S <- rep(S, each = 2)[-2 * length(S)]
  cdf = (1 - S)
  time[time == Inf] <- max(time)
  return(data.frame(time, cdf))
}

extract_curve_df <- function(mod){
  if(inherits(mod, "icfit")){
    if(length(mod$strata)==1){
      return(data.frame(strata = 1, extract_lines(mod)))
    } else {
      curve_strata <- 
        data.frame(strata = map(names(mod$strata), \(y) rep(y, times = mod$strata[y])) |> list_c()) |>
        mutate(index = row_number()) %>% 
        split(~.$strata) |> 
        map(\(i) extract_lines(mod, range = i$index)) |> 
        list_rbind(names_to = "strata") |>
        mutate(strata = str_remove(strata, "^.*=")) |>
        tibble()
    }
  } else if (inherits(mod, "ic_npList")){
    strata_names <- setNames(names(mod$scurves), names(mod$scurves))
    curve_strata <- map(strata_names, \(x){
      data.frame(intervals = mod$scurves[[x]]$Tbull_ints,
                 cdf =  1-mod$scurves[[x]]$S_curves$baseline) |>
        pivot_longer(c(intervals.lower, intervals.upper), names_to = "interval", values_to = "time") |>
        select(time, cdf) |>
        mutate(cdf = lag(cdf)) |>
        slice(-1) |>
        mutate(time = ifelse(row_number()==1, 0, time))
    }) |>
      list_rbind(names_to = "strata")
  }
  return(curve_strata)
}

# LOL had to adjust the imputation function to use SP models
imputeCens_np <- function (fit, newdata = NULL, imputeType = "fullSample", samples = 1) {
  if (is.null(newdata)) 
    newdata <- fit$getRawData()
  yMat <- icenReg:::expandY(fit$formula, newdata, fit)
  p1 <- icenReg:::getFitEsts(fit, newdata, q = as.numeric(yMat[, 1]))
  p2 <- icenReg:::getFitEsts(fit, newdata, q = as.numeric(yMat[, 2]))
  p2 <- ifelse(is.nan(p2), 1, p2)
  ans <- matrix(nrow = length(p1), ncol = samples)
  storage.mode(ans) <- "double"
  
  if (imputeType == "fixedParSample") {
    isBayes <- is(fit, "bayes_fit")
    if (isBayes) {
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    for (i in 1:samples) {
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[, 1]
      theseImputes[isLow] <- yMat[isLow, 1]
      isHi <- theseImputes > yMat[, 2]
      theseImputes[isHi] <- yMat[isHi, 2]
      rownames(ans) <- rownames(newdata)
      ans <- fastMatrixInsert(theseImputes, ans, colNum = i)
    }
    if (isBayes) 
      setSamplablePars(fit, orgCoefs)
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  
  if (imputeType == "fullSample") {
    isSP <- is(fit, "sp_fit")
    isBayes <- is(fit, "bayes_fit")
    for (i in 1:samples) {
      orgCoefs <- icenReg:::getSamplablePars(fit)
      if (isBayes) {
        sampledCoefs <- icenReg:::sampBayesPar(fit)
      }
      else if (!isSP) {
        coefVar <- icenReg:::getSamplableVar(fit)
        sampledCoefs <- icenReg:::sampPars(orgCoefs, coefVar)
      }
      else {
        sampledCoefs <- icenReg:::getBSParSample(fit)
      }
      icenReg:::setSamplablePars(fit, sampledCoefs)
      p1 <- icenReg:::getFitEsts(fit, newdata, q = as.numeric(yMat[, 
                                                                   1]))
      p2 <- icenReg:::getFitEsts(fit, newdata, q = as.numeric(yMat[, 
                                                                   2]))
      p2 <- ifelse(is.nan(p2), 1, p2)
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- icenReg:::getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[, 1]
      theseImputes[isLow] <- yMat[isLow, 1]
      isHi <- theseImputes > yMat[, 2]
      theseImputes[isHi] <- yMat[isHi, 2]
      icenReg:::fastMatrixInsert(theseImputes, ans, colNum = i)
      icenReg:::setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  stop("imputeType type not recognized.")
}

# trim_weights() is from the pewmethods package which has no CRAN release so I
# just nabbed it here. See https://github.com/pewresearch/pewmethods/blob/master/R/trim_weights.R
trim_weights <- function (weight, lower_quantile = 0.01, upper_quantile = 0.99, minval = NULL, maxval = NULL, strict = FALSE) {
  design <- svydesign(ids = ~1, data = data.frame(weight = weight), weights = ~weight)
  lower <- quantile(weights(design), lower_quantile)
  upper <- quantile(weights(design), upper_quantile)
  
  if (!is.null(minval)) {
    lower <- ifelse(lower >= minval, lower, minval)
  }
  
  if (!is.null(maxval)) {
    upper <- ifelse(upper <= maxval, upper, maxval)
  }
  
  trimmed <- trimWeights(design, lower = lower, upper = upper, strict = strict)
  return(weights(trimmed))
}



