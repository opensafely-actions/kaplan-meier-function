
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Get disclosure-safe Kaplan-Meier estimates.
# The function requires an origin date, an event date, and a censoring date, which are converted into a (time , indicator) pair that is passed to `survival::Surv`
# Estimates are stratified by the `exposure` variable, and additionally by any `subgroups`
# Counts are rounded to midpoint values defined by `count_min`.
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----

library('here')
library('glue')
library('tidyverse')
library('survival')


## import local functions ----

source(here::here("analysis", "time-rounding.R"))

## parse command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  df_input <- "output/extract.arrow"
  dir_output <- "output/km_estimates/"
  exposure <- c("sex")
  subgroups <- c("age_group")
  origin_date <- "first_vax_date"
  event_date <- "second_vax_date"
  censor_date <- character() # "censor_date"
  weight <- character()
  min_count <- as.integer("6")
  method <- "constant"
  max_fup <- as.numeric("365")
  smooth <- as.logical("FALSE")
  smooth_df <- as.integer("4")
  concise <- as.logical("TRUE")
  plot <- as.logical("FALSE")
  contrast <- as.logical("TRUE")
  filename_suffix <- as.character("")
} else {

  library("optparse")

  option_list <- list(
    make_option("--df_input", type = "character",
                help = "[default: Must be specified] character. The input dataset .arrow filename. feather/arrow format is enforced to ensure date types are preserved.",
                metavar = "df_input.arrow"),
    make_option("--dir_output", type = "character",
                help = "[default: must be specified] character. The output directory. All requested output files (eg 'estimates.arrow', 'contrasts.arrow') will be placed in this directory. See also: 'filename_suffix' argument.",
                metavar = "/output/"),
    make_option("--exposure", type = "character", default = character(),
                help = "[default: NULL] character. The name of an exposure variable in the input dataset. Must be binary or not given. All outputs will be stratified by this variable. This could be an exposure in the usual sense, or it could (mis)used to show different types of events (as long as the censoring structure is the same). If not specified, no stratification will occur.",
                metavar = "exposure_varname"),
    make_option("--subgroups", type = "character", default = character(),
                help = "[default: NULL] The name of a subgroup variable or list of variable names. If a subgroup variable is used, analyses will be stratified as exposure * ( subgroup1, subgroup2, ...). If not specified, no stratification will occur.",
                metavar = "subgroup_varname"),
    make_option("--origin_date", type = "character",
                help = "[default: must be specified] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') in the input dataset that represents the start of follow-up.",
                metavar = "origin_varname"),
    make_option("--event_date", type = "character",
                help = "[default: must be specified] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') in the input dataset that represents the event date.",
                metavar = "event_varname"),
    make_option("--censor_date", type = "character", default = character(),
                help = "[default: NULL] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') that represents the censoring date. If not specified, then no censoring occurs except at `max_fup` time.",
                metavar = "censor_varname"),
    make_option("--weight", type = "character", default = character(),
                help = "[default: NULL] The name of a numeric variable that represents balancing / sampling weights. If not specified, then no weighting occurs.",
                metavar = "censor_varname"),
    make_option("--min_count", type = "integer", default = 6,
                help = "[default: %default] integer. The minimum permissable event and censor counts for each 'step' in the KM curve. This ensures that at least `min_count` events occur at each event time.",
                metavar = "min_count"),
    make_option("--method", type = "character", default = "constant",
                help = "[default: %default] character. The interpolation method after rounding. The 'constant' method leaves the event times unchanged after rounding, making the KM curve have bigger, fewer steps. The 'linear' method linearly interpolates between rounded events times (then rounds to the nearest day), so that the steps appear more natural.",
                metavar = "method"),
    make_option("--max_fup", type = "numeric", default = Inf,
                help = "[default: %default] numeric. The maximum follow-up time after the origin date. If event variables are dates, then this will be days.",
                metavar = "max_fup"),
    make_option("--smooth", type = "logical", default = FALSE,
                help = "[default: %default] logical. Should Kaplan-Meier estimates be smoothed on the log cumulative hazard scale (TRUE) or not (FALSE). ",
                metavar = "TRUE/FALSE"),
    make_option("--smooth_df", type = "logical", default = 4,
                help = "[default: %default]. integer. Degrees of freedom to use for the smoother. Unused if smooth=FALSE.",
                metavar = "TRUE/FALSE"),
    make_option("--concise", type = "logical", default = TRUE,
                help = "[default: %default] logical. Should the outputted table only report core variables (defined here as exposure, subgroups, time, number at risk, cumulative number of events, cumulative incidence, and confidence limits) (TRUE) or should it report everything (FALSE)?",
                metavar = "TRUE/FALSE"),
    make_option("--plot", type = "logical", default = FALSE,
                help = "[default: %default] logical. Should Kaplan-Meier plots be created in the output directory? These are fairly basic plots for sense-checking purposes.",
                metavar = "TRUE/FALSE"),
    make_option("--contrast", type = "logical", default = TRUE,
                help = "[default: %default] logical. Should Kaplan-Meier curves for a binary exposure be compared to estimate risk difference, risk ratio, and survival ratio? Ignored if exposure is not supplied.",
                metavar = "TRUE/FALSE"),
    make_option("--filename_suffix", type = "character", default = "",
                help = "[default: %default] character. This will be appended to the end of all outputted files. This is useful if you want to re-run the KM action across different arguments, but put outputs from all actions in the same directory.",
                metavar = "TRUE/FALSE")
  )

  opt_parser <- OptionParser(
    usage = "kaplan-meier-function:[version] [options]",
    option_list = option_list
  )

  opt <- parse_args(opt_parser)
  list2env(opt, .GlobalEnv)
}

# Use quasiquotation for passing exposure and subgroup stratification variables
# around the place
# use `syms()` instead of `sym()` even though it's possible to pull
# only one exposure or subgroup variable is from the args (without hacking!)
# this ensures that if `exposure` or `subgroups` is not used,
# the quasiquotation still works inside ggplot, transmute, etc

exposure_syms <- syms(exposure)
subgroup_syms <- syms(subgroups)

# Create output directory ----

dir_output <- here::here(dir_output)
fs::dir_create(dir_output)


# Import and process person-level data  ----

## Import ----
data_patients <-
  arrow::read_feather(here::here(df_input))

## Derive variables ----

if(length(censor_date)==0) {
  # if censor date is not specified, then create a censor_date variable in the dataset, taking value `as.Date(Inf)`
  data_patients$censor_date <- as.Date(Inf)
  censor_date <- "censor_date"
}

is.weighted <- length(weight)>0
if(!is.weighted) {
  # if weight is not specified, then create a weight variable in the dataset, taking value `1L`
  data_patients$.weight <- 1L
  weight <- ".weight"
}

data_tte <-
  data_patients |>
  transmute(
    patient_id,
    !!!exposure_syms,
    !!!subgroup_syms,
    .weight = .data[[weight]],
    event_date = as.Date(.data[[event_date]]),
    origin_date = as.Date(.data[[origin_date]]),
    censor_date = pmin(
      as.Date(.data[[censor_date]]),
      origin_date + max_fup,
      na.rm=TRUE
    ),
    event_time = tte(origin_date, event_date, censor_date, na.censor=FALSE),
    event_indicator = censor_indicator(event_date, censor_date),
  )

if(max_fup==Inf) max_fup <- max(data_tte$event_time)+1

## tests ----

if(length(exposure)>0){
  stopifnot("exposure variable must be binary or have two levels" = (length(unique(data_patients[[exposure]])) == 2))
}

stopifnot("censoring dates must be non-missing" = all(!is.na(data_tte$censor_date)))

stopifnot("origin dates must be non-missing" = all(!is.na(data_tte$origin_date)))

times_count <- table(cut(data_tte$event_time, c(-Inf, 0, 1, Inf), right=FALSE, labels= c("<0", "0", ">0")), useNA="ifany")

if(!identical(as.integer(times_count), c(0L, 0L, nrow(data_tte)))) {
  print(times_count)
  stop("all event times must be strictly positive")
}

# Calculate max follow-up time available in the data ----
# and print to log file

if(length(exposure)>0){
  max_time_data <-
    data_tte |>
    group_by(!!!exposure_syms) |>
    summarise(
      max_fup_time = max(event_time),
      max_event_time = max(event_time[event_indicator])
    )

  cat("maximum follow-up time in exposure levels [", paste0(max_time_data[[exposure]], collapse=", "), "] is [", paste0(max_time_data$max_fup_time, collapse= ", "), "]", "\n")
  cat("maximum event time in exposure levels [", paste0(max_time_data[[exposure]], collapse=", "), "] is [", paste0(max_time_data$max_event_time, collapse= ", "), "]", "\n")
} else {
  max_time_data <-
    data_tte |>
    summarise(
      max_fup_time = max(event_time),
      max_event_time = max(event_time[event_indicator])
    )
  cat("maximum follow-up time is [", paste0(max_time_data$max_fup_time, collapse= ", "), "]", "\\n")
  cat("maximum event time is [", paste0(max_time_data$max_event_time, collapse= ", "), "]", "\\n")
}


# Calculate KM estimates ------

## Run `survfit` across each level of exposure and subgroup ----
## do this independently rather than using stratification or covariates
## because it makes variable name handling easier

# for each exposure level and subgroup level, pass data through `survival::survfit` to get KM table
data_surv <-
  data_tte |>
  dplyr::group_by(!!!subgroup_syms, !!!exposure_syms) |>
  tidyr::nest() |>
  dplyr::mutate(
    surv_obj_tidy = purrr::map(data, ~ {
      survival::survfit(
        survival::Surv(event_time, event_indicator) ~ 1,
        data = .x,
        conf.type="log-log",
        weight = .weight
      ) |>
      broom::tidy() |>
      tidyr::complete(
        time = seq_len(max_fup), # fill in 1 row for each day of follow up
        fill = list(n.event = 0L, n.censor = 0L) # fill in zero events on those days
      ) |>
      tidyr::fill(n.risk, .direction = c("up"))
    }),
  ) |>
  select(-data) |>
  tidyr::unnest(surv_obj_tidy)


## Round the count values in the survival data ----
# round event times such that no event time has fewer than `min_count` events
# recalculate KM estimates based on these rounded event times

round_km <- function(.data, min_count=0, method="constant") {

  # min_count == 0 means no rounding.
  # precision is 0 rather than 1 because if using weighting, then
  # there may be non-integer counts
  if(min_count==0){
    rounded_data <-
      .data |>
      mutate(
        N = max(n.risk, na.rm = TRUE),
        # rounded to `min_count - (min_count/2)`
        cml.event = cumsum(n.event),
        cml.censor = cumsum(n.censor),
        cml.eventcensor = cml.event + cml.censor,
      )
  } else {
    rounded_data <-
      .data |>
      mutate(
        N = max(n.risk, na.rm = TRUE),
        # rounded to `min_count - (min_count/2)`
        cml.event = round_cmlcount(cumsum(n.event), time, min_count, method, integer.counts = !is.weighted),
        cml.censor = round_cmlcount(cumsum(n.censor), time, min_count, method, integer.counts = !is.weighted),
        cml.eventcensor = cml.event + cml.censor,
        # re-derive counts from cumulative data
        n.event = diff(c(0, cml.event)),
        n.censor = diff(c(0, cml.censor)),
        n.risk = roundmid_any(N, min_count) - lag(cml.eventcensor, 1, 0)
      )
  }

  rounded_data1 <-
    rounded_data |>
    mutate(
      cml.nrisk = cumsum(n.risk),
      cml.rate = cml.event / cml.nrisk,
      inc = n.event / n.risk,
      # KM estimate for event of interest, combining censored and competing events as censored
      summand = (1 / (n.risk - n.event)) - (1 / n.risk), # = n.event / ((n.risk - n.event) * n.risk) but re-written to prevent integer overflow
      surv = cumprod(1 - n.event / n.risk),
      # standard errors on survival scale
      surv.se = surv * sqrt(cumsum(summand)),# greenwood's formula
      # surv.low = surv + qnorm(0.025)*surv.se,
      # surv.high = surv + qnorm(0.975)*surv.se,
      ## standard errors on log scale
      surv.ln.se = surv.se / surv,
      # surv.low = exp(log(surv) + qnorm(0.025)*surv.ln.se),
      # surv.high = exp(log(surv) + qnorm(0.975)*surv.ln.se),
      ## standard errors on complementary log-log scale
      surv.cll = log(-log(surv)), # this is equivalent to the log cumulative hazard
      surv.cll.se = if_else(surv==1, 0, sqrt((1 / log(surv)^2) * cumsum(summand))), # assume SE is zero until there are events -- makes plotting easier
      surv.low = exp(-exp(surv.cll + qnorm(0.975) * surv.cll.se)),
      surv.high = exp(-exp(surv.cll + qnorm(0.025) * surv.cll.se)),
      #cumulative incidence (= complement of survival)
      cmlinc = 1 - surv,
      cmlinc.se = surv.se,
      cmlinc.ln.se = surv.ln.se,
      cmlinc.low = 1 - surv.high,
      cmlinc.high = 1 - surv.low,
      # restricted mean survival time.
      # https://doi.org/10.1186/1471-2288-13-152
      rmst = cumsum(surv), # this only works if one row per day using fill_times! otherwise need cumsum(surv*int)
      rmst.se = sqrt(((2* cumsum(time*surv)) - (rmst^2))/n.risk), # this only works if one row per day using fill_times! otherwise need sqrt(((2* cumsum(time*interval*surv)) - (rmst^2))/n.risk)
      rmst.low = rmst + (qnorm(0.025) * rmst.se),
      rmst.high = rmst + (qnorm(0.975) * rmst.se),
    ) |>
    # filter(
    #   !(n.event==0 & n.censor==0 & !fill_times) # remove times where there are no events (unless all possible event times are requested with fill_times)
    # ) |>
    select(
      !!!subgroup_syms,
      !!!exposure_syms,
      time,
      cml.nrisk, cml.event, cml.censor,
      n.risk, n.event, n.censor,
      inc,
      #surv, surv.se, surv.low, surv.high,
      cml.rate,
      cmlinc, cmlinc.se, cmlinc.low, cmlinc.high,
      rmst, rmst.se, rmst.low, rmst.high,
    )

  rounded_data1
}

data_surv_unrounded <- round_km(data_surv, 0L)
data_surv_rounded <- round_km(data_surv, min_count, method=method)

## Select a smaller set of variable ----
## if requested via concise argument

data_surv_rounded_output <-
  if(concise){
    data_surv_rounded |>
    select(
      !!!subgroup_syms,
      !!!exposure_syms,
      time,
      n.risk, n.censor, n.event,
      cml.event,
      cmlinc, cmlinc.low, cmlinc.high
    )
  } else {
    data_surv_rounded_output <- data_surv_rounded
  }


## output to disk ----
## include both arrow and csvv formats here - if you don't want one of them,
## don't include it in the `output:` slot in the action

## write arrow to disk
arrow::write_feather(data_surv_rounded_output, fs::path(dir_output, glue("estimates{filename_suffix}.arrow")))
## write csv to disk
write_csv(data_surv_rounded_output, fs::path(dir_output, glue("estimates{filename_suffix}.csv")))


# Smoothing via parametric survival ----
## if requested via `smooth` argument

smooth_km <- function(.data, smooth_df){
  # smooth the KM curve on the complementary log-log scale (ie, smooth the log cumulative hazard)
  # using rtpm2 package
  library('rstpm2')
  data_smoothed <-
    .data |>
    dplyr::group_by(!!!subgroup_syms, !!!exposure_syms) |>
    tidyr::nest() |>
    dplyr::mutate(
      surv_obj = purrr::map(data, ~ {
        stpm2(
          survival::Surv(event_time, event_indicator) ~ 1,
          data = .x,
          df = smooth_df,
          weights = .weight
        )
      }),
      surv_smooth = purrr::map(surv_obj, ~ {
        new_data <- data.frame(event_time=seq_len(max_fup))
        surv_predict <- predict(
          .x,
          newdata = new_data,
          type = "surv",
          level = 0.95,
          se.fit = TRUE
        )
        # not yet available as "rmst currently only for single value"
        # rmst_predict <- predict(
        #   .x,
        #   newdata=new_data,
        #   type="rmst",
        #   level=0.95,
        #   se.fit=TRUE
        # )
        hazard_predict <- predict(
          .x,
          newdata = new_data,
          type = "hazard",
          level = 0.95,
          se.fit = TRUE
        )
        tibble(
          time = seq_len(max_fup),
          surv = surv_predict$Estimate,
          surv.low = surv_predict$lower,
          surv.high = surv_predict$upper,
          inc = n.event / n.risk,
          cml.nrisk = cumsum(n.risk),
          cml.rate = cml.event / cml.n.risk,
          cmlinc = 1 - surv,
          cmlinc.low = 1 - surv.high,
          cmlinc.high = 1 - surv.low,
          rmst = cumsum(surv),
          rmst.low = cumsum(surv.low),
          rmst.high = cumsum(surv.high),
          hazard = hazard_predict$Estimate,
          hazard.low = hazard_predict$lower,
          hazard.high = hazard_predict$upper,
        )
      }),
    ) |>
    dplyr::select(!!!subgroup_syms, !!!exposure_syms, surv_smooth) |>
    tidyr::unnest(surv_smooth)

  data_smoothed
}


if(smooth){
  data_surv_smoothed <- smooth_km(data_tte, smooth_df)

  ## Select a smaller set of variable
  ## if requested via concise argument

  data_surv_smoothed_output <-
    if(concise){
      data_surv_smoothed |>
      select(
        !!!subgroup_syms,
        !!!exposure_syms,
        time,
        n.risk, n.censor, n.event,
        cml.event, inc,
        cmlinc, cmlinc.low, cmlinc.high
      )
    } else {
      data_surv_smoothed
    }

  ## output to disk ----
  ## include both arrow and csv formats here - if you don't want one of them,
  ## don't include it in the `output:` slot in the action

  ## write arrow to disk
  arrow::write_feather(data_surv_smoothed_output, fs::path(dir_output, glue("estimates{filename_suffix}.arrow")))
  ## write csv to disk
  write_csv(data_surv_smoothed_output, fs::path(dir_output, glue("estimates{filename_suffix}.csv")))
}

# Plot KM curves ----
## if requested via `plot` argument

km_plot <- function(.data) {

  data_with_time0 <-
    .data |>
    mutate(
      "{exposure}" := as.factor(!!!exposure_syms),
      lagtime = lag(time, 1, 0), # assumes the time-origin is zero
    ) %>%
    group_modify(
      ~ add_row(
        .x,
        time = 0, # assumes time origin is zero
        lagtime = 0,
        cmlinc = 0,
        cmlinc.low = 0,
        cmlinc.high = 0,
        .before = 0
      )
    )
  ggplot_init <- if(length(exposure)==0L){
    ggplot(data_with_time0)
  } else {
    exposure_sym <- sym(exposure)
    ggplot(data_with_time0, aes(group = !!exposure_sym, colour = !!exposure_sym, fill = !!exposure_sym))
  }
    ggplot_init +
    geom_step(aes(x = time, y = cmlinc), direction = "vh") +
    geom_step(aes(x = time, y = cmlinc), direction = "vh", linetype = "dashed", alpha = 0.5) +
    geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high), alpha = 0.1, colour = "transparent") +
    facet_grid(rows = vars(!!!subgroup_syms)) +
    scale_color_brewer(type = "qual", palette = "Set1", na.value = "grey") +
    scale_fill_brewer(type = "qual", palette = "Set1", guide = "none", na.value = "grey") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(xlim = c(0, NA)) +
    labs(
      x = "Time",
      y = "Cumulative Incidence",
      colour = NULL,
      title = NULL
    ) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(colour = "black"),
      panel.grid.minor.x = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(.05, .95),
      legend.justification = c(0, 1),
    )
}

## output to disk ----
if(plot){
  if(smooth){
    km_plot <- km_plot(data_surv_smoothed)
  } else{
    km_plot <- km_plot(data_surv_rounded)
  }
  ggsave(filename = fs::path(dir_output, glue("plot{filename_suffix}.png")), km_plot, width = 20, height = 20, units = "cm")
}


# Calculate contrasts between exposure groups ----
# if requested via 'contrast` argument

contrast_km <- function(.data) {

  .data |>
    filter(
      time != 0
    ) |>
    ungroup() |>
    mutate(
      # convert exposure variable to a 0/1 binary variable, with 0 the reference level
      "{exposure}" := as.integer(as.factor(!!!exposure_syms)) - 1L
    ) |>
    pivot_wider(
      id_cols = all_of(c(subgroups, "time")),
      #names_glue = "{.value}_{exposure}", #but this doesn't work because scoping, quosures, something something
      names_from = c(!!!exposure_syms),
      names_sep = "_",
      values_from = c(
        n.risk, n.event, n.censor,
        cml.nrisk, cml.event, cml.rate,
        cmlinc, cmlinc.se, cmlinc.low, cmlinc.high#,
        #rmst, rmst.low, rmst.high
      )
    ) |>
    mutate(
      n.nonevent_0 = n.risk_0 - n.event_0,
      n.nonevent_1 = n.risk_1 - n.event_1,

      ## cumulative quantities using information during time [0,t] (not just [t])

      # cumulative incidence rate ratio
      cmlirr = cml.rate_1 / cml.rate_0,
      cmlirr.ln.se = sqrt((1 / cml.event_0) + (1 / cml.event_1)),
      cmlirr.ll = exp(log(cmlirr) + qnorm(0.025) * cmlirr.ln.se),
      cmlirr.ul = exp(log(cmlirr) + qnorm(0.975) * cmlirr.ln.se),

      # survival ratio, standard error, and confidence limits
      sr = (1-cmlinc_1) / (1-cmlinc_0),
      sr.ln.se = (cmlinc.se_0 / (1-cmlinc_0)) + (cmlinc.se_1 / (1-cmlinc_1)),
      sr.ll = exp(log(sr) + qnorm(0.025) * sr.ln.se),
      sr.ul = exp(log(sr) + qnorm(0.975) * sr.ln.se),

      # risk ratio, standard error, and confidence limits, using delta method
      rr = cmlinc_1 / cmlinc_0,
      # cirr.ln = log(cirr),
      rr.ln.se = sqrt((cmlinc.se_1 / cmlinc_1)^2 + (cmlinc.se_0 / cmlinc_0)^2),
      rr.ll = exp(log(rr) + qnorm(0.025) * rr.ln.se),
      rr.ul = exp(log(rr) + qnorm(0.975) * rr.ln.se),

      # risk difference, standard error and confidence limits, using delta method
      rd = cmlinc_1 - cmlinc_0,
      rd.se = sqrt((cmlinc.se_0^2) + (cmlinc.se_1^2)),
      rd.ll = rd + qnorm(0.025) * rd.se,
      rd.ul = rd + qnorm(0.975) * rd.se,
    ) |>
    select(
      # remove rows relating to individual curves
      -ends_with("0"),
      -ends_with("1"),
      -ends_with(".se"),
    )
}

if((length(exposure)>0) & contrast){
  data_contrasts_rounded <- contrast_km(data_surv_rounded)

  ## output to disk ----
  arrow::write_feather(data_contrasts_rounded, fs::path(dir_output, glue("contrasts{filename_suffix}.arrow")))
  write_csv(data_contrasts_rounded, fs::path(dir_output, glue("contrasts{filename_suffix}.csv")))
}

