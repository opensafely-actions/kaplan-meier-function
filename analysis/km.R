
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
library("optparse")

## parse command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  df_input <- "output/input.feather"
  dir_output <- "output/km_estimates/"
  exposure <- c("sex")
  subgroups <- c("previous_covid_test")
  origin_date <- "start_date"
  event_date <- "death_date"
  censor_date <- "dereg_date"
  min_count <- as.integer("6")
  method <- "linear"
  max_fup <- as.numeric("200")
  fill_times <- as.logical("TRUE")
  smooth <- as.logical("FALSE")
  smooth_df <- as.integer("4")
  plot <- as.logical("FALSE")
} else {

  option_list <- list(
    make_option("--df_input", type = "character", default = NULL,
                help = "Input dataset .feather filename [default %default]. feather format is enforced to ensure date types are preserved.",
                metavar = "filename.feather"),
    make_option("--dir_output", type = "character", default = NULL,
                help = "Output directory [default %default].",
                metavar = "output"),
    make_option("--exposure", type = "character", default = NULL,
                help = "Exposure variable name in the input dataset [default %default]. All outputs will be stratified by this variable.",
                metavar = "exposure_varname"),
    make_option("--subgroups", type = "character", default = NULL,
                help = "Subgroup variable name or list of variable names [default %default]. If subgroups are used, analyses will be stratified as exposure * ( subgroup1, subgroup2, ...). If NULL, no stratification will occur.",
                metavar = "subgroup_varnames"),
    make_option("--origin_date", type = "character", default = NULL,
                help = "Time-origin variable name in the input dataset [default %default]. Should refer to a date variable, or a character of the form YYYY-MM-DD.",
                metavar = "origin_varname"),
    make_option("--event_date", type = "character", default = NULL,
                help = "Event variable name in the input dataset [default %default]. Should refer to a date variable, or a character of the form YYYY-MM-DD.",
                metavar = "event_varname"),
    make_option("--censor_date", type = "character", default = NULL,
                help = "Censor variable name in the input dataset [default %default]. Should refer to a date variable, or a character of the form YYYY-MM-DD.",
                metavar = "censor_varname"),
    make_option("--min_count", type = "integer", default = 6,
                help = "The minimum permissable event and censor counts for each 'step' in the KM curve [default %default]. This ensures that at least `min_count` events occur at each event time.",
                metavar = "min_count"),
    make_option("--method", type = "character", default = "constant",
                help = "Interpolation method after rounding [default %default]. The 'constant' method leaves the event times unchanged after rounding, making the KM curve have bigger, fewer steps. The 'linear' method linearly interpolates between rounded events times (then rounds to the nearest day), so that the steps appear more natural.",
                metavar = "method"),
    make_option("--max_fup", type = "numeric", default = Inf,
                help = "The maximum time. If event variables are dates, then this will be days. [default %default]. ",
                metavar = "max_fup"),
    make_option("--fill_times", type = "logical", default = TRUE,
                help = "Should Kaplan-Meier estimates be provided for all possible event times (TRUE) or just observed event times (FALSE) [default %default]. ",
                metavar = "TRUE/FALSE"),
    make_option("--smooth", type = "logical", default = FALSE,
                help = "Should Kaplan-Meier estimates be smoothed on the log cumulative hazard scale (TRUE) or not (FALSE) [default %default]. ",
                metavar = "TRUE/FALSE"),
    make_option("--smooth_df", type = "logical", default = 4,
                help = "Degrees of freedom to use for the smoother [default %default]. Unused if smooth=FALSE.",
                metavar = "smooth_df"),
    make_option("--plot", type = "logical", default = TRUE,
                help = "Should Kaplan-Meier plots be created in the output folder? [default %default]. These are fairly basic plots for sense-checking purposes.",
                metavar = "TRUE/FALSE")
  )

  opt_parser <- OptionParser(usage = "km:[version] [options]", option_list = option_list)
  opt <- parse_args(opt_parser)

  df_input <- opt$df_input
  dir_output <- opt$dir_output
  exposure <- opt$exposure
  subgroups <- opt$subgroups
  origin_date <- opt$origin_date
  event_date <- opt$event_date
  censor_date <- opt$censor_date
  min_count <- opt$min_count
  method <- opt$method
  max_fup <- opt$max_fup
  fill_times <- opt$fill_times
  smooth <- opt$smooth
  smooth_df <- opt$smooth_df
  plot <- opt$plot
}

exposure_sym <- sym(exposure)
subgroup_syms <- syms(subgroups)


# create output directories ----

dir_output <- here::here(dir_output)
fs::dir_create(dir_output)

# survival functions -----

censor <- function(event_date, censor_date, na.censor=TRUE){
  # censors event_date to on or before censor_date
  # if na.censor = TRUE then returns NA if event_date>censor_date, otherwise returns min(event_date, censor_date)
  if (na.censor)
    dplyr::if_else(event_date>censor_date, as.Date(NA_character_), as.Date(event_date))
  else
    dplyr::if_else(event_date>censor_date, as.Date(censor_date), as.Date(event_date))
}

censor_indicator <- function(event_date, censor_date){
  # returns FALSE if event_date is censored by censor_date, or if event_date is NA. Otherwise TRUE
  # if censor_date is the same day as event_date, it is assumed that event_date occurs first (ie, the event happened)

  stopifnot("all censoring dates must be non-missing" =  all(!is.na(censor_date)))
  !((event_date>censor_date) | is.na(event_date))
}

tte <- function(origin_date, event_date, censor_date, na.censor=FALSE){
  # returns time-to-event date or time to censor date, which is earlier
  if (na.censor) {
    censor_date <- dplyr::if_else(censor_date>event_date, as.Date(NA), censor_date)
  }
  as.numeric(pmin(event_date-origin_date, censor_date-origin_date, na.rm=TRUE))
}

ceiling_any <- function(x, to=1){
  # round to nearest 100 millionth to avoid floating point errors
  #ceiling(plyr::round_any(x/to, 1/100000000))*to
  x - (x-1)%%to + (to-1)
}

floor_any <- function(x, to=1){
  x - x%%to
}

roundmid_any <- function(x, to=1){
  # like ceiling_any, but centers on (integer) midpoint of the rounding points
  ceiling(x/to)*to - (floor(to/2)*(x!=0))
}



round_cmlcount <- function(x, time, min_count, method="linear", integer.times=TRUE) {

  # take a vector of cumulative counts and round them according to...
  stopifnot("x must be non-descreasing" = all(diff(x)>=0))
  stopifnot("x must be integer" = all(x %% 1 ==0))

  # round events such that the are no fewer than min_count events per step
  # steps are then shifted by ` - floor(min_count/2)` to remove bias
  if(method=="constant") {
    rounded_counts <- roundmid_any(x, min_count)
  }

  # as above, but then linearly-interpolate event times between rounded steps
  # this will also linearly interpolate event if _true_ counts are safe but "steppy" -- can we avoid this by not over-interpolating?
  if(method=="linear") {
    x_ceiling <- ceiling_any(x, min_count)
    x_mid <- roundmid_any(x, min_count)
#    naturally_steppy <- which((x - x_mid) == 0)
    x_rle <- rle(x_ceiling)

    # get index locations of step increases
    steptime <- c(0,time[cumsum(x_rle$lengths)])

    # get cumulative count at each step
    stepheight <- c(0,x_rle$values)

    rounded_counts <- approx(x=steptime, y=stepheight, xout = time, method=method)$y
    if(integer.times) rounded_counts <- floor(rounded_counts)
  }

  return (rounded_counts)
}


# import and process person-level data  ----

## Import ----
data_patients <-
  arrow::read_feather(here::here(df_input)) %>%
  mutate(across(.cols = ends_with("_date"), ~ as.Date(.x)))

## Derive time to event (tte) variables ----
data_tte <-
  data_patients %>%
  transmute(
    patient_id,
    !!exposure_sym,
    !!!subgroup_syms,
    event_date = as.Date(.data[[event_date]]),
    origin_date = as.Date(.data[[origin_date]]),
    censor_date = pmin(as.Date(.data[[censor_date]]), origin_date + max_fup, na.rm=TRUE),
    event_time = tte(origin_date, event_date, censor_date, na.censor=FALSE),
    event_indicator = censor_indicator(event_date, censor_date),
  )

if(max_fup==Inf) max_fup <- max(data_tte$event_time)+1

## tests ----

stopifnot("censoring dates must be non-missing" = all(!is.na(data_tte$censor_date)))

stopifnot("origin dates must be non-missing" = all(!is.na(data_tte$origin_date)))

times_count <- table(cut(data_tte$event_time, c(-Inf, 0, 1, Inf), right=FALSE, labels= c("<0", "0", ">0")), useNA="ifany")
if(!identical(as.integer(times_count), c(0L, 0L, nrow(data_tte)))) {
  print(times_count)
  stop("all event times must be strictly positive")
}


# Get KM estimates ------

for (subgroup_i in subgroups) {

  #subgroup_i = "previous_covid_test"

  # for each exposure level and subgroup level, pass data through `survival::Surv` to get KM table
  data_surv <-
    data_tte %>%
    dplyr::mutate(
      .subgroup = .data[[subgroup_i]]
    ) %>%
    dplyr::group_by(.subgroup, !!exposure_sym) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      surv_obj = purrr::map(data, ~ {
        survival::survfit(survival::Surv(event_time, event_indicator) ~ 1, data = .x, conf.type="log-log")
      }),
      surv_obj_tidy = purrr::map(surv_obj, ~ {
        tidied <- broom::tidy(.x)
        if(fill_times){ # return survival table for each day of follow up
          tidied <- tidied %>%
          tidyr::complete(
            time = seq_len(max_fup), # fill in 1 row for each day of follow up
            fill = list(n.event = 0, n.censor = 0) # fill in zero events on those days
          ) %>%
          tidyr::fill(n.risk, .direction = c("up")) # fill in n.risk on each zero-event day
        }
      }),
    ) %>%
    dplyr::select(.subgroup, !!exposure_sym, surv_obj_tidy) %>%
    tidyr::unnest(surv_obj_tidy)

  # round event times such that no event time has fewer than `min_count` events
  # recalculate KM estimates based on these rounded event times
  round_km <- function(.data, min_count) {
    .data %>%
      mutate(
        N = max(n.risk, na.rm = TRUE),

        # rounded to `min_count - (min_count/2)`
        cml.event = round_cmlcount(cumsum(n.event), time, min_count),
        cml.censor = round_cmlcount(cumsum(n.censor), time, min_count),
        cml.eventcensor = cml.event + cml.censor,
        n.event = diff(c(0, cml.event)),
        n.censor = diff(c(0, cml.censor)),
        n.risk = roundmid_any(N, min_count) - lag(cml.eventcensor, 1, 0),

        # KM estimate for event of interest, combining censored and competing events as censored
        summand = (1 / (n.risk - n.event)) - (1 / n.risk), # = n.event / ((n.risk - n.event) * n.risk) but re-written to prevent integer overflow
        surv = cumprod(1 - n.event / n.risk),

        # standard errors on survival scale
        surv.se = surv * sqrt(cumsum(summand)), # greenwood's formula
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

        #risk (= complement of survival)
        risk = 1 - surv,
        risk.se = surv.se,
        risk.ln.se = surv.ln.se,
        risk.low = 1 - surv.high,
        risk.high = 1 - surv.low
      ) %>%
      filter(
        !(n.event==0 & n.censor==0 & !fill_times) # remove times where there are no events (unless all possible event times are requested with fill_times)
      ) %>%
      mutate(
        lagtime = lag(time, 1, 0), # assumes the time-origin is zero
        interval = time - lagtime,
      ) %>%
      transmute(
        .subgroup_var = subgroup_i,
        .subgroup,
        !!exposure_sym,
        time, lagtime, interval,
        cml.event, cml.censor,
        n.risk, n.event, n.censor,
        surv, surv.se, surv.low, surv.high,
        risk, risk.se, risk.low, risk.high
      )
  }

  #data_surv_unrounded <- round_km(data_surv, 1)
  data_surv_rounded <- round_km(data_surv, min_count)
  ## write to disk
  arrow::write_feather(data_surv_rounded, fs::path(dir_output, glue("km_estimates_{subgroup_i}.feather")))

  if(smooth){

    # smooth the KM curve on the complementary log-log scale (ie, smooth the log cumulative hazard)
    # using rtpm2 pacakge

    library('rstpm2')

    data_surv_smoothed <-
      data_tte %>%
      dplyr::mutate(
        .subgroup = .data[[subgroup_i]]
      ) %>%
      dplyr::group_by(.subgroup, !!exposure_sym) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        surv_obj = purrr::map(data, ~ {
          stpm2(survival::Surv(event_time, event_indicator) ~ 1, data = .x, df=smooth_df)
        }),

        surv_smooth = purrr::map(surv_obj, ~ {

          new_data <- data.frame(event_time=seq_len(max_fup))

          surv_predict <- predict(
            .x,
            newdata=new_data,
            type="surv",
            level=0.95,
            se.fit=TRUE
          )


          hazard_predict <- predict(
            .x,
            newdata=new_data,
            type="hazard",
            level=0.95,
            se.fit=TRUE
          )

          tibble(
            time = seq_len(max_fup),
            lagtime = lag(time, 1, 0), # assumes the time-origin is zero
            interval = time - lagtime,
            surv = surv_predict$Estimate,
            surv.low = surv_predict$lower,
            surv.high = surv_predict$upper,
            risk = 1 - surv,
            risk.low = 1 - surv.high,
            risk.high = 1 - surv.low,
            hazard = hazard_predict$Estimate,
            hazard.low = hazard_predict$lower,
            hazard.high = hazard_predict$upper,
          )
        }),
      ) %>%
      dplyr::select(.subgroup, !!exposure_sym, surv_smooth) %>%
      tidyr::unnest(surv_smooth)

      ## write to disk
      arrow::write_feather(data_surv_smoothed, fs::path(dir_output, glue("km_estimates_{subgroup_i}.feather")))
  }

  if(plot){

    km_plot <- function(.data) {
      .data %>%
        group_modify(
          ~ add_row(
            .x,
            time = 0, # assumes time origin is zero
            lagtime = 0,
            surv = 1,
            surv.low = 1,
            surv.high = 1,
            risk = 0,
            risk.low = 0,
            risk.high = 0,
            .before = 0
          )
        ) %>%
        ggplot(aes(group = !!exposure_sym, colour = !!exposure_sym, fill = !!exposure_sym)) +
        geom_step(aes(x = time, y = risk), direction = "vh") +
        geom_step(aes(x = time, y = risk), direction = "vh", linetype = "dashed", alpha = 0.5) +
        geom_rect(aes(xmin = lagtime, xmax = time, ymin = risk.low, ymax = risk.high), alpha = 0.1, colour = "transparent") +
        facet_grid(rows = vars(.subgroup)) +
        scale_color_brewer(type = "qual", palette = "Set1", na.value = "grey") +
        scale_fill_brewer(type = "qual", palette = "Set1", guide = "none", na.value = "grey") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
        coord_cartesian(xlim = c(0, NA)) +
        labs(
          x = "Days since origin",
          y = "Kaplan-Meier estimate",
          colour = NULL,
          title = NULL
        ) +
        theme_minimal() +
        theme(
          axis.line.x = element_line(colour = "black"),
          panel.grid.minor.x = element_blank(),
          legend.position = c(.05, .95),
          legend.justification = c(0, 1),
        )
    }

    km_plot_rounded <- km_plot(data_surv_rounded)
    ggsave(filename = fs::path(dir_output, glue("km_plot_rounded_{subgroup_i}.png")), km_plot_rounded, width = 20, height = 20, units = "cm")

    if(smooth){
      km_plot_smoothed <- km_plot(data_surv_smoothed)
      ggsave(filename = fs::path(dir_output, glue("km_plot_smoothed_{subgroup_i}.png")), km_plot_smoothed, width = 20, height = 20, units = "cm")
    }
  }
}
