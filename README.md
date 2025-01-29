

<!-- README.md is generated from the README.qmd file. 
Please edit that file and re-render in R `....` -->

# kaplan-meier-function

This is the code and configuration for `kaplan-meier-function`, an R
reusable action for the OpenSAFELY framework.

This action uses the Kaplan-Meier estimator to calculate the cumulative
incidence of an event over time, from a given origin date, possibly
stratified by an exposure variable and/or a subgroup variable. Event
times may be censored.

There is EXPERIMENTAL support for statistical disclosure control in the
action, with further assurance tests and methodological explanations
pending. We recommend comparing rounded KM estimates to the unrounded KM
estimates (prior to release) to ensure that the rounded estimates are a
sensible approximation.

Additional arguments provide control over smoothing, maximum follow-up
duration, plotting, and other options – a full list of available
arguments is given below.

## Usage

This action should be specified in the `project.yaml` file of your
OpenSAFELY study repository, where you should replace `[version]` with
[the latest
tag](https://github.com/opensafely-actions/kaplan-meier-function/tags),
e.g., `v0.0.1`. Note that no space is allowed between
`kaplan-meier-function:` and `[version]`.

``` yaml

my_kaplan_meier_function_action:
  run: kaplan-meier-function:v0.0.1
    --df_input=cohort.csv
    --df_output=output/km_estimates.csv
    ...# more arguments here
  outputs:
    highly_sensitive:
      output: output/km_estimates/.*feather

...
```

The arguments to the action are specified using the flags style, i.e.,
`--argname=argvalue`. The available arguments to this action are as
follows:

    Usage: km:[version] [options]


    Options:
        --df_input=FILENAME.FEATHER
            Input dataset .feather filename [default: NULL]. feather format is enforced to ensure date types are preserved.

        --dir_output=OUTPUT
            Output directory [default: NULL].

        --exposure=EXPOSURE_VARNAME
            Exposure variable name in the input dataset [default: NULL]. All outputs will be stratified by this variable. This could be an exposure in the usual sense, or it could (mis)used to show different types of events (as long as the censoring structure is the same)

        --subgroups=SUBGROUP_VARNAMES
            Subgroup variable name or list of variable names [default: NULL]. If a subgroup variable is used, analyses will be stratified as exposure * ( subgroup1, subgroup2, ...). If NULL, no stratification will occur.

        --origin_date=ORIGIN_VARNAME
            Time-origin variable name in the input dataset [default: NULL]. Should refer to a date variable, or a character of the form YYYY-MM-DD.

        --event_date=EVENT_VARNAME
            Event variable name in the input dataset [default: NULL]. Should refer to a date variable, or a character of the form YYYY-MM-DD.

        --censor_date=CENSOR_VARNAME
            Censor variable name in the input dataset [default: NULL]. Should refer to a date variable, or a character of the form YYYY-MM-DD.

        --min_count=MIN_COUNT
            The minimum permissable event and censor counts for each 'step' in the KM curve [default: 6]. This ensures that at least `min_count` events occur at each event time.

        --method=METHOD
            The interpolation method after rounding [default: constant]. The 'constant' method leaves the event times unchanged after rounding, making the KM curve have bigger, fewer steps. The 'linear' method linearly interpolates between rounded events times (then rounds to the nearest day), so that the steps appear more natural.

        --max_fup=MAX_FUP
            The maximum follow-up time after the origin date. If event variables are dates, then this will be days. [default: Inf]. 

        --smooth=TRUE/FALSE
            Should Kaplan-Meier estimates be smoothed on the log cumulative hazard scale (TRUE) or not (FALSE) [default: FALSE]. 

        --smooth_df=SMOOTH_DF
            Degrees of freedom to use for the smoother [default: TRUE]. Unused if smooth=FALSE.

        --concise=TRUE/FALSE
            Should the outputted table only report core variables (defined here as exposure, subgroups, time, number at risk, cumulative number of events, cumulative incidence, and confidence limits) (TRUE) or should it report everything (FALSE)? [default: TRUE].

        --plot=TRUE/FALSE
            Should Kaplan-Meier plots be created in the output folder? [default: TRUE]. These are fairly basic plots for sense-checking purposes.

        -h, --help
            Show this help message and exit

For a more complete example of the arguments this action uses, see [this
repo’s `project.yaml` file](./project.yaml).

The main script for this action is [`analyis/km.R`](./analysis/km.R).

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for
electronic health records research in the NHS, with a focus on public
accountability and research quality.

Read more at [OpenSAFELY.org](https://opensafely.org).

# Licences

As standard, research projects have a MIT license.
