

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

The arguments to the action are specified using the flags style, i.e.,
`--argname=argvalue`. The available arguments to this action are as
follows:

    Usage: kaplan-meier-function:[version] [options]


    Options:
        --df_input=DF_INPUT.FEATHER
            [default: NULL] character. The input dataset .feather filename. feather format is enforced to ensure date types are preserved.  Must be specified.

        --dir_output=/OUTPUT/
            [default: NULL] character. The output directory. Must be specified.

        --exposure=EXPOSURE_VARNAME
            [default: character(0)] character. The name of an exposure variable in the input dataset. All outputs will be stratified by this variable. This could be an exposure in the usual sense, or it could (mis)used to show different types of events (as long as the censoring structure is the same). If not specified, no stratification will occur.

        --subgroups=SUBGROUP_VARNAME
            character. The name of a subgroup variable or list of variable names. If a subgroup variable is used, analyses will be stratified as exposure * ( subgroup1, subgroup2, ...). If not specified, no stratification will occur.

        --origin_date=ORIGIN_VARNAME
            The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') in the input dataset that represents the start of follow-up. Must be specified.

        --event_date=EVENT_VARNAME
            The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') in the input dataset that represents the event date. Must be specified.

        --censor_date=CENSOR_VARNAME
            [default: character(0)] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') that represents the censoring date. If not specified, then no censoring occurs except at `max_fup` time.

        --min_count=MIN_COUNT
            [default: 6] integer. The minimum permissable event and censor counts for each 'step' in the KM curve. This ensures that at least `min_count` events occur at each event time.

        --method=METHOD
            [default: constant] character. The interpolation method after rounding. The 'constant' method leaves the event times unchanged after rounding, making the KM curve have bigger, fewer steps. The 'linear' method linearly interpolates between rounded events times (then rounds to the nearest day), so that the steps appear more natural.

        --max_fup=MAX_FUP
            [default: Inf] numeric. The maximum follow-up time after the origin date. If event variables are dates, then this will be days.

        --smooth=TRUE/FALSE
            [default: FALSE] logical. Should Kaplan-Meier estimates be smoothed on the log cumulative hazard scale (TRUE) or not (FALSE). 

        --smooth_df=TRUE/FALSE
            [default: TRUE]. integer. Degrees of freedom to use for the smoother. Unused if smooth=FALSE.

        --concise=TRUE/FALSE
            [default: TRUE] logical. Should the outputted table only report core variables (defined here as exposure, subgroups, time, number at risk, cumulative number of events, cumulative incidence, and confidence limits) (TRUE) or should it report everything (FALSE)?

        --plot=TRUE/FALSE
            [default: TRUE] logical. Should Kaplan-Meier plots be created in the output folder? These are fairly basic plots for sense-checking purposes.

        -h, --help
            Show this help message and exit

For example, the relevant action in the `project.yaml` file might look
something like this:

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
