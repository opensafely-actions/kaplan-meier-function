# kaplan-meier-function

This is the code and configuration for `kaplan-meier-function`, an R reusable action for the OpenSAFELY framework.

Given a dataset containing start dates, event dates, censoring dates, exposure groups, the action returns Kaplan-Meier estimates for these characteristics. 
Additional arguments provide control over estimation within sub-groups, the level of rounding, smoothing, maximum follow-up duration, and other elements.

There is EXPERIMENTAL support for statistical disclosure control in this repo, with further assurance and methodological explanations pending. 
We recommend comparing rounded KM estimates to the unrounded KM estimates (prior to release) to ensure that the rounded curves are a sensible approximation

Further details about this reusable action to follow.


## Usage

The arguments/options to the action are specified using the flags style
(i.e., `--argname=argvalue`), the arguments and how to use them can be
viewed in the [`km.R`](./analysis/km.R) script.

This action can be specified in the `project.yaml` with its options at
their default values as follows, where you should replace `[version]`
with the latest tag from
[here](https://github.com/opensafely-actions/kaplan-meier-function/tags),
e.g., `v0.0.1`. Note that no space is allowed between
`kaplan-meier-function:` and `[version]`.

``` yaml

my_kaplan_meier_function_action:
  run: kaplan-meier-function:v0.0.1
    --df_input=cohort.csv
    --df_output=output/km_estimates.csv
    ...
  outputs:
    highly_sensitive:
      output: output/km_estimates.csv

...
```

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic
health records research in the NHS, with a focus on public accountability and
research quality.

Read more at [OpenSAFELY.org](https://opensafely.org).

# Licences

As standard, research projects have a MIT license.

