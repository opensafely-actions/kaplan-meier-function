from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv  # NOQA



start_date = "2020-03-01"
end_date = "2020-12-01"

study = StudyDefinition(
    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 1,
        "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
        "float": {"distribution": "normal", "mean": 25, "stddev": 5},
    },
    
    population = patients.satisfying(
      f"""
        registered
        AND NOT dead
      """,
      registered=patients.registered_as_of(
        "start_date",
      ),
      dead=patients.died_from_any_cause(
        on_or_before="start_date",
        returning="binary_flag",
      ),
    ),
   
    start_date = patients.fixed_value(start_date),
  
    end_date = patients.fixed_value(end_date),
    
    age=patients.age_as_of( 
      "start_date",
      return_expectations={
        "int": {"distribution": "normal", "mean": 60, "stddev": 10},
        "incidence": 1
      }
    ),
    
    age_group = patients.categorised_as(
      {
        "Unknown": "DEFAULT",
        "under 65": "age<65",
        "65+": "age >= 65",
      },
      return_expectations={
        "rate": "universal",
        "category": {"ratios": {"Unknown": 0.02, "under 65": 0.6, "65+": 0.38}},
      },
    ),
    
    sex=patients.sex(
      return_expectations={
        "rate": "universal",
        "category": {"ratios": {"M": 0.49, "F": 0.51}},
        "incidence": 1,
      }
    ),
  
    region=patients.registered_practice_as_of(
      "start_date",
      returning="nuts1_region_name",
      return_expectations={
      "rate": "universal",
        "category": {
          "ratios": {
            "North East": 0.1,
            "North West": 0.1,
            "Yorkshire and The Humber": 0.2,
            "East Midlands": 0.1,
            "West Midlands": 0.1,
            "East": 0.1,
            "London": 0.1,
            "South East": 0.1,
            "South West": 0.1
            #"" : 0.01
          },
        },
      },
    ),
    
    previous_covid_test=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="any",
      on_or_before="start_date",
      returning="binary_flag",
      restrict_to_earliest_specimen_date=False,
      return_expectations={
        "incidence": 0.5,
      }
    ),
  
    death_date=patients.died_from_any_cause(
      returning="date_of_death",
      date_format="YYYY-MM-DD",
      return_expectations={
        "incidence": 0.02,
      }
    ),
    
    dereg_date=patients.date_deregistered_from_all_supported_practices(
      on_or_after="start_date",
      date_format="YYYY-MM-DD",
      return_expectations={
        "incidence": 0.02,
      }
    ),
)
