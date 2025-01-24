
from ehrql import (
    case,
    create_dataset,
    days,
    when,
    maximum_of
)

from ehrql.tables.tpp import (
  patients,
  practice_registrations,
  vaccinations,
  ons_deaths
)


# initialise dataset

dataset = create_dataset()

dataset.configure_dummy_data(population_size=1000)


# get first covid vaccination for each patient

first_covid_vaccination = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .where(vaccinations.date.is_on_or_after("2020-03-01"))
  .sort_by(vaccinations.date)
  .first_for_patient()
)

first_vax_date = first_covid_vaccination.date


# get second covid vaccination for each patient

second_covid_vaccination = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .where(vaccinations.date.is_after(first_vax_date))
  .sort_by(vaccinations.date)
  .first_for_patient()
)

second_vax_date = second_covid_vaccination.date


# check whether patient was alive and registered on first vax date

registered_patients = practice_registrations.for_patient_on(first_vax_date)

alive = (ons_deaths.date>=first_vax_date) | ons_deaths.date.is_null()


# define dataset poppulation

dataset.define_population(
  registered_patients.exists_for_patient()
  & alive
  & first_covid_vaccination.exists_for_patient()
  & (patients.age_on(first_vax_date) >=16)
)


# grouping variables

dataset.sex = patients.sex

dataset.age = patients.age_on(first_vax_date)

dataset.age_group = case(
  when(dataset.age < 50).then("aged <50"),
  when(dataset.age < 70).then("aged 50-69"),
  when(dataset.age >= 70).then("aged 70+"),
  otherwise="unknown",
)


# start of follow up variable

dataset.first_vax_date = first_vax_date


# outcome event variable

dataset.second_vax_date = second_vax_date


# censoring event variables

dataset.deregistration_date = registered_patients.end_date

dataset.death_date = ons_deaths.date

dataset.censoring_date = minimum_of(dataset.deregistration_date, dataset.death_date)
