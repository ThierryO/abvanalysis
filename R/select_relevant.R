#' Select the relevant observation of a species
#'
#' Relevant locations have at least two observations in different years. Relevant periods have average numbers of at least 5% of the most important period.
#' @param observation.species The output from \code{\link{read_observation_species}}
#' @inheritParams calculate_weight
#' @importFrom n2kanalysis select_factor_count_strictly_positive select_factor_threshold
#' @importFrom n2khelper check_dataframe_variable
#' @importFrom dplyr %>% select_ group_by_ mutate_ inner_join summarise_ ungroup transmute_
#' @export
select_relevant <- function(observation, observation.species){
  check_dataframe_variable(
    df = observation,
    variable = c("ObservationID", "LocationID", "Stratum", "Year", "Period", "Weight"),
    name = "observation"
  )
  check_dataframe_variable(
    df = observation.species,
    variable = c("ObservationID", "Count"),
    name = "observation.species"
  )

  observation.species <- observation %>%
    select_(
      ~ObservationID, ~LocationID, ~Stratum, ~Year, ~Period, ~Weight
    ) %>%
    group_by_(~Stratum) %>%
    mutate_(
      Sampled = ~ n_distinct(LocationID)
    ) %>%
    inner_join(observation.species, by = "ObservationID") %>%
    mutate_(
      Count = ~ifelse(is.na(Count), 0, Count)
    ) %>%
    as.data.frame()
  observation.species <- select_factor_count_strictly_positive(
    observation = observation.species,
    variable = "Period",
    threshold = 0.15,
    relative = TRUE
  )
  if (nrow(observation.species) == 0) {
    return(NULL)
  }
  observation.species <- select_factor_count_strictly_positive(
    observation = observation.species,
    variable = c("LocationID", "Year"),
    threshold = 2
  )
  if (nrow(observation.species) == 0) {
    return(NULL)
  }

  weights <- observation.species %>%
    group_by_(~Stratum, ~Weight, ~Sampled) %>%
    summarise_(Observed = ~n_distinct(LocationID)) %>%
    ungroup() %>%
    transmute_(
      ~Stratum,
      Weight = ~ Weight * Observed / Sampled, Weight = ~Weight / sum(Weight)
    )

  return(
    list(
      ObservationSpecies = observation.species %>%
        select_(~ObservationID, ~Count),
      Weights = weights
    )
  )
}
