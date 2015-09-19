#' Select the relevant observation of a species
#'
#' Relevant locations have at least two observations in different years. Relevant periods have average numbers of at least 5% of the most important period.
#' @param observation.species The output from \code{\link{read_observation_species}}
#' @inheritParams calculate_weight
#' @importFrom n2kanalysis select_factor_count_strictly_positive select_factor_threshold
#' @importFrom n2khelper check_dataframe_variable
#' @export
select_relevant <- function(observation, observation.species){
  check_dataframe_variable(
    df = observation,
    variable = c("ObservationID", "LocationID", "Cycle", "Year", "Period", "Stratum", "Weight"),
    name = "observation"
  )
  check_dataframe_variable(
    df = observation.species,
    variable = c("ObservationID", "Count"),
    name = "observation.species"
  )

  observation.species <- merge(
    observation.species,
    observation[, c("ObservationID", "LocationID", "Cycle", "Year", "Period", "Stratum", "Weight")],
    all.y = TRUE
  )
  observation.species$Count[is.na(observation.species$Count)] <- 0
  sampled <- unique(observation.species[, c("Stratum", "Year", "LocationID")])
  sampled.year <- as.data.frame(table(Stratum = sampled$Stratum, Year = sampled$Year))
  sampled <- unique(observation.species[, c("Stratum", "Cycle", "LocationID")])
  sampled.cycle <- as.data.frame(table(Stratum = sampled$Stratum, Cycle = sampled$Cycle))
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

  sampled <- unique(observation.species[, c("Stratum", "Year", "LocationID")])
  observed.year <- as.data.frame(table(Stratum = sampled$Stratum, Year = sampled$Year))
  sampled <- unique(observation.species[, c("Stratum", "Cycle", "LocationID")])
  observed.cycle <- as.data.frame(table(Stratum = sampled$Stratum, Cycle = sampled$Cycle))
  selected.strata <- unique(observation.species[, c("Stratum", "Weight")])
  weight.year <- merge(sampled.year, observed.year, by = c("Stratum", "Year"))
  weight.year <- merge(weight.year, selected.strata)
  weight.year$WeightYear <- weight.year$Weight * weight.year$Freq.y / weight.year$Freq.x
  total <- aggregate(WeightYear ~ Year, data = weight.year, FUN = sum)
  weight.year <- merge(weight.year, total, by = "Year")
  weight.year$Weight <- weight.year$WeightYear.x / weight.year$WeightYear.y

  weight.cycle <- merge(sampled.cycle, observed.cycle, by = c("Stratum", "Cycle"))
  weight.cycle <- merge(weight.cycle, selected.strata)
  weight.cycle$WeightCycle <- weight.cycle$Weight * weight.cycle$Freq.y / weight.cycle$Freq.x
  total <- aggregate(WeightCycle ~ Cycle, data = weight.cycle, FUN = sum)
  weight.cycle <- merge(weight.cycle, total, by = "Cycle")
  weight.cycle$Weight <- weight.cycle$WeightCycle.x / weight.cycle$WeightCycle.y

  return(
    list(
      ObservationSpecies = observation.species[, c("ObservationID", "Count")],
      WeightYear = weight.year[weight.year$Weight > 0, c("Year", "Stratum", "Weight")],
      WeightCycle = weight.cycle[weight.year$Cycle > 0, c("Cycle", "Stratum", "Weight")]
    )
  )
}
