#' Prepare all datasets and a to do list of models
#' @inheritParams prepare_analysis_dataset
#' @inheritParams prepare_dataset
#' @export
#' @importFrom n2khelper check_single_strictly_positive_integer check_path list_files_git git_sha
#' @importFrom n2kanalysis status
#' @importFrom plyr d_ply
#' @importFrom dplyr inner_join
#' @importFrom assertthat assert_that is.count
prepare_analysis <- function(
  raw.connection, analysis.path = ".", min.observation = 1
){
  assert_that(is.count(min.observation))
  path <- check_path(
    paste0(analysis.path, "/"),
    type = "directory",
    error = FALSE
  )
  if (is.logical(path)) {
    dir.create(path = analysis.path, recursive = TRUE)
    path <- check_path(paste0(analysis.path, "/"), type = "directory")
  }
  analysis.path <- path

  observation <- read_delim_git(
    file = "observation.txt",
    connection = raw.connection
  )
  check_dataframe_variable(
    df = observation,
    variable = c(
      "ObservationID", "DatasourceID", "LocationID", "SubLocationID", "Year",
      "Period", "Stratum"
    ),
    name = "observation.txt"
  )

  locationgrouplocation <- read_delim_git(
    file = "locationgrouplocation.txt",
    connection = raw.connection
  )
  check_dataframe_variable(
    df = locationgrouplocation,
    variable = c("LocationID", "LocationGroupID"),
    name = "locationgrouplocation.txt"
  )

  observation <- inner_join(
    observation,
    locationgrouplocation,
    by = "LocationID"
  )
  rm(locationgrouplocation)

  message("Prepare analysis per species")
  utils::flush.console()
  rawdata.files <- list_files_git(
    connection = raw.connection,
    pattern = "^[0-9]*\\.txt$"
  )
  if (length(rawdata.files) == 0) {
    warning("Nothing to do")
    return(invisible(NULL))
  }
  analysis <- do.call(rbind, lapply(
    rawdata.files,
    prepare_analysis_dataset,
    min.observation = min.observation,
    observation = observation,
    raw.connection = raw.connection,
    analysis.path = analysis.path
  ))
  if (nrow(analysis) == 0) {
    return(invisible(NULL))
  }
  current.status <- status(analysis.path)
  analysis <- merge(
    analysis,
    current.status[, c("FileFingerprint", "Status")]
  )

  message("\nPrepare model comparison")
  utils::flush.console()
  d_ply(
    .data = analysis,
    .variables = c("LocationGroupID", "SpeciesGroupID"),
    .fun = prepare_analysis_comparison,
    raw.connection = raw.connection,
    analysis.path = analysis.path
  )

  analysis <- analysis[analysis$Covariate %in% c("fYear", "fCycle"), ]
  species.id <- read_delim_git(
    file = "species.txt",
    connection = raw.connection
  )
  check_dataframe_variable(
    df = species.id,
    variable = c("SpeciesGroupID", "SpeciesID"),
    name = "species.txt"
  )
  analysis <- merge(analysis, species.id)
  analysis$SpeciesGroupID <- NULL

  speciesgroupid <- read_delim_git(
    file = "speciesgroup.txt",
    connection = raw.connection
  )
  check_dataframe_variable(
    df = speciesgroupid,
    variable = c("SpeciesGroupID", "SpeciesID"),
    name = "species.txt"
  )
  analysis <- merge(analysis, speciesgroupid)

  message("Prepare Composite indices")
  utils::flush.console()
  d_ply(
    .data = analysis,
    .variables = c("LocationGroupID", "SpeciesGroupID", "Covariate"),
    .fun = prepare_analysis_composite,
    raw.connection = raw.connection,
    analysis.path = analysis.path
  )

  return(invisible(NULL))
}
