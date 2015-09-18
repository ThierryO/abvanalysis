#' Read the counts and save them to git
#' @param this.species a dataframe with the ExternalCode and the SpeciesGroupID
#' @param observation a data.frame with observations
#' @param first.year oldest year of the data
#' @param last.year most recent year of the data
#' @inheritParams prepare_dataset
#' @export
#' @importFrom n2khelper check_dataframe_variable odbc_get_id  check_single_strictly_positive_integer odbc_get_multi_id get_sha1
#' @importFrom n2kanalysis get_analysis_version
#' @importFrom RODBC sqlQuery
prepare_dataset_species_observation <- function(
  this.species, observation, result.channel, source.channel, raw.connection, scheme.id, first.year, last.year
){
  check_dataframe_variable(
    df = this.species,
    variable = c("ExternalCode", "SpeciesGroupID"),
    name = "this.species"
  )
  if(nrow(this.species) != 1){
    stop("this.species must contain exactly one row")
  }
  check_dataframe_variable(
    df = observation,
    variable = c("ObservationID", "SubLocationID"),
    name = "observation"
  )
  scheme.id <- check_single_strictly_positive_integer(scheme.id, name = "scheme.id")
  
  import.date <- Sys.time()
  observation.species <- read_observation_species(
    species.id = this.species$ExternalCode,
    source.channel = source.channel
  )
  if (nrow(observation.species) == 0) {
    observation.species <- NULL
  } else {
    observation.species <- select_relevant(observation, observation.species)
  }
  
  filename <- paste0(this.species$SpeciesGroupID, ".txt")
  
  if (is.null(observation.species)) {
    observation.species.sha <- NA
    status.id <- odbc_get_multi_id(
      data = data.frame(Description = "insufficient data"),
      id.field = "ID", 
      merge.field = c("Description"), 
      table = "AnalysisStatus",
      channel = result.channel, 
      create = TRUE,
      select = TRUE
    )$ID
  } else {
    observation.species <- observation.species[
      order(observation.species$ObservationID),
      c("ObservationID", "Count")
    ]
    observation.species.sha <- write_delim_git(
      x = observation.species, 
      file = filename, 
      connection = raw.connection
    )
    status.id <- odbc_get_multi_id(
      data = data.frame(Description = "waiting"),
      id.field = "ID", 
      merge.field = c("Description"), 
      table = "AnalysisStatus",
      channel = result.channel, 
      create = TRUE,
      select = TRUE
    )$ID
  }
  
  import.id <- odbc_get_multi_id(
    data = data.frame(Description = "import"),
    id.field = "ID", 
    merge.field = "Description", 
    table = "ModelType",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )$ID
  
  model.set <- data.frame(
    ModelTypeID = import.id,
    FirstYear = first.year,
    LastYear = last.year,
    Duration = last.year - first.year + 1
  )
  model.set.id <- odbc_get_multi_id(
    data = model.set,
    id.field = "ID", 
    merge.field = c("ModelTypeID", "FirstYear", "LastYear", "Duration"), 
    table = "ModelSet",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )$ID
  
  version <- get_analysis_version(sessionInfo())
  
  r.package <- odbc_get_multi_id(
    data = version@RPackage[, c("Description", "Version")],
    id.field = "ID",
    merge.field = c("Description", "Version"),
    table = "RPackage",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )
  r.package <- merge(r.package, version@RPackage)[, c("Fingerprint", "ID")]
  colnames(r.package)[colnames(r.package) == "Fingerprint"] <- "RPackage"
  colnames(r.package)[colnames(r.package) == "ID"] <- "RPackageID"
  analysis.package <- merge(version@AnalysisVersionRPackage, r.package)[, c("AnalysisVersion", "RPackageID")]

  analysis.version <- odbc_get_multi_id(
    data = data.frame(Description = version@AnalysisVersion$Fingerprint),
    id.field = "ID",
    merge.field = "Description",
    table = "AnalysisVersion",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )
  colnames(analysis.version)[colnames(analysis.version) == "ID"] <- "AnalysisVersionID"
  colnames(analysis.version)[colnames(analysis.version) == "Description"] <- "AnalysisVersion"
  analysis.package <- merge(analysis.version, analysis.package)[, c("AnalysisVersionID", "RPackageID")]
  odbc_get_multi_id(
    data = analysis.package,
    id.field = "ID",
    merge.field = c("AnalysisVersionID", "RPackageID"),
    table = "AnalysisVersionRPackage",
    channel = result.channel, 
    create = TRUE,
    select = FALSE
  )
  
  version.id <- analysis.version$AnalysisVersionID
  
  sql <- paste0("
    SELECT
      ID, FileName, PathName, Fingerprint
    FROM
      Dataset
    WHERE
      PathName = '", raw.connection@LocalPath, "' AND
      Filename IN ('observation.txt', 'species.txt', 'speciesgroup.txt', 'locationgrouplocation.txt') AND
      Obsolete = 0
  ")
  location.ds <- sqlQuery(channel = result.channel, query = sql, stringsAsFactors = FALSE)
  fingerprint <- get_sha1(sort(c(location.ds$Fingerprint, observation.species.sha)))
  
  analysis <- data.frame(
    ModelSetID = model.set.id,
    LocationGroupID = odbc_get_id(
      table = "LocationGroup",
      variable = c("Description", "SchemeID"),
      value = c("Vlaanderen", scheme.id),
      channel = result.channel
    ),
    SpeciesGroupID = this.species$SpeciesGroupID,
    AnalysisVersionID = version.id,
    AnalysisDate = import.date,
    StatusID = status.id,
    Fingerprint = fingerprint
  )
  analysis.id <- odbc_get_multi_id(
    data = analysis,
    id.field = "ID", 
    merge.field = c("ModelSetID", "LocationGroupID", "SpeciesGroupID", "AnalysisVersionID", "Fingerprint"),
    table = "Analysis",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )$ID
  
  if (!is.na(observation.species.sha)) {
    dataset <- data.frame(
      FileName = filename,
      PathName = raw.connection@LocalPath,
      Fingerprint = observation.species.sha,
      ImportDate = import.date,
      Obsolete = FALSE
    )
    location.ds <- rbind(
      location.ds,
      odbc_get_multi_id(
        data = dataset,
        id.field = "ID", 
        merge.field = c("FileName", "PathName", "Fingerprint"), 
        table = "Dataset", 
        channel = result.channel, 
        create = TRUE,
        select = TRUE
      )
    )
  }
  analysis.dataset <- data.frame(
    AnalysisID = analysis.id,
    DatasetID = location.ds$ID
  )
  odbc_get_multi_id(
    data = analysis.dataset,
    id.field = "ID", 
    merge.field = c("AnalysisID", "DatasetID"), 
    table = "AnalysisDataset", 
    channel = result.channel, 
    create = TRUE,
    select = FALSE
  )
  return(analysis[, c("SpeciesGroupID", "Fingerprint")])
}
