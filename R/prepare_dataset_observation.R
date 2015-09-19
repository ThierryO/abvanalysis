#' Read the observations and save them to git and the results database
#' @inheritParams prepare_dataset
#' @export
#' @importFrom assertthat assert_that is.count
#' @importFrom n2khelper read_delim_git odbc_get_multi_id connect_result
#' @importFrom dplyr %>% select_ mutate_ inner_join
#' @importFrom n2kanalysis mark_obsolete_dataset
#' @importFrom RODBC odbcClose
prepare_dataset_observation <- function(
  source.channel, result.channel, raw.connection, attribute.connection, scheme.id
){
  assert_that(is.count(scheme.id))

  import.date <- Sys.time()
  observation <- read_observation(
    source.channel = source.channel,
    result.channel = result.channel
  )
  stratum <- read_delim_git(
    file = "habitat.txt",
    connection = attribute.connection
  )
  check_dataframe_variable(
    df = stratum,
    variable = c("ExternalCode", "Stratum"),
    name = "habitat.txt"
  )
  observation <- stratum %>%
    select_(~ExternalCode, ~Stratum) %>%
    mutate_(Stratum = ~ factor(Stratum)) %>%
    group_by_(~Stratum) %>%
    mutate_(Weight = ~n_distinct(ExternalCode)) %>%
    inner_join(observation, by = "ExternalCode")
  data.field.type <- odbc_get_multi_id(
    data = data.frame(Description = "Observation"),
    id.field = "ID",
    merge.field = "Description",
    table = "DatafieldType",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  data.field <- unique(observation[, "DatasourceID", drop = FALSE])
  data.field$TableName <- "tblWaarnemingPunt"
  data.field$PrimaryKey <- "WRPT_ID"
  data.field$TypeID <- data.field.type$ID
  odbc_get_multi_id(
    data = data.field,
    id.field = "ID",
    merge.field = c("DatasourceID", "TypeID"),
    table = "Datafield",
    channel = result.channel,
    create = TRUE,
    select = FALSE
  )

  # TO DO
  #   - check for observations < 2007
  #   - check for multiple observations per sublocation - time combination
  #   - check for observations outside the 1/3 - 15/7 range/

  # store the locations
  result.datasource <-  result_datasource_id(result.channel = result.channel)
  data.field.type <- odbc_get_multi_id(
    data = data.frame(Description = "Location"),
    id.field = "ID",
    merge.field = "Description",
    table = "DatafieldType",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  data.field <- data.frame(
    DatasourceID = c(
      result.datasource,
      unique(observation$DatasourceID)
    ),
    TableName = c("Location", "tblUTM1"),
    PrimaryKey = c("ID", "UTM1_CDE"),
    TypeID = data.field.type$ID
  )
  data.field <- odbc_get_multi_id(
    data = data.field,
    id.field = "ID",
    merge.field = c("DatasourceID", "TypeID"),
    table = "Datafield",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  data.field <- data.field[data.field$DatasourceID != result.datasource, ]
  colnames(data.field) <- gsub("^ID$", "DatafieldID", colnames(data.field))

  main.location <- unique(observation[, c("ExternalCode", "DatasourceID")])
  main.location <- merge(main.location, data.field[, c("DatasourceID", "DatafieldID")])
  main.location$Description <- main.location$ExternalCode

  database.id <- odbc_get_multi_id(
    data = main.location[, c("Description", "ExternalCode", "DatafieldID")],
    id.field = "ID",
    merge.field = c("ExternalCode", "DatafieldID"),
    table = "Location",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  database.id <- merge(database.id, data.field[, c("DatasourceID", "DatafieldID")])
  observation <- merge(
    observation,
    database.id
  )
  colnames(observation) <- gsub("^ID$", "LocationID", colnames(observation))
  observation$DatafieldID <- NULL

  location.group <- data.frame(
    Description = "Vlaanderen",
    SchemeID = scheme.id
  )
  location.group.id <- odbc_get_multi_id(
    data = location.group,
    id.field = "ID",
    merge.field = c("SchemeID", "Description"),
    table = "LocationGroup",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )$ID
  location.group.location <- data.frame(
    LocationGroupID = location.group.id,
    LocationID = database.id$ID
  )
  location.group.location <- location.group.location[
    order(location.group.location$LocationGroupID, location.group.location$LocationID),
  ]
  odbc_get_multi_id(
    data = location.group.location,
    id.field = "ID",
    merge.field = c("LocationGroupID", "LocationID"),
    table = "LocationGroupLocation",
    channel = result.channel,
    create = TRUE,
    select = FALSE
  )

  data.field.type <- odbc_get_multi_id(
    data = data.frame(Description = "SubLocation"),
    id.field = "ID",
    merge.field = "Description",
    table = "DatafieldType",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  data.field <- data.frame(
    DatasourceID = c(result.datasource, unique(observation[, "DatasourceID"])),
    TableName = c("Location", "tblWaarnemingPunt"),
    PrimaryKey = c("ID", "WRPT_PTN"),
    TypeID = data.field.type$ID
  )
  data.field <- odbc_get_multi_id(
    data = data.field,
    id.field = "ID",
    merge.field = c("DatasourceID", "TypeID"),
    table = "Datafield",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  data.field <- data.field[data.field$DatasourceID != result.datasource, ]
  colnames(data.field) <- gsub("^ID$", "DatafieldID", colnames(data.field))
  sub.location <- unique(
    observation[, c("LocationID", "ExternalCode", "SubExternalCode", "DatasourceID")]
  )
  observation$ExternalCode <- NULL
  sub.location <- merge(sub.location, data.field)
  sub.location$DatasourceID <- NULL
  sub.location$TypeID <- NULL
  colnames(sub.location) <- gsub("^LocationID$", "ParentLocationID", colnames(sub.location))
  sub.location$Description <- paste(sub.location$ExternalCode, sub.location$SubExternalCode)
  sub.location$ExternalCode <- NULL
  colnames(sub.location) <- gsub(
    "^SubExternalCode$",
    "ExternalCode",
    colnames(sub.location)
  )
  database.id <- odbc_get_multi_id(
    data = sub.location,
    id.field = "ID",
    merge.field = c("ParentLocationID", "ExternalCode", "DatafieldID"),
    table = "Location",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  database.id <- merge(database.id, data.field)
  database.id$DatafieldID <- NULL
  database.id$TypeID <- NULL
  observation <- merge(
    observation,
    database.id,
    by.x = c("DatasourceID", "LocationID", "SubExternalCode"),
    by.y = c("DatasourceID", "ParentLocationID", "ExternalCode")
  )
  colnames(observation) <- gsub("^ID$", "SubLocationID", colnames(observation))

  observation <- observation[
    order(
      observation$Stratum,
      observation$LocationID,
      observation$SubLocationID,
      observation$Year,
      observation$Period
    ),
    c("DatasourceID", "ObservationID", "Stratum", "LocationID", "SubLocationID", "Year", "Period", "Weight")
  ]

  location.group.location.sha <- write_delim_git(
    x = location.group.location,
    file = "locationgrouplocation.txt",
    connection = raw.connection
  )
  observation.sha <- write_delim_git(
    x = observation[, c("DatasourceID", "ObservationID", "Stratum", "LocationID", "SubLocationID", "Year", "Period")],
    file = "observation.txt",
    connection = raw.connection
  )

  dataset <- data.frame(
    FileName = c("observation.txt",  "locationgrouplocation.txt"),
    PathName = raw.connection@LocalPath,
    Fingerprint = c(observation.sha, location.group.location.sha),
    ImportDate = import.date,
    Obsolete = FALSE
  )
  database.id <- odbc_get_multi_id(
    data = dataset,
    id.field = "ID",
    merge.field = c("FileName", "PathName", "Fingerprint"),
    table = "Dataset",
    channel = result.channel,
    create = TRUE
  )
  mark_obsolete_dataset(channel = result.channel)

  observation$Cycle <- (observation$Year - min(observation$Year)) %/% 3
  return(observation)
}
