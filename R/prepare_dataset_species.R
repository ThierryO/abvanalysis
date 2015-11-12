#' Read the species groups and save them to git and the results database
#' @inheritParams prepare_dataset
#' @export
#' @importFrom n2khelper check_single_strictly_positive_integer get_nbn_key_multi get_nbn_name odbc_get_multi_id connect_result
#' @importFrom n2kanalysis mark_obsolete_dataset
#' @importFrom reshape2 dcast
#' @importFrom dplyr %>% mutate_
prepare_dataset_species <- function(
  source.channel, result.channel, raw.connection, scheme.id
){
  scheme.id <- check_single_strictly_positive_integer(scheme.id, name = "scheme.id")
  
  import.date <- Sys.time()
  species.list <- read_specieslist(
    source.channel = source.channel, 
    result.channel = result.channel
  )
  
  species <- get_nbn_key_multi(species.list$Species, orders = c("la", "nl", "en"))
  species$Description <- species$DutchName

  data.field.type <- data.frame(Description = "SourceSpecies")
  data.field.type <- odbc_get_multi_id(
    data = data.field.type,
    id.field = "ID",
    merge.field = "Description",
    table = "DatafieldType",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  data.field <- unique(species[, c("DatasourceID", "TableName", "PrimaryKey")])
  data.field$TypeID <- data.field.type$ID
  data.field <- odbc_get_multi_id(
    data = data.field,
    id.field = "ID",
    merge.field = c("DatasourceID", "TypeID"),
    table = "Datafield",
    channel = result.channel,
    create = TRUE,
    select = TRUE
  )
  colnames(data.field) <- gsub("^ID$", "DatafieldID", colnames(data.field))
  species <- merge(species, data.field)
  
  database.id <- odbc_get_multi_id(
    data = species[, c("ExternalCode", "DatafieldID", "Description")],
    id.field = "ID", merge.field = c("ExternalCode", "DatafieldID"), 
    table = "Sourcespecies",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )
  species <- merge(database.id, species)
  colnames(species) <- gsub("^ID$", "SourcespeciesID", colnames(species))
  species$Description <- NULL
  

  nbn.species <- dcast(
    get_nbn_name(species$NBNKey), 
    formula = NBNKey ~ Language, 
    value.var = "Name"
  )
  colnames(nbn.species) <- gsub("^en$", "EnglishName", colnames(nbn.species))
  colnames(nbn.species) <- gsub("^la$", "ScientificName", colnames(nbn.species))
  colnames(nbn.species) <- gsub("^nl$", "DutchName", colnames(nbn.species))
  colnames(nbn.species) <- gsub("^fr$", "FrenchName", colnames(nbn.species))
  
  database.id <- odbc_get_multi_id(
    data = nbn.species,
    id.field = "ID", merge.field = "NBNKey", 
    table = "Species",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )
  colnames(database.id) <- gsub("^ID$", "SpeciesID", colnames(database.id))
  nbn.species <- merge(nbn.species, database.id)
  
  species <- merge(species, database.id)
  odbc_get_multi_id(
    data = species[, c("SpeciesID", "SourcespeciesID")],
    id.field = "ID", merge.field = c("SpeciesID", "SourcespeciesID"), 
    table = "SpeciesSourcespecies",
    channel = result.channel, 
    create = TRUE,
    select = FALSE
  )
  
  
  #define each species as a species group
  nbn.species$Description <- nbn.species$DutchName
  if (anyNA(nbn.species$Description)) {
    warning(
      "Species without DutchName in NBN:\n",
      paste(
        nbn.species$NBNKey[is.na(nbn.species$Description)], 
        nbn.species$ScientificName[is.na(nbn.species$Description)], 
        collapse = "\n"
      )
    )
    nbn.species$Description[is.na(nbn.species$Description)] <- nbn.species$ScientificName[is.na(nbn.species$Description)]
  }
  if (anyNA(nbn.species$Description)) {
    stop("At least one species with neither DutchName and ScienticName. NA is not acceptable as Speciesgroup Description")
  }
  species.group <- data.frame(
    Description = nbn.species$Description,
    SchemeID = scheme.id
  )
  database.id <- odbc_get_multi_id(
    data = species.group,
    id.field = "ID", 
    merge.field = c("Description", "SchemeID"), 
    table = "SpeciesGroup",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )
  species.group <- merge(database.id, species.group)
  colnames(species.group) <- gsub("^ID", "SpeciesGroupID", colnames(species.group))
  
  species.group.species <- merge(
    species.group[, c("Description", "SpeciesGroupID")],
    nbn.species[, c("Description", "SpeciesID")]
  )
  species.group.species$Description <- NULL
  odbc_get_multi_id(
    data = species.group.species,
    id.field = "ID", 
    merge.field = c("SpeciesGroupID", "SpeciesID"), 
    table = "SpeciesGroupSpecies",
    channel = result.channel, 
    create = TRUE,
    select = FALSE
  )
  
  output <- merge(
    species.group.species,
    species[, c("SpeciesID", "ExternalCode", "DatasourceID", "DatafieldID")]
  )
  
  # Define other speciesgroups
  species.group <- data.frame(
    Description = unique(species.list$Speciesgroup$Description),
    SchemeID = scheme.id
  )
  database.id <- odbc_get_multi_id(
    data = species.group,
    id.field = "ID", 
    merge.field = c("Description", "SchemeID"), 
    table = "SpeciesGroup",
    channel = result.channel, 
    create = TRUE,
    select = TRUE
  )
  species.group <- merge(database.id, species.group)
  colnames(species.group) <- gsub("^ID", "SpeciesGroupID", colnames(species.group))
  
  species.group.species <- merge(
    species.group[, c("Description", "SpeciesGroupID")],
    species.list$Speciesgroup
  )
  species.group.species$Description <- NULL
  species.group.species <- merge(
    species.group.species,
    species[, c("ExternalCode", "SpeciesID")]
  )
  species.group.species$ExternalCode <- NULL
  odbc_get_multi_id(
    data = species.group.species,
    id.field = "ID", 
    merge.field = c("SpeciesGroupID", "SpeciesID"), 
    table = "SpeciesGroupSpecies",
    channel = result.channel, 
    create = TRUE,
    select = FALSE
  )
  
  species.group.species <- species.group.species[
    order(species.group.species$SpeciesGroupID, species.group.species$SpeciesID),
    c("SpeciesGroupID", "SpeciesID")
  ]
  
  output$Seed <- sample(.Machine$integer.max, nrow(output))
  
  species.sha <- write_delim_git(
    x = output[
      order(output$SpeciesGroupID, output$SpeciesID), 
      c("SpeciesGroupID", "SpeciesID", "Seed")
    ], 
    file = "species.txt", 
    connection = raw.connection
  )
  species.group.sha <- write_delim_git(
    x = species.group.species, file = "speciesgroup.txt", connection = raw.connection
  )
  dataset <- data.frame(
    FileName = c("species.txt", "speciesgroup.txt"),
    PathName = raw.connection@LocalPath,
    Fingerprint = c(species.sha, species.group.sha),
    ImportDate = import.date,
    Obsolete = FALSE
  )
  odbc_get_multi_id(
    data = dataset,
    id.field = "ID", 
    merge.field = c("FileName", "PathName", "Fingerprint"), 
    table = "Dataset", 
    channel = result.channel, 
    create = TRUE,
    select = FALSE
  )  
  mark_obsolete_dataset(channel = result.channel)
  
  return(output)
}
