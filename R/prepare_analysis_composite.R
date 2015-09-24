#' Prepare the analysis files for the composite indices
#' @inheritParams prepare_analysis_dataset
#' @inheritParams prepare_dataset
#' @param dataset a dataframe with the metadata of the analysis
#' @export
#' @importFrom n2kanalysis n2k_composite
#' @importFrom assertthat assert_that is.count
#' @importFrom dplyr %>% add_rownames mutate_ select_ arrange_ filter_
prepare_analysis_composite <- function(dataset, raw.connection, analysis.path){
  analysis.path <- check_path(paste0(analysis.path, "/"), type = "directory")
  check_dataframe_variable(
    df = dataset,
    variable = c(
      "LocationGroupID", "SpeciesGroupID", "Covariate", "FileFingerprint",
      "StatusFingerprint", "Seed"
    ),
    name = "dataset"
  )

  metadata <- read_delim_git("metadata.txt", connection = raw.connection)
  scheme.id <- metadata$Value[metadata$Key == "SchemeID"]
  assert_that(is.count(scheme.id))
  first.year <- metadata$Value[metadata$Key == "FirstImportedYear"]
  assert_that(is.count(first.year))
  last.year <- metadata$Value[metadata$Key == "LastImportedYear"]
  assert_that(is.count(last.year))
  dataset$ParentAnalysis <- dataset$FileFingerprint
  dataset$ParentStatusFingerprint <- dataset$StatusFingerprint
  dataset$ParentStatus <- dataset$Status

  extractor <- switch(
    dataset$Covariate[1],
    "fCycle" = function(model){
      if (!require(dplyr)) {
        stop("'dplyr' is required")
      }
      if (!require(INLA)) {
        stop("'INLA' is required")
      }
      parameter <- model$summary.lincomb.derived %>%
        add_rownames("Value") %>%
        mutate_(Variance = ~sd ^ 2) %>%
        select_(~Value, Estimate = ~mean, ~Variance) %>%
        arrange_(~Value)
      reference <- parameter %>%
        head(1)
      parameter %>%
        mutate_(Estimate = ~Estimate - reference$Estimate)
    },
    "cYear" = function(model){
      if (!require(dplyr)) {
        stop("'dplyr' is required")
      }
      if (!require(INLA)) {
        stop("'INLA' is required")
      }
      parameter <- model$summary.lincomb.derived %>%
        add_rownames("Value") %>%
        mutate_(Variance = ~sd ^ 2) %>%
        select_(~Value, Estimate = ~mean, ~Variance)
      reference <- parameter %>%
        filter_(~Value != "Trend") %>%
        mutate_(Year = ~as.numeric(Value)) %>%
        arrange_(~Year) %>%
        head(1)
      parameter %>%
        mutate_(
          Estimate = ~Estimate - ifelse(Value == "Trend", 0, reference$Estimate)
        )
    },
    "fYear" = function(model){
      if (!require(dplyr)) {
        stop("'dplyr' is required")
      }
      if (!require(INLA)) {
        stop("'INLA' is required")
      }
      parameter <- model$summary.lincomb.derived %>%
        add_rownames("Value") %>%
        mutate_(Variance = ~sd ^ 2) %>%
        select_(~Value, Estimate = ~mean, ~Variance)
      reference <- parameter %>%
        mutate_(Year = ~as.numeric(Value)) %>%
        arrange_(~Year) %>%
        head(1)
      parameter %>%
        mutate_(Estimate = ~Estimate - reference$Estimate)
    }
  )

  analysis <- n2k_composite(
    parent.status = dataset[
      ,
      c("ParentAnalysis", "ParentStatusFingerprint", "ParentStatus")
    ],
    seed = min(dataset$Seed),
    scheme.id = scheme.id,
    species.group.id = dataset$SpeciesGroupID[1],
    location.group.id = dataset$LocationGroupID[1],
    model.type = paste("composite index:", dataset$Covariate[1]),
    formula = paste("~", dataset$Covariate[1]),
    first.imported.year = first.year,
    last.imported.year = last.year,
    extractor = extractor,
    analysis.date = max(dataset$AnalysisDate)
  )
  file.fingerprint <- get_file_fingerprint(analysis)
  filename <- paste0(analysis.path, file.fingerprint, ".rda")
  if (!file.exists(filename)) {
    save(analysis, file = filename)
  }
}
