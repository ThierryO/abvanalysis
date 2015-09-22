#' Create the analysis dataset based on the available raw data
#'
#' This analysis fits an unweighted model but adds the stratum effect. The
#'    indices per year of the different strata are combined with a linear
#'    combination into a single index per year.
#' @return A data.frame with the species id number of rows in the analysis
#'    dataset, number of precenses in the analysis datset and SHA-1 of the
#'    analysis dataset or NULL if not enough data.
#' @importFrom n2khelper check_path check_dataframe_variable git_recent
#' @importFrom assertthat assert_that is.count
#' @importFrom plyr ddply
#' @importFrom dplyr rename_
#' @importClassesFrom n2kanalysis n2kInlaNbinomial
#' @importFrom n2kanalysis n2k_inla_nbinomial
#' @importMethodsFrom n2kanalysis get_file_fingerprint get_status_fingerprint get_seed
#' @export
#' @param rawdata.file The file with the counts per visit
#' @param observation the dataframe with the visits and location group
#' @param analysis.path the path to store the rda files for the analysis
#' @param min.observation The minimum number of positive observations
#'    (Count > 0)
#' @inheritParams prepare_dataset
prepare_analysis_dataset <- function(
  rawdata.file,
  observation,
  raw.connection,
  analysis.path,
  min.observation
){
  assert_that(is.count(min.observation))
  metadata <- read_delim_git("metadata.txt", connection = raw.connection)
  scheme.id <- metadata$Value[metadata$Key == "SchemeID"]
  assert_that(is.count(scheme.id))
  first.year <- metadata$Value[metadata$Key == "FirstImportedYear"]
  assert_that(is.count(first.year))
  last.year <- metadata$Value[metadata$Key == "LastImportedYear"]
  assert_that(is.count(last.year))

  parent <- read_delim_git(file = "parent.txt", connection = raw.connection)
  check_dataframe_variable(
    df = parent,
    variable = c("SpeciesGroupID", "Fingerprint"),
    name = "parent.txt"
  )

  analysis.path <- check_path(paste0(analysis.path, "/"), type = "directory")
  check_dataframe_variable(
    df = observation,
    variable = c(
      "ObservationID", "DatasourceID", "LocationID", "SubLocationID", "Year",
      "Period", "Stratum", "LocationGroupID"
    ),
    name = "observation"
  )

  species.observation <- read_delim_git(
    file = rawdata.file,
    connection = raw.connection
  )
  check_dataframe_variable(
    df = species.observation,
    variable = c("ObservationID", "Count"),
    name = rawdata.file
  )
  weight.year <- read_delim_git(
    file = paste0("weight_year_", rawdata.file),
    connection = raw.connection
  )
  check_dataframe_variable(
    df = weight.year,
    variable = c("Year", "Stratum", "Weight"),
    name = rawdata.file
  )
  weight.year  <- rename_(weight.year, WeightYear = ~Weight)
  weight.cycle <- read_delim_git(
    file = paste0("weight_cycle_", rawdata.file),
    connection = raw.connection
  )
  check_dataframe_variable(
    df = weight.cycle,
    variable = c("Cycle", "Stratum", "Weight"),
    name = rawdata.file
  )
  weight.cycle  <- rename_(weight.cycle, WeightCycle = ~Weight)
  analysis.date <- git_recent(
    file = rawdata.file,
    connection = raw.connection
  )$Date
  speciesgroup.id <- as.integer(gsub("\\.txt$", "", rawdata.file))
  species.id <- read_delim_git("species.txt", connection = raw.connection)
  seed <- species.id$Seed[species.id$SpeciesGroupID == speciesgroup.id]
  parent <- parent$Fingerprint[parent$SpeciesGroupID == speciesgroup.id]

  rawdata <- merge(species.observation, observation)

  message(speciesgroup.id, " ", appendLF = FALSE)
  utils::flush.console()
  analysis <- ddply(
    .data = rawdata,
    .variables = "LocationGroupID",
    .fun = function(dataset){
      locationgroup.id <- dataset$LocationGroupID[1]

      if (sum(dataset$Count > 0) < min.observation) {
        analysis <- n2k_inla_nbinomial(
          scheme.id = scheme.id,
          species.group.id = speciesgroup.id,
          location.group.id = locationgroup.id,
          analysis.date = analysis.date,
          model.type =
            "inla nbinomial: fYear + Stratum + Period + Location + SubLocation",
          formula = "Count ~ 1",
          first.imported.year = first.year,
          last.imported.year = last.year,
          data = dataset,
          status = "insufficient data",
          parent = parent,
          seed = seed
        )
        filename <- paste0(
          analysis.path,
          get_file_fingerprint(analysis),
          ".rda"
        )
        if (!file.exists(filename)) {
          save(analysis, file = filename)
        }
        return(NULL)
      }

      dataset$fStratum <- factor(dataset$Stratum)
      multi.stratum <- length(levels(dataset$fStratum)) > 1
      if (multi.stratum) {
        design.variable <- "fStratum"
        dataset <- merge(dataset, weight.year, by = c("Stratum", "Year"))
      } else {
        design.variable <- character(0)
      }
      design <- character(0)

      dataset$fYear <- factor(dataset$Year)
      if (length(levels(dataset$fYear)) == 1) {
        if (multi.stratum) {
          trend <- "fYear * fStratum"
        } else {
          trend <- "fYear"
        }
        trend.variable <- "fYear"
      } else {
        dataset$cYear <- dataset$Year - max(dataset$Year)
        if (multi.stratum) {
          trend <- c("fYear * fStratum", "0 + fStratum + cYear:fStratum")
        } else {
          trend <- c("fYear", "cYear")
        }
        trend.variable <- c("fYear", "cYear")
        cycle.label <- seq(min(dataset$Year), max(dataset$Year), by = 3)
        if (length(cycle.label) > 1) {
          cycle.label <- paste(cycle.label, cycle.label + 2, sep = " - ")
          dataset$Cycle <- (dataset$Year - min(dataset$Year)) %/% 3
          dataset$fCycle <- factor(dataset$Cycle, labels = cycle.label)

          if (multi.stratum) {
            trend <- c(trend, "fCycle * fStratum")
            dataset <- merge(dataset, weight.cycle, by = c("Stratum", "Cycle"))
          } else {
            trend <- c(trend, "fCycle")
          }
          trend.variable <- c(trend.variable, "fCycle")
        }
      }

      dataset$fLocation <- factor(dataset$LocationID)
      if (length(levels(dataset$fLocation)) > 1) {
        design <- c(design, "f(fLocation, model = \"iid\")")
        design.variable <- c(design.variable, "fLocation")
      }

      dataset$fSubLocation <- factor(dataset$SubLocationID)
      if (length(levels(dataset$fSubLocation)) > 1) {
        design <- c(design, "f(fSubLocation, model = \"iid\")")
        design.variable <- c(design.variable, "fSubLocation")
      }

      if (multi.stratum) {
        weight.formula <- paste("~", trend)
      }

      dataset$fPeriod <- factor(dataset$Period)
      if (length(levels(dataset$fPeriod)) > 1) {
        design <- c(design, "fPeriod")
        design.variable <- c(design.variable, "fPeriod")
        if (multi.stratum) {
          weight.formula <- paste(weight.formula, "fPeriod", sep = "+")
        }
      }
      design <- paste(design, collapse = " + ")
      covariates <- paste(trend, design, sep = " + ")

      fingerprint <- do.call(rbind, lapply(seq_along(covariates), function(i){
        if (multi.stratum) {
          if (trend.variable[i] == "cYear") {
            lc <- get_linear_lincomb(
              dataset = dataset,
              time.var = trend.variable[i],
              formula = as.formula(weight.formula[i]),
              stratum.var = "fStratum",
              weight = "WeightYear"
            )
          } else {
            lc <- get_nonlinear_lincomb(
              dataset = dataset,
              time.var = trend.variable[i],
              formula = as.formula(weight.formula[i]),
              stratum.var = "fStratum",
              weight = ifelse(
                trend.variable[i] == "fYear",
                "WeightYear",
                "WeightCycle"
              )
            )
          }
        } else {
          lc <- NULL
        }
        data <- dataset[
          ,
          c(
            "ObservationID", "DatasourceID", "Count",
            trend.variable[i], design.variable
          )
        ]
        model.type <- paste(
          "inla nbinomial:",
          trend.variable[i], "+ Period + Location + SubLocation"
        )
        formula <- paste("Count ~", covariates[i])
        analysis <- n2k_inla_nbinomial(
          scheme.id = scheme.id,
          species.group.id = speciesgroup.id,
          location.group.id = locationgroup.id,
          analysis.date = analysis.date,
          seed = seed,
          model.type = model.type,
          formula = formula,
          first.imported.year = first.year,
          last.imported.year = last.year,
          data = data,
          parent = parent,
          lin.comb = lc
        )
        file.fingerprint <- get_file_fingerprint(analysis)
        filename <- paste0(analysis.path, file.fingerprint, ".rda")
        if (!file.exists(filename)) {
          save(analysis, file = filename)
        }
        data.frame(
          ModelType = model.type,
          Covariate = trend.variable[i],
          FileFingerprint = file.fingerprint,
          StatusFingerprint = get_status_fingerprint(analysis),
          Seed = get_seed(analysis),
          stringsAsFactors = FALSE
        )
      }))
      return(
        cbind(
          LocationGroupID = locationgroup.id,
          AnalysisDate = analysis.date,
          fingerprint
        )
      )
    }
  )
  if (nrow(analysis) > 0) {
    analysis$SpeciesGroupID <- speciesgroup.id
  }
  return(analysis)
}
