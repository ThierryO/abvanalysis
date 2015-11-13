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
#' @importFrom dplyr rename_ mutate_ arrange_ %>% inner_join select_ distinct_ left_join group_by_ filter_ select_
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
  speciesgroup.id <- as.integer(gsub("\\.txt$", "", rawdata.file))
  message(speciesgroup.id, " ", appendLF = FALSE)
  utils::flush.console()

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
  weights <- read_delim_git(
    file = paste0("weights_", rawdata.file),
    connection = raw.connection
  )
  check_dataframe_variable(
    df = weights,
    variable = c("Stratum", "Weight"),
    name = rawdata.file
  )
  weights <- mutate_(weights, fStratum = ~factor(Stratum)) %>%
    arrange_(~fStratum)
  analysis.date <- git_recent(
    file = rawdata.file,
    connection = raw.connection
  )$Date

  species.id <- read_delim_git("species.txt", connection = raw.connection)
  seed <- species.id$Seed[species.id$SpeciesGroupID == speciesgroup.id]
  parent <- parent$Fingerprint[parent$SpeciesGroupID == speciesgroup.id]

  rawdata <- left_join(observation, species.observation, by = "ObservationID") %>%
    group_by_(~Period) %>%
    mutate_(Relevant = ~any(!is.na(Count))) %>%
    filter_(~Relevant) %>%
    group_by_(~LocationID) %>%
    mutate_(Relevant = ~any(!is.na(Count))) %>%
    filter_(~Relevant) %>%
    select_(~-Relevant) %>%
    mutate_(Count = ~ifelse(is.na(Count), 0, Count)) %>%
    as.data.frame()

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
        stratum.weight <- model.matrix(~0 + fStratum:Weight, data = weights) %>%
          colSums() %>%
          t()
        colnames(stratum.weight) <- gsub(
          ":Weight$",
          "",
          colnames(stratum.weight)
        )
        design.variable <- "fStratum"
      } else {
        design.variable <- character(0)
      }
      design <- character(0)

      dataset$fYear <- factor(dataset$Year)
      if (length(levels(dataset$fYear)) == 1) {
        if (multi.stratum) {
          trend <- "0 + fStratum"
        } else {
          trend <- "1"
        }
        trend.variable <- "fYear"
      } else {
        dataset$cYear <- dataset$Year - max(dataset$Year)
        if (multi.stratum) {
          dataset$cStratum <- as.integer(dataset$fStratum)
          design.variable <- list(
            c(design.variable, "cYear"),
            design.variable
          )
          trend <- c(
            "0 + fStratum + fYear +
            f(cYear, model = \"rw1\", replicate = as.integer(fStratum))",
            "0 + fStratum + cYear:fStratum"
          )
        } else {
          trend <- c("fYear", "cYear")
        }
        trend.variable <- c("fYear", "cYear")
        cycle.label <- seq(min(dataset$Year), max(dataset$Year), by = 3)
        if (length(cycle.label) > 1) {
          cycle.label <- paste(cycle.label, cycle.label + 2, sep = " - ")
          dataset$cCycle <- (dataset$Year - min(dataset$Year)) %/% 3
          dataset$fCycle <- factor(dataset$cCycle, labels = cycle.label)

          if (multi.stratum) {
            design.variable <- c(
              design.variable,
              list(c(design.variable[[2]], "cCycle"))
            )
            trend <- c(
              trend,
              "0 + fStratum + fCycle +
              f(cCycle, model = \"rw1\", replicate = as.integer(fStratum))"
            )
          } else {
            trend <- c(trend, "fCycle")
          }
          trend.variable <- c(trend.variable, "fCycle")
        }
      }

      dataset$fLocation <- factor(dataset$LocationID)
      if (length(levels(dataset$fLocation)) > 1) {
        design <- c(design, "f(fLocation, model = \"iid\")")
        design.variable <- lapply(design.variable, c, "fLocation")
      }

      dataset$fSubLocation <- factor(dataset$SubLocationID)
      if (length(levels(dataset$fSubLocation)) > 1) {
        design <- c(design, "f(fSubLocation, model = \"iid\")")
        design.variable <- lapply(design.variable, c, "fSubLocation")
      }

      if (multi.stratum) {
        weight.formula <- paste("~", trend)
      }

      dataset$fPeriod <- factor(dataset$Period)
      if (length(levels(dataset$fPeriod)) > 1) {
        design <- c(design, "fPeriod")
        design.variable <- lapply(design.variable, c, "fPeriod")
        if (multi.stratum) {
          weight.formula <- paste(weight.formula, "fPeriod", sep = "+")
        }
      }
      design <- paste(design, collapse = " + ")
      covariates <- paste(trend, design, sep = " + ")

      fingerprint <- do.call(rbind, lapply(seq_along(covariates), function(i){
        if (multi.stratum) {
          if (trend.variable[i] == "cYear") {
            cyear <- c(unique(dataset$cYear), 1)
            time.fixed <- outer(cyear, stratum.weight[1, ])
            rownames(time.fixed) <- c(head(cyear, -1), "Trend")
            colnames(time.fixed) <- paste0(colnames(time.fixed), ":cYear")
            lc <- cbind(
              rbind(
                stratum.weight[rep(1, nrow(time.fixed) - 1), ],
                0
              ),
              time.fixed
            )
          } else {
            time.fixed <- dataset %>%
              select_(trend.variable[i]) %>%
              distinct_() %>%
              model.matrix(object = as.formula(paste("~", trend.variable[i])))
            rownames(time.fixed) <- levels(dataset[, trend.variable[i]])
            random <- matrix(
              weights$Weight,
              ncol = nrow(time.fixed) * nrow(weights),
              nrow = nrow(time.fixed)
            ) * do.call(cbind, rep(list(diag(nrow(time.fixed))), nrow(weights)))
            colnames(random) <- rep(
              gsub("^f", "c", trend.variable[i]),
              ncol(random)
            )
            lc <- cbind(
              stratum.weight[rep(1, nrow(time.fixed)), ],
              time.fixed[, -1],
              random
            )
          }
        } else {
          lc <- NULL
        }

        data <- dataset %>%
          select_(
            ~ObservationID, ~DatasourceID, ~Count,
            trend.variable[i], .dots = design.variable[[i]]
          ) %>%
          arrange_(.dots = design.variable[[i]], trend.variable[i])
        model.type <- paste(
          "inla nbinomial:",
          trend.variable[i], "+ Stratum + Period + Location + SubLocation"
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
