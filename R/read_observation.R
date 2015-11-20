#' Return the set of all visited locations
#' @export
#' @inheritParams prepare_dataset
#' @importFrom RODBC sqlQuery
#' @importFrom dplyr %>% mutate_ filter_ group_by_ summarise_ semi_join arrange_ select_
#' @importFrom lubridate floor_date year
#' @importFrom n2khelper cut_date
read_observation <- function(source.channel, result.channel){
  check_dbtable_variable(
    table = "tblWaarneming",
    variable = c(
      "WRNG_ID", "WRNG_DTE_BGN", "WRNG_UTM1_CDE", "WRNG_UCWT_CDE",
      "WRNG_WGST_CDE"
    ),
    channel = source.channel
  )
  check_dbtable_variable(
    table = "tblWaarnemingPunt",
    variable = c("WRPT_ID", "WRPT_PTN", "WRPT_WRNG_ID", "WRPT_BZT"),
    channel = source.channel
  )

  latest.year <- as.integer(format(Sys.time(), "%Y")) - 1
  if (Sys.time() < as.POSIXct(format(Sys.time(), "%Y-03-01"))) {
    latest.year <- latest.year - 1
  }

  sql <- paste0("
    SELECT
      WRPT_ID AS ObservationID,
      WRNG_DTE_BGN AS Timestamp,
      WRNG_UTM1_CDE AS ExternalCode,
      WRPT_PTN AS SubExternalCode,
      WRNG_USR_CRE AS Username
    FROM
        tblWaarneming
      INNER JOIN
        tblWaarnemingPunt
      ON
        tblWaarneming.WRNG_ID = tblWaarnemingPunt.WRPT_WRNG_ID
    WHERE
      WRNG_WRNG_ID IS NULL AND
      WRNG_DTE_BGN < '", latest.year, "-07-16' AND
      WRPT_BZT = 1 AND
      WRNG_UCWT_CDE IN ('ABV', 'LSABV', 'IJK') AND
      WRNG_WGST_CDE <> 'NV'
    ORDER BY
      WRNG_UTM1_CDE, WRNG_DTE_BGN
  ")
  observation <- sqlQuery(
    channel = source.channel,
    query = sql,
    stringsAsFactors = FALSE
  ) %>%
    mutate_(
      Date = ~floor_date(Timestamp, unit = "day"),
      Period = ~cut_date(
        x = Date,
        dm = c("1-3", "16-4", "1-6", "16-7"),
        include.last = FALSE
      ),
      Year = ~year(Date),
      DatasourceID = ~datasource_id(result.channel = result.channel)
    ) %>%
    filter_(~!is.na(Period))
  repeat.visit <- observation %>%
    group_by_(~ExternalCode) %>%
    summarise_(nYear = ~n_distinct(Year)) %>%
    filter_(~nYear > 1)
  observation <- observation %>%
    semi_join(repeat.visit, by = "ExternalCode") %>%
    arrange_(~ObservationID, ~SubExternalCode) %>%
    select_(
      ~DatasourceID, ~ObservationID, ~ExternalCode, ~SubExternalCode, ~Year,
      ~Period
    )
  levels(observation$Period) <- 1:3
  return(observation)
}
