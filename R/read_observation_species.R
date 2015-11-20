#' Read the observations of a given species
#' @export
#' @importFrom RODBC sqlQuery
#' @importFrom n2khelper check_single_strictly_positive_integer check_dbtable_variable
#' @param species.id The external id of the species
#' @inheritParams prepare_dataset
read_observation_species <- function(species.id, source.channel){
  species.id <- check_single_strictly_positive_integer(species.id, name = "species.id")
  check_dbtable_variable(
    table = "tblWaarneming",
    variable = c("WRNG_ID", "WRNG_UCWT_CDE", "WRNG_WGST_CDE"),
    channel = source.channel
  )
  check_dbtable_variable(
    table = "tblWaarnemingPunt",
    variable = c("WRPT_PTN", "WRPT_WRNG_ID", "WRPT_BZT"),
    channel = source.channel
  )
  check_dbtable_variable(
    table = "tblWaarnemingMeting",
    variable = c("WRME_WRNG_ID", "WRME_PNT", "WRME_ANT", "WRME_SPEC_CDE"),
    channel = source.channel
  )

  sql <- paste0("
    SELECT
      visit.ObservationID AS ObservationID,
      observed.Count AS Count
    FROM
        (
          SELECT
            WRPT_ID AS ObservationID,
            WRNG_ID AS ExternalCode,
            WRPT_PTN AS SubExternalCode
          FROM
              tblWaarneming
            INNER JOIN
              tblWaarnemingPunt
            ON
              tblWaarneming.WRNG_ID = tblWaarnemingPunt.WRPT_WRNG_ID
          WHERE
            WRNG_WRNG_ID IS NULL AND
            WRPT_BZT = 1 AND
            WRNG_UCWT_CDE IN ('ABV', 'LSABV', 'IJK') AND
            WRNG_WGST_CDE <> 'NV'
        ) AS visit
      INNER JOIN
        (
          SELECT
            WRME_WRNG_ID AS ExternalCode,
            Left([WRME_PNT],3) AS SubExternalCode,
            Sum(WRME_ANT) AS Count
          FROM
            tblWaarnemingMeting
          WHERE
            WRME_SPEC_CDE = ", species.id, "
          GROUP BY
            WRME_WRNG_ID,
            Left([WRME_PNT],3)
        ) AS observed
      ON
        visit.ExternalCode = observed.ExternalCode AND
        visit.SubExternalCode = observed.SubExternalCode
  ")

  observed.count <- sqlQuery(
    channel = source.channel,
    query = sql,
    stringsAsFactors = FALSE
  )
  return(observed.count)
}
