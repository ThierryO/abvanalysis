#' Calculate the matrix for linear combinations of strata
#' @inheritParams get_nonlinear_lincomb
#' @export
#' @importFrom n2khelper check_dataframe_variable
#' @importFrom dplyr %>% group_by_ summarise_ inner_join mutate_ select_
get_linear_lincomb <- function(
  dataset,
  time.var,
  stratum.var = "fStratum",
  formula
){
  check_dataframe_variable(
    df = dataset,
    variable = c(time.var, stratum.var, "Weight", "fPeriod")
  )

  available.weight <- dataset %>%
    group_by_(stratum.var) %>%
    summarise_(Weight = ~mean(Weight)) %>%
    mutate_(
      Weight = ~Weight / sum(Weight),
      fPeriod = ~sort(unique(dataset$fPeriod))[1]
    )
  trend.weight <- available.weight %>%
    merge(
      matrix(1, dimnames = list(1, time.var))
    )
  available.weight <- available.weight %>%
    merge(
      unique(dataset[, time.var, drop = FALSE])
    )
  mm <- available.weight %>%
    model.matrix(object = formula) * available.weight$Weight
  old.names <- colnames(mm)
  mm <- data.frame(available.weight %>% select_(ID = time.var), mm)
  weights <- aggregate(. ~ ID, mm, FUN = sum)
  colnames(weights)[-1] <- old.names
  mm.trend <- trend.weight %>%
    model.matrix(object = formula) * trend.weight$Weight
  mm.trend[, -grep(":cYear", colnames(mm.trend))] <- 0
  mm.trend <- data.frame(ID = "Trend", mm.trend)
  weights.trend <- aggregate(. ~ ID, mm.trend, FUN = sum)
  colnames(weights.trend)[-1] <- old.names
  weights <- rbind(weights, weights.trend)
  rownames(weights) <- weights$ID
  as.matrix(weights[, -1])
}
