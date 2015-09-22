#' Calculate the matrix for linear combinations of strata
#' @param dataset the raw dataset
#' @param time.var the name of the time variable
#' @param stratum.var the name of the stratum variable
#' @param formula the formula for the model.matrix
#' @param weight the name of the weight variable
#' @export
#' @importFrom n2khelper check_dataframe_variable
#' @importFrom dplyr %>% select_ distinct group_by_ mutate_
#' @importFrom assertthat assert_that is.string
get_nonlinear_lincomb <- function(
  dataset,
  time.var,
  stratum.var = "fStratum",
  formula,
  weight
){
  assert_that(is.string(time.var))
  assert_that(is.string(stratum.var))
  assert_that(is.string(weight))
  check_dataframe_variable(
    df = dataset,
    variable = c(time.var, stratum.var, weight, "fPeriod")
  )
  available.weight <- dataset %>%
    select_(time.var, stratum.var, weight) %>%
    distinct() %>%
    group_by_(time.var) %>%
    mutate_(
      Weight = weight,
      Weight = ~Weight / sum(Weight),
      fPeriod = ~sort(unique(dataset$fPeriod))[1]
    )
  mm <- available.weight %>%
    model.matrix(object = formula) * available.weight$Weight
  old.names <- colnames(mm)
  mm <- data.frame(available.weight %>% select_(ID = time.var), mm)
  weights <- aggregate(. ~ ID, mm, FUN = sum)
  colnames(weights)[-1] <- old.names
  rownames(weights) <- weights$ID
  as.matrix(weights[, -1])
}
