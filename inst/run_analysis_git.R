library(n2kanalysis)
result.channel <- n2khelper::connect_result()
datasource.id <- abvanalysis::result_datasource_id(
  result.channel = result.channel
)
result <- get_result(
  x = "~/analysis/abv", #nolint
  datasource.id = datasource.id,
  keep.fingerprint = FALSE,
  n.cluster = parallel::detectCores() - 1
)
save(result, file = "~/analysis/output/abv.rda") #nolint

library(n2kanalysis)
result.channel <- n2khelper::connect_result()
load(file = "~/analysis/output/abv.rda") #nolint
import_result(result = result, result.channel = result.channel)
