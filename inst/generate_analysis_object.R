library(abvanalysis)
raw.connection <- n2khelper::git_connection(
  repo.path = "~/n2k/https/rawdata", #nolint
  local.path = "abv",
  username = username,
  password = password,
  commit.user = "abvanalysis",
  commit.email = "bmk@inbo.be"
)
prepare_analysis(
  raw.connection = raw.connection,
  analysis.path = "~/analysis/abv" #nolint
)
