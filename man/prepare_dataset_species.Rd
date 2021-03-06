% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/prepare_dataset_species.R
\name{prepare_dataset_species}
\alias{prepare_dataset_species}
\title{Read the species groups and save them to git and the results database}
\usage{
prepare_dataset_species(source.channel, result.channel, raw.connection,
  scheme.id)
}
\arguments{
\item{source.channel}{An open ODBC connection to the source database}

\item{result.channel}{An open RODBC connection to the results database}

\item{raw.connection}{a git-connection object to write the output to}

\item{scheme.id}{The id of the scheme}
}
\description{
Read the species groups and save them to git and the results database
}

