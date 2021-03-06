% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insert_table.R
\name{build_insert_query}
\alias{build_insert_query}
\title{Build query to insert R data frame into SQL table}
\usage{
build_insert_query(values, table, .con)
}
\arguments{
\item{values}{`data.frame` of values to write to SQL database}

\item{table}{Name of target SQL table, as character}

\item{.con}{[\code{DBIConnection}]:A DBI connection object obtained from \code{DBI::dbConnect()}.}
}
\description{
Build query to insert R data frame into SQL table
}
