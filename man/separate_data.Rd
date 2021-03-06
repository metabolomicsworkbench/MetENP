% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate_data.R
\name{separate_data}
\alias{separate_data}
\title{Separates user data into metadata and metabolomics data with metabolite measurements}
\usage{
separate_data(example, typedata, sample_col)
}
\arguments{
\item{example}{user data file}

\item{typedata}{type of data desired from the user data uploaded: metadata or metabolomics data}

\item{sample_col}{TRUE or FALSE, TRUE if sample names are in the first column, FALSE if metabolites names are in the
first column.}
}
\description{
Separates user data into metadata and metabolomics data with metabolite measurements
}
\examples{
separate_data('example.txt', "metadata",FALSE)
}
