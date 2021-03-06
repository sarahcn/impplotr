% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_functions.R
\name{maskedSummary}
\alias{maskedSummary}
\title{Summarize masked variant metrics}
\usage{
maskedSummary(maf_thresh = 0.05, mets_all, out_dir, start_chr, end_chr)
}
\arguments{
\item{maf_thresh}{Minor allele frequency threshold for grouping masked variants}

\item{mets_all}{Dataframe of imputation metrics created by \code{\link{readMetrics}}}

\item{out_dir}{Output directory}

\item{start_chr}{Start chromosome (1:23)}

\item{end_chr}{End chromosome (1:23)}
}
\value{
Writes a csv file of masked metric summaries
}
\description{
Summarize masked variant metrics
}
