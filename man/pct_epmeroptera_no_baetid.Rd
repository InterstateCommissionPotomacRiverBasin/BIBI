% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composition_metrics.R
\name{pct_epmeroptera_no_baetid}
\alias{pct_epmeroptera_no_baetid}
\title{Percentage of Ephemeropteran Individuals Minus Baetid}
\usage{
pct_epmeroptera_no_baetid(order.wide, family.wide)
}
\arguments{
\item{order.wide}{= Taxonomic counts aggregated at the order level
classification in a wide data format. Use the wide function to
prepare the data.}
}
\value{
The percentage of individuals identified as Ephemeropterans
(Order: Ephemeroptera) minus the percentage of baetid (Family: Baetidae)
 in the sample.  This metric typically decreases with degradation.
}
\description{
Percentage of Ephemeropteran Individuals Minus Baetid
}
