\name{dataca20}
\alias{dataca20}
\docType{data}
\title{
Calcium Content In Soil Samples.
}
\description{
This data set contains the calcium content measured in soil samples taken from the 0-20cm layer at 178 locations within a certain study area divided in three sub-areas. The elevation at each location was also recorded.See \code{geoR} package for details.
}
\usage{data("dataca20")}
\format{
  A data frame with 178 observations on the following 3 variables.
  \describe{
    \item{\code{east}}{X Coordinate.}
    \item{\code{north}}{Y coordinate.}
    \item{\code{calcont}}{Calcium content measured in \eqn{mmol_c/dm^3}.}
    \item{\code{altitude}}{A vector with the elevation of each sampling location,in meters.}
    \item{\code{area}}{A factor indicating the sub area to which the locations belongs.}
  }
}

\references{
Oliveira, M. C. N. (2003). M{\'e}todos de estimac{\~a}o de par{\^a}metros em modelos geoestat{\'i}sticos com diferentes estruturas de covari{\^a}ncias: uma aplicac{\~a}o ao teor de c{\'a}lcio no solo. Ph.D. thesis, ESALQ/USP/Brasil.
}


