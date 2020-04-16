\name{dnsrprioroba}
\alias{dnsrprioroba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective prior density for the NSR model
}
\description{
It calculates the density function \eqn{\pi(\phi)} (up to a proportionality constant) for the NSR model using the based reference, Jeffreys' rule and Jeffreys' independent  priors. In this context \eqn{\phi} corresponds to the range parameter.
}
\usage{
dnsrprioroba(x,trend="cte",prior="reference",coords.col=1:2,
kappa=0.5,cov.model="exponential",data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The \eqn{\phi} quantil value.
}

 \item{trend}{
Builds the trend matrix in accordance to a specification of the mean provided by the user. See \code{DETAILS} below.}


 \item{prior}{
Objective prior densities avaiable for the TSR model: ( \code{reference}: Reference based, \code{jef.rul}: Jeffreys' rule, \code{jef.ind}: Jeffreys' independent).
}

\item{coords.col}{
A vector with the column numbers corresponding to the spatial coordinates.
}

\item{kappa}{
Shape parameter of the covariance function (fixed).
}


\item{cov.model}{ Covariance functions available for the TSR
model. \code{matern}: Matern, \code{pow.exp}: power exponential, \code{exponential}:exponential, \code{cauchy}: Cauchy, \code{spherical}: Spherical.}

\item{data}{
Data set with 2D spatial coordinates, the response and optional covariates.
}

}

\details{
Denote as \eqn{\bold{c}=(c_{1},c_{2})} the coordinates of a spatial location. \code{trend} defines the design matrix as:

\itemize{
\item \code{0} (zero,without design matrix) Only valid for the Independent Jeffreys' prior
\item \code{"cte"}, the design matrix is such that mean function \eqn{\mu(\bold{c})=\mu} is constant over the region.
\item \code{"1st"}, the design matrix is such that mean function becames a first order polynomial on the coordinates:
\deqn{\mu(\bold(c))=\beta_0+ \beta_1c_1+\beta_2c_2}
\item \code{"2nd"}, the design matrix is such that mean function \eqn{\mu(\bold{c})=\mu} becames a second order polynomial on the coordinates:
\deqn{\mu(\bold(c))=\beta_0+ \beta_1c_1+\beta_2c_2 + \beta_3c_{1}^2+ \beta_4c_{2}^2+ \beta_5c_1c_2}
\item \code{~model} a model specification to include covariates (external trend) in the model.
}

}
\value{
Prior density of x=\eqn{\phi}
}
\references{
Berger, J.O,  De Oliveira, V. and Sans\'{o}, B. (2001). Objective Bayesian Analysis of Spatially Correlated Data. Journal of the American Statistical Association., 96, 1361 -- 1374.
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dtsrposoba}},\code{\link{dtsrprioroba}},\code{\link{dnsrposoba}}
}
\examples{

data(dataelev)## data using by Berger et. al (2001)

######### Using reference prior ###########
dnsrprioroba(x=20,kappa=0.3,cov.model="matern",data=dataelev)


######### Using jef.rule prior###########
dnsrprioroba(x=20,prior="jef.rul",kappa=0.3,cov.model="matern",
data=dataelev)

######### Using  jef.ind prior ###########
dnsrprioroba(x=20,prior="jef.ind",trend=0,
kappa=0.3,cov.model="matern",data=dataelev)
}
