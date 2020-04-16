\name{dtsrprioroba}
\alias{dtsrprioroba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective prior density for the TSR model
}
\description{
It calculates the density function \eqn{\pi(\phi,\nu)} (up to a proportionality constant) for the TSR model using the based reference, Jeffreys' rule and Jeffreys' independent  priors. In this context \eqn{\phi} corresponds to the range parameter and \eqn{\nu} to the degrees of freedom.
}
\usage{
dtsrprioroba(x,trend="cte",prior="reference",coords.col=1:2,
kappa=0.5,cov.model="exponential",data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector with the quanties \eqn{(\phi,\nu)}
}

 \item{trend}{
Builds the trend matrix in accordance to a specification of the mean provided by the user. See \code{DETAILS} below.}


 \item{prior}{
Objective prior densities avaiable for the TSR model: ( \code{reference}: Reference based, \code{jef.rul}: Jeffreys' rule, \code{jef.ind}: Jeffreys' independent)
}

\item{coords.col}{
A vector with the column numbers corresponding to the spatial coordinates.
}

\item{kappa}{
Shape parameter of the covariance function (fixed)
}


\item{cov.model}{ Covariance functions available for the TSR
model. \code{matern}: Matern, \code{pow.exp}: power exponential, \code{exponential}:exponential, \code{cauchy}: Cauchy, \code{spherical}: Spherical}

\item{data}{
Data set with 2D spatial coordinates, the response and optional covariates
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
Density of x=(\eqn{\phi,\nu})
}
\references{
Ordonez, J. A., Prates, M. O., Matos, L. A., and Lachos, V. H. (2020). Objective Bayesian analysis for spatial Student-t regression models. arXiv preprint arXiv:2004.04341.
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dtsrposoba}},\code{\link{dnsrprioroba}},\code{\link{dnsrposoba}}
}
\examples{
data(dataca20)

######### Using reference prior and a constant trend###########
dtsrprioroba(x=c(6,100),kappa=0.3,cov.model="matern",data=dataca20)


######### Using jef.rule prior and 1st trend###########
dtsrprioroba(x=c(6,100),prior="jef.rul",trend=~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20)

######### Using  jef.ind prior ###########
dtsrprioroba(x=c(6,100),prior="jef.ind",trend=0,
kappa=0.3,cov.model="matern",data=dataca20)
}
