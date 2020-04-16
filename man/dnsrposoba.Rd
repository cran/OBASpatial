\name{dnsrposoba}
\alias{dnsrposoba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective posterior density for the NSR model
}
\description{
It calculates the density function \eqn{\pi(\phi)} (up to a proportionality constant) for the TSR model using the based reference, Jeffreys' rule, Jeffreys' independent and vague priors. In this context \eqn{\phi} corresponds to the range parameter.
}
\usage{
dnsrposoba(x,formula,prior="reference",coords.col=1:2,
kappa=0.5,cov.model="exponential",data,asigma=2.1,intphi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The \eqn{\phi} quantil value.
}

 \item{formula}{
A valid formula for a linear regression model.
}


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

\item{asigma}{
Value of \eqn{a} for  vague prior.
}

\item{intphi}{
An interval for \eqn{\phi} used for vague prior.
}

}

\details{
The posterior distribution is computed for this priors under the improper family \eqn{\frac{\pi(\phi)}{(\sigma^2)^a}}. For the vague prior, it was considered the structure where a priori, \eqn{\phi} folows an uniform distribution on the interval \code{intphi}.

For the Jeffreys independent prior, this family of priors generates improper posterior distribution when intercept is considered for the mean function.

}
\value{
Posterior density of x=\eqn{\phi}.
}
\references{
Berger, J.O,  De Oliveira, V. and Sans\'{o}, B. (2001). Objective Bayesian Analysis of Spatially Correlated Data. Journal of the American Statistical Association., 96, 1361 -- 1374.
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dtsrposoba}},\code{\link{dtsrprioroba}},\code{\link{dnsrprioroba}}
}
\examples{


data(dataelev)

######### Using reference prior ###########
dnsrposoba(x=5,prior="reference",formula=elevation~1,
kappa=1,cov.model="matern",data=dataelev)

######### Using Jeffreys' rule prior ###########
dnsrposoba(x=5,prior="jef.rul",formula=elevation~1,
kappa=1,cov.model="matern",data=dataelev)

######### Using vague independent prior ###########
dnsrposoba(x=5,prior="vague",formula=elevation~1,
kappa=0.3,cov.model="matern",data=dataelev,intphi=c(0.1,10))

}
