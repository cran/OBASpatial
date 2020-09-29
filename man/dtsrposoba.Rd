\name{dtsrposoba}
\alias{dtsrposoba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective posterior density for the TSR model
}
\description{
It calculates the density function \eqn{\pi(\phi,\nu)} (up to a proportionality constant) for the TSR model using the based reference, Jeffreys' rule, Jeffreys' independent and vague priors. In this context \eqn{\phi} corresponds to the range parameter and \eqn{\nu} to the degrees of freedom.
}
\usage{
dtsrposoba(x,formula,prior="reference",coords.col=1:2,
kappa=0.5,cov.model="exponential",data,asigma=2.1,intphi,intnu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector with the quanties \eqn{(\phi,\nu)}. For the vague prior x must be a three dimension vector \eqn{(\phi,\nu,\lambda)} with \eqn{\lambda} a number in the interval \eqn{(0.02,0.5)}. See \code{DETAILS} below.
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

\item{intnu}{
An interval for \eqn{\nu} used for vague prior.
}


}

\details{
The posterior distribution is computed for this priors under the improper family \eqn{\frac{\pi(\phi,\nu)}{(\sigma^2)^a}}. For the vague prior, it was considered the structure \eqn{\pi(\phi,\nu,\lambda)=\phi(\phi)\pi(\nu|\lambda)\pi(\lambda)} where a priori, \eqn{\phi} follows an uniform distribution on the interval \code{intphi}, \eqn{\nu|\lambda~ Texp(\lambda,A)} with A the interval given by the argument \code{intnu} and \eqn{\lambda~unif(0.02,0.5)}.

For the Jeffreys independent prior, this family of priors generates improper posterior distribution when intercept is considered for the mean function.

}
\value{
Posterior density of x=(\eqn{\phi,\nu}) for the reference based, Jeffreys' rule and Jeffreys' independent priors. For the vague the result is the posterior density of x=(\eqn{\phi,\nu,\lambda})
}
\references{
Ordonez, J.A, M.O. Prattes, L.A. Matos, and V.H. Lachos (2020+). Objective Bayesian analysis for spatial Student-t regression models (Submitted).
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dnsrposoba}},\code{\link{dtsrprioroba}},\code{\link{dnsrprioroba}}
}
\examples{
data(dataca20)

######### Using reference prior ###########
dtsrposoba(x=c(5,11),prior="reference",formula=calcont~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20)

######### Using Jeffreys' rule prior ###########
dtsrposoba(x=c(5,11),prior="jef.rul",formula=calcont~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20)


######### Using Jeffreys' independent prior ###########
dtsrposoba(x=c(5,11),prior="jef.ind",formula=calcont~altitude+area
,kappa=0.3,cov.model="matern",data=dataca20)

######### Using vague independent prior ###########
dtsrposoba(x=c(5,11,.3),prior="vague",formula=calcont~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20,intphi=c(0.1,10),
intnu=c(4.1,30))
}
