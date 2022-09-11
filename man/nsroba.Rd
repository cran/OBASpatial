\name{nsroba}
\alias{nsroba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian estimation for the NSR model.
}
\description{
This function performs Bayesian estimation of \eqn{\theta=(\bold{\beta},\sigma^2,\phi)} for the NSR model using the based reference, Jeffreys' rule ,Jeffreys' independent and vague  priors.
}
\usage{
nsroba(formula, method="median",
prior = "reference",coords.col = 1:2,kappa = 0.5,
cov.model = "matern", data,asigma=2.1, intphi = "default",
ini.pars, burn=500, iter=5000, thin=10,
cprop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A valid formula for a linear regression model.
}


  \item{method}{
  Method to estimate (\eqn{\bold{beta},\sigma,\phi}). The methods availables are \code{"mean"},\code{"median"} and \code{"mode"}.
}

  \item{prior}{
Objective prior densities avaiable for the TSR model: ( \code{reference}: Reference based, \code{jef.rul}: Jeffreys' rule, \code{jef.ind}: Jeffreys' independent, \code{vague}, Vague).
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
Value of \eqn{a} for  the vague prior.
}

\item{intphi}{
An interval for \eqn{\phi} used for the uniform proposal. See \code{DETAILS} below.
}

  \item{ini.pars}{
Initial values for \eqn{(\sigma^2,\phi)} in that order.
}
  \item{burn}{
Number of observations considered in the burning process.
}
  \item{iter}{
Number of iterations for the sampling procedure.
}
  \item{thin}{
Number of observations considered in the thin process.
}

  \item{cprop}{
A constant related to the acceptance probability
(Default = NULL indicates that cprop is computed as the interval length of intphi). See
\code{DETAILS} below.
}

}
\details{
For the \code{"unif"} proposal, it was considered the structure where a priori, \eqn{\phi} follows an uniform distribution on the interval \code{intphi}. By default, this interval is computed using the empirical range of data as well as the constant \code{cprop}.

For the Jeffreys independent prior, the sampling procedure generates improper posterior distribution when intercept is considered for the mean function.
}
\value{

\item{$dist}{
Joint sample (matrix object) obtaining for (\eqn{\bold{beta},\sigma^2,\phi}).
}
\item{$betaF}{
Sample obtained for \eqn{\bold{beta}}.
}

\item{$sigmaF}{
Sample obtained for \eqn{\sigma^2}.
}

\item{$phiF}{
Sample obtained for \eqn{\phi}.
}

\item{$coords}{
Spatial data coordinates.
}

\item{$kappa}{
Shape parameter of the covariance function.
}

\item{$X}{
Design matrix of the model.
}

\item{$type}{
Covariance function of the model.
}

\item{$theta}{
Bayesian estimator of (\eqn{\bold{beta},\sigma,\phi}).
}

\item{$y}{
Response variable.
}

\item{$prior}{
Prior density considered.
}
}
\references{
Berger, J.O,  De Oliveira, V. and Sanso, B. (2001). Objective Bayesian Analysis of Spatially Correlated Data. Journal of the American Statistical Association., 96, 1361 -- 1374.
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dnsrposoba}},\code{\link{dtsrprioroba}},\code{\link{dnsrprioroba},\link{tsroba}}
}

\examples{

\dontshow{
set.seed(25)
data(dataelev)


######covariance matern: kappa=0.5
res=nsroba(elevation~1, kappa = 0.5, cov.model = "matern", data=dataelev,
ini.pars=c(10,390),iter=2,burn=0,thin=1)

summary(res)

}


\donttest{

set.seed(25)
data(dataelev)


######covariance matern: kappa=0.5
res=nsroba(elevation~1, kappa = 0.5, cov.model = "matern", data=dataelev,
ini.pars=c(10,390))

summary(res)

######covariance matern: kappa=1
res1=nsroba(elevation~1, kappa = 1, cov.model = "matern", data=dataelev,
ini.pars=c(10,390))

summary(res1)

######covariance matern: kappa=1.5
res2=nsroba(elevation~1, kappa = 1.5, cov.model = "matern", data=dataelev,
ini.pars=c(10,390))

summary(res2)



}
}

