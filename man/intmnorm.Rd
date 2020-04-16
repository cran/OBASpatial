\name{intmnorm}
\alias{intmnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal posterior density for a model.
}
\description{
It calculates the marginal density density for a model \eqn{M} (up to a proportionality constant) for the NSR model using the based reference, Jeffreys' rule, Jeffreys' independent and vague  priors.  In this context \eqn{\phi} corresponds to the range parameter.
}
\usage{
intmnorm(formula,prior="reference",coords.col=1:2,kappa=0.5,
cov.model="exponential",data,asigma=2.1,intphi,maxEval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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

\item{maxEval}{
Maximum number of iterations for the integral computation.
}


}

\details{
Let \eqn{m_k} a parametric model with parameter vector \eqn{\theta_k}. Under the TSR model and the prior density proposal:

\deqn{\frac{\pi(\phi)}{(\sigma^2)^a}}

we have that the marginal density is given by:

\deqn{\int L(\theta_{m_k})\pi(m_k)dm_k}

This quantity can be useful as a criteria for model selection. The computation of \eqn{m_k} could be compute demanding depending on the number of iterations in \code{maxEval}.
}
\value{
Marginal density of the model \eqn{m_k} for the reference based, Jeffreys' rule, Jeffreys' independent and vague priors.
}
\references{
Berger, J.O,  De Oliveira, V. and Sans\'{o}, B. (2001). Objective Bayesian Analysis of Spatially Correlated Data. Journal of the American Statistical Association., 96, 1361 -- 1374.
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dnsrposoba}},\code{\link{dtsrprioroba}},\code{\link{dnsrprioroba}}
}
\examples{

\dontshow{
data(dataca20)

set.seed(25)
data(dataelev)## data using by Berger et. al (2001)

######### Using reference prior ###########
m1=intmnorm(prior="reference",formula=elevation~1,
kappa=0.5,cov.model="matern",data=dataelev,maxEval=1)

log(m1)


######### Using reference prior kappa=1 ###########
m2=intmnorm(prior="reference",formula=elevation~1,
kappa=1,cov.model="matern",data=dataelev,maxEval=1)
log(m2)

######### Using reference prior kappa=1.5 ###########
m3=intmnorm(prior="reference",formula=elevation~1
,kappa=1.5,cov.model="matern",data=dataelev,maxEval=1)
log(m3)

tot=m1+m2+m3

########posterior probabilities: higher probability:
#########prior="reference", kappa=1
p1=m1/tot
p2=m2/tot
p3=m3/tot


}


\donttest{

data(dataca20)

set.seed(25)
data(dataelev)## data using by Berger et. al (2001)

######### Using reference prior ###########
m1=intmnorm(prior="reference",formula=elevation~1,
kappa=0.5,cov.model="matern",data=dataelev,maxEval=1000)

log(m1)


######### Using reference prior kappa=1 ###########
m2=intmnorm(prior="reference",formula=elevation~1,
kappa=1,cov.model="matern",data=dataelev,maxEval=1000)
log(m2)

######### Using reference prior kappa=1.5 ###########
m3=intmnorm(prior="reference",formula=elevation~1
,kappa=1.5,cov.model="matern",data=dataelev,maxEval=1000)
log(m3)

tot=m1+m2+m3

########posterior probabilities: higher probability:
#########prior="reference", kappa=1
p1=m1/tot
p2=m2/tot
p3=m3/tot
}
}
