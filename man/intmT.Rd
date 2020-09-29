\name{intmT}
\alias{intmT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal posterior density for a model.
}
\description{
It calculates the marginal density density for a model \eqn{M} (up to a proportionality constant) for the TSR model using the based reference, Jeffreys' rule, Jeffreys' independent and vague  priors.  In this context \eqn{\phi} corresponds to the range parameter and \eqn{\nu} to the degrees of freedom.
}
\usage{
intmT(formula,prior="reference",coords.col=1:2,kappa=0.5,
cov.model="exponential",data,asigma,intphi="default",intnu=c(4.1,Inf),maxEval)
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

\item{intnu}{
An interval for \eqn{\nu} used for vague prior.
}

\item{maxEval}{
Maximum number of iterations for the integral computation.
}


}

\details{
Let \eqn{m_k} a parametric model with parameter vector \eqn{\theta_k}. Under the TSR model and the prior density proposal:

\deqn{\frac{\pi(\phi,\nu)}{(\sigma^2)^a}}

we have that the marginal density is given by:

\deqn{\int L(\theta_{m_k})\pi(m_k)dm_k}

This quantity can be useful as a criteria for model selection. The computation of \eqn{m_k} could be compute demanding depending on the number of iterations in \code{maxEval}.
}
\value{
Marginal density of the model \eqn{m_k} for the reference based, Jeffreys' rule, Jeffreys' independent and vague priors.
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

\dontshow{
set.seed(25)
data(dataca20)

######### Using reference prior ###########
m1=intmT(prior="reference",formula=calcont~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20,maxEval=1)

}

\donttest{

set.seed(25)
data(dataca20)



######### Using reference prior ###########
m1=intmT(prior="reference",formula=calcont~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000)

######### Using Jeffreys' rule prior ###########
m1j=intmT(prior="jef.rul",formula=calcont~altitude+area,
kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000)


######### Using Jeffreys' independent prior ###########
m1ji=intmT(prior="jef.ind",formula=calcont~altitude+area
,kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000)

m1v=intmT(prior="vague",formula=calcont~altitude+area
,kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000,intphi="default")




tot=m1+m1j+m1ji+m1v

########posterior probabilities: higher probability:
#########prior="reference", kappa=0.3


p1=m1/tot
pj=m1j/tot
pji=m1ji/tot
pv=m1v/tot


}
}
