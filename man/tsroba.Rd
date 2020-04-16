\name{tsroba}
\alias{tsroba}

\title{
Bayesian estimation for the TSR model.
}
\description{
This function performs Bayesian estimation of \eqn{\theta=(\bold{\beta},\sigma^2,\phi)} for the TSR model using the based reference, Jeffreys' rule ,Jeffreys' independent and vague  priors.
}
\usage{
tsroba(formula, proposal="unif", method="median",candpar,
prior = "reference",coords.col = 1:2,kappa = 0.5,
cov.model = "matern", data,asigma=2.1, intphi = "default",
intnu="default",ini.pars,burn=500, iter=5000,thin=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A valid formula for a linear regression model.
}
  \item{proposal}{
  Proposal for \eqn{\phi}, the range parameter. The proposals availables are \code{"unif"} and \code{"gamma"}. See \code{DETAILS} below.
}
 \item{method}{
  Method to estimate (\eqn{\bold{beta},\sigma,\phi,\nu}). The methods availables are \code{"mean"},\code{"median"} and \code{"mode"}.
}
  \item{candpar}{
Scale parameter when the proposal is \code{"gamma"}, unnecessary when proposal=\code{"uniform"}
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
An interval for \eqn{\phi} used for the uniform proposal. See \code{DETAILS} below.
}

\item{intnu}{
An interval for \eqn{\nu} used for the uniform proposal. See \code{DETAILS} below.
}

  \item{ini.pars}{
Initial values for \eqn{(\sigma^2,\phi,\nu)} in that order.
}
  \item{burn}{
Number of observations considered in burning process.
}
  \item{iter}{
Number of iterations for the sampling procedure.
}
  \item{thin}{
Number of observations considered in thin process.
}
}
\details{
For the prior proposal, it was considered the structure \eqn{\pi(\phi,\nu,\lambda)=\phi(\phi)\pi(\nu|\lambda)\pi(\lambda)}. The \code{"unif"} denotes the proposal where a priori, \eqn{\phi} follows an uniform distribution on the interval \code{intphi}, by default, this interval is computed using the empirical range of data.\eqn{\nu|\lambda~ Texp(\lambda,A)} with A the interval given by the argument \code{intnu} and \eqn{\lambda~unif(0.02,0.5)}

For the Jeffreys independent prior, the sampling procedure generates improper posterior distribution when intercept is considered for the mean function.
}
\value{

\item{dist}{
Joint sample (matrix object) obtaining for (\eqn{\bold{beta},\sigma^2,\phi}).
}
\item{betaF}{
Sample obtained for \eqn{\bold{beta}}.
}

\item{sigmaF}{
Sample obtained for \eqn{\sigma^2}.
}

\item{phiF}{
Sample obtained for \eqn{\phi}.
}

\item{nuF}{
Sample obtained for \eqn{\phi}.
}


\item{coords}{
Spatial data coordinates.
}

\item{kappa}{
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
Ordonez, J. A., Prates, M. O., Matos, L. A., and Lachos, V. H. (2020). Objective Bayesian analysis for spatial Student-t regression models. arXiv preprint arXiv:2004.04341.
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{dnsrposoba}},\code{\link{dtsrprioroba}},\code{\link{dnsrprioroba}},\code{\link{tsroba}}
}
\examples{

\dontshow{

set.seed(25)
data(dataca20)
d1=dataca20[1:158,]

xpred=model.matrix(calcont~altitude+area,data=dataca20[159:178,])
xobs=model.matrix(calcont~altitude+area,data=dataca20[1:158,])
coordspred=dataca20[159:178,1:2]

######covariance matern: kappa=0.3 prior:reference
res=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
           ini.pars=c(10,3,10),iter=2,burn=0,thin=1)

}


\donttest{

set.seed(25)
data(dataca20)
d1=dataca20[1:158,]

xpred=model.matrix(calcont~altitude+area,data=dataca20[159:178,])
xobs=model.matrix(calcont~altitude+area,data=dataca20[1:158,])
coordspred=dataca20[159:178,1:2]

######covariance matern: kappa=0.3 prior:reference
res=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
           ini.pars=c(10,3,10))

summary(res)

######covariance matern: kappa=0.3 prior:jef.rul
res1=tsroba(calcont~altitude+area, kappa = 0.3,
            data=d1,prior="jef.rul",ini.pars=c(10,3,10))

summary(res1)

######covariance matern: kappa=0.3 prior:jef.ind
res2=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,prior="jef.ind",
ini.pars=c(10,3,10))

summary(res2)

######covariance matern: kappa=0.3 prior:vague
res3=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,prior="vague",
ini.pars=c(10,3,10))

summary(res3)

####obtaining posterior probabilities
###(just comparing priors with kappa=0.3).
###the real aplication (see Ordonez et.al) consider kappa=0.3,0.5,0.7.

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

####posterior probabilities#####
p1=m1/tot
pj=m1j/tot
pji=m1ji/tot
pv=m1v/tot


##########MSPE#######################################

pme=tsrobapred(res,xpred=xpred,coordspred=coordspred)
pme1=tsrobapred(res1,xpred=xpred,coordspred=coordspred)
pme2=tsrobapred(res2,xpred=xpred,coordspred=coordspred)
pme3=tsrobapred(res3,xpred=xpred,coordspred=coordspred)

mse=mean((pme-dataca20$calcont[159:178])^2)
mse1=mean((pme1-dataca20$calcont[159:178])^2)
mse2=mean((pme2-dataca20$calcont[159:178])^2)
mse3=mean((pme3-dataca20$calcont[159:178])^2)

}
}
