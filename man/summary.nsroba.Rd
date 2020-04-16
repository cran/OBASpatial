\name{summary.nsroba}
\alias{summary.nsroba}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Summary of a nsroba object
Summarize the out put of the normal spatial regression
}
\description{
  \code{summary} method for class "nsroba".
}

\usage{
\method{summary}{nsroba}(object,\dots)
}

\arguments{
    \item{object}{
%%     ~~Describe \code{est} here~~
An object of the class "nsroba" (see \code{\link{nsroba}} function).
}
\item{\dots}{Additional arguments.}
   }

\value{
 \item{mean.str}{
Estimates for the mean structure parameters \eqn{\bold{beta}}.
}
 \item{var.str}{
Estimates for the variance structure parameters \eqn{\sigma^2, \phi}.
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
}
\references{
Berger, J.O,  De Oliveira, V. and Sans\'{o}, B. (2001). Objective Bayesian Analysis of Spatially Correlated Data. Journal of the American Statistical Association., 96, 1361 -- 1374.
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
resa=nsroba(elevation~1, kappa = 0.5, cov.model = "matern", data=dataelev,
ini.pars=c(10,3),iter=2,burn=0,thin=1)
summary(resa)

}


\donttest{

set.seed(25)
data(dataelev)


######covariance matern: kappa=0.5
res=nsroba(elevation~1, kappa = 0.5, cov.model = "matern", data=dataelev,
ini.pars=c(10,3))

summary(res)



}
}
