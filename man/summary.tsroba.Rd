\name{summary.tsroba}
\alias{summary.tsroba}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Summary of a tsroba object
Summarize the out put of the stundent-t spatial regression
}
\description{
  \code{summary} method for class "tsroba".
}

\usage{
\method{summary}{tsroba}(object, \dots)
}

\arguments{
    \item{object}{
%%     ~~Describe \code{est} here~~
An object of the class "tsroba" (see \code{\link{tsroba}} function).
}

\item{\dots}{Additional arguments.}
   }

\value{
 \item{mean.str}{
Estimates for the mean structure parameters \eqn{\bold{beta}}.
}
 \item{var.str}{
Estimates for the variance structure parameters \eqn{\sigma^2, \phi, \nu}.
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
Sample obtained for \eqn{\nu}.
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

######covariance matern: kappa=0.3 prior:reference
rest=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
           ini.pars=c(10,3,10),iter=2,burn=0,thin=1)

summary(rest)

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
}
}
