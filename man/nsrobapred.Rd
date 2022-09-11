\name{nsrobapred1}
\alias{nsrobapred1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction under Normal Objective Bayesian Analysis (OBA).
}
\description{
This function uses the sampling distribution of parameters obtained from the function \code{tsroba} to predict values at unknown locations.
}
\usage{
nsrobapred1(xpred, coordspred, obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xpred}{
Values of the X design matrix for prediction coordinates.
}
  \item{coordspred}{
Points coordinates to be predicted.
}
  \item{obj}{
%%     ~~Describe \code{est} here~~
object of the class "nsroba" (see \code{\link{nsroba}} function).
}
}
\details{
This function predicts using the sampling distribution of parameters obtained from the function \code{nsroba} and the conditional normal distribution of the predicted values given the data.
}
\value{
This function returns a vector with the predicted values at the specified locations.
}
\references{
Berger, J.O,  De Oliveira, V. and Sanso, B. (2001). Objective Bayesian Analysis of Spatially Correlated Data. Journal of the American Statistical Association., 96, 1361 -- 1374.

Diggle, P. and P. Ribeiro (2007).Model-Based Geostatistics.  Springer Series in Statistics.

}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}


\seealso{
\code{\link{nsroba}},\code{\link{tsrobapred}}
}
\examples{


\dontshow{
set.seed(25)
data(dataelev)
d1=dataelev[1:42,]

reselev=nsroba(elevation~1, kappa = 0.5, cov.model = "matern", data=d1,
ini.pars=c(10,3),intphi=c(0.8,10),iter=2,burn=0,thin=1)

datapred1=dataelev[43:52,]
coordspred1=datapred1[,1:2]
nsrobapred1(obj=reselev,coordspred=coordspred1,xpred=rep(1,10))
}


\donttest{
set.seed(25)
data(dataelev)
d1=dataelev[1:42,]

reselev=nsroba(elevation~1, kappa = 0.5, cov.model = "matern", data=d1,
ini.pars=c(10,3),intphi=c(0.8,10))

datapred1=dataelev[43:52,]
coordspred1=datapred1[,1:2]
nsrobapred1(obj=reselev,coordspred=coordspred1,xpred=rep(1,10))
}
}
