\name{tsrobapred}
\alias{tsrobapred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction under Student-t Objective Bayesian Analysis (OBA).
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function uses the sampling distribution of parameters obtained from the function \code{tsroba} to predict values at unknown locations.
}
\usage{
tsrobapred(obj,xpred,coordspred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{obj}{
%%     ~~Describe \code{est} here~~
object of the class "tsroba" (see \code{\link{tsroba}} function).
}
  \item{xpred}{
%%     ~~Describe \code{xpred} here~~
Values of the X design matrix for prediction coordinates.
}
  \item{coordspred}{
%%     ~~Describe \code{coordspred} here~~
Points coordinates to be predicted.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function predicts using the sampling distribution of parameters obtained from the function \code{tsroba} and the conditional Student-t distribution of the predicted values given the data.
}
\value{
This function returns a vector with the predicted values at the specified locations.
}
\references{
Diggle, P. and P. Ribeiro (2007).Model-Based Geostatistics.  Springer Series in Statistics.

Ordonez, J.A, M.O. Prattes, L.A. Matos, and V.H. Lachos (2020+). Objective Bayesian analysis for spatial Student-t regression models. (Submitted)
}
\author{
Jose A. Ordonez, Marcos O. Prates, Larissa A. Matos, Victor H. Lachos.
}

\seealso{
\code{\link{tsroba}},\code{\link{nsrobapred}}
}
\examples{
\dontshow{
set.seed(25)
data(dataca20)
d1=dataca20[1:158,]

######covariance matern: kappa=0.3 prior:reference
res=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
ini.pars=c(10,3,10),iter=2,thin=1,burn=0)

datapred=dataca20[159:178,]
formula=calcont~altitude+area
xpred=model.matrix(formula,data=datapred)

tsrobapred(res,xpred=xpred,coordspred=dataca20[159:178,1:2])
}
\donttest{
set.seed(25)
data(dataca20)
d1=dataca20[1:158,]

######covariance matern: kappa=0.3 prior:reference
res=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
ini.pars=c(10,3,10),iter=50,thin=1,burn=5)

datapred=dataca20[159:178,]
formula=calcont~altitude+area
xpred=model.matrix(formula,data=datapred)

tsrobapred(res,xpred=xpred,coordspred=dataca20[159:178,1:2])
}

}
