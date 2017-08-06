\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,kaps-method}
\title{Predict new values using the fitted object "kaps".}
\description{
	This function provides the predicted subgroups and its test statistic.
}
\usage{
\S4method{predict}{kaps}(object, newdata, type = c("predict","kaps"))
}
\arguments{
  \item{object}{object from \code{kaps}.}
  \item{newdata}{An optional argument in which the name of predicted object is located. If omitted, the dataset used in the model fitting is utilized.}
  \item{type}{a type of prediction. If the type is "predict", it predicts subgroups based on the fitted model. If the type is "kaps", it returns the overall and its worst-pair test statistic.}
}
\seealso{
  \code{\link{kaps}}
}
\keyword{methods}
