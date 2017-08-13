\name{kaps-package}
\alias{kaps-package}
\docType{package}
\title{
	K-adaptive partitioning for survival data.
}
\description{
This package provides some routines to conduct a K-adaptive partitioning algorithm, which divides the dataset into K heterogeneous subgroups based on the information from a prognostic factor.
}
\details{
\tabular{ll}{
Package: \tab kaps\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2014-11-01\cr
License: \tab GPL-3\cr
LazyLoad: \tab no\cr
}
This package contains some routines to conduct a \emph{K}-adaptive partitioning for survival data (kaps) algorithm. A function \code{kaps()} is an implementation version of our algorithm which provides minimax-based partitioning rule. 
}
\author{
Soo-Heang Eo <hanansh@korea.ac.kr> \cr
Seung-Mo Hong <smhong28@gmail.com> \cr
HyungJun Cho <hj4cho@korea.ac.kr>
}

\references{
Eo, S. H., Kang, H. J., Hong, S. M., and Cho, H. (2013). K-adaptive partitioning for survival data, with an application to cancer staging. \emph{arXiv preprint arXiv:1306.4615}.
}
\seealso{
    \code{\link{kaps}} 
}
\keyword{ package }