

#' Class \code{"kaps"}
#' 
#' A S4 class for \emph{K}-adaptive partitioning for survival data (kaps).
#' 
#' 
#' @name kaps-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("kaps")}. The most important slot is \code{groupID}, which is a
#' vector consisting of the information about classified subgroups.
#' @keywords classes
#' @examples
#' 
#' showClass("kaps")
#' 
NULL





#' K-adaptive partitioning for survival data.
#' 
#' This package provides some routines to conduct a K-adaptive partitioning
#' algorithm, which divides the dataset into K heterogeneous subgroups based on
#' the information from a prognostic factor.
#' 
#' \tabular{ll}{ Package: \tab kaps\cr Type: \tab Package\cr Version: \tab
#' 1.0.2\cr Date: \tab 2014-11-01\cr License: \tab GPL-3\cr LazyLoad: \tab
#' no\cr } This package contains some routines to conduct a \emph{K}-adaptive
#' partitioning for survival data (kaps) algorithm. A function \code{kaps()} is
#' an implementation version of our algorithm which provides minimax-based
#' partitioning rule.
#' 
#' @name kaps-package
#' @docType package
#' @author Soo-Heang Eo <hanansh@@korea.ac.kr> \cr Seung-Mo Hong
#' <smhong28@@gmail.com> \cr HyungJun Cho <hj4cho@@korea.ac.kr>
#' @seealso \code{\link{kaps}}
#' @references Eo, S. H., Kang, H. J., Hong, S. M., and Cho, H. (2013).
#' K-adaptive partitioning for survival data, with an application to cancer
#' staging. \emph{arXiv preprint arXiv:1306.4615}.
#' @keywords package
NULL





#' Visualize an object "kaps"
#' 
#' \code{plot} method for "kaps" with extended facilities. It provides four
#' panels consisting of a scatter plot, a Kaplan-Meier survival curve, an
#' overall p-values, and a plot with the worst-pair p-values against K.
#' 
#' This function generates four plots. The top left panel is the scatterplot of
#' survival times against the selected prognostic factor with the line fitted
#' by local censored regression using \code{locfit}. The top right panel is a
#' Kaplan-Meier survival curve for the subgroups selected with the optimal
#' \emph{K}. At the bottom are displayed the plots of the overall and
#' worst-pair p-values against K. The dotted lines indicate thresholds for
#' significance 0.05.The outputs for a specific K can also be printed out with
#' the argument K.
#' 
#' For the sake of the Kaplan-Meier curve with estimated subgroups, in
#' addition, the function \code{km.curve} is provided.
#' 
#' @name plot
#' @aliases plot,kaps-method plot,kaps,missing-method
#' @docType methods
#' @param x an object from \code{kaps}
#' @param y the "y" argument is not used in the plot-method for "kaps" object.
#' @param K a scalar object that plots the Kaplan-Meier survival curves for the
#' K. If missing, it works with selected K in the model fitting.
#' @param \dots other arguments to the \code{\link[=graphics]{plot.default}}
#' function can be passed here.
#' @seealso \code{\link{kaps}} \cr \code{\link{km.curve}}
#' @keywords methods
NULL





#' Predict new values using the fitted object "kaps".
#' 
#' This function provides the predicted subgroups and its test statistic.
#' 
#' 
#' @name predict
#' @aliases predict predict,kaps-method
#' @docType methods
#' @param object object from \code{kaps}.
#' @param newdata An optional argument in which the name of predicted object is
#' located. If omitted, the dataset used in the model fitting is utilized.
#' @param type a type of prediction. If the type is "predict", it predicts
#' subgroups based on the fitted model. If the type is "kaps", it returns the
#' overall and its worst-pair test statistic.
#' @seealso \code{\link{kaps}}
#' @keywords methods
NULL





#' Print an object "kaps" with specific information about K
#' 
#' It functions like \code{show} but the only difference is the output with the
#' specific information about K.
#' 
#' 
#' @name print
#' @aliases print,kaps-method
#' @docType methods
#' @param x an object from \code{kaps}
#' @param K a scalar object to determine the number of subgroups K. If missing,
#' the estimated subgroup K is selected.
#' @seealso \code{\link{kaps}}
#' @keywords methods
NULL





#' Show an object "kaps"
#' 
#' It returns the outputs of the object "kaps" consisting of three parts. The
#' first part displays the model formula with a dataset and the selected number
#' for K. Next, the information regarding the selection of an optimal set of
#' cut-off points is provided. Lastly, the p-values of pairwise two-sample test
#' comparisons among all the pairs of subgroups are provided.
#' 
#' 
#' @name show
#' @aliases show show,kaps-method
#' @docType methods
#' @param object object from \code{kaps}.
#' @seealso \code{\link{kaps}}
#' @keywords methods
NULL





#' Summarize an object "kaps"
#' 
#' This function provides the tabloid information with survival median, 1-, 3-,
#' and 5 years actual survival time for each partition.
#' 
#' 
#' @name summary
#' @aliases summary summary,kaps-method
#' @docType methods
#' @param object object with the class \code{kaps}.
#' @param K scalar object to determine the number of subgroups K. If missing,
#' the estimated subgroup K is selected.
#' @seealso \code{\link{kaps}}
#' @keywords methods
NULL





#' toy example
#' 
#' toy dataset
#' 
#' 
#' @name toy
#' @docType data
#' @return \item{meta}{covariate variable that describes the number of
#' metastatic lymph nodes} \item{status}{censoring status} \item{time}{time to
#' event}
#' @keywords datasets
NULL



