#' K-adaptive partitioning for survival data
#' 
#' Conduct \emph{K}-adaptive partitioning algorithm for survival data
#' 
#' This function provides routines to conduct KAPS algorithm which is designed
#' to classify cut-off values by the minimax-based rule.
#' 
#' @param formula Formula object with a response on the left hand side of the
#' '~' operator, and the covariate terms on the right side. The response has to
#' be a survival object with survival time and censoring status in the
#' \link[survival:Surv]{Surv} function. For more details, see
#' \link[Formula:Formula]{Formula} page.
#' @param data data frame with variables used in formula. It needs at least
#' three variables including survival time, censoring status, and a covariate.
#' Multivariate covariates can be supported with "+" sign.
#' @param K number of subgroups used in the model fitting. The default value is
#' 2:4 which means finding optimal subgroups ranging from 2 to 4.
#' @param type Select a type of algorithm in order to find optimal number of
#' subgroups. Two options are provided: \code{perm} and \code{NULL}. The
#' \code{perm} chooses subgroups using permutation procedures, while the
#' \code{NULL} passes a optimal selection algorithm.
#' @param mindat the minimum number of observations at each subgroup. The
#' default value is 5\% of observations.
#' @param list() a list of tuning parameters with the class, "kapsOptions". For
#' more details, see \link{kaps.control}.
#' @return The function returns an object with class "kaps" with the following
#' slots.  \item{list("call")}{evaluated function call}\item{:}{evaluated
#' function call} \item{list("formula")}{formula to be used in the model
#' fitting}\item{:}{formula to be used in the model fitting}
#' \item{list("data")}{data to be used in the model fitting}\item{:}{data to be
#' used in the model fitting} \item{list("groupID")}{information about the
#' subgroup classification}\item{:}{information about the subgroup
#' classification} \item{list("index")}{an index for the optimal subgroup among
#' the candidate K}\item{:}{an index for the optimal subgroup among the
#' candidate K} \item{list("X")}{test statistic with the worst pair of
#' subgroups for the split set s}\item{:}{test statistic with the worst pair of
#' subgroups for the split set s} \item{list("Z")}{the overall test statistic
#' with K subgroups using the split set s}\item{:}{the overall test statistic
#' with K subgroups using the split set s} \item{list("pair")}{selected pair of
#' subgroups}\item{:}{selected pair of subgroups}
#' \item{list("split.var")}{selected covariate in the model fitting
#' }\item{:}{selected covariate in the model fitting }
#' \item{list("split.pt")}{selected set of cut-off points}\item{:}{selected set
#' of cut-off points} \item{list("mindat")}{minimum number of observations at a
#' subgroup}\item{:}{minimum number of observations at a subgroup}
#' \item{list("test.stat")}{Bonferroni corrected p-value matrix. The first row
#' means overall p-values and the second one denotes p-values of the worst-pair
#' against K. The column in the matrix describes the order of
#' K.}\item{:}{Bonferroni corrected p-value matrix. The first row means overall
#' p-values and the second one denotes p-values of the worst-pair against K.
#' The column in the matrix describes the order of K.}
#' \item{list("over.stat.sample")}{adjusted overall test statistic by
#' Bootstrapping}\item{:}{adjusted overall test statistic by Bootstrapping}
#' \item{list("pair.stat.sample")}{adjusted worst-pair test statistic by
#' Bootstrapping}\item{:}{adjusted worst-pair test statistic by Bootstrapping}
#' \item{list("groups")}{candidate K used in the argument}\item{:}{candidate K
#' used in the argument} \item{list("results")}{list object about the results
#' of each candidate K}\item{:}{list object about the results of each candidate
#' K} \item{list("Options")}{tuning parameters}\item{:}{tuning parameters}
#' @author Soo-Heang Eo \email{eo.sooheang@@gmail.com} \cr Seung-Mo Hong
#' \email{smhong28@@gmail.com} \cr HyungJun Cho \email{hj4cho@@korea.ac.kr} \cr
#' @seealso \code{\link{show}}, \code{\link{plot}}, \code{\link{predict}},
#' \code{\link{print}} and \code{\link{summary}} for the convenient use of
#' kaps() \cr \code{\link{kaps.control}} to control kaps() more detail \cr
#' \code{\link{count.mindat}} to calculate minimum subgroup sample size
#' @references Eo, S. H., Kang, H. J., Hong, S. M., and Cho, H. (2013).
#' K-adaptive partitioning for survival data, with an application to cancer
#' staging. \emph{arXiv preprint arXiv:1306.4615}.
#' @keywords kaps
#' @examples
#' 
#'   \dontrun{
#'     data(toy)
#'     f <- Surv(time, status) ~ meta
#'     # Fit kaps algorithm without cross-validation.
#'     # It means the step to finding optimal K is not entered.
#'     fit1 <- kaps(f, data = toy, K = 3)
#' 
#'     # show the object of kaps (it contains apss S4 class)
#'     fit1
#' 
#'     # plot Kaplan-Meire estimates
#'     plot(fit1)
#' 
#'     # Fit kaps algorithm for selection optimal number of subgropus.
#'     fit2 <- kaps(f, data = toy, K= 2:4) 
#'     fit2
#' 
#'     # plot outputs with subgroup selection
#'     require(locfit) # for scatterplot smoothing
#'     plot(fit2)
#' 
#'     print(fit2,K=2)
#'     summary(fit2)
#'     summary(fit2,K=2)
#' 
#'     # require(party)
#'     # fit4 <- ctree(f, data = toy)
#'   }
#' 
kaps <- function(formula, data, K = 2:4, mindat, type = c("perm", "NULL"), ...){

  ##########################
  ##### pre-processing step
  if (missing(mindat)) mindat <- floor(nrow(data) * 0.05)
  if (!Formula::is.Formula(formula)) {
    formula <- Formula::as.Formula(formula)
  }

  # for minor options used in kaps algorithm
  minors <- kaps.control(...)
  if (any(K == 1)) 
    stop("the number of subgroups (K) have to be greater than subjects.")
  
  n <- nrow(data) # total number of observations concerning with time 
  rownames(data) <- 1:n
  
  if (n == 0L) stop("0 (non-NA) cases.")
  if (length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
  if (length(K) > 10) stop("the maximum number of subgroups (K) is too large.")
  
  call <- match.call() 
  
	######################################################
	###	Finding optimal subgroups
	type <- match.arg(type)

	if (length(K) == 1 & type == "NULL") {
		result <- kaps.fit(
		  formula = formula, 
		  data = data, 
		  K = K, 
		  mindat = mindat, 
		  minors = minors
		  )
		test.stat2 <- adj.test(result)
		result@call <- call
		result@Z <- test.stat2[1] # overall statistic 
		result@X <- test.stat2[2] # test statistic for the worst pair
		return(result)
	} else if (type == "test") {
		cat("Now, finding optimal number of subgroups (K) by test estimates. \n")
		test.stat <- kaps.test(formula, data, K, minors)
	} else if (type == "perm") {
	
		cat("Now, finding optimal number of subgroups (K) by KAPS-permutation. \n")
		fit <- lapply(K, kaps.fit, formula = formula, data = data, mindat = mindat, minors = minors)
		# 1: overall test statistic
		# 2: overall p-value
		# 3: worst pair test statistic
		# 4: worst pair p-value
		test.stat <- sapply(fit, kaps.perm, permute = TRUE)
		# 1: overall p-value
		# 2: worst pair p-value
		# choose worst pairwise permutation test

		test.stat2 <- sapply(fit, adj.test)
		test.stat2 <- as.matrix(test.stat2)

		index <- 1:ncol(test.stat)
		index.pair <- test.stat[2,] <= 0.05

		if (all(index.pair == FALSE)) {
			# No significant subgroups		
			result <- fit[[1]]
			result@index <- as.integer(0)
			result@groups <- K
			attr(result@groups,"names") <- paste("K<",K,sep="")
		} else {
			index <- index[sum(index.pair)]
			cat("Now, selecting a set of cut-off points...\n")

			result <- fit[[index]]
			result@index <- as.integer(index)
			result@groups <- K
			attr(result@groups,"names") <- paste("K=",K,sep="")
		}

		result@Z <- as.vector(test.stat2[1, ]) #overall statistic 
		result@X <- as.vector(test.stat2[2, ]) #optimal pair p-value
		result@results <- fit # results objects
		result@test.stat <- test.stat # Bonferroni corrected p-values
		result@call <- call
		result@Options <- minors
		return(result)
	
	} else if (type == "NULL") {
		index = 1
	}

	######################################################
	# Obtain log-rank statistics at K
	### parallel computing in order to find optimal k subgroups
  cat("Now, selecting a set of cut-off points...\n")
  fit <- lapply(K, kaps.fit, formula = formula, data = data, mindat= mindat, minors = minors)
  test.stat2 <- sapply(fit, adj.test)
  test.stat2 <- as.matrix(test.stat2)

  result <- fit[[index]]
  result@index <- as.integer(index)
  result@groups <- K
  attr(result@groups,"names") <- paste("K=",K,sep="")
  result@Z <- as.vector(test.stat2[1, ]) # overall test statistic 
  result@X <- as.vector(test.stat2[2, ]) # test statistic for the worst pair
  result@results <- fit # result objects
  result@call <- call
  result@Options <- minors
  return(result)
}

############################################################
# Fit KAPS with test estimates approach 
############################################################
kaps.test <- function(formula, data, K = 2:5, minors = kaps.control()){

	#options(warn = -1)
	if (any(K == 1)) stop("the number of subgroups (K) is greater than 1.")

	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n

  if (n == 0L) stop("0 (non-NA) cases.")
	if (length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	#if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")

	#call = match.call()

	# kaps by test estimates approach
	fold.over.stat <- matrix(0, nrow = 1, ncol = length(K)) 
	fold.pair.stat <- matrix(0, nrow = 1, ncol = length(K))
	result <- matrix(0, nrow = length(K), ncol = 4)
	colnames(result) <- c("over_pval", "over_stat", "pair_pval", "pair_stat")
	rownames(result) <- paste("K=",K,sep="")

	## CHECK ME: reduce computing time by parallel computing
	index <- sample(1:n, floor(n*0.7))
	learning <- data[index,, drop = FALSE]
	test <- data[-index,, drop = FALSE]
	mindat <- floor(nrow(learning) * 0.05)
	
	fit <- lapply(K, kaps.fit, formula = formula, data = learning, mindat = mindat, minors = minors)
	test.stat <- sapply(fit, kaps.perm, newdata = test)			
	rownames(test.stat) <- c("perm_overall_pval", "perm_min_pair_pval")
	colnames(test.stat) <- paste("K=",K,sep="")
	print(round(test.stat,3))

	fold.over.stat[1,] <- test.stat[1,]
	fold.pair.stat[1,] <- test.stat[2,]

	### output
	result[,1] <- fold.over.stat
	#result[,2] <- apply(fold.over.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	result[,3] <- fold.pair.stat
	#result[,4] <- apply(fold.pair.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	#result <- as.data.frame(result)
	return(result)
}
