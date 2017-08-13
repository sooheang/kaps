kaps.fit <- function(formula, data, K, mindat, minors) {
  #######################
  ## pre-processing step
  n <- nrow(data)  
  rownames(data) <- 1:n

  if (n == 0L) stop("0 (non-NA) cases.")
  if (!is.data.frame(data)) 
    data <- as.data.frame(data)
  
  ## Set the object of result
  result <- new("kaps")
  result@formula <- formula
  
  ##################################
  ##### Model Fitting
  X <- Formula::model.part(formula, data = data, rhs = 1, drop = FALSE)
  f <- update(formula, . ~ subgroups)
  
  ## treat pre-determined split points
  if (is.null(attr(minors@pre.pt,"names"))) { 
    pt.set <- lapply(X, function(x) sort(unique(x)) )
  } else {
    pre.name <- colnames(X) != attr(minors@pre.pt, "names")
    pt.set <- lapply(X[,pre.name, drop = FALSE], function(x) sort(unique(x)) )
    pt.set <- c(pt.set, minors@pre.pt)
    X <- X[,attr(pt.set, "names"), drop = FALSE]
	}
  
  ##################################
  ## treat pre-determined ranges
  scope <- minors@scope
	if (!is.null(attr(scope, "names"))) {
		scope.vars <- which(names(pt.set) == attr(scope, "names"))
		pt.set.name <- names(pt.set)
	}
	
	result@X <- 0
	pt.set <- lapply(pt.set, function(x,upper, lower) x[x <= upper & x >= lower], upper = minors@upper.limit, lower = minors@lower.limit)

	####################################
	# NEED FORTRAN ROUTINE UPDATE
	# 1. group.sel() goes to FORTRAN 90 routine
	# 2. gorup.sel() with loop goes to FORTRAN 90 routine 

	for (i in 1: length(pt.set)){
		if (!is.null(attr(scope, "names"))){
			if (i %in% scope.vars){
				rngs <- scope[[pt.set.name[i]]]  
				pt.set[[i]] <- pt.set[[i]][ pt.set[[i]] >= rngs[1] & pt.set[[i]] <= rngs[2] ]
			}
		}

		ord <- order(X[,i])
		data <- data[ord,]
		X <- X[ord, , drop = FALSE]
		test.where <- group.sel(x.vec = X[,i], 
			pt = pt.set[[i]], 
			K = K, 
			mindat = mindat, 
			data = data, 
			f = f, 
			minors = minors)
		x.test <- test.where$test
		
		##### Result
		if (x.test[1,test.where$index] >= result@X) {
			index <- test.where$index
			result@X <- x.test[1,index] # pairwise test statistic X
			result@pair <- x.test[3,index] # pair for selected
			result@groupID <- test.where$where
			result@split.pt <- sapply(unique(result@groupID), function(x,y,where) max(y[where == x]), y = X[,i], where = result@groupID) 
			result@split.pt <- result@split.pt[-length(result@split.pt)]
			result@split.var <- colnames(X[,i,drop =FALSE])
			result@data <- data
		}
	}

	####################################
	### results
	result@groups <- K
	attr(result@groups,"names") <- paste("G=", K, sep="")
	result@mindat <- mindat
	result@Options <- minors
	return(result)
}
