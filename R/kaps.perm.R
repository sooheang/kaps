############################################################
# Fit KAPS with Permutation test by test estimates approach
# This algorithm select best subgroups by permutation test
# in which the data is randomly divided into 
# learning and validation set
kaps.perm <- function(fit, newdata, permute = TRUE){
  #if(!inherits(fit, "kaps")) stop("This function requires the object of kaps class as a main object.")
  if(missing(newdata)){
    dat.perm = fit@data
    dat.perm$subgroups <- factor(fit@groupID)
  } else{
    if(!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
    dat.perm <-	predict(fit, newdata)
    dat.perm$subgroups <- factor(dat.perm$Group)
  }
  res.perm = matrix(1, nrow = 1, ncol = 2)
  
  if(length(levels(dat.perm$subgroups)) == 1){
    return(res.perm)
  } else{
    Ku = fit@groups
    Kul <- length(Ku)
    
    
    formula = fit@formula
    minors = fit@Options
    
    formula1 = update(formula, ~ subgroups)
    formula1 <- as.formula(formula1)
    
    # Set storages
    perm.over.pval = perm.pair.stat = perm.pair.pval = c()
    res.perm = matrix(0, nrow = 1, ncol = 2)
    #res.perm <- as.data.frame(res.perm)
    #colnames(res.perm) <- c("overall.stat", "overall.pval", "pair.stat", "pair.pval")
    #colnames(res.perm) <- c( "perm_overall_pval", "perm_min_pair_pval")
    ######################################################################
    # GOD by permutation test with asymptotic distribution
    # require(coin)
    # cat("Now, permutation kaps is working. Please, wait a minute.^^\n")
    # 1. overall permuatation test statistic
    if(permute){
      tmp <- logrank_test(formula = formula1, data = dat.perm, distribution = approximate(B = minors@N.perm))
      
      # Naive approach
      #tmp <- clinfun::permlogrank(formula1, dat.perm)
      #perm.over.stat <- statistic(tmp)
      perm.over.pval <- switch(minors@correct,
                               Bf = pvalue(tmp) * Kul * (Kul-1) * 0.5,
                               Adj.Bf = pvalue(tmp) * Kul,
                               None = pvalue(tmp)
      )
      perm.over.pval <- ifelse(perm.over.pval > 1, 1, perm.over.pval)
      #perm.over.pval <- pvalue(tmp) * Ku * (Ku-1) * 0.5
      #names(perm.over.stat) <- "overall_statistic"
      #surv_test(formula = Surv(Y,cens) ~ Group, data = dat.perm, distribution = exact())
    } else{
      tmp <- survdiff(formula = formula, data = dat.perm, rho = minors@rho)
      perm.over.pval <- 1 - pchisq(tmp$chisq, df = Ku - 1)
      perm.over.pval <- switch(minors@correct,
                               Bf = perm.over.pval * Kul * (Kul-1) * 0.5,
                               Adj.Bf = perm.over.pval * Kul,
                               None = perm.over.pval
      )
      perm.over.pval <- ifelse(perm.over.pval > 1, 1, perm.over.pval)
    }
    
    
    # 2. pairwise permutation test statistic
    #pair = combn(Ku, 2)
    pair = combnat(Ku, 2)
    if(minors@shortcut) pair <- pair[, !duplicated(pair[1,]), drop = FALSE]
    
    for(i in 1:ncol(pair)){
      dat.perm.tmp = dat.perm[dat.perm$subgroups %in% pair[,i],]
      #if(length(levels(dat.perm.tmp$Group)) <= 1) {
      #	perm.pair.stat <- c(perm.pair.stat, 100)
      #	perm.pair.pval <- c(perm.pair.pval, 0)
      #} else if(nrow(dat.perm.tmp) == 0){
      #	perm.pair.stat <- c(perm.pair.stat, 100)
      #	perm.pair.pval <- c(perm.pair.pval, 0)
      #} else if(all(table(dat.perm.tmp$Group) >= fit@mindat)){
      #	perm.pair.stat <- c(perm.pair.stat, 100)
      #	perm.pair.pval <- c(perm.pair.pval, 0)
      #}	else {
      #	tmp <- surv_test(formula = formula1, data = dat.perm.tmp, distribution = approximate(B = minors@N.perm))
      #	perm.pair.stat <- c(perm.pair.stat, statistic(tmp))
      #	perm.pair.pval <- c(perm.pair.pval, pvalue(tmp))
      #}
      if(permute){
        tmp <- try(logrank_test(
          formula = formula1, 
          data = dat.perm.tmp, 
          distribution = approximate(B = minors@N.perm)
        ), 
        silent = TRUE
        )
        
        if(class(tmp) != "try-error"){
          perm.pair.stat <- c(perm.pair.stat, statistic(tmp))
          pair.pval <- switch(minors@correct,
                              Bf = pvalue(tmp) * Kul * (Kul-1) * 0.5,
                              Adj.Bf = pvalue(tmp) * Kul,
                              None = pvalue(tmp)
          )
          perm.pair.pval <- c(perm.pair.pval, pair.pval)
        } 
      } else{
        tmp <- try(survdiff(formula = formula1, data = dat.perm.tmp, rho = minors@rho), silent = TRUE)
        
        if(class(tmp) != "try-error"){
          perm.pair.stat <- c(perm.pair.stat, tmp$chisq)
          pair.pval <- 1 - pchisq(tmp$chisq, 1)
          pair.pval <- switch(minors@correct,
                              Bf = pair.pval * Kul * (Kul-1) * 0.5,
                              Adj.Bf = pair.pval * Kul,
                              None = pair.pval
          )
          perm.pair.pval <- c(perm.pair.pval, pair.pval)
        } 
      }
    }
    perm.pair.pval <- ifelse(perm.pair.pval > 1, 1, perm.pair.pval)
    #names(perm.pair.stat) <- apply(pair, 2, function(x) paste(x, collapse = "&"))
    
    ##########
    # results
    min.ind = which.min(perm.pair.stat)
    
    #res.perm[1,1] <- perm.over.stat
    res.perm[1,1] <- perm.over.pval
    #res.perm[1,3] <- perm.pair.stat[min.ind]
    res.perm[1,2] <- ifelse(length(min.ind) == 0, 1.00, perm.pair.pval[min.ind])
    return(res.perm)
  }
}
