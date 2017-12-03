#' Control tuning parameters for "kaps" object
#' 
#' Tuning parameters for an object "kaps"
#' 
#' 
#' @param pre.pt list parameter that treats pre-specified split candidates. Use
#' the option as list(var = split points), i.e., x = 1:100
#' @param scope list parameter that treats pre-determined split range. Use the
#' option as list(var = ranges), i.e., x = c(1,100)
#' @param rho scalar parameter that controls the type of logrank test. See
#' \link[=survdiff]{survdiff}.
#' @param V numeric parameter that determines the number of folds in the
#' cross-validation subgroup selection.
#' @param ncl integer parameter that determines the number of cores to improve
#' computing power
#' @param lower.limit numeric parameter that treats pre-determined overall
#' lower bound. Default is 0.
#' @param upper.limit numeric parameter that treats pre-determined overall
#' upper bound. Default is 100.
#' @param shortcut logical parameter. If shortcut = TRUE, we skip the
#' off-diagonal matrix in pairwise-comparison to reduce computational cost. The
#' default value is TRUE.
#' @param N.perm numeric parameter that gives the number of permutation samples
#' used in the kaps algorithm. The default value is 9999.
#' @param N.boot numeric parameter that gives the number of Bootstrap samples
#' used in the bootstrap and permuting kaps algorithm. The default value is
#' 200.
#' @param alpha numeric parameter that provides a significant level in the
#' process of Bootstrap and permuting algorithm.
#' @param splits character parameter that determines the kind of pairwise test.
#' Default is logrank test. At this stage, the option \code{exact} is not
#' working.
#' @param correct character parameter to select the criteria for the multiple
#' comparison in the simple permuting kaps algorithm.
#' @param p.adjust.methods character parameter to select the criteria for the
#' multiple comparison.
#' @seealso \code{\link{kaps}}
#' @keywords kaps
kaps.control <- function(pre.pt = list(), scope = list(), rho = 0, V = 5, 
  ncl = 1, lower.limit = 0, upper.limit = 100, shortcut = TRUE, N.perm = 9999, 
  N.boot = 200, alpha = 0.05, splits = c("logrank", "exact"), 
  correct = c("Adj.Bf", "Bf", "None"), 
  p.adjust.methods = c("none", "holm", "hochberg", "hommel", 
                       "bonferroni",  "BH", "BY", "fdr")){
  p.adjust.methods <- match.arg(p.adjust.methods)
  splits <- match.arg(splits)
  correct <- match.arg(correct)
  res <- new("kapsOptions")
  res@V <- V
  res@lower.limit <- lower.limit
  res@upper.limit <- upper.limit
  res@N.perm <- N.perm
  res@N.boot <- N.boot
  res@alpha <- alpha
  res@pre.pt <- pre.pt
  res@scope <- scope
  res@rho <- rho
  res@shortcut <- shortcut
  res@splits <- splits
  res@ncl <- as.integer(ncl)
  res@correct <- correct
  res@p.adjust.methods <- p.adjust.methods
  return(res)
}
