%% kaps Version History

\name{NEWS}
\title{kaps News}
\encoding{UTF-8}

% --------------------------------------------------------------------------- %
\section{Changes in version 1.1.4 2017-12-02}{
  \itemize{
    \item Correct a bug: src/*.o file format recognising on Windows  
  }
}

% --------------------------------------------------------------------------- %
\section{Changes in version 1.1.3 2017-08-13}{
  \itemize{
    \item Update appVeyor.yml for Windows CI.
  }
}

% --------------------------------------------------------------------------- %
\section{Changes in version 1.1.2 2017-08-05}{
  \itemize{
    \item Add registerDynamicSymbols.c
    \item Aad unit.test rountine based on the testthat package.
    \item Update survival function routine using coin::logrank_test().
    \item Bug fixed and performance improved
  }
}
% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.3 2015-08-22}{
  \itemize{
    \item Update the NAMESPACE to avoid masking problem
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.2 2014-11-01}{
  \itemize{
    \item Imrpove man pages 
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.0 2014-03-17}{
  \itemize{
    \item Release an official version
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.26 2013-12-15}{
  \itemize{
    \item Bug fixed in simple permuting algorithm
    \item Adjust a minimum sample size rule to improve end-cut preference
    \item Add validation methods based on the Fisher's combined probability test 
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.25 2013-11-20}{
  \itemize{
    \item Add correction methods for simple permuting kaps()
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.24 2013-11-18}{
  \itemize{
    \item add simple permuting KAPS algorithm
    \item add update() easily to update the criteria for the selection of optimal number of subgroups
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.23 2013-11-09}{
  \itemize{
    \item add new decision rule, random permuting kaps algorithm, \code{type = "random"}, to find optimal subgroups
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.22 2013-11-07}{
  \itemize{
    \item add new parameter, \code{alpha}, to provide the standard of significance in the process of Bootstrap and permutations
    }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.21 2013-11-05}{
  \itemize{
    \item add new decision rules called kaps with simple permutation to determine the number of optimal subgroups
    \item add new comparison rules, i.e., mean (\code{mean}), median (\code{median}), censoring-related trimmed mean (\code{trim}), and test-based approah (\code{test}), for the bootstrap and permuting kaps \code{kaps.conrtrol(boot.sel)}
    }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.20 2013-11-03}{
  \itemize{
    \item bugfix in function the selection of optimal subgroups \code{kaps.boot()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.19 2013-10-30}{
  \itemize{
    \item bugfix in the parallel Bootstrap algorithm \code{kaps.boot()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.18 2013-10-20}{
  \itemize{
    \item consider new routine when subgroups is not significant
    }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.17 2013-10-03}{
  \itemize{
    \item added support for the \pkg{foreach} and \pkg{doMC} packages to utilize parallel computing procedure
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.15 2013-10-01}{
  \itemize{
    \item removed routines for the tree-structured model, \code{lrtree()}
    \item log-rank tree, \code{lrtree()} moves a new package \pkg{SurvTree}
    \item added new permutation based rules for the identification of optimal subgroups
  }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.9.14 2013-08-31}{
  \itemize{
    \item changed S4 class from \code{apss} to \code{kaps}
    \item added new decision rules to find the number of optimal subgroups based on the permutation test
  }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 0.7.0 2012-07-15}{
  \itemize{
    \item published the beta version of the kaps algorithm on CRAN
  }
}



