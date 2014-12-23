\name{ppc}
\alias{ppc}
\alias{post.bott}
\title{
  Data to illustrate the posterior predictive checks for the data
  \code{\link{human}}. \code{ppc} and \code{\link{human}} are used to
  illustrate model selection and parameter inference in an ABC framework
  (see the package's vignette for more details).
}
\description{
  \code{data(ppc)} loads in the data frame \code{post.bott}, which
  contains the summary statistics calculated from data simulated a
  posteriori under the bottleneck model (see \code{data(human)} and the
  package's vignette for more details).
}
\usage{data(ppc)}
\format{
  The \code{post.bott} data frame contains the following columns:
  \describe{
    \item{\code{pi}}{
      The mean nucleotide diversity over 50 loci.
    }
    \item{\code{TajD.m}}{
      The mean of Tajima's D statistic over 50 loci.
    }
    \item{\code{TajD.v}}{
      The variance of Tajima's D statistic over 50 loci.
    }
  }
  Each row represents a simulation. 1000 simulations
  were performed under the bottleneck model.
  }
\seealso{\code{\link{human}}}
\keyword{datasets}
% Converted by Sd2Rd version 1.15.
