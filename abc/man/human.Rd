\name{human}
\alias{human}
\alias{stat.voight}
\alias{par.italy.sim}
\alias{stat.3pops.sim}
\alias{models}
\title{
A set of R objects containing observed data from three human
populations, and simulated data under three different demographic
models. The data set is used to illustrate model selection and parameter
inference in an ABC framework (see the package's vignette for more details).
}
\description{
  \code{data(human)} loads in four R objects: \code{stat.voight} is a
  data frame with 3 rows and 3 columns and contains the observed summary
  statistics for three human populations, \code{stat.3pops.sim} is also a
  data frame with 150,000 rows and 3 columns and contains the simulated
  summary statistics, \code{models} is a vector of character strings of
  length 150,000 and contains the model indices, \code{par.italy.sim} is a
  data frame with 50,000 rows and 4 columns and contains the parameter
  values that were used to simulate data under a population bottleneck
  model. The corresponding summary statistics can be subsetted from the
  \code{stat.3pops.sim} object as \code{subset(stat.3pops.sim,
  subset=models=="bott")}.
}
\usage{data(human)}
\format{
  The \code{stat.voight} data frame contains the following columns:
  \describe{
    \item{\code{pi}}{
      The mean nucleotide diversity over 50 loci in 3 human
      populations, Hausa, Italian, and Chinese.
    }
    \item{\code{TajD.m}}{
      The mean of Tajima's D statistic over 50 loci in 3 human
      populations, Hausa, Italian, and Chinese.
    }
    \item{\code{TajD.v}}{
      The variance of Tajima's D statistic over 50 loci in 3 human
      populations, Hausa, Italian, and Chinese.
    }
  }
  Each row represents a simulation. Under each model 50,000 simulations
  were performed. Row names indicate the type of demographic model.
  
  The \code{stat.3pops.sim} data frame contains the following columns:
  \describe{
    \item{\code{pi}}{
      The mean of nucleotide diversity over 50 simulated loci under 3
      demographic scenarios: constant size population, population
      bottleneck, and population expansion.
    }
    \item{\code{TajD.m}}{
      The mean of Tajima's D statistic over 50 simulated loci under 3
      demographic scenarios: constant size population, population
      bottleneck, and population expansion.
    }
    \item{\code{TajD.v}}{
      The variance of Tajima's D statistic over 50 simulated loci under
      3 demographic scenarios: constant size population, population
      bottleneck, and population expansion.
    }
  }
  Each row represents a simulation. Under each model 50,000 simulations
  were performed. Row names indicate the type of demographic model.

  The \code{par.italy.sim} data frame contains the following columns:
  \describe{
    \item{\code{Ne}}{
      The effective population size.
    }
    \item{\code{a}}{
      The intensity of the bottleneck (i.e. the ratio of the population sizes
      before and during the bottleneck).
    }
    \item{\code{duration}}{
      The duration of the bottleneck.
    }
    \item{\code{start}}{
      The start of the bottleneck.
    }
  }
  Each row represents a simulation. 
  
  \code{models} contains the names of the demographic models.  
}

\details{
  Data is provided to estimate the posterior probabilities of classical
  demographic scenarios in three human populations: Hausa, Italian, and
  Chinese. These three populations represent the three continents:
  Africa, Europe, Asia, respectively. \code{par.italy.sim} may then used
  to estimate the ancestral population size of the European population
  assuming a bottleneck model.
  
  It is generally believed that African human populations are expanding,
  while human populations from outside of Africa have gone through a
  population bottleneck. Tajima's D statistic has been classically used
  to detect changes in historical population size. A negative Tajima's D
  signifies an excess of low frequency polymorphisms, indicating
  population size expansion. While a positive Tajima's D indicates low
  levels of both low and high frequency polymorphisms, thus a sign of a
  population bottleneck. In constant size populations, Tajima's D is
  expected to be zero.

  With the help of the \code{human} data one can reach these expected
  conclusions for the three human population samples, in accordance with
  the conclusions of Voight et al. (2005) (where the observed statistics
  was taken from), but using ABC.
}
\source{
  The observed statistics were taken from Voight et al. 2005 (Table
  1.). Also, the same input parameters were used as in Voight et al. 2005 to
  simulate data under the three demographic models. Simulations were
  performed using the software ms and the summary statistics were
  calculated using sample_stats (Hudson 1983).
}
\references{
  B. F. Voight, A. M. Adams, L. A. Frisse, Y. Qian, R. R. Hudson and
  A. Di Rienzo (2005) Interrogating multiple aspects of variation in a
  full resequencing data set to infer human population size
  changes. PNAS \bold{102}, 18508-18513.

  Hudson, R. R. (2002) Generating samples under a Wright-Fisher neutral
  model of genetic variation. Bioinformatics \bold{18} 337-338.
}
\seealso{\code{\link{ppc}}}
\keyword{datasets}
% Converted by Sd2Rd version 1.15.
