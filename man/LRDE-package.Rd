\name{LRDE-package}
\docType{package}
\alias{LRDE}
\alias{LRDE-package}

\title{LRDE: Differential Expression Analysis for Long-Read RNA-Seq Data}

\description{
The \pkg{LRDE} package provides statistical methods for differential
expression analysis of long-read RNA sequencing (RNA-Seq) data using a
hurdle negative binomial generalized linear model (hurdle-NB GLM).

It implements procedures for:
\itemize{
  \item Estimation of sample-specific size factors for normalization
  \item Modeling zero inflation via group-specific expression probabilities
  \item Gene-wise (tag-wise) dispersion estimation
  \item Statistical testing for differential expression
}

These methods are designed to address key challenges in long-read RNA-Seq data,
including limited sample sizes and excess zero counts (dropout events).
}

\details{
The main functions in this package include:
\itemize{
  \item \code{\link{prepareDGE}}: Prepare count data for analysis. Converts supported
  input types (\code{matrix}, \code{data.frame}, \code{DGEList},
  \code{DESeqDataSet}, and \code{SummarizedExperiment}) into a standardized format.

  \item \code{\link{sizeFactorsEst}}: Estimate sample-specific size factors for normalization.

  \item \code{\link{tagwiseEst}}: Estimate gene-specific (tag-wise) dispersion parameters
  for a hurdle negative binomial model using prior information from bin-level estimates.

  \item \code{\link{hurdle.LRT}}: Perform gene-wise likelihood ratio tests (LRT) for differential expression.

  \item \code{\link{hurdle.Wald.Test}}: Perform gene-wise Wald tests for differential expression.
}

\bold{Typical workflow:}
\enumerate{
  \item Prepare data using \code{prepareDGE}
  \item Normalize counts with \code{sizeFactorsEst}
  \item Estimate tag-wise dispersions using \code{tagwiseEst}
  \item Perform differential expression testing with
        \code{hurdle.LRT} or \code{hurdle.Wald.Test}
}
}

\author{
Ziyang Liu \email{jacobleo773@gmail.com}
}

\seealso{
\code{\link{prepareDGE}},
\code{\link{sizeFactorsEst}},
\code{\link{tagwiseEst}},
\code{\link{hurdle.LRT}},
\code{\link{hurdle.Wald.Test}}
}

\examples{
# Load the package
library(LRDE)

# Simulate count data
set.seed(123)
mat <- matrix(rnbinom(300, size = 5, mu = 5), nrow = 50)
grp <- factor(c("A", "A", "A", "B", "B", "B"))

# Prepare data
y <- prepareDGE(mat, grp)

# Normalize counts
y <- sizeFactorsEst(y)

# Estimate dispersions
y <- tagwiseEst(y)

# Differential expression testing
y <- hurdle.Wald.Test(y)
y <- hurdle.LRT(y)

# Access results
head(y$lrt_stats)
head(y$p.values)
}

\keyword{package}
