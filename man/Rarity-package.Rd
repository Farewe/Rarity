\name{Rarity-package}
\alias{Rarity-package}
\alias{Rarity}
\docType{package}
\title{
  A package to calculate rarity indices for species and assemblages of species for conservation studies
}
\description{
  This package allows calculation of rarity weights for species and indices of rarity for assemblages of species according to different methods (Leroy et al. 2012, In Press).
  
  The methods developped in this package are based on occurrence (presence-absence) data. Species occurrence is transformed in rarity weights according to various methods. On the basis of the calculated rarity weights, the Index of Relative Rarity can be calculated for assemblages of species (see Leroy et al. 2012 and Leroy et al. in press)
}
\details{
  This package works with two important functions.
  
  First, the function \code{rWeights} is designed to calculate rarity weights according to different weighting function. This package implements the flexible weighting function integrating the rarity cutoff point (Leroy et al. 2012, in press).
  
  Second, the function \code{Irr} is designed to calculate the Index of Relative Rarity for assemblages of species. The Index of Relative Rarity is the average of rarity weights of species in an assemblage, normalized between 0 and 1  (Leroy et al. 2012, in press).
}
\author{
  Boris Leroy
  
  Maintainer: Boris Leroy <leroy.boris@gmail.com>
}
\references{
  Leroy B., Petillon J., Gallon R., Canard A., & Ysnel F. (2012) Improving occurrence-based rarity metrics in conservation studies by including multiple rarity cut-off points. \emph{Insect Conservation and Diversity}, \bold{5}, 159-168.
  
  Leroy B., Canard A., & Ysnel F. In Press. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}, \bold{19}, 794-803.
}

\keyword{ package }

\examples{
  # Example dataset: occurrences of spider species at two different scales
  data(spid.occ)
  head(spid.occ)
  
  # Calculation of rarity weights
  rarity.weights <- rWeights(occData = spid.occ, rCutoff = "Gaston")
  
  # Generation of a random assemblage matrix
  assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 50, replace = TRUE),
                              assemblage.2 = sample(c(0, 1), 50, replace = TRUE),
                              assemblage.3 = sample(c(0, 1), 50, replace = TRUE),
                              assemblage.4 = sample(c(0, 1), 50, replace = TRUE),
                              assemblage.5 = sample(c(0, 1), 50, replace = TRUE))
  # Random attribution of species names
  rownames(assemblages.matrix) <- sample(rownames(spid.occ), 50, replace = FALSE)
  head(assemblages.matrix)
  
  # Calculation of rarity indices of assemblages
  Irr(assemblages = assemblages.matrix, W = rarity.weights)
}