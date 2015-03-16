\name{Irr}
\alias{Irr}
\alias{Irr,data.frame,data.frame-method}
\alias{Irr,data.frame,matrix-method}
\alias{Irr,data.frame,vector-method}
\alias{Irr,matrix,data.frame-method}
\alias{Irr,matrix,matrix-method}
\alias{Irr,matrix,vector-method}
\alias{Irr,vector,data.frame-method}
\alias{Irr,vector,matrix-method}
\alias{Irr,vector,vector-method}
\title{
Index of Relative Rarity
}
\description{
Calculate the Index of Relative Rarity for a single or a set of assemblages of species.
}
\usage{
Irr(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assemblages}{
vector, matrix or data.frame. Data of the assemblage(s) of species for which the Index will be calculated. If matrix or data.frame, species are in rows and assemblages are in columns (see details).
}
  \item{W}{
vector, matrix or data.frame. Rarity weights of the species. If matrix or data.frame, species are in rows (see details). 
}
  \item{abundance}{
logical. If \code{TRUE}, then abundance values in \code{assemblages} are taken into account to calculate the Index of Relative Rarity.
}
  \item{Wmin}{
decimal or vector. Minimum weight(s)
}
  \item{Wmax}{
decimal or vector. Maximum weight(s)
}
}
\details{
\code{assemblages} can contain a single (vector) or several (matrix or data.frame) assemblages. If \code{assemblages} is a vector, values must be named with species names (e.g. with \code{names(assemblages) <- sp.names}). If \code{assemblages} is a matrix or a data.frame, then rows must be named with species names.

All the species of \code{assemblages} must be contained in \code{W}

\code{W} can contain a single vector of rarity weights, or a matrix or a data.frame of rarity weights (e.g. to calculate different Indices of Relative Rarity for different scales). If \code{W} is a vector, values must be named with species names. If \code{W} is a matrix or a data.frame, then species are in rows and each column contains rarity weights. Rows must be named with species names.

If the input of \code{rWeights} is directly entered, then columns \code{Q} (occurrence) and \code{R} (rarity status) will be dropped.
}
\value{
A data.frame containing the Indices of Relative Rarity and the species richness of the assemblage(s).
}
\references{
Leroy B., Petillon J., Gallon R., Canard A., & Ysnel F. (2012) Improving occurrence-based rarity metrics in conservation studies by including multiple rarity cut-off points. \emph{Insect Conservation and Diversity}, \bold{5}, 159-168.

Leroy B., Canard A., & Ysnel F. In Press. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}, \bold{19}, 794-803.

}
\author{
Boris Leroy leroy.boris@gmail.com
}

\seealso{
\code{\link{rWeights}}
}
\examples{
# Input rarity weights
data(spid.occ)

# Example of a single scale dataset
regional.occ <- spid.occ$occurMA
names(regional.occ) <- rownames(spid.occ)
head(regional.occ)

# Preparation of rarity weights
rarity.weights <- rWeights(regional.occ)

# Generation of an assemblage matrix
assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.2 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.3 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.4 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.5 = sample(c(0, 1), 50, replace = TRUE))
# Random attribution of names to the sampled species
rownames(assemblages.matrix) <- sample(names(regional.occ), 50, replace = FALSE)
head(assemblages.matrix)

# Calculation of Irr
Irr(assemblages.matrix, rarity.weights)


# Example of a multi scale dataset
rarity.weights <- rWeights(spid.occ, extended = TRUE)
head(rarity.weights)

# Generation of an assemblage matrix
assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.2 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.3 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.4 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.5 = sample(c(0, 1), 50, replace = TRUE))
rownames(assemblages.matrix) <- sample(names(regional.occ), 50, replace = FALSE)
head(assemblages.matrix)

# Calculation of Irr
Irr(assemblages.matrix, rarity.weights)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
