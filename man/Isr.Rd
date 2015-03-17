\name{Isr}
\alias{Isr}
\alias{Isr,data.frame,data.frame-method}
\alias{Isr,data.frame,matrix-method}
\alias{Isr,data.frame,vector-method}
\alias{Isr,matrix,data.frame-method}
\alias{Isr,matrix,matrix-method}
\alias{Isr,matrix,vector-method}
\alias{Isr,vector,data.frame-method}
\alias{Isr,vector,matrix-method}
\alias{Isr,vector,vector-method}
\title{
Index of Summed Rarity
}
\description{
Calculates the Index of Relative Rarity for a single or a set of assemblages of species. The index of summed rarity consists in summing the rarity weights of species in the assemblage, with the possibility of normalising between 0 (the assemblage has only minimum-weight species) and 1 (value of the complete pool of species (see details)).

}
\usage{
Isr(assemblages, W, abundance = FALSE, Wmin = min(W), normalise = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assemblages}{
\code{vector}, \code{matrix} or \code{data.frame}. Data of the assemblage(s) of species for which the Index will be calculated. If \code{matrix} or \code{data.frame}, species are in rows and assemblages are in columns (see details).
}
  \item{W}{
\code{vector}, \code{matrix} or \code{data.frame}. Rarity weights of the species. If \code{matrix} or \code{data.frame}, species are in rows (see details). 
}
  \item{abundance}{
logical. If \code{TRUE}, then abundance values in \code{assemblages} are taken into account to calculate the Index of Relative Rarity.
}
  \item{Wmin}{
decimal or \code{vector}. Minimum weight(s)
}
  \item{normalise}{
logical. If \code{TRUE}, then the Isr will be rescaled between 0 and 1 (see formula in details). Cannot be used if \code{abundance = TRUE}.
}
}
\details{
\code{assemblages} can contain a single (vector) or several (matrix or data.frame) assemblages. If \code{assemblages} is a vector, values must be named with species names (e.g. with \code{names(assemblages) <- sp.names}). If \code{assemblages} is a matrix or a data.frame, then rows must be named with species names.

All the species of \code{assemblages} must be contained in \code{W}

\code{W} can contain a single vector of rarity weights, or a matrix or a data.frame of rarity weights (e.g. to calculate different Indices of Summed Rarity for different scales). If \code{W} is a vector, values must be named with species names. If \code{W} is a matrix or a data.frame, then species are in rows and each column contains rarity weights. Rows must be named with species names.

If the input of \code{rWeights} is directly entered, then columns \code{Q} (occurrence) and \code{R} (rarity status) will be dropped.

Formula implemented in this function:
\itemize{
\item{Occurrence-based Isr, not normalised between 0 and 1:
\deqn{\sum w_i}{\sum wi}
where \emph{wi} is the weight of the \emph{i}th species in the assemblage
}
\item{Abundance-based Isr, not normalised between 0 and 1:
\deqn{\sum a_i \times w_i}{\sum ai wi}
where \emph{ai} and \emph{wi} are respectively the abundance and weight of the \emph{i}th species in the assemblage
}
\item{Occurrence-based Isr, normalised between 0 and 1:
\deqn{\frac{\sum(w_i - w_{min})}{\sum(W_j - w_{min})}}{\sum(wi - wmin) / \sum(Wj - wmin)}
where \emph{wi} is the weight of the \emph{i}th species in the assemblage, \emph{wmin} the minimum weight and \emph{Wj} the weight of the \emph{j}th species in the regional species pool.
}
}
}
\section{Warning}{
  Column names of \code{W} with names like "Q", "R", "Q1" "Q2" "Q3" and so on, "R1", "R2" "R3" and so on will be ignored.
}
\value{
A data.frame containing the Indices of Relative Rarity and the species richness of the assemblage(s).
}
\references{
Leroy B., Petillon J., Gallon R., Canard A., & Ysnel F. (2012) Improving occurrence-based rarity metrics in conservation studies by including multiple rarity cut-off points. \emph{Insect Conservation and Diversity}, \bold{5}, 159-168.

Leroy B., Canard A., & Ysnel F. 2013. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}, \bold{19}, 794-803.

}
\author{
Boris Leroy leroy.boris@gmail.com
}
\seealso{
\code{\link{rWeights}}, \code{\link{Irr}}
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

# Calculation of Isr
Isr(assemblages.matrix, rarity.weights)
Isr(assemblages.matrix, rarity.weights, normalise = TRUE)


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

# Calculation of Isr
Isr(assemblages.matrix, rarity.weights)
Isr(assemblages.matrix, rarity.weights, normalise = TRUE)
}
