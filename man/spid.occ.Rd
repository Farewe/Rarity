\name{spid.occ}
\alias{spid.occ}
\docType{data}
\title{
Occurrences of spiders of Western France at two different scales.
}
\description{
\code{spid.occ} is a \code{data.frame} providing the occurrences of the 708 spider species of Western France at 2 different scales.
}
\usage{data(spid.occ)}
\format{
  A data frame with 708 rows (species) and 2 columns.
  \describe{
    \item{\code{occurMA}}{a numeric vector. Occurrence of spider species at the regional scale (Western France).}
    \item{\code{occurWP}}{a numeric vector. Occurrence of spider species at the western Palearctic scale.}
  }
}
\details{
The regional occurrence is obtained from our spider database of Western France. It was defined as the number of grid cells (16 x 16 cells grid, approximately 683 km^2 per cell) in which the species were found. The western Palearctic occurrence of species was defined as the number of "biogeopolitical" units in which species were found (Canard, 2005; taxonomy and occurrence of species updated in 2011). This estimation comes from the catalogue Islands/groups of islands were separated from continental countries for this estimation. See Leroy et al. (2013 - refs below) for a detailed description of this dataset. Note that the regional occurrence of species may have very slightly changed since the publication of the article because of small additions or modifications to our regional database.
}
\source{
Leroy B., Canard A. & Ysnel F. \emph{2013}. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}, \bold{19}, 794-803.
}
\references{
Canard, A. (2005) Catalogue of spider species from Europe and the Mediterranean basin. \emph{Revue Arachnologique}, \bold{15}, 1-408.

Petillon, J., Courtial, C., Canard, A. & Ysnel, F. (2007a) First assessment of spider rarity in Western France. \emph{Revista iberica de aracnologia}, \bold{15}, 105-113.

Ysnel, F., Petillon, J., Gerard, E. & Canard, A. (2008) Assessing the conservation value of the spider fauna across the West Palearctic area. \emph{Journal of Arachnology}, \bold{36}, 457-463.

Leroy B., Canard A., & Ysnel F. 2013. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}, \bold{19}, 794-803.

}
\examples{
data(spid.occ)
head(spid.occ)
# Congruency between regional and western Palearctic occurrences
corPlot(spid.occ, method = "pearson", xlab = "Regional occurrence", 
        ylab = "Western Palearctic occurrence") # Shows the moderate congruency between regional and
                                                # biogeographic estimations of occurrence
# Calculation of multi-scale rarity weights      
rarity.weights <- rWeights(spid.occ, extended = TRUE)
head(rarity.weights)
}
\keyword{datasets}
