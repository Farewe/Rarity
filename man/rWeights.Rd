\name{rWeights}
\alias{rWeights}
\alias{rWeights-method}
\alias{rWeights,data.frame-method}
\alias{rWeights,matrix-method}
\alias{rWeights,vector-method}
\title{
Calculation of rarity weights at a single or at multiple scales
}
\description{
Calculate rarity weights for a single scale or for multiple scales on the basis of the selected weighting function(s).
}
\usage{
rWeights(occData, Qmax = max(occData), Qmin = min(occData), 
         wMethods = "W", rCutoff = "Gaston", normalized = T, 
         assemblages, extended = F, rounding = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{occData}{
vector, matrix or data.frame. Occurrence data for a single scale (vector) or several scales (matrix or data.frame). 
}
  \item{Qmax}{
integer. Maximum occurrence (see details). By default, the maximum occurrence of the dataset is used (i.e., maximum occurrence among the provided set of species), however it can be changed to another value, e.g. to provide the number of possible sites. 
}
  \item{Qmin}{
integer. Minimum occurrence (see details). By default, the minimum occurrence of the dataset is used (i.e., minimum occurrence among the provided set of species).
}
  \item{wMethods}{
\code{W}, \code{invQ} or \code{oldW}. Chosen function to weight species occurrence (see details)
}
  \item{rCutoff}{
a decimal or a vector of values between 0 and 1, or \code{Gaston} or \code{Leroy}. Indicates the rarity cutoff(s) or the method to use to calculate the rarity cutoff(s). The rarity cut-off is by default calculated as a percentage of the maximum occurrence  (see details)
}
  \item{normalized}{
TRUE or FALSE. If TRUE, then weights are normalized between 0 and 1.
}
  \item{assemblages}{
matrix or data.frame. Set of assemblages of species to calculate the rarity cutoff point(s) with the \code{Leroy}  method (optional)
}
  \item{extended}{
TRUE or FALSE. Useful in case of multiple scales only. If TRUE, then weights will be given for every input scale in addition to multiscale weights. If FALSE, then only multiscale weights will be provided. 
}
  \item{rounding}{
An integer or FALSE. If an integer is provided, then the values of weights will be rounded according to this value. If FALSE, weights will not be rounded.
}
}
\details{
NA values are not (yet) allowed (remove them from the dataset before running the function).



To calculate single-scale weights, simply provide a vector with species occurrences.
To calculate multiscale rarity weights, provide either a matrix or a data.frame where species are in rows, and each column provides occurrence for a particular scale.



The minimum and maximum weights can be set manually, or automatically calculated with the default parameters. Defaults parameters : if \code{occData} is a vector, \code{Qmin = min(Q)} and \code{Qmax = max(Q)}. If \code{occData} is a matrix or a data.frame, \code{Qmin = apply(occData, 2, min)} and \code{Qmax = apply(occData, 2, max)}



Three weighting methods are available (more will become available later):

\enumerate{
\item{\code{W}: This is the method described in Leroy et al. (2013). We recommend using this method for both single and multiscale weight calculations.
\deqn{\exp(-(\frac{Q_{i} - Q_{min}}{r_j \times Q_{max} - Q_{min}}\times0.97 + 1.05)^2)}{exp(-(((Qi - Qmin)/(r * Qmax - Qmin)) * 0.97 + 1.05)^2)}
where \emph{Qi} is the occurrence of species \emph{i}, Qmin and Qmax are respectively the minimum and maximum occurrences in the species pool and \emph{r} is the choosen rarit cut-off point (as a percentage of occurrence).
}
\item{\code{invQ}: This is the inverse of occurrence
\deqn{\frac{1}{Q_i}}{1/Qi}
where \emph{Qi} is the occurrence of the \emph{i}th species. The inverse of the occurrence should be avoided as a weighting procedure because it cannot be adjusted to the considered species pool, and it does not attribute 0 weights to common species (see discussion in Leroy et al. (2012)).
}
\item{\code{oldW}: This is the original method described in Leroy et al. (2012). As this method was improved in Leroy et al. (2013), we recommend to rather use \code{W}. Formula:
\deqn{\exp(-(\frac{Q_i}{Q_{min}} n + 1)^2)}{exp(-((Qi/Qmin) * n + 1)^2)} 
where \emph{Qi} is the occurrence of species \emph{i}, Qmin is the minimum occurrence in the species pool, and \emph{n} is and adjustment coefficient numerically approximated to fit the choosen rarity cut-off point.
}
}



For methods \code{W} and \code{oldW}, a rarity cutoff point is required. The rarity cutoff point can either be entered manually (a single value for a single scale, a vector of values for multiple scales), or the methods of \code{Gaston} or \code{Leroy} can be used (see references):

 - \code{Gaston} method: the rarity cutoff point is the first quartile of species occurrences, i.e. rare species are the 25 percent species with the lowest occurrence.
 
 - \code{Leroy} method: the rarity cutoff point is the occurrence at which the average proportion of rare species in local assemblages is 25 percent. This method requires \code{assemblages} to calculate the average proportion of rare species in  assemblages.
}
\value{
A data.frame containing the results : species occurrences, rarity statuses, rarity weights and the used rarity cut-offs.

-  If \code{occData} is a vector (single scale weights): A data.frame with 4 columns : \code{Q} (species occurrence), \code{R} (species rarity status), \code{W}, (species rarity weights), \code{cut.off} (rarity cut-off used for weight calculation)

-	If \code{occData} is matrix or a data.frame (multiscale rarity weights): A data.frame with \emph{n} columns \code{Q} (species occurrences), \emph{n} columns \code{R} (species rarity statuses), one (if \code{extended = F}) or \emph{n} + 1 (if \code{extended = T}) columns \code{W} (species rarity weights) where \emph{n} is the number of scales (number of columns of \code{occData}), \emph{n} columns \code{cut.off} (rarity cut-offs used for weight calculation).

By default, weights are rounded to 3 digits, which should be sufficient in most cases. Another number of digits can also be chosen; or simply changing \code{rounding} to \code{FALSE} will remove the rounding.
}
\references{
Leroy B., Petillon J., Gallon R., Canard A., & Ysnel F. (2012) Improving occurrence-based rarity metrics in conservation studies by including multiple rarity cut-off points. \emph{Insect Conservation and Diversity}, \bold{5}, 159-168.

Leroy B., Canard A., & Ysnel F. 2013. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}, \bold{19}, 794-803.

}
\author{
Boris Leroy leroy.boris@gmail.com
}

\seealso{
\code{\link{Irr}}, \code{\link{Isr}}
}
\examples{
# 1. Single scale rarity weights
data(spid.occ)
head(spid.occ)

regional.occ <- spid.occ$occurMA
names(regional.occ) <- rownames(spid.occ)
head(regional.occ)

# Calculation of rarity weights at a single scale:
rWeights(regional.occ, rCutoff = "Gaston")
rWeights(regional.occ, rCutoff = 0.1)
rWeights(regional.occ, wMethods = "invQ")
rWeights(regional.occ, wMethods = c("W", "invQ"))

# Calculation of rarity weights with the method of Leroy
# Creating a fictive assemblage matrix of 5 assemblages
# Warning: this is to provide an example of how the function works!
# The correct use of this method requires a matrix of actually sampled species.
assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 708, replace = TRUE),
                            assemblage.2 = sample(c(0, 1), 708, replace = TRUE),
                            assemblage.3 = sample(c(0, 1), 708, replace = TRUE),
                            assemblage.4 = sample(c(0, 1), 708, replace = TRUE),
                            assemblage.5 = sample(c(0, 1), 708, replace = TRUE))
rownames(assemblages.matrix) <- names(regional.occ) # Rownames of assemblages.matrix must 
                                                 # correspond to rownames in occurrences
head(assemblages.matrix)
                                          
rWeights(regional.occ, wMethods = "W", rCutoff = "Leroy", assemblages = assemblages.matrix)

# 2. Multiscale rarity weights
data(spid.occ)
head(spid.occ)

rWeights(spid.occ, wMethods = "W", rCutoff = "Gaston")
rWeights(spid.occ, wMethods = "W", rCutoff = "Gaston", extended = TRUE)
rWeights(spid.occ, wMethods = c("W", "invQ"), rCutoff = "Gaston", extended = TRUE)
rWeights(spid.occ, wMethods = c("W", "invQ"), rCutoff = "Leroy", 
         assemblages = assemblages.matrix, extended = TRUE) # Provided that you have 
                                             # created "assemblages.matrix" as above
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
