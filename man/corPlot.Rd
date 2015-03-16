\name{corPlot}
\alias{corPlot}
\title{
Correlation plots
}
\description{
Produces a scatter plot with the distribution of points according to the tested variables in the lower left triangle and the correlation values in the upper left triangle.
}
\usage{
corPlot(df, method = "spearman", digits = 2, na.action = "keep", ties.method = "average",
        title = "", xlab = "variable.name", ylab = "variable.name", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
a matrix or a data.frame. The values of variables (e.g., indices) to be compared are in columns. 
}
  \item{method}{
\code{pearson}, \code{spearman} or \code{kendall}. This is the method to be used to produce the plot, according to \code{\link[stats:cor.test]{cor.test}}. See details.
}
  \item{digits}{
Number of digits to round the correlation values on the plot.
}
  \item{na.action}{
for controlling the treatment of NAs in \code{spearman} or \code{kendall} plots. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed; if "keep" they are kept with rank NA. See \code{\link[base:rank]{rank}}.
}
  \item{ties.method}{
\code{average}, \code{first}, \code{random}, \code{max}, or \code{min}; a character string specifying how ties are treated in \code{spearman} or \code{kendall} plots. See \code{\link[base:rank]{rank}} for details.
}
  \item{title}{
Title of the plot.
}
  \item{xlab}{
a character string for labelling x axes. \code{variable.name} (default value) will produce automatic labelling according to column names of \code{df}. Otherwise, either a single string can be entered, or a vector of strings of length equal to the number of columns of \code{df}.
}
  \item{ylab}{
a character string for labelling y axes. \code{variable.name} (default value) will produce automatic labelling according to column names of \code{df}. Otherwise, either a single string can be entered, or a vector of strings of length equal to the number of columns of \code{df}.
}
  \item{...}{
Further arguments to be passed to the individual plots. See \code{\link[graphics:plot]{plot}} and \code{\link[graphics:par]{par}}
}
}
\details{
The lower half shows the scatter plots of values or ranks of variables. The upper half shows the corresponding correlation coefficients (significativity: 0 '***' 0.001 '**' 0.01 '*' 0.05 '-' 0.1 ' ' 1). The diagonal shows the considered variables and the number of individuals available for each.

If the chosen method is \code{pearson}, then the actual values of the variables will be plotted. If the chosen method is a rank-based method, \code{spearman} or \code{kendall}, then the ranks will be plotted.

}
\section{Warning}{
A high number of variables will likely result in a slow generation of plots and a poor readability.
Above 10 variables, the readability is greatly reduced.
}
\author{
Boris Leroy leroy.boris@gmail.com
}

\seealso{
\code{\link[arm:corrplot]{corrplot}} in package \code{arm}
}
\examples{
# Comparisons of species occurrences estimated from 2 different scales
data(spid.occ)
corPlot(spid.occ, method = "pearson")

# Another example:
# Correlation between different variables measured on the same individuals
data(iris)
corPlot(iris[, 1:4], method = "pearson")
corPlot(iris[, 1:4], method = "spearman")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
