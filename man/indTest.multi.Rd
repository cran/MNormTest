% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MNormTest.R
\name{indTest.multi}
\alias{indTest.multi}
\title{Multivariate Normal Independence Test}
\usage{
indTest.multi(data, subdim = FALSE, alpha = 0.05, verbose = TRUE)
}
\arguments{
\item{data}{The data matrix which is a matrix or data frame. Each column represents a random variable.}

\item{subdim}{The dimensions of submatrices. The default is FALSE, which means the independence of all components of the random vector will be tested.}

\item{alpha}{The significance level. Default is 0.05.}

\item{verbose}{A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.}
}
\value{
An object of class "testResult", which is a list with the following elements:
\item{Conclusion}{The conclusion of the test.}
\item{Stat}{A data frame containing the statistics, p value and critical value.}
\item{SampMean}{The sample mean.}
\item{SampA}{The sample deviation.}
\item{SampAii}{The sample deviation of submatrices.}
\item{Df}{The degree of freedom.}
\item{b}{The Modified factor of the statistic.}
}
\description{
Test whether a set of multivariate normal random vectors are independent. The null hypothesis is "H0: The random vectors are independent of each other".
}
\examples{
data(iris)
chart <- iris[, 1:4]
# carry out the test
test1 <- indTest.multi(chart)
test2 <- indTest.multi(chart, subdim = c(2, 1, 1))
test3 <- indTest.multi(chart, verbose = FALSE)
# get the elements
test1$Stat
test1$SampMean
test2$SampAii
}
\references{
Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.92-94.
}
\author{
Xifeng Zhang
}
