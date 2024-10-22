#' @title Single Mean Vector Hypothesis Testing
#' @description Test whether the mean vector of a single multivariate normal population is equal to a certain value when the population covariance matrix is known or unknown. The null hypothesis is that "H0: mu = mu0".
#' @author Xifeng Zhang
#' @param data The data matrix which is a matrix or data frame.
#' @param mu0 The mean vector when the null hypothesis is true.
#' @param Sigma0 The population covariance matrix. Default is FALSE which means the covariance matrix is unknown.
#' @param alpha The significance level. Default is 0.05.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.66-68.
#' @return An object of class "testResult", which is a list with the following elements:
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMean}{The sample mean.}
#' \item{SampA}{The sample deviation.}
#' \item{Df}{The degree of freedom.}
#' @export
#'
#' @examples
#' data(iris)
#' X <- iris[, 1:4]
#' mu0 <- c(5.8, 3.0, 4.3, 1.3)
#' # carry out the test
#' test1 <- meanTest.single(X, mu0)
#' test2 <- meanTest.single(X, mu0, Sigma0 = diag(1, 4))
#' test3 <- meanTest.single(X, mu0, verbose = FALSE)
#' # get the elements
#' test1$Stat
#' test1$SampMean
#' test1$SampA
#' test1$Df
meanTest.single <- function(data, mu0, Sigma0 = FALSE, alpha = 0.05, verbose = TRUE) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- apply(data, 2, mean)

  if (isFALSE(Sigma0)) {
    A <- (n - 1) * cov(data)
    T2 <- (n - 1) * n * t(X.bar - mu0) %*% solve(A) %*% (X.bar - mu0)
    F <- (n - p) / ((n - 1) * p) * T2

    p.value <- min(pf(F, p, n - p), 1 - pf(F, p, n - p))
    critical.value <- qf(1 - alpha, p, n - p)
    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

    statistics <- data.frame(
      "Value" = c(T2, F),
      "p value" = c(" ", p.value),
      "Critical Value" = c(" ", critical.value)
    )
    rownames(statistics) <- c("Hotelling T2", "F")

    dof <- data.frame("df" = c(paste("T2(", p, ",", n - 1, ")"), paste("F(", p, ",", n - p, ")")))
    rownames(dof) <- c("Hotelling T2", "F")

    testResult <- list(Conclusion = reject, Stat = statistics, SampMean = X.bar, SampA = A, Df = dof)

    if (verbose) {
      return(testResult)
    } else {
      message("H0: mu = (", toString(mu0), ")' when Sigma is unknown\n")
      return(testResult)
    }
  } else {
    T0 <- n * t(X.bar - mu0) %*% solve(Sigma0) %*% (X.bar - mu0)
    p.value <- min(pchisq(T0, p), 1 - pchisq(T0, p))
    critical.value <- qchisq(1 - alpha, p)
    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

    statistics <- data.frame(
      "Value" = c(T0),
      "p value" = c(p.value),
      "Critical Value" = c(critical.value)
    )
    rownames(statistics) <- c("Chi2")

    testResult <- list(Conclusion = reject, Stat = statistics, SampMean = X.bar, SampA = Sigma0, Df = paste("Chi2(", p, ")"))

    if (verbose) {
      return(testResult)
    } else {
      message("H0: mu = (", toString(mu0), ")' when Sigma is known\n")
      return(testResult)
    }
  }
}

#' @title Two Mean Vectors Hypothesis Testing
#' @description Test whether the mean vectors of two multivariate normal populations are equal when the covariance matrices are equal or unequal. The null hypothesis is that "H0: mu1 = mu2".
#' @author Xifeng Zhang
#' @param data1 A matrix or data frame of group 1.
#' @param data2 A matrix or data frame of group 2.
#' @param alpha The significance level. Default is 0.05.
#' @param equal A boolean value. Default is TRUE. If TRUE, the covariance matrix is equal. If FALSE, the covariance matrix is not equal.
#' @param method A string value. Default is "None". When equal is FALSE, you must choose a method in "Coupled" or "Transformed". Choose "Coupled" when the sample size of two groups is equal. Choose "Transformed" when the sample size of two groups is not equal.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.76-80.
#' @return An object of class "testResult", which is a list with the following elements:
#' Return when the param equal is TRUE.
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMean1}{The sample mean of group 1.}
#' \item{SampMean2}{The sample mean of group 2.}
#' \item{SampA1}{The sample deviation of group 1.}
#' \item{SampA2}{The sample deviation of group 2.}
#' \item{MixSampA}{The mixed sample deviation.}
#' \item{Df}{The degree of freedom.}
#' Return when the param equal is FALSE and method is "Coupled".
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMeanC}{The sample mean of coupled data.}
#' \item{SampAC}{The sample deviation of coupled data.}
#' \item{Df}{The degree of freedom.}
#' \item{dataC}{The coupled data.}
#' Return when the param equal is FALSE and method is "Transformed".
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMeanT}{The sample mean of transformed data.}
#' \item{SampAT}{The sample deviation of transformed data.}
#' \item{Df}{The degree of freedom.}
#' \item{dataT}{The transformed data. Return when the param equal is FALSE and method is "Transformed".}
#' @export
#'
#' @examples
#' data(iris)
#' X <- iris[1:50, 1:4]
#' Y <- iris[51:100, 1:4]
#' # carry out the test
#' test1 <- meanTest.two(X, Y)
#' test2 <- meanTest.two(X, Y, verbose = TRUE)
#' test3 <- meanTest.two(X, Y, equal = FALSE, method = "Coupled")
#' test4 <- meanTest.two(X, Y, equal = FALSE, method = "Transformed")
#' # get the elements
#' test1$Stat
#' test1$SampMean1
#' test3$SampMeanC
#' test4$dataT
meanTest.two <- function(data1, data2, alpha = 0.05, equal = TRUE, method = c("None", "Coupled", "Transformed"), verbose = TRUE) {
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p <- ncol(data1)

  if (equal) {
    data1 <- as.matrix(data1)
    data2 <- as.matrix(data2)

    X.bar <- apply(data1, 2, mean)
    A1 <- (n1 - 1) * cov(data1)
    Y.bar <- apply(data2, 2, mean)
    A2 <- (n2 - 1) * cov(data2)
    A <- (A1 + A2) / (n1 + n2 - 2)

    T2 <- (n1 * n2 / (n1 + n2)) * t(X.bar - Y.bar) %*% solve(A) %*% (X.bar - Y.bar)
    F <- (n1 + n2 - p - 1) / ((n1 + n2 - 2) * p) * T2

    critical.Value <- qf(1 - alpha, p, n1 + n2 - p - 1)
    p.value <- min(pf(F, p, n1 + n2 - p - 1), 1 - pf(F, p, n1 + n2 - p - 1))
    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

    statistics <- data.frame(
      "Value" = c(T2, F),
      "p value" = c(" ", p.value),
      "Critical Value" = c(" ", critical.Value)
    )
    rownames(statistics) <- c("Hotelling T2", "F")

    dof <- data.frame("df" = c(paste("T2(", p, n1 + n2 - 2, ")"), paste("F(", p, n1 + n2 - p - 1, ")")))
    rownames(dof) <- c("Hotelling T2", "F")

    testResult <- list(Conclusion = reject, Stat = statistics, SampMean1 = X.bar, SampMean2 = Y.bar, SampA1 = A1, SampA2 = A2, Df = dof, MixSampA = A)

    if (verbose) {
      return(testResult)
    } else {
      message("H0: mu1 = mu2 with unknown but equal covariance matrix\n")
      return(testResult)
    }
  } else {
    if (method == "None") {
      stop("Please choose a method in 'Coupled' or 'Transformed'")
    } else if (method == "Coupled" && n1 == n2) {
      dataCoup <- data1 - data2
      dataC <- as.matrix(dataCoup)

      Z.bar <- apply(dataC, 2, mean)
      Z.A <- (n1 - 1) * cov(dataC)

      T2 <- (n1 - 1) * t(Z.bar) %*% solve(Z.A) %*% Z.bar
      F <- (n1 - p) / ((n1 - 1) * p) * T2

      critical.Value <- qf(1 - alpha, p, n1 - p)
      p.value <- min(pf(F, p, n1 - p), 1 - pf(F, p, n1 - p))
      reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

      statistics <- data.frame(
        "Value" = c(T2, F),
        "p value" = c(" ", p.value),
        "Critical Value" = c(" ", critical.Value)
      )
      rownames(statistics) <- c("Hotelling T2", "F")

      dof <- data.frame("df" = c(paste("T2(", p, ",", n1 - 1, ")"), paste("F(", p, ",", n1 - p, ")")))
      rownames(dof) <- c("Hotelling T2", "F")

      testResult <- list(Conclusion = reject, Stat = statistics, SampMeanC = Z.bar, SampAC = Z.A, Df = dof, dataC = dataCoup)

      if (verbose) {
        return(testResult)
      } else {
        message("H0: mu1 = mu2 with unknown and different covariance matrix\n")
        return(testResult)
      }
    } else if (method == "Coupled" && n1 != n2) {
      stop("The sample size of two groups should be equal when using Coupled method!")
    } else if (method == "Transformed") {
      if (n1 > n2) {
        temp <- data1
        data1 <- data2
        data2 <- temp
        rm(temp)
        n1 <- nrow(data1)
        n2 <- nrow(data2)
      } else {
        n1 <- nrow(data1)
        n2 <- nrow(data2)
      }
      dataTranf <- data1 - sqrt(n1 / n2) * data2 + 1 / sqrt(n1 * n2) * apply(data2[1:n1, ], 2, sum) - 1 / n2 * apply(data2, 2, sum)

      data.z <- as.matrix(dataTranf)

      Z.bar <- apply(data.z, 2, mean)
      Z.A <- (n1 - 1) * cov(data.z)

      T2 <- (n1 - 1) * t(Z.bar) %*% solve(Z.A) %*% Z.bar
      F <- (n1 - p) / ((n1 - 1) * p) * T2

      critical.Value <- qf(1 - alpha, p, n1 - p)
      p.value <- min(pf(F, p, n1 - p), 1 - pf(F, p, n1 - p))
      reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

      statistics <- data.frame(
        "Value" = c(T2, F),
        "p value" = c(" ", p.value),
        "Critical Value" = c(" ", critical.Value)
      )
      rownames(statistics) <- c("Hotelling T2", "F")

      dof <- data.frame(
        "df" = c(paste("T2(", p, ",", n1 - 1, ")"), paste("F(", p, ",", n1 - p, ")"))
      )
      rownames(dof) <- c("Hotelling T2", "F")

      testResult <- list(Conclusion = reject, Stat = statistics, SampMeanT = Z.bar, SampAT = Z.A, Df = dof, dataT = dataTranf)

      if (verbose) {
        return(testResult)
      } else {
        message("H0: mu1 = mu2, with unknown but different covariance matrix\n")
        return(testResult)
      }
    }
  }
}

#' @title Multiple Mean Vectors Hypothesis Testing
#' @description Test whether the mean vectors of multiple multivariate normal populations are all equal when the covariance matrices are equal. Suppose we have k populations, the null hypothesis is that "H0: mu1 = mu2 = ... = muk". There are two approximations (Bartlett's chi2 and Rao's F) to compute the p-value and the critical value. The realized value of the Wilks Lambda statistic and its degrees of freedom are also provided. If you want to perform an exact test, consult the Wilks Lambda statistic quantile table yourself, depending on the realized value of the statistic and its degrees of freedom.
#' @author Xifeng Zhang
#' @param X The data matrix which is a matrix or data frame.
#' @param label A vector of group labels.
#' @param alpha The significance level. Default is 0.05.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.80-83.
#' @return An object of class "testResult", which is a list with the following elements:
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMeanT}{The sample mean.}
#' \item{SampMeanWithin}{The sample mean of each group.}
#' \item{SdTotal}{The total sample deviation.}
#' \item{SdBetween}{The sample deviation between group.}
#' \item{SdWithin}{The sample deviation of each group.}
#' \item{SdWithinT}{The sample deviation within group.}
#' \item{Df}{The degree of freedom.}
#' \item{sampleSize}{The sample size of each group.}
#' @export
#'
#' @examples
#' data(iris)
#' chart <- iris[, 1:4]
#' species <- iris[, 5]
#' # carry out the test
#' test1 <- meanTest.multi(chart, species)
#' test2 <- meanTest.multi(chart, species, verbose = FALSE)
#' # get the elements
#' test1$Stat
#' test1$SampMeanT
#' test1$sampleSize
meanTest.multi <- function(X, label, alpha = 0.05, verbose = TRUE) {
  data <- cbind(X, label)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(unique(label))

  X.bar <- apply(X, 2, mean)
  T <- (n - 1) * cov(X)
  data.split <- split(data[, 1:p], data[, p + 1])
  Within.mean <- lapply(data.split, function(x) {
    apply(x, 2, mean)
  })
  At <- lapply(data.split, function(x) {
    (nrow(x) - 1) * cov(x)
  })
  A <- Reduce("+", At)
  nt <- sapply(data.split, nrow)

  B <- T - A

  Lambda <- det(A) / det(T)
  df.Lambda <- c("df1" = p, "df2" = n - k, "df3" = k - 1)

  Bartlett.chi2 <- -(df.Lambda[2] + df.Lambda[3] - (df.Lambda[1] + df.Lambda[3] + 1) / 2) * log(Lambda)
  df.chi2 <- df.Lambda[1] * df.Lambda[3]
  p.value.chi2 <- min(pchisq(Bartlett.chi2, df.chi2), 1 - pchisq(Bartlett.chi2, df.chi2))
  critical.value.chi2 <- qchisq(1 - alpha, df.chi2)

  t <- df.Lambda[2] + df.Lambda[3] - (df.Lambda[1] + df.Lambda[3] + 1) / 2
  s <- sqrt((df.Lambda[1]^2 * df.Lambda[3]^2 - 4) / (df.Lambda[1]^2 + df.Lambda[3]^2 - 5))
  lambda <- (df.Lambda[1] * df.Lambda[3] - 2) / 4
  Rao.F <- (1 - Lambda^(1 / s)) * (t * s - 2 * lambda) / (Lambda^(1 / s)) * (df.Lambda[1] * df.Lambda[3])
  df.F <- c("df1" = df.Lambda[1] * df.Lambda[3], "df2" = round(t * s - 2 * lambda))

  p.value.F <- min(pf(Rao.F, df.F[1], df.F[2]), 1 - pf(Rao.F, df.F[1], df.F[2]))
  critical.value.F <- qf(1 - alpha, df.F[1], df.F[2])

  reject <- paste("Bartlett's Chi-square: ", ifelse(p.value.chi2 < alpha, "Reject", "Not Reject"), "; Rao's F: ", ifelse(p.value.F < alpha, "Reject", "Not Reject"))

  statistics <- data.frame(
    "Value" = c(Lambda, Bartlett.chi2, Rao.F),
    "p value" = c(" ", p.value.chi2, p.value.F),
    "Critical Value" = c(" ", critical.value.chi2, critical.value.F)
  )
  rownames(statistics) <- c("Wilks Lambda", "Bartlett's Chi2", "Rao's F")

  dof <- data.frame(
    "df" = c(paste("Lambda(", df.Lambda[1], ",", df.Lambda[2], ",", df.Lambda[3], ")"), paste("Chi2(", df.chi2, ")"), paste("F(", df.F[1], ",", df.F[2], ")"))
  )
  rownames(dof) <- c("Wilks Lambda", "Bartlett's Chi2", "Rao's F")

  testResult <- list(Conclusion = reject, Stat = statistics, SampMeanT = X.bar, SampMeanWithin = Within.mean, SdTotal = T, SdBetween = B, SdWithin = At, SdWithinT = A, Df = dof, sampleSize = nt)

  if (verbose) {
    return(testResult)
  } else {
    message("H0: mu1 = mu2 = ... = muk, k = ", k, "\n")
    return(testResult)
  }
}

#' @title Single Covariance Matrix Hypothesis Testing
#' @description Test whether the covariance matrix is equal to a certain value. The null hypothesis is "H0: Sigma = Sigma0" or "H0: Sigma = sigma^2 * Sigma0".
#' @author Xifeng Zhang
#' @param data The data matrix which is a matrix or data frame.
#' @param Sigma0 The covariance matrix when the null hypothesis is true.
#' @param ball A boolean value. Default is FALSE. If FALSE, test whether the covariance matrix is Sigma0 (known), which means the null hypothesis is "H0: Sigma = Sigma0". If TRUE and the Sigma0 is a unit matrix, the Mauchly's ball test will be performed. If TRUE but Sigma0 (known) is not a unit matrix, the covariance matrix will be tested to see if it is sigma^2*Sigma0 (sigma^2 is unknown), which means the null hypothesis is "H0: Sigma = sigma^2 * Sigma0".
#' @param alpha The significance level. Default is 0.05.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @importFrom Rmpfr mpfr
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.83-88.
#' @return An object of class "testResult", which is a list with the following elements:
#' Return when ball is FALSE.
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMean}{The sample mean.}
#' \item{SampA}{The sample deviation.}
#' \item{Df}{The degree of freedom.}
#' Return when ball is TRUE
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMean}{The sample mean.}
#' \item{SampA}{The sample deviation.}
#' \item{sigma.hat}{The estimation of sigma^2.}
#' \item{Df}{The degree of freedom.}
#' @export
#'
#' @examples
#' data(iris)
#' X <- iris[, 1:4]
#' # carry out the test
#' test1 <- covTest.single(X, diag(1, 4))
#' test2 <- covTest.single(X, diag(1, 4), ball = TRUE)
#' test3 <- covTest.single(X, diag(2, 4), ball = TRUE)
#' test4 <- covTest.single(X, diag(1, 4), verbose = FALSE)
#' # get the elements
#' test1$Stat
#' test2$Df
#' test3$sigma.hat
covTest.single <- function(data, Sigma0, ball = FALSE, alpha = 0.05, verbose = TRUE) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- apply(data, 2, mean)
  A <- (n - 1) * cov(data)

  if (isFALSE(ball)) {
    lambda <- exp(-1 / 2 * sum(diag(A %*% solve(Sigma0)))) * Rmpfr::mpfr(det(A %*% solve(Sigma0)), 2048)^(n / 2) * exp(n * p / 2) * Rmpfr::mpfr(n, 2048)^(-n * p / 2)
    chi2 <- -2 * log(lambda)

    lambda <- as.numeric(lambda)
    chi2 <- as.numeric(chi2)

    df.chi2 <- p * (p + 1) / 2

    p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
    critical.value <- qchisq(1 - alpha, df.chi2)

    statistics <- data.frame(
      "Value" = c(lambda, chi2),
      "p value" = c(" ", p.value),
      "Critical Value" = c(" ", critical.value)
    )
    rownames(statistics) <- c("Likelihood Ratio", "Chi2")

    dof <- paste("Chi2(", df.chi2, ")")

    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

    testResult <- list(Conclusion = reject, Stat = statistics, SampMean = X.bar, SampA = A, Df = dof)

    if (verbose) {
      return(testResult)
    } else {
      message("H0: Sigma = Sigma0\n")
      return(testResult)
    }
  } else {
    sigma.hat <- sum(diag(solve(Sigma0) %*% A)) / (n * p)
    lambda <- (Rmpfr::mpfr(det(solve(Sigma0) %*% A), 2048)^(n / 2)) / (Rmpfr::mpfr(sum(diag(solve(Sigma0)) %*% A) / p, 2048)^(n * p / 2))
    W <- lambda^(2 / n)
    chi2 <- -((n - 1) - (2 * p^2 + p + 2) / (6 * p)) * log(W)

    lambda <- as.numeric(lambda)
    W <- as.numeric(W)
    chi2 <- as.numeric(chi2)

    df.chi2 <- p * (p + 1) / 2 - 1

    p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
    critical.value <- qchisq(1 - alpha, df.chi2)
    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

    statistics <- data.frame(
      "Value" = c(lambda, W, chi2),
      "p value" = c(" ", " ", p.value),
      "Critical Value" = c(" ", " ", critical.value)
    )
    rownames(statistics) <- c("Likelihood Ratio", "W", "Chi2")

    dof <- paste("Chi2(", df.chi2, ")")

    testResult <- list(Conclusion = reject, Stat = statistics, SampMean = X.bar, SampA = A, Df = dof, sigma.hat = sigma.hat)

    if (verbose) {
      return(testResult)
    } else {
      message("H0: Sigma = sigma^2 * Sigma0\n")
      return(testResult)
    }
  }
}

#' @title Multiple Covariance Matrix Hypothesis Testing
#' @description Test whether the covariance matrices of multiple multivariate normal populations are all equal. Suppose we have k populations, the null hypothesis is "H0: Sigma1 = Sigma2 = ... = Sigmak".
#' @author Xifeng Zhang
#' @param X The data matrix which is a matrix or data frame.
#' @param label A vector of group labels.
#' @param alpha The significance level. Default is 0.05.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @importFrom Rmpfr mpfr
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.88-89.
#' @return An object of class "testResult", which is a list with the following elements:
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMeanT}{The sample mean.}
#' \item{SampMeanWithin}{The sample mean of each group.}
#' \item{SdWithin}{The sample deviation of each group.}
#' \item{SdWithinT}{The sample deviation within group.}
#' \item{Df}{The degree of freedom.}
#' \item{sampleSize}{The sample size of each group.}
#' \item{d}{The Modified factor of the statistic.}
#' @export
#'
#' @examples
#' data(iris)
#' chart <- iris[, 1:4]
#' species <- iris[, 5]
#' # carry out the test
#' test1 <- covTest.multi(chart, species)
#' test2 <- covTest.multi(chart, species, verbose = FALSE)
#' # get the elements
#' test1$Stat
#' test1$SampMeanT
#' test1$sampleSize
covTest.multi <- function(X, label, alpha = 0.05, verbose = TRUE) {
  data <- data.frame(cbind(X, label))
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(unique(label))

  X.bar <- apply(X, 2, mean)

  data.split <- split(data[, 1:p], data[, p + 1])
  nt <- sapply(data.split, nrow)
  Within.mean <- lapply(data.split, function(x) {
    apply(x, 2, mean)
  })
  At <- lapply(data.split, function(x) {
    (nrow(x) - 1) * cov(x)
  })
  A <- Reduce("+", At)

  St <- lapply(data.split, function(x) {
    cov(x)
  })
  St.root <- lapply(data.split, function(x) {
    Rmpfr::mpfr(det(cov(x)), 2048)^(-(nrow(x) - 1) / 2)
  })
  St.sum <- lapply(data.split, function(x) {
    (nrow(x) - 1) * log(Rmpfr::mpfr(det(cov(x)), 2048))
  })
  St.star <- lapply(data.split, function(x) {
    cov(x) * (nrow(x) - 1) / nrow(x)
  })
  St.star.root <- lapply(data.split, function(x) {
    Rmpfr::mpfr(det(cov(x) * (nrow(x) - 1) / nrow(x)), 2048)^(-nrow(x) / 2)
  })

  lambda <- Rmpfr::mpfr(det(A / nrow(data)), 2048)^(-n / 2) / Rmpfr::mpfr(Reduce("*", St.star.root), 2048)
  lambda.star <- Rmpfr::mpfr(det(A / (n - k)), 2048)^(-(n - k) / 2) / Rmpfr::mpfr(Reduce("*", St.root), 2048)
  M <- -2 * log(lambda.star)
  M.star <- (n - k) * log(det(A / (n - k))) - Rmpfr::mpfr(Reduce("+", St.sum), 2048)

  lambda <- as.numeric(lambda)
  lambda.star <- as.numeric(lambda.star)
  M <- as.numeric(M)
  M.star <- as.numeric(M.star)

  if (identical(nt, rep(1, k)) == TRUE) {
    d <- (2 * p^2 + 3 * p - 1) * (k + 1) / (6 * (p + 1) * (n - k))
  } else {
    d <- (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (k - 1)) * (sum(1 / (nt - 1)) - 1 / (n - k))
  }

  chi2 <- (1 - d) * M
  df.chi2 <- p * (p + 1) * (k - 1) / 2
  p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
  critical.value <- qchisq(1 - alpha, df.chi2)
  reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

  statistics <- data.frame(
    "Value" = c(lambda, lambda.star, M.star, chi2),
    "p value" = c(" ", " ", " ", p.value),
    "Critical Value" = c(" ", "", " ", critical.value)
  )
  rownames(statistics) <- c("Likelihood Ratio", "Likelihood Ratio (Adjusted)", "M", "Chi2")

  dof <- paste("Chi2(", df.chi2, ")")

  testResult <- list(Conclusion = reject, Stat = statistics, SampMeanT = X.bar, SampMeanWithin = Within.mean, SdWithin = At, SdWithinT = A, Df = dof, sampleSize = nt, d = d)

  if (verbose) {
    return(testResult)
  } else {
    message("H0: Sigma1 = Sigma2 = ... = Sigmak, k = ", k, "\n")
    return(testResult)
  }
}

#' @title Mean and Covariance Matrix Hypothesis Testing (Simultaneously)
#' @description Test whether the mean vectors and covariance matrices of multiple multivariate normal populations are all equal simultaneously. Suppose we have k populations, the null hypothesis is "H0: mu1 = mu2 = ... = muk and Sigma1 = Sigma2 = ... = Sigmak".
#' @author Xifeng Zhang
#' @param X The data matrix which is a matrix or data frame.
#' @param label A vector of group labels.
#' @param alpha The significance level. Default is 0.05.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @importFrom Rmpfr mpfr
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.90-91.
#' @return An object of class "testResult", which is a list with the following elements:
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMeanT}{The sample mean.}
#' \item{SampMeanWithin}{The sample mean of each group.}
#' \item{SdTotal}{The total sample deviation.}
#' \item{SdWithin}{The sample deviation of each group.}
#' \item{SdWithinT}{The sample deviation within group.}
#' \item{Df}{The degree of freedom.}
#' \item{sampleSize}{The sample size of each group.}
#' \item{d}{The Modified factor of the statistic.}
#' @export
#'
#' @examples
#' data(iris)
#' chart <- iris[, 1:4]
#' species <- iris[, 5]
#' # carry out the test
#' test1 <- meancov.Test(chart, species)
#' test2 <- meancov.Test(chart, species, verbose = FALSE)
#' # get the elements
#' test1$Stat
#' test1$SampMeanT
meancov.Test <- function(X, label, alpha = 0.05, verbose = TRUE) {
  data <- data.frame(cbind(X, label))
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(unique(label))

  X.bar <- apply(X, 2, mean)

  data.split <- split(data[, 1:p], data[, p + 1])
  nt <- sapply(data.split, nrow)
  Within.mean <- lapply(data.split, function(x) {
    apply(x, 2, mean)
  })
  At <- lapply(data.split, function(x) {
    (nrow(x) - 1) * cov(x)
  })
  A <- Reduce("+", At)
  T <- (n - 1) * cov(X)

  lambda <- Reduce("*", lapply(data.split, function(x) {
    Rmpfr::mpfr(det((nrow(x) - 1) * cov(x)), 2048)^(nrow(x) / 2)
  })) / Rmpfr::mpfr(det(T), 2048)^(n / 2) * Rmpfr::mpfr(n, 2048)^(n * p / 2) / Rmpfr::mpfr(Reduce("*", lapply(data.split, function(x) {
    Rmpfr::mpfr(nrow(x), 2048)^(nrow(x) * p / 2)
  })), 2048)

  lambda.star <- Reduce("*", lapply(data.split, function(x) {
    Rmpfr::mpfr(det((nrow(x) - 1) * cov(x)), 2048)^((nrow(x) - 1) / 2)
  })) / Rmpfr::mpfr(det(T), 2048)^((n - k) / 2) * Rmpfr::mpfr((n - k), 2048)^((n - k) * p / 2) / Rmpfr::mpfr(Reduce("*", lapply(data.split, function(x) {
    Rmpfr::mpfr((nrow(x) - 1), 2048)^((nrow(x) - 1) * p / 2)
  })), 2048)

  d <- (Reduce("+", lapply(data.split, function(x) {
    sum(1 / (nrow(x) - 1))
  })) - sum(1 / (n - k))) * ((2 * p^2 + 3 * p - 1) / (6 * (p + 3) * (k - 1))) - (p - k + 2) / ((n - k) * (p + 3))

  M <- -2 * log(lambda.star)
  chi2 <- (1 - d) * M

  df.chi2 <- p * (p + 3) * (k - 1) / 2

  lambda <- as.numeric(lambda)
  lambda.star <- as.numeric(lambda.star)
  M <- as.numeric(M)
  chi2 <- as.numeric(chi2)

  p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
  critical.value <- qchisq(1 - alpha, df.chi2)
  reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

  statistics <- data.frame(
    "Value" = c(lambda, lambda.star, M, chi2),
    "p value" = c(" ", " ", " ", p.value),
    "Critical Value" = c(" ", " ", " ", critical.value)
  )
  rownames(statistics) <- c("Likelihood Ratio", "Likelihood Ratio (Modified)", "M", "Chi2")

  dof <- paste("Chi2(", df.chi2, ")")

  testResult <- list(Conclusion = reject, Stat = statistics, SampMeanT = X.bar, SampMeanWithin = Within.mean, SdTotal = T, SdWithin = At, SdWithinT = A, Df = dof, sampleSize = nt, d = d)

  if (verbose) {
    return(testResult)
  } else {
    message("H0: mu1 = mu2 = ... = muk, Sigma1 = Sigma2 = ... = Sigmak, k = ", k, "\n")
    return(testResult)
  }
}

#' @title Multivariate Normal Independence Test
#' @description Test whether a set of multivariate normal random vectors are independent. The null hypothesis is "H0: The random vectors are independent of each other".
#' @author Xifeng Zhang
#' @param data The data matrix which is a matrix or data frame. Each column represents a random variable.
#' @param subdim The dimensions of submatrices. The default is FALSE, which means the independence of all components of the random vector will be tested.
#' @param alpha The significance level. Default is 0.05.
#' @param verbose A boolean value. Default is TRUE. If TRUE, the null hypothesis will be displayed. If FALSE, the test will be carried out silently.
#' @import stats
#' @import utils
#' @references Huixuan, Gao. Applied Multivariate Statistical Analysis. Peking University Press, 2005: pp.92-94.
#' @return An object of class "testResult", which is a list with the following elements:
#' \item{Conclusion}{The conclusion of the test.}
#' \item{Stat}{A data frame containing the statistics, p value and critical value.}
#' \item{SampMean}{The sample mean.}
#' \item{SampA}{The sample deviation.}
#' \item{SampAii}{The sample deviation of submatrices.}
#' \item{Df}{The degree of freedom.}
#' \item{b}{The Modified factor of the statistic.}
#' @export
#'
#' @examples
#' data(iris)
#' chart <- iris[, 1:4]
#' # carry out the test
#' test1 <- indTest.multi(chart)
#' test2 <- indTest.multi(chart, subdim = c(2, 1, 1))
#' test3 <- indTest.multi(chart, verbose = FALSE)
#' # get the elements
#' test1$Stat
#' test1$SampMean
#' test2$SampAii
indTest.multi <- function(data, subdim = FALSE, alpha = 0.05, verbose = TRUE) {
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- round(apply(data, 2, mean), 4)
  A <- (n - 1) * cov(data)
  Aii <- list()

  if (isFALSE(subdim)) {
    subdim <- rep(1, p)
  } else {
    subdim <- subdim
  }

  subcol <- c(0, cumsum(subdim))

  for (i in seq_along(subdim)) {
    datai <- data[, (subcol[i] + 1):subcol[i + 1]]
    Aii <- append(Aii, list((n - 1) * cov(as.matrix(datai))))
  }

  V <- (det(A) / Reduce("*", lapply(Aii, det)))
  lambda <- V^(n / 2)

  b <- n - 3 / 2 - (p^3 - sum(subdim^3)) / (3 * (p^2 - sum(subdim^2)))
  chi2 <- -b * log(V)

  df.chi2 <- 1 / 2 * (p * (p + 1) - sum(subdim * (subdim + 1)))

  p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
  critical.value <- qchisq(1 - alpha, df.chi2)
  reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

  statistics <- data.frame(
    "Value" = c(V, lambda, chi2),
    "p value" = c(" ", " ", p.value),
    "Critical Value" = c(" ", " ", critical.value)
  )
  rownames(statistics) <- c("V", "Likelihood Ratio", "Chi2")

  dof <- paste("Chi2(", df.chi2, ")")

  testResult <- list(Conclusion = reject, Stat = statistics, SampMean = X.bar, SampA = A, SampAii = Aii, Df = dof, b = b)

  if (verbose) {
    return(testResult)
  } else {
    message("H0: The random vectors are independent of each other\n")
    return(testResult)
  }
}
