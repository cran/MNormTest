
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MNormTest

<!-- badges: start -->

[![R-CMD-check](https://github.com/Astringency/MNormTest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Astringency/MNormTest/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Astringency/MNormTest/graph/badge.svg)](https://app.codecov.io/gh/Astringency/MNormTest)
<!-- badges: end -->

## Code of Conduct

Please note that the MNormTest project is released with a [Contributor
Code of
Conduct](https://astringency.github.io/MNormTest/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## 介绍 (Introduction)

MNormTest的目标是提供一套用于多元正态分布参数假设检验的函数，包括检验单个均值向量、两个均值向量、多个均值向量、单个协方差矩阵、多个协方差矩阵、同时检验均值和协方差矩阵，以及检验多元正态随机向量的独立性.
方法原理参考了高惠璇（2005,
<ISBN:9787301078587）的《应用多元统计分析>》.

The objective of MNormTest is to provide a set of functions for
hypothesis testing of the parameters of multivariate normal
distributions, including the testing of a single mean vector, two mean
vectors, multiple mean vectors, a single covariance matrix, multiple
covariance matrices, a mean and a covariance matrix simultaneously, and
the testing of independence of multivariate normal random vectors. The
methods are based on the book “Applied Multivariate Statistical
Analysis” by Huixuan Gao (2005, <ISBN:9787301078587>).

## 下载 (Installation)

- CRAN:

首先，您可以从[CRAN](https://CRAN.R-project.org)安装MNormTest的发布版本:

Firstly, you can install the released version of MNormTest from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MNormTest")
```

- GitHub:

其次，您可以从[GitHub](https://github.com/)安装MNormTest的开发版本:

What’s more, You can install the development version of MNormTest from
[GitHub](https://github.com/) with:

``` r
install.packages("pak")
pak::pak("Astringency/MNormTest")
```

## 加载 (Load)

在使用MNormTest包之前，您需要先加载它:

Before using MNormTest, you need to load it:

``` r
library(MNormTest)
```

## 使用 (Usage)

MNormTest提供了7个函数：

MNormTest provides 7 functions:

`meanTest.single`, `meanTest.two`, `meanTest.multi`, `covTest.single`,
`covTest.multi`, `meancov.Test`, `ind.Test.multi`.

### 均值向量检验 (Mean Vector Test)

`meanTest.single`被用于单个多元正态总体（总体协方差阵已知和未知时）的均值向量的检验.详细内容可参照函数的说明文档.

`meanTest.single` is used for testing the mean vector of a single
multivariate normal population (when the population covariance matrix is
known and unknown). For more details, please refer to the documentation
of the function.

``` r
?meanTest.single
# or ??meanTest.single
```

`meanTest.two`被用于两个多元正态总体（总体协方差阵未知相等或不等）的均值向量的检验.
在总体协方差阵不相等时，函数提供了两种处理方法将其化为单样本问题：构造配对数据（两样本量相同）`Coupled`和构造修正后的数据（两样本量不同）`Tranformed`，详细内容可参照函数的说明文档.

`meanTest.two` is used for testing the mean vectors of two multivariate
normal populations (when the population covariance matrices are unknown
and equal or unequal). When the population covariance matrices are
unequal, the function provides two methods to reduce it to a single
sample problem: constructing paired data (two samples of the same size)
`Coupled` and constructing transformed data (two samples of different
sizes) `Tranformed`. For more details, please refer to the documentation
of the function.

``` r
?meanTest.two
# or ??meanTest.two
```

`meanTest.multi`被用于多个多元正态总体（总体协方差阵已知且相等）的均值向量的检验，采用Wilks
$\Lambda$统计量进行检验.
但由于分布的特殊性，函数提供了两种近方法（Bartlett’s $\chi^2$和Rao’s
$F$）的检验结果，并给出了Wilks
$\Lambda$统计量的实现值及其自由度（需指定`full = TRUE`），若您想进行精确检验，则请根据统计量的实现值和自由度自行查找临界值.
详细内容可参照函数的说明文档.

`meanTest.multi` is used for testing the mean vectors of multiple
multivariate normal populations (when the population covariance matrices
are known and equal), using Wilks $\Lambda$ statistic for testing.
However, due to the special distribution, the function provides the test
results of two approximate methods (Bartlett’s $\chi^2$ and Rao’s $F$)
and gives the observed value and degrees of freedom of Wilks $\Lambda$
statistic (need to specify `full = TRUE`), if you want to conduct an
exact test, please find the critical value according to the observed
value and degrees of freedom of the statistic. For more details, please
refer to the documentation of the function.

``` r
?meanTest.multi
# or ??meanTest.multi
```

### 协方差矩阵检验 (Covariance Matrix Test)

`covTest.single`被用于单个多元正态总体（总体均值未知时）的协方差矩阵的检验.
包含了2种情况：

`covTest.single` is used for testing the covariance matrix of a single
multivariate normal population (when the population mean is unknown).
The function contains two cases:

$$
H_0: \Sigma = \Sigma_0, \quad H_1: \Sigma \neq \Sigma_0 (\Sigma_0 \text{ is known})
$$

$$
H_0: \Sigma = \sigma^2 \Sigma_0, \quad H_1: \Sigma \neq \sigma^2 \Sigma_0 (\sigma^2 \text{ is unknown}, \Sigma_0 \text{ is known})
$$

在第二种情况中，若$\Sigma_0 = I$，则检验通常称为球形检验.
您可以通过设置不同的函数参数来选择检验的类型，详细内容可参照函数的说明文档.

In the second case, if $\Sigma_0 = I$, the test is usually called a
sphericity test. You can choose the type of test by setting different
parameters of the function. For more details, please refer to the
documentation of the function.

``` r
?covTest.single
# or ??covTest.single
```

`covTest.multi`被用于多个多元正态总体（总体均值未知时）的协方差矩阵的检验.
详细内容可参照函数的说明文档.

`covTest.multi` is used for testing the covariance matrices of multiple
multivariate normal populations (when the population mean is unknown).
For more details, please refer to the documentation of the function.

``` r
?covTest.multi
# or ??covTest.multi
```

### 均值和协方差矩阵的同时检验 (Mean and Covariance Matrix Test)

`meancov.Test`被用于多个多元正态总体的均值和协方差矩阵的同时检验.
详细内容可参照函数的说明文档.

`meancov.Test` is used for testing the mean and covariance matrices of
multiple multivariate normal populations simultaneously. For more
details, please refer to the documentation of the function.

``` r
?meancov.Test
# or ??meancov.Test
```

### 多元正态随机向量的独立性检验 (Independence Test of Multivariate Normal Random Vectors)

`ind.Test.multi`被用于多元正态随机向量的独立性检验.
详细内容可参照函数的说明文档.

`ind.Test.multi` is used for testing the independence of multivariate
normal random vectors. For more details, please refer to the
documentation of the function.

``` r
?ind.Test.multi
# or ??ind.Test.multi
```

## 例 (Example)

下面是几个简单的例子，展示了如何使用MNormTest包进行多元正态分布参数的假设检验，我们使用著名的鸢尾花数据集`iris`进行多总体协方阵的检验.
更多示例可以参考每个函数的说明文档.

Here are some simple examples that show how to use MNormTest to test the
parameters of multivariate normal distributions. We use the famous iris
dataset `iris` to test the covariance matrices of multiple populations.
More examples can be found in the documentation of each function.

我们首先来看一下数据集`iris`的前几行:

First, let’s take a look at the first few rows of the dataset `iris`:

``` r
data(iris)
head(iris)
```

该数据集包含了150个样品，每个样品包含5个变量：4个特征（`Sepal.Length`,
`Sepal.Width`, `Petal.Length`,
`Petal.Width`）和1个类别（`Species`）.为了方便后面的检验，我们将二者分开:

The dataset contains 150 samples, each sample contains 5 variables: 4
features (`Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width`)
and 1 category (`Species`). For the convenience of the following tests,
we separate them:

``` r
chart <- iris[, 1:4]
species <- iris[, 5]
```

鸢尾花的类别有三种：`setosa`, `versicolor`, `virginica`.
这便是我们的三个总体（这里假设三个总体均服从四元正态分布），我们现在想要检验它们的协方差矩阵是否相等.
检验假设如下：

There are three categories of iris: `setosa`, `versicolor`, `virginica`.
These are our three populations (assuming that the three populations all
follow a four-dimensional normal distribution), and we now want to test
whether their covariance matrices are equal. The hypotheses are as
follows:

$$
H_0: \Sigma_1 = \Sigma_2 = \Sigma_3, \quad H_1: \text{not all } \Sigma_i \text{ are equal}
$$

给定显著性水平$\alpha = 0.05$，我们可以使用`covTest.multi`进行检验：

Given the significance level $\alpha = 0.05$, we can use `covTest.multi`
for testing:

``` r
test.iris <- covTest.multi(chart, species)
```

如果您想查看检验的原假设，可以设置 `verbose = FALSE`:

If you want to see the null hypothesis of the test, you can set
`verbose = FALSE`:

``` r
test.iris <- covTest.multi(chart, species, verbose = FALSE)
```

计算得$\chi^2$统计量为140.94，$p = 0.0000 < 0.05 = \alpha$（$p$值很小，在`R`的运算中直接算成0了），故拒绝原假设，认为在显著性水平$\alpha = 0.05$下，三种鸢尾花的协方差矩阵不相等.

The calculated $\chi^2$ statistic is 140.94,
$p = 0.0000 < 0.05 = \alpha$ ($p$ is very small, and it is directly
calculated as 0 in `R`), so we reject the null hypothesis and conclude
that the covariance matrices of the three categories of iris are not
equal at the significance level $\alpha = 0.05$.

此外，您可以从`test.iris`中提取感兴趣的信息，如检验统计量的实现值、$p$值、临界值和自由度:

In addition, you can extract the information of interest from
`test.iris`, such as the observed value of the test statistic, $p$
value, critical value, and degrees of freedom:

``` r
test.iris$Stat
test.iris$Df
```

更多示例可以参考每个函数的说明文档.

More examples can be found in the documentation of each function.
