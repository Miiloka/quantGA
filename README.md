# quantGA

`quantGA` provides a quantile regression estimator based on a genetic algorithm. The package relies on `Rcpp` and `RcppArmadillo` for the heavy computations performed in C++.

## Prerequisites

- **R** (version 4.0 or later recommended)
- The **Rcpp** and **RcppArmadillo** packages (as well as a working C++ toolchain)

Install the development version of `quantGA` by simply using `devtools::install_github("https://github.com/Miiloka/quantGA.git")`

## Basic usage

Below is a simplified example adapted from `.ignore/tests.R` showing how to generate a small data set and estimate a median regression model:

```r
set.seed(45)

n <- 5000
x1 <- rnorm(n)
x2 <- runif(n, -1, 1)
beta_true <- c(2, -1, 0.5)
errors <- rt(n, df = 3) # heavy-tailed

y <- as.numeric(cbind(1, x1, x2) %*% beta_true + errors)

df <- data.frame(y = y, x1 = x1, x2 = x2)

result <- GA.QR(
  formula = y ~ x1 + x2,
  data    = df,
  tau     = 0.5
)

print(result)
```

`result` is an object of class `quantileGA` which can be printed, plotted or summarised with the usual S3 methods provided in this package.

## Simulation

For a more complete Monte Carlo simulation illustrating the use of `GA.QR` at several quantile levels, see the script [`Simus.R`](.ignore/Simus.R) located in the `.ignore` folder.
