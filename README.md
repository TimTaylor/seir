
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Exploration of deSolve with compiled code

<!-- badges: start -->
<!-- badges: end -->

## Example

Note, at present, the implementations assumes:

- a closed population (no births, deaths or ageing)
- no time dependent parameters or interventions.

``` r
library(seir)

# demography
contact_data <- suppressMessages(
    socialmixr::contact_matrix(
        socialmixr::polymod,
        countries = "United Kingdom",
        age.limits = c(0,20,65),
        symmetric = TRUE
    )
)

contact_matrix <- t(contact_data$matrix)
demography <- contact_data$demography$population
n <- nrow(contact_matrix)

# initial conditions
init <- matrix(
    c(S = 0.999999, E = 0, I = 0.000001, R = 0),
    nrow = n,
    ncol = 4L,
    byrow = TRUE
)

# epi parameters
alpha <- 1/3
beta <- 1.5/7
gamma <- 1/7

# time stuff
time_end <- 600
increment <- 1

# R implementation
out_r <- seir_r(alpha, beta, gamma, contact_matrix, demography, init, time_end)

# C implementation
out_c <- seir_c(alpha, beta, gamma, contact_matrix, demography, init, time_end)

# compare outputs
infective_index <- (2*n + 2):(3*n + 1)
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
op <- par(mai=rep(0.5, 4))
matplot(out_r[,infective_index], ylab = "", xlab = "", type = "l", main = "R")
matplot(out_c[,infective_index], ylab = "", xlab = "", type = "l", main = "C")
par(mai=c(0,0,0,0))
plot.new()
legend(
    title = "age category",
    x="center", ncol=n, legend=rownames(contact_matrix), fill = seq_len(n)
)
```

![](man/figures/README-example-1.png)<!-- -->

``` r
par(op)
microbenchmark::microbenchmark(
    R = seir_r(alpha, beta, gamma, contact_matrix, demography, init, time_end),
    C = seir_c(alpha, beta, gamma, contact_matrix, demography, init, time_end)
)
#> Unit: microseconds
#>  expr      min        lq       mean     median         uq       max neval cld
#>     R 9799.274 10117.782 10732.0709 10199.3280 10433.3535 16997.538   100  a 
#>     C  741.053   794.137   867.8647   859.2965   938.5335  1173.579   100   b
```
