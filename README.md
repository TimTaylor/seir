
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Exploration of deSolve with compiled code

<!-- badges: start -->
<!-- badges: end -->

## Example

Note, at present, the implementations assumes:

- a closed population (no births, deaths or ageing);
- single (age-dependent) vaccination;
- cumulative contact reduction interventions;
- no other time dependent parameters.

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
init_i <- 1e-4
init <- matrix(
    c(S = 1-init_i, E = 0, I = init_i, R = 0),
    nrow = n,
    ncol = 4L,
    byrow = TRUE
)

# epi parameters
infectious_period <- 5
preinfectious_period <- 1.5
r0 <- 1.3

alpha <- 1 / preinfectious_period
beta <- r0 / infectious_period
gamma <- 1 / infectious_period

# time stuff
time_end <- 600
increment <- 1

is <- c(60,80)
ie <- c(240,140)
ii <- cbind(c(0.3,0.01,0.01), c(0.01,0.3,0.01))

# R implementation
out_r <- seir_r(
    alpha, beta, gamma,
    contact_matrix, demography,
    init, time_end,
    intervention_start_times = c(60,80), intervention_end_times = c(240,140),
    intervention = cbind(c(0.3,0.01,0.01), c(0.01,0.3,0.01))
)

# C implementation
out_c <- seir_c(
    alpha, beta, gamma,
    contact_matrix, demography,
    init, time_end,
    intervention_start_times = c(60,80), intervention_end_times = c(240,140),
    intervention = cbind(c(0.3,0.01,0.01), c(0.01,0.3,0.01))
)

# compare outputs
s_index <- (0*n + 2):(1*n + 1)
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
op <- par(mai=rep(0.5, 4))

matplot(out_r[1:600, s_index] - out_r[2:601, s_index], ylab = "", xlab = "", type = "l", main = "R ODE")
matplot(out_c[1:600, s_index] - out_c[2:601, s_index], ylab = "", xlab = "", type = "l", main = "C ODE")
par(mai=c(0,0,0,0))
plot.new()
legend(
    title = "age category", x="center",
    legend=rownames(contact_matrix),
    ncol=n, col = seq_len(n), lty = seq_len(n)
)
```

![](man/figures/README-example-1.png)<!-- -->

``` r
par(op)
microbenchmark::microbenchmark(
    R = seir_r(
        alpha, beta, gamma,
        contact_matrix, demography,
        init, time_end,
        intervention_start_times = c(60,80), intervention_end_times = c(240,140),
        intervention = cbind(c(0.3,0.01,0.01), c(0.01,0.3,0.01))
    ),
    C = seir_c(
        alpha, beta, gamma,
        contact_matrix, demography,
        init, time_end,
        intervention_start_times = c(60,80), intervention_end_times = c(240,140),
        intervention = cbind(c(0.3,0.01,0.01), c(0.01,0.3,0.01))
    )
)
#> Unit: microseconds
#>  expr       min       lq      mean    median        uq       max neval cld
#>     R 22044.896 22980.14 23903.176 23255.611 24957.828 30075.917   100  a 
#>     C   930.805   992.98  1102.562  1084.883  1212.781  1553.091   100   b
```
