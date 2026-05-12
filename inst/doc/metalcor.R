## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
                      collapse = TRUE,
                      comment = "#>",
                      fig.width = 6, fig.height = 4, out.width = "100%"
                  )

## -----------------------------------------------------------------------------
# number of studies
S <- 3
# number of loci
m <- 10000
# this is the true covariance matrix of this simulation
R <- matrix(
    c( 1.0, 0.2, 0.0,
      0.2, 1.0, 0.1,
      0.0, 0.1, 1.0 ),
    nrow = 3, ncol = 3
)

# simulate the correlated Z scores first
# use Cholesky decomposition to simulate via affine transform
L <- chol( R )
Z <- matrix( rnorm( m * 3 ), nrow = m, ncol = 3 ) %*% L

# construct some trivial study summary statistics with the above Z scores
library(tibble)
library(dplyr)
study1 <- tibble(
    id = 1 : m,
    chr = ceiling( id / m * 10 ),
    pos = id,
    n = 1000,
    z = Z[, 1],
    se = runif( m, 0, 1 / sqrt( n ) ),
    beta = z * se
)
# this copies id/chr/pos between studies, adds study-specific info.
# note that, besides the different Z scores, there's also different sample sizes per study,
# and standard errors scale inversely proportional to the square root of study size
study2 <- study1 %>%
    select( id, chr, pos) %>%
    mutate(
        n = 2000,
        z = Z[, 2],
        se = runif( m, 0, 1 / sqrt( n ) ),
        beta = z * se
    )
study3 <- study1 %>%
    select( id, chr, pos) %>%
    mutate(
        n = 10000,
        z = Z[, 3],
        se = runif( m, 0, 1 / sqrt( n ) ),
        beta = z * se
    )
# for realism, randomly subset studies so only 90% of loci survive, different sets per study
study1 <- study1[ sort( sample.int( m, m * 0.9) ), ]
study2 <- study2[ sort( sample.int( m, m * 0.9) ), ]
study3 <- study3[ sort( sample.int( m, m * 0.9) ), ]

## -----------------------------------------------------------------------------
library(metalcor)

# gather studies into a list
studies <- list( study1, study2, study3 )
# pass to metalcor!
out <- metalcor( studies )
# you can quickly inspect the return value, it is a list with two named elements:
# - assoc: the meta-analyzed association table
# - R: the estimated covariance matrix between studies
out

## -----------------------------------------------------------------------------
mean( out$assoc$z^2 )

## -----------------------------------------------------------------------------
# default median estimate
estimate_R( Z )
# mean estimate
estimate_R( Z, median = FALSE )
# standard sample covariance estimate is similar but not identical
cov( Z, use = 'pairwise.complete.obs' )

## -----------------------------------------------------------------------------
n <- 10000
rho <- 0.33

# draw random values from this distribution
x <- rprodcor( n, rho )
x <- sort( x )
# this is the theoretical density function at those points
y <- dprodcor( x, rho )

# compare empirical and theoretical densities
hist( x, breaks = 100, freq = FALSE )
lines( x, y, col = 'red' )

# compare empirical and theoretical CDF now
yc <- pprodcor( x, rho )
plot( ecdf( x ) )
lines( x, yc, lty = 2, col = 'red' )

## -----------------------------------------------------------------------------
x <- ( (-200) : 200 ) / 100
plot( x, dprodcor( x, -1 ), type ='l', col = 'red', ylab = 'Density' )
lines( x, dprodcor( x, -0.5 ), col = 'orange' )
lines( x, dprodcor( x, 0 ), col = 'black' )
lines( x, dprodcor( x, 0.5 ), col = 'green' )
lines( x, dprodcor( x, 1 ), col = 'blue' )
legend(
    'topright',
    legend = c( -1, -0.5, 0, 0.5, 1 ),
    col = c( 'red', 'orange', 'black', 'green', 'blue' ),
    lty = 1,
    title = expression( rho )
)

## -----------------------------------------------------------------------------
plot( x, pprodcor( x, -1 ), type ='l', col = 'red', ylim = c(0,1), ylab = 'CDF' )
lines( x, pprodcor( x, -0.5 ), col = 'orange' )
lines( x, pprodcor( x, 0 ), col = 'black' )
lines( x, pprodcor( x, 0.5 ), col = 'green' )
lines( x, pprodcor( x, 1 ), col = 'blue' )
legend(
    'bottomright',
    legend = c( -1, -0.5, 0, 0.5, 1 ),
    col = c( 'red', 'orange', 'black', 'green', 'blue' ),
    lty = 1,
    title = expression( rho )
)

## -----------------------------------------------------------------------------
p <- ( 100 : 900 ) / 1000
plot( p, qprodcor( p, -1 ), type ='l', col = 'red', ylim = c(-2,2), ylab = 'Quantile' )
lines( p, qprodcor( p, -0.5 ), col = 'orange' )
lines( p, qprodcor( p, 0 ), col = 'black' )
lines( p, qprodcor( p, 0.5 ), col = 'green' )
lines( p, qprodcor( p, 1 ), col = 'blue' )
legend(
    'bottomright',
    legend = c( -1, -0.5, 0, 0.5, 1 ),
    col = c( 'red', 'orange', 'black', 'green', 'blue' ),
    lty = 1,
    title = expression( rho )
)

## -----------------------------------------------------------------------------
x <- ( (-500) : 500 ) / 1000
plot( x, rho_from_median( x ), type ='l', ylab = 'Rho', xlab = 'Median' )
# coarse approximation slope is given by chisq median
em <- qchisq( 0.5, 1 ) # [1] 0.4549364
abline( 0, 1/em, lty = 2, col = 'gray' )

