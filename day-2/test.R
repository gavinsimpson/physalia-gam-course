test_fun <- function(n = 200, sx = 0.3, sz = 0.4, scale = 0.1,
                     seed = NULL) {
    require("tibble")
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        runif(1)
    }
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    }
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    xz <- runif(2 * n)
    take <- seq_len(n)
    x <- xz[take]
    z <- xz[-take]
    f <- 1.2 * exp(-(x - 0.2)^2 / sx^2 - (z - 0.3)^2 / sz^2) +
      08. * exp(-(x - 0.7)^2 / sx^2 - (z - 0.8)^2 / sz^2)
    tibble(y = f + rnorm(n) * scale, x = x, z = z)
}

test_df <- test_fun(seed = 42)

# fit a model to these data (y is response) x and z are spatial covariates
