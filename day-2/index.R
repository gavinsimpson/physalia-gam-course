## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = 'svg', echo = TRUE, message = FALSE, warning = FALSE,
                      fig.height = 6, fig.width = 1.777777 * 6)

library('here')
library('mgcv')
library('gratia')
library('ggplot2')
library('purrr')
library('mvnfast')
library("tibble")
library('patchwork')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('dplyr')
library('gganimate')

## plot defaultsFriday 14th
theme_set(theme_bw(base_size = 16, base_family = 'Fira Sans'))

## constants
anim_width <- 1000
anim_height <- anim_width / 1.77777777
anim_dev <- 'png'
anim_res <- 200


## ----hadcrut-temp-example, echo = FALSE---------------------------------------
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))

## Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_line() + 
    geom_point() +
    labs(x = 'Year', y = expression(Temeprature ~ degree * C))
gtemp_plt


## ----hadcrut-temp-example, echo = FALSE---------------------------------------
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))

## Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_line() + 
    geom_point() +
    labs(x = 'Year', y = expression(Temeprature ~ degree * C))
gtemp_plt


## ----hadcrut-temp-polynomial, echo = FALSE------------------------------------
p <- c(1,3,8,15)
N <- 300
newd <- with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N)))
polyFun <- function(i, data = data) {
    lm(Temperature ~ poly(Year, degree = i), data = data)
}
mods <- lapply(p, polyFun, data = gtemp)
pred <- vapply(mods, predict, numeric(N), newdata = newd)
colnames(pred) <- p
newd <- cbind(newd, pred)
polyDat <- gather(newd, Degree, Fitted, - Year)
polyDat <- mutate(polyDat, Degree = ordered(Degree, levels = p))
gtemp_plt +
  geom_line(data = polyDat,
            mapping = aes(x = Year, y = Fitted, colour = Degree),
            size = 1.5, alpha = 0.9) +
    scale_color_brewer(name = "Degree", palette = "PuOr") +
    theme(legend.position = "right")


## ----read-hadcrut, echo = TRUE------------------------------------------------
library('readr')
library('dplyr')
URL <-  "https://bit.ly/hadcrutv4"
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))


## ----show-hadcrut, echo = TRUE------------------------------------------------
gtemp


## ----hadcrutemp-fitted-gam, echo = TRUE, results = 'hide'---------------------
library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)


## ----hadcrutemp-fitted-gam, echo = FALSE--------------------------------------
library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)


## ----hadcrtemp-plot-gam, echo = FALSE-----------------------------------------
N <- 300
newd <- as_tibble(with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N))))
pred <- as_tibble(as.data.frame(predict(m, newdata = newd, se.fit = TRUE,
                                        unconditional = TRUE)))
pred <- bind_cols(newd, pred) %>%
    mutate(upr = fit + 2 * se.fit, lwr = fit - 2*se.fit)

ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_point() +
    geom_ribbon(data = pred,
                mapping = aes(ymin = lwr, ymax = upr, x = Year), alpha = 0.4, inherit.aes = FALSE,
                fill = "#fdb338") +
    geom_line(data = pred,
              mapping = aes(y = fit, x = Year), inherit.aes = FALSE, size = 1, colour = "#025196") +
    labs(x = 'Year', y = expression(Temeprature ~ degree*C))


## ----hadcrtemp-plot-gam, echo = FALSE-----------------------------------------
N <- 300
newd <- as_tibble(with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N))))
pred <- as_tibble(as.data.frame(predict(m, newdata = newd, se.fit = TRUE,
                                        unconditional = TRUE)))
pred <- bind_cols(newd, pred) %>%
    mutate(upr = fit + 2 * se.fit, lwr = fit - 2*se.fit)

ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_point() +
    geom_ribbon(data = pred,
                mapping = aes(ymin = lwr, ymax = upr, x = Year), alpha = 0.4, inherit.aes = FALSE,
                fill = "#fdb338") +
    geom_line(data = pred,
              mapping = aes(y = fit, x = Year), inherit.aes = FALSE, size = 1, colour = "#025196") +
    labs(x = 'Year', y = expression(Temeprature ~ degree*C))


## ----smooth-fun-animation, results = FALSE, echo = FALSE----------------------
f <- function(x) {
    x^11 * (10 * (1 - x))^6 + ((10 * (10 * x)^3) * (1 - x)^10)
}

draw_beta <- function(n, k, mu = 1, sigma = 1) {
    rmvn(n = n, mu = rep(mu, k), sigma = diag(rep(sigma, k)))
}

weight_basis <- function(bf, x, n = 1, k, ...) {
    beta <- draw_beta(n = n, k = k, ...)
    out <- sweep(bf, 2L, beta, '*')
    colnames(out) <- paste0('f', seq_along(beta))
    out <- as_tibble(out)
    out <- add_column(out, x = x)
    out <- pivot_longer(out, -x, names_to = 'bf', values_to = 'y')
    out
}

random_bases <- function(bf, x, draws = 10, k, ...) {
    out <- rerun(draws, weight_basis(bf, x = x, k = k, ...))
    out <- bind_rows(out)
    out <- add_column(out, draw = rep(seq_len(draws), each = length(x) * k),
                      .before = 1L)
    class(out) <- c("random_bases", class(out))
    out
}

plot.random_bases <- function(x, facet = FALSE) {
    plt <- ggplot(x, aes(x = x, y = y, colour = bf)) +
        geom_line(lwd = 1, alpha = 0.75) +
        guides(colour = FALSE)
    if (facet) {
        plt + facet_wrap(~ draw)
    }
    plt
}

normalize <- function(x) {
    rx <- range(x)
    z <- (x - rx[1]) / (rx[2] - rx[1])
    z
}

set.seed(1)
N <- 500
data <- tibble(x     = runif(N),
               ytrue = f(x),
               ycent = ytrue - mean(ytrue),
               yobs  = ycent + rnorm(N, sd = 0.5))

k <- 10
knots <- with(data, list(x = seq(min(x), max(x), length = k)))
sm <- smoothCon(s(x, k = k, bs = "cr"), data = data, knots = knots)[[1]]$X
colnames(sm) <- levs <- paste0("f", seq_len(k))
basis <- pivot_longer(cbind(sm, data), -(x:yobs), names_to = 'bf')
basis

set.seed(2)
bfuns <- random_bases(sm, data$x, draws = 20, k = k)

smooth <- bfuns %>%
    group_by(draw, x) %>%
    summarise(spline = sum(y)) %>%
    ungroup()

p1 <- ggplot(smooth) +
    geom_line(data = smooth, aes(x = x, y = spline), lwd = 1.5) +
    labs(y = 'f(x)', x = 'x') +
    theme_minimal(base_size = 16, base_family = 'Fira Sans')

smooth_funs <- animate(
    p1 + transition_states(draw, transition_length = 4, state_length = 2) + 
    ease_aes('cubic-in-out'),
    nframes = 200, height = anim_height, width = anim_width, res = anim_res, dev = anim_dev)

anim_save('resources/spline-anim.gif', smooth_funs)


## ----basis-functions, fig.height=6, fig.width = 1.777777*6, echo = FALSE------
ggplot(basis,
       aes(x = x, y = value, colour = bf)) +
    geom_line(lwd = 2, alpha = 0.5) +
    guides(colour = FALSE) +
    labs(x = 'x', y = 'b(x)') +
    theme_minimal(base_size = 20, base_family = 'Fira Sans')


## ----basis-function-animation, results = 'hide', echo = FALSE-----------------
bfun_plt <- plot(bfuns) +
    geom_line(data = smooth, aes(x = x, y = spline),
              inherit.aes = FALSE, lwd = 1.5) +
    labs(x = 'x', y = 'f(x)') +
    theme_minimal(base_size = 14, base_family = 'Fira Sans')

bfun_anim <- animate(
    bfun_plt + transition_states(draw, transition_length = 4, state_length = 2) + 
    ease_aes('cubic-in-out'),
    nframes = 200, height = anim_height, width = anim_width, res = anim_res, dev = anim_dev)

anim_save('resources/basis-fun-anim.gif', bfun_anim)


## ----example-data-figure, fig.height=6, fig.width = 1.777777*6, echo = FALSE----
data_plt <- ggplot(data, aes(x = x, y = ycent)) +
    geom_line(col = 'goldenrod', lwd = 2) +
    geom_point(aes(y = yobs), alpha = 0.2, size = 3) +
    labs(x = 'x', y = 'f(x)') +
    theme_minimal(base_size = 20, base_family = 'Fira Sans')
data_plt


## ----basis-functions-anim, results = "hide", echo = FALSE---------------------
sm2 <- smoothCon(s(x, k = k, bs = "cr"), data = data, knots = knots)[[1]]$X
beta <- coef(lm(ycent ~ sm2 - 1, data = data))
wtbasis <- sweep(sm2, 2L, beta, FUN = "*")
colnames(wtbasis) <- colnames(sm2) <- paste0("F", seq_len(k))
## create stacked unweighted and weighted basis
basis <- as_tibble(rbind(sm2, wtbasis)) %>%
    add_column(x = rep(data$x, times = 2),
               type = rep(c('unweighted', 'weighted'), each = nrow(sm2)),
               .before = 1L)
##data <- cbind(data, fitted = rowSums(scbasis))
wtbasis <- as_tibble(rbind(sm2, wtbasis)) %>%
    add_column(x      = rep(data$x, times = 2),
               fitted = rowSums(.),
               type   = rep(c('unweighted', 'weighted'), each = nrow(sm2))) %>%
    pivot_longer(-(x:type), names_to = 'bf')
basis <- pivot_longer(basis, -(x:type), names_to = 'bf')

p3 <- ggplot(data, aes(x = x, y = ycent)) +
    geom_point(aes(y = yobs), alpha = 0.2) +
    geom_line(data = basis,
              mapping = aes(x = x, y = value, colour = bf),
              lwd = 1, alpha = 0.5) +
    geom_line(data = wtbasis,
              mapping = aes(x = x, y = fitted), lwd = 1, colour = 'black', alpha = 0.75) +
    guides(colour = FALSE) +
    labs(y = 'f(x)', x = 'x') +
    theme_minimal(base_size = 16, base_family = 'Fira Sans')

crs_fit <- animate(p3 + transition_states(type, transition_length = 4, state_length = 2) + 
                   ease_aes('cubic-in-out'),
                   nframes = 100, height = anim_height, width = anim_width, res = anim_res,
                   dev = anim_dev)

anim_save('./resources/gam-crs-animation.gif', crs_fit)


## ----hadcrut-temp-penalty, echo = FALSE---------------------------------------
K <- 40
lambda <- c(10000, 1, 0.01, 0.00001)
N <- 300
newd <- with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N)))
fits <- lapply(lambda, function(lambda) gam(Temperature ~ s(Year, k = K, sp = lambda), data = gtemp))
pred <- vapply(fits, predict, numeric(N), newdata = newd)
op <- options(scipen = 100)
colnames(pred) <- lambda
newd <- cbind(newd, pred)
lambdaDat <- gather(newd, Lambda, Fitted, - Year)
lambdaDat <- transform(lambdaDat, Lambda = factor(paste("lambda ==", as.character(Lambda)),
                                                  levels = paste("lambda ==", as.character(lambda))))

gtemp_plt + geom_line(data = lambdaDat, mapping = aes(x = Year, y = Fitted, group = Lambda),
                      size = 1, colour = "#e66101") +
    facet_wrap( ~ Lambda, ncol = 2, labeller = label_parsed)
options(op)


## ----birds-1, echo = TRUE-----------------------------------------------------
library("here"); library("readr"); library("dplyr")
larks <-  read_csv(here("data", "larks.csv"),
                   col_types = "ccdddd")

larks <- larks %>%
  mutate(crestlark = factor(crestlark),
         linnet = factor(linnet),
         e = x / 1000,
         n = y / 1000)
head(larks)


## ----birds-2, fig.width = 5, fig.height = 6, echo = FALSE---------------------
ggplot(larks, aes(x = e, y = n, colour = crestlark)) + geom_point(size = 0.5) +
  coord_fixed() + scale_colour_discrete(na.value = '#bbbbbb33') +
  labs(x = NULL, y = NULL)


## ----birds-gam-1, echo = TRUE-------------------------------------------------
crest <- gam(crestlark ~ s(e, n, k = 100),
             data = larks,
             family = binomial,
             method = 'REML')


## ----birds-gam-2, echo = TRUE-------------------------------------------------
summary(crest)


## ----munge-larks, echo = TRUE-------------------------------------------------
larks2 <- larks %>%
  mutate(crestlark = as.numeric(as.character(crestlark)),  # to numeric
         linnet = as.numeric(as.character(linnet)),
         tet_n = rep(1, nrow(larks)),                      # counter for how many grid cells (tetrads) we sum
         N = rep(1, nrow(larks)),                          # number of observations, 1 per row currently
         N = if_else(is.na(crestlark), NA_real_, N)) %>%   # set N to NA if observation taken
  group_by(QUADRICULA) %>%                                 # group by the larger grid square
  summarise(across(c(N, crestlark, linnet, tet_n, e, n),
                   ~ sum(., na.rm = TRUE))) %>%            # sum all needed variables
  mutate(e = e / tet_n, n = n / tet_n)                     # rescale to get avg E,N coords for QUADRICULA

## fit binomial GAM
crest2 <- gam(cbind(crestlark, N - crestlark) ~ s(e, n, k = 100),
              data = larks2, family = binomial, method = 'REML')


## ----crest-3, echo = TRUE-----------------------------------------------------
crest3 <- gam(cbind(crestlark, N - crestlark) ~
                s(e, n, k = 100),
              data = larks2, family = quasibinomial,
              method = 'REML')
crest3$scale # should be == 1


## ----gam-check-aggregated-lark, echo = TRUE, fig.width = 4.5, fig.height = 4----
ggplot(tibble(Fitted = fitted(crest2),
              Resid = resid(crest2)),
       aes(Fitted, Resid)) + geom_point()


## ----co2-example-1, echo = TRUE-----------------------------------------------
library("gratia")
south <- read_csv(here("data", "south_pole.csv"), col_types = "ddd")
south


## ----co2-example-2, echo = TRUE, fig.align = "center", out.width = "95%"------
ggplot(south, aes(x = c.month, y = co2)) + geom_line()


## ----co2-example-3, echo = TRUE-----------------------------------------------
m_co2 <- gam(co2 ~ s(c.month, k = 300, bs = "cr"), data = south, method = 'REML')
summary(m_co2)


## ----co2-example-4, echo = TRUE-----------------------------------------------
new_df <- with(south, tibble(c.month = 1:(nrow(south) + 36)))
fv <- fitted_values(m_co2, data = new_df, scale = "response")
fv


## ----co2-example-5, echo = TRUE, fig.height = 5-------------------------------
ggplot(fv, aes(x = c.month, y = fitted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(data = south, aes(c.month, co2), col = 'red') +
    geom_line(alpha = 0.4)


## ----co2-example-6, echo = TRUE-----------------------------------------------
m2_co2 <- gam(co2 ~ s(month, bs = "cc") + s(c.month, bs = "cr", k = 300),
              data = south, method = 'REML',
              knots = list(month = c(0.5, 12.5)))


## ----co2-example-7, echo = TRUE-----------------------------------------------
summary(m2_co2)


## ----co2-example-8, echo = TRUE-----------------------------------------------
nr <- nrow(south)
new_df <- with(south,
               tibble(c.month = 1:(nr + 36),
                      month = rep(seq_len(12), length = nr + 36)))
fv2 <- fitted_values(m2_co2, data = new_df, scale = "response")
fv2


## ----co2-example-9, echo = TRUE, fig.height = 5-------------------------------
ggplot(fv2, aes(x = c.month, y = fitted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(data = south, aes(c.month, co2), col = 'red') +
    geom_line(alpha = 0.4)


## -----------------------------------------------------------------------------
library("dplyr")
data(mcycle, package = "MASS")
mcycle <- as_tibble(mcycle)
mcycle


## ----plot-mcycle, eval = FALSE------------------------------------------------
## plt <- mcycle %>%
##   ggplot(aes(x = times, y = accel)) +
##     geom_point() +
##     labs(x = "Milliseconds after impact",
##          y = expression(italic(g)))
## plt


## ----plot-mcycle, fig.width = 6, fig.height = 4, echo = FALSE-----------------
plt <- mcycle %>%
  ggplot(aes(x = times, y = accel)) +
    geom_point() +
    labs(x = "Milliseconds after impact",
         y = expression(italic(g)))
plt


## -----------------------------------------------------------------------------
new_df <- with(mcycle, tibble(times = evenly(times, n = 100)))
bfun <- basis(s(times), data = new_df)
bfun


## ----draw-tprs-basis, fig.show = "hide"---------------------------------------
draw(bfun)


## ----draw-tprs-basis, out.width = "95%", echo = FALSE-------------------------
draw(bfun)


## ----draw-tprs-basis-facetted, fig.show = "hide"------------------------------
draw(bfun) + facet_wrap(~ bf)


## ----draw-tprs-basis-facetted, out.width = "90%", echo = FALSE----------------
draw(bfun) + facet_wrap(~ bf)


## ---- out.width = "90%", fig.align = "center"---------------------------------
bfun <- basis(s(times, bs = "cr"), data = new_df)
draw(bfun) + facet_wrap(~ bf)


## ---- echo = FALSE, out.width = "95%", fig.align = "center"-------------------
K <- 7
knots <- with(mcycle, list(times = evenly(times, n = K)))
bfun <- basis(s(times, bs = "cr", k = K), data = new_df, knots = knots)
draw(bfun) + facet_wrap(~ bf) +
  geom_rug(data = as.data.frame(knots), aes(x = times), sides = "b", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, alpha = 0.5)


## -----------------------------------------------------------------------------
## how many basis functions?
K <- 7

## create the knots list
knots <- with(mcycle,
              list(times = evenly(times, n = K)))

## provide `knots` to functions
bfun <- basis(s(times, bs = "cr", k = K), data = new_df, knots = knots)

model <- gam(accel ~ s(times, bs = "cr", k = K),
             data = mcycle, method = "REML",
             knots = knots) # <- specify knots to the model


## ----cc-basis-default, results = FALSE----------------------------------------
month_df <- with(south, tibble(month = evenly(month, n = 100)))
bfun <- basis(s(month, bs = "cc"), data = month_df)


## ----draw-crs-basis-facetted, fig.align = "center", out.width = "95%"---------
draw(bfun) + facet_wrap(~ bf)


## ----plot-cc-basis, eval = FALSE----------------------------------------------
## knots <- list(month = c(0.5, 12.5))
## bfun <- basis(s(month, bs = "cc"), data = month_df, knots = knots)
## draw(bfun) + facet_wrap(~ bf)


## ----plot-cc-basis, echo = FALSE----------------------------------------------
knots <- list(month = c(0.5, 12.5))
bfun <- basis(s(month, bs = "cc"), data = month_df, knots = knots)
draw(bfun) + facet_wrap(~ bf)


## ----whole-basis-proces, echo = FALSE, fig.height = 4, fig.width = 1.777777 * 6----
K <- 13
df <- data.frame(x = seq(0, 1, length = 200))
knots <- data.frame(x = seq(0, 1, length.out = 11))
bs <- basis(s(x, bs = "ps", k = K), data = df,
    knots = list(x = seq(-3, 13) / 10))

# let's weight the basis functions (simulating model coefs)
set.seed(1)
betas <- data.frame(bf = factor(seq_len(K)), beta = rnorm(K))

unwtd_bs_plt <- bs |>
    draw() +
    geom_vline(aes(xintercept = x), data = knots, linetype = "dotted",
        alpha = 0.5)

# we need to merge the weights for each basis function with the basis object
bs <- bs |>
    left_join(betas, by = join_by("bf" == "bf")) |>
    mutate(value_w = value * beta)

# weighted basis
wtd_bs_plt <- bs |>
    ggplot(aes(x = x, y = value_w, colour = bf, group = bf)) +
    geom_line(show.legend = FALSE) +
    geom_vline(aes(xintercept = x), data = knots, linetype = "dotted",
        alpha = 0.5) +
    labs(y = expression(f(x)), x = "x")

# now we want to sum the weighted basis functions for each value of `x`
spl <- bs |>
    group_by(x) |>
    summarise(spline = sum(value_w))

take <- c(83, 115)
pts <- bs |>
    group_by(bf) |>
    slice(take)

# now plot
bs_plt <- bs |>
    ggplot(aes(x = x, y = value_w, colour = bf, group = bf)) +
    geom_line(show.legend = FALSE) +
    geom_line(aes(x = x, y = spline), data = spl, linewidth = 1.25,
              inherit.aes = FALSE) +
    geom_vline(aes(xintercept = x), data = knots, linetype = "dotted",
        alpha = 0.5) +
    geom_vline(xintercept = c(df$x[take]), linetype = "dashed",
        alpha = 1) +
    geom_point(data = pts, aes(x = x, y = value_w, colour = bf, group = bf),
        size = 2, show.legend = FALSE) +
    geom_point(data = slice(spl, take), aes(x = x, y = spline),
        size = 3, colour = "red", inherit.aes = FALSE) +
    labs(y = expression(f(x)), x = "x")

unwtd_bs_plt + wtd_bs_plt + bs_plt + plot_layout(ncol = 3)


## ----whole-basis-proces-2-model-----------------------------------------------
dat <- data_sim("eg1", seed = 4)
m <- gam(y ~ s(x0) + s(x1) + s(x2, bs = "bs") + s(x3),
         data = dat, method = "REML")


## ----whole-basis-proces-2-model-draw, fig.height = 5, fig.width = 1.777777 * 6----
draw(m) + plot_layout(ncol = 4)


## ----whole-basis-proces-2, echo = FALSE, fig.height = 6, fig.width = 1.777777 * 6----
# data to evaluate the basis at
# using the CRAN version of {gratia}, we need `m`
ds <- data_slice(m, x2 = evenly(x2, n = 200))
# from 0.9.0 (or current GitHub version) you can do
# ds <- data_slice(dat, x2 = evenly(x2, n = 200))

# generate a tidy representation of the fitted basis functions
x2_bs <- basis(m, term = "s(x2)", data = ds)

# compute values of the spline by summing basis functions at each x2
x2_spl <- x2_bs |>
    group_by(x2) |>
    summarise(spline = sum(value))

# evaluate the spline at the same values as we evaluated the basis functions
x2_sm <- smooth_estimates(m, "s(x2)", data = ds) |>
    add_confint()

take <- c(65, 175)
pts <- x2_bs |>
    group_by(bf) |>
    slice(take)

# now plot
x2_bs |>
    ggplot(aes(x = x2, y = value, colour = bf, group = bf)) +
    geom_line(show.legend = FALSE) +
    geom_ribbon(aes(x = x2, ymin = lower_ci, ymax = upper_ci),
                data = x2_sm,
                inherit.aes = FALSE, alpha = 0.2) +
    geom_line(aes(x = x2, y = est), data = x2_sm,
              linewidth = 1.5, inherit.aes = FALSE) +
    geom_vline(xintercept = c(ds$x2[take]), linetype = "dashed",
        alpha = 1) +
    geom_point(data = pts, aes(x = x2, y = value, colour = bf, group = bf),
        size = 2, show.legend = FALSE) +
    geom_point(data = slice(x2_sm, take), aes(x = x2, y = est),
        size = 3, colour = "red", inherit.aes = FALSE) +
    labs(y = expression(f(x2)), x = "x2")


## ---- out.width = "85%", fig.align = "center"---------------------------------
new_df <- with(mcycle, tibble(times = evenly(times, n = 100)))
bfun <- basis(s(times), data = new_df, constraints = TRUE)
draw(bfun) + facet_wrap(~ bf)


## -----------------------------------------------------------------------------
m <-  gam(accel ~ s(times), dat = mcycle, method = "REML")
S <- penalty(m, smooth = "s(times)")
S


## ---- echo = FALSE, dev = "png", out.width = "100%", dpi = 300----------------
library("patchwork")
p1 <- draw(bfun) + facet_wrap(~ bf)
p2 <- draw(S)
p1 + p2 + plot_layout(ncol = 2)


## -----------------------------------------------------------------------------
m <-  gam(accel ~ s(times, bs = "cr"), data = mcycle, method = "REML")
S <- penalty(m, smooth = "s(times)")
# draw(S)


## ---- echo = FALSE, dev = "png", out.width = "100%", dpi = 300----------------
library("patchwork")
bfun_cc <- basis(s(times, bs = "cr"), data = new_df, constraints = TRUE)
p1 <- draw(bfun_cc) + facet_wrap(~ bf)
p2 <- draw(S)
p1 + p2 + plot_layout(ncol = 2)


## ----draw-mcycle, out.width = "90%", fig.align = "center"---------------------
m <- gam(accel ~ s(times), data = mcycle, method = "REML")
draw(m)


## ----draw-four-fun-sim, fig.show = "hide"-------------------------------------
df <- data_sim("eg1", n = 1000, seed = 42)
df
m_sim <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
             data = df, method = "REML")
draw(m_sim)


## ----draw-four-fun-sim-plot, echo = FALSE, out.width = "95%"------------------
draw(m_sim)


## ----draw-mcycle-options, fig.show = "hide"-----------------------------------
draw(m_sim,
     residuals = TRUE,           # add partial residuals
     overall_uncertainty = TRUE, # include uncertainty due to constant
     unconditional = TRUE,       # correct for smoothness selection
     rug = FALSE)                # turn off rug plot


## ----draw-mcycle-options-plot, out.width = "95%", echo = FALSE----------------
draw(m_sim, residuals = TRUE, overall_uncertainty = TRUE,
     unconditional = TRUE, rug = FALSE)

