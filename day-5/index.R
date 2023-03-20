## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = FALSE, dev = "svg", echo = TRUE, message = FALSE,
                      warning = FALSE,
                      fig.height = 6, fig.width = 1.777777 * 6)

library("gridGraphics")
library("here")
library("mgcv")
library("qgam")
library("gratia")
library("ggplot2")
library("forcats")
library("purrr")
library("mvnfast")
library("tibble")
library("patchwork")
library("tidyr")
library("knitr")
library("viridis")
library("readr")
library("dplyr")
library("sf")

## plot defaults
theme_set(theme_bw(base_size = 16, base_family = "Fira Sans"))



## ----xaringan-tile-view, echo = FALSE, eval = TRUE----------------------------
# If you don't have xaringanExtra or can't get it installed, just change
#   eval = TRUE above to eval = FALSE, you don't need it to use the slidedeck
#
# To install it, you ned the remotes pkg installed then run the next line
#   to install xaringanExtra from GitHUb
#
# remotes::install_github("gadenbuie/xaringanExtra")
#
# then this code chunk will work
xaringanExtra::use_tile_view()


## ----correlated-data-eg, fig.show = "hide"------------------------------------
set.seed(321)
n <- 100
time <- 1:n
xt <- time/n
Y <- (1280 * xt^4) * (1- xt)^4
y <- as.numeric(Y + arima.sim(list(ar = 0.3713),
                              n = n))
df <- tibble(y = y, time = time, f = Y)

# plot
plt <- ggplot(df, aes(x = time, y = y)) +
  geom_point() +
  geom_line(aes(y = f),
            col = "steelblue", lwd = 2)
plt


## ----correlated-data-eg, echo = FALSE, fig.width = 6, fig.height = 6----------
set.seed(321)
n <- 100
time <- 1:n
xt <- time/n
Y <- (1280 * xt^4) * (1- xt)^4
y <- as.numeric(Y + arima.sim(list(ar = 0.3713),
                              n = n))
df <- tibble(y = y, time = time, f = Y)

# plot
plt <- ggplot(df, aes(x = time, y = y)) +
  geom_point() +
  geom_line(aes(y = f),
            col = "steelblue", lwd = 2)
plt


## ----fit-correlated-data-eg, fig.show = "hide"--------------------------------
# standard fit
m_reml <- gam(y ~ s(time, k = 20), data = df,
              method = "REML")
# use GCV
m_gcv <- gam(y ~ s(time, k = 20), data = df)

# fitted values
fv_reml <- fitted_values(m_reml)
fv_gcv <- fitted_values(m_gcv)

# plot
plt + geom_line(data = fv_reml,
               aes(x = time, y = fitted),
               col = "red") +
  geom_line(data = fv_gcv,
            aes(x = time, y = fitted),
            col = "darkgreen")


## ----fit-correlated-data-eg, echo = FALSE, fig.width = 6, fig.height = 6------
# standard fit
m_reml <- gam(y ~ s(time, k = 20), data = df,
              method = "REML")
# use GCV
m_gcv <- gam(y ~ s(time, k = 20), data = df)

# fitted values
fv_reml <- fitted_values(m_reml)
fv_gcv <- fitted_values(m_gcv)

# plot
plt + geom_line(data = fv_reml,
               aes(x = time, y = fitted),
               col = "red") +
  geom_line(data = fv_gcv,
            aes(x = time, y = fitted),
            col = "darkgreen")


## ----fit-correlated-data-gamm, fig.show = "hide"------------------------------
# standard fit
m_ar1 <- gamm(y ~ s(time, k = 20), data = df,
              correlation = corAR1(form = ~ 1),
              method = "REML")

# fitted values
fv_ar1 <- fitted_values(m_ar1$gam)

# plot
plt +
  geom_ribbon(data = fv_ar1,
              aes(ymin = lower, ymax = upper,
                  y = NULL),
              alpha = 0.2, fill = "hotpink") +
  geom_line(data = fv_ar1,
            aes(x = time, y = fitted),
            col = "hotpink", lwd = 1.5)


## ----fit-correlated-data-gamm, echo = FALSE, fig.width = 6, fig.height = 6----
# standard fit
m_ar1 <- gamm(y ~ s(time, k = 20), data = df,
              correlation = corAR1(form = ~ 1),
              method = "REML")

# fitted values
fv_ar1 <- fitted_values(m_ar1$gam)

# plot
plt +
  geom_ribbon(data = fv_ar1,
              aes(ymin = lower, ymax = upper,
                  y = NULL),
              alpha = 0.2, fill = "hotpink") +
  geom_line(data = fv_ar1,
            aes(x = time, y = fitted),
            col = "hotpink", lwd = 1.5)


## ----gaussian-distributions-plt, echo = FALSE---------------------------------
x <- seq(8, -8, length = 500)
df <- data.frame(density = c(dnorm(x, 0, 1), dnorm(x, 0, 2), dnorm(x, 2, 1), dnorm(x, -2, 1)),
                 x = rep(x, 4),
                 distribution = factor(rep(c("mean = 0; var = 1", "mean = 0; var = 4",
                                             "mean = 2; var = 1", "mean = -2; var = 1"), each = 500),
                                       levels = c("mean = 0; var = 1", "mean = 0; var = 4",
                                                  "mean = 2; var = 1", "mean = -2; var = 1")))
plt1 <- ggplot(subset(df, distribution %in% c("mean = 0; var = 1", "mean = 0; var = 4")),
               aes(x = x, y = density, colour = distribution)) +
    geom_line(size = 1) + theme(legend.position = "top") +
    guides(col = guide_legend(title = "Distribution", nrow = 2, title.position = "left")) +
    labs(x = "x", y = "Probability density")

plt2 <- ggplot(subset(df, distribution %in% c("mean = 2; var = 1", "mean = -2; var = 1")),
               aes(x = x, y = density, colour = distribution)) +
    geom_line(size = 1) + theme(legend.position = "top") +
    guides(col = guide_legend(title = "Distribution", nrow = 2, title.position = "left")) +
    labs(x = "x", y = "Probability density")

plt <- plt1 + plt2
plt


## ----eval = FALSE-------------------------------------------------------------
## gam(list(accel ~ s(times, k = 20, bs = "ad"),
##                ~ s(times, k = 10)),
##          data = mcycle,
##          method = "REML", # <== IIRC REML is only option for these LSS
##          family = gaulss())


## ----qgam-example, results = "hide"-------------------------------------------
## Load the mcycle data
data(mcycle, package = "MASS")

## Fit QGAM using an adaptive smoother
m_q <- qgam(accel ~ s(times, k = 20, bs = "ad"),
            data = mcycle, # <== no family either
            qu = 0.8) #<<


## ----qgam-example-plt, echo = FALSE, out.width = "70%", fig.align = "center"----
## predict to visualise the estimated quantile
## new data to predict at
new_df <- with(mcycle,
               tibble(times = seq(round(min(times)), round(max(times)),
                                  length.out = 200)))
pred <- predict(m_q, newdata = new_df, se.fit = TRUE) %>%
  as.data.frame() %>%
  as_tibble() %>%
  setNames(c("est", "se")) %>%
  add_confint() %>%
  bind_cols(new_df)

## plot
pred %>%
  ggplot(aes(x = times)) +
    geom_point(data = mcycle, aes(x = times, y = accel)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
    geom_line(aes(y = est))

