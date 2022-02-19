## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = "svg", echo = TRUE, message = FALSE,
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


## ----load-shrimp--------------------------------------------------------------
shrimp <- read.csv(here("data", "trawl_nl.csv"))


## ----shrimp-richness----------------------------------------------------------
m_rich <- gam(richness ~ s(year),
              family = poisson,
              method = "REML",
              data = shrimp)


## ----richness-violin, fig.height=5, fig.width=5, echo=FALSE-------------------
ggplot(shrimp) +
  geom_violin(aes(x = richness, y = factor(year))) +
    labs(x = "Number of species", y = "Year")


## ----draw-richness-gam, out.width = "90%", fig.align = "center"---------------
draw(m_rich)


## ----biom-space-time-plot, fig.height=8, fig.width=15, echo=FALSE, dev="png", dpi = 300----
coast <- read_sf(here("data", "nl_coast.shp"))
ggplot(shrimp) +
  geom_point(aes(x = long, y = lat, size = shrimp), alpha = 0.5) +
  geom_sf(data = coast) +
  facet_wrap(~year, ncol = 5)


## ----fit-shrimp-space-time----------------------------------------------------
m_spt <- gam(shrimp ~ te(x, y, year, d = c(2,1),
                         bs = c("tp", "cr"), k = c(20, 5)),
             data = shrimp,
             family = tw(),
             method = "REML")


## ----predict-newdata----------------------------------------------------------
new_year <- with(shrimp, tibble(year = seq(min(year), max(year), length.out = 100)))
pred <- predict(m_rich, newdata = new_year, se.fit = TRUE, type = "link")
pred <- bind_cols(new_year, as_tibble(as.data.frame(pred)))
pred


## ----predict-newdata-resp-----------------------------------------------------
ilink <- inv_link(m_rich)                         # inverse link function
crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE) # or just `crit <- 2`
pred <- mutate(pred, richness = ilink(fit),
               lwr = ilink(fit - (crit * se.fit)), # lower...
               upr = ilink(fit + (crit * se.fit))) # upper credible interval
pred


## ----plot-predictions-richness, fig.height = 4--------------------------------
ggplot(pred, aes(x = year)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(aes(y = richness)) + labs(y = "Species richness", x = NULL)


## ----spt-example-predict------------------------------------------------------
sp_new <- with(shrimp, expand.grid(x = seq_min_max(x, n = 100), y = seq_min_max(y, n = 100),
                                   year = unique(year)))
sp_pred <- predict(m_spt, newdata = sp_new, se.fit = TRUE) # link scale is default
sp_pred <- bind_cols(as_tibble(sp_new), as_tibble(as.data.frame(sp_pred)))
sp_pred


## ----spt-example-response-scale-----------------------------------------------
ilink <- inv_link(m_spt)
too_far <- exclude.too.far(sp_pred$x, sp_pred$y, shrimp$x, shrimp$y, dist = 0.1)
sp_pred <- sp_pred %>% mutate(biomass = ilink(fit),
                              biomass = case_when(too_far ~ NA_real_,
                                                  TRUE ~ biomass))
sp_pred


## ----spt-example-plot, fig.height = 5.5, dev="png", dpi = 300-----------------
ggplot(sp_pred, aes(x = x, y = y, fill = biomass)) + geom_raster() +
    scale_fill_viridis_c(option = "plasma") + facet_wrap(~ year, ncol = 5) + coord_equal()


## ----show-m-spt---------------------------------------------------------------
m_spt


## ----shrimp-ti-model----------------------------------------------------------
m_ti <- gam(shrimp ~ ti(x, y, year, d = c(2, 1), bs = c("tp", "cr"), k = c(20, 5)) +
                s(x, y, bs = "tp", k = 20) +
                s(year, bs = "cr", k = 5),
            data = shrimp, family = tw, method = "REML")


## ----summary-spt-ti-----------------------------------------------------------
smooths(m_ti)


## ----pred-data-ti-model-------------------------------------------------------
ti_new <- with(shrimp, expand.grid(x = mean(x), y = mean(y), year = seq_min_max(year, n = 100)))

ti_pred <- predict(m_ti, newdata = ti_new, se.fit = TRUE,
                   exclude = c("ti(x,y,year)", "s(x,y)")) #<<

ti_pred <- bind_cols(as_tibble(ti_new), as_tibble(as.data.frame(ti_pred))) %>%
    mutate(biomass = ilink(fit),
           lwr = ilink(fit - (crit * se.fit)),
           upr = ilink(fit + (crit * se.fit)))


## ----pred-data-ti-model-terms, results = "hide"-------------------------------
predict(m_ti, newdata = ti_new, se.fit = TRUE, terms = "s(year)")


## ----plot-ti-marginal-trend, fig.height = 5-----------------------------------
ggplot(ti_pred, aes(x = year)) + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(aes(y = biomass)) + labs(y = "Biomass", x = NULL)


## ----predict-via-fitted-values, out.width = "70%", fig.align = "center"-------
ti_pred2 <- fitted_values(m_ti, data = ti_new,
                          scale = "response",
                          exclude = c("ti(x,y,year)", "s(x,y)")) #<<

ggplot(ti_pred2, aes(x = year)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(aes(y = fitted)) + labs(y = "Biomass", x = NULL)


## ----echo = FALSE, out.width = "80%"------------------------------------------
knitr::include_graphics(here("day-3/resources",
                             "miller-bayesian-gam-interpretation-fig.svg"))


## ----richness-coefs-----------------------------------------------------------
sm_year <- get_smooth(m_rich, "s(year)") # extract the smooth object from model
idx <- gratia:::smooth_coefs(sm_year)    # indices of the coefs for this smooth
idx

beta <- coef(m_rich)                     # vector of model parameters
beta[idx]                                # coefs for this smooth


## ----richness-vcov, results = "hide"------------------------------------------
Vb <- vcov(m_rich) # default is the bayesian covariance matrix
Vb


## ----richness-vcov-print, echo = FALSE----------------------------------------
op <- options(width = 170)
Vb
options(op)


## ----richness-xp-matrix-------------------------------------------------------
new_year <- with(shrimp, tibble(year = seq_min_max(year, n = 100)))
Xp <- predict(m_rich, newdata = new_year, type = "lpmatrix")
dim(Xp)


## ----richness-reduce-xp-------------------------------------------------------
Xp <- Xp[, idx, drop = FALSE]
dim(Xp)


## ----richness-simulate-params-------------------------------------------------
set.seed(42)
beta_sim <- rmvn(n = 20, beta[idx], Vb[idx, idx, drop = FALSE])
dim(beta_sim)


## ----richness-posterior-draws, fig.height = 5, fig.show = "hide"--------------
sm_draws <- Xp %*% t(beta_sim)
dim(sm_draws)
matplot(sm_draws, type = "l")


## ----richness-posterior-draws, fig.height = 5, fig.width = 5, echo = FALSE, results = "hide"----
sm_draws <- Xp %*% t(beta_sim)
dim(sm_draws)
matplot(sm_draws, type = "l")


## ----plot-posterior-smooths, fig.height = 5-----------------------------------
sm_post <- smooth_samples(m_rich, "s(year)", n = 20, seed = 42)
draw(sm_post)


## ----posterior-sim-model------------------------------------------------------
beta <- coef(m_rich)   # vector of model parameters
Vb <- vcov(m_rich)     # default is the bayesian covariance matrix
Xp <- predict(m_rich, type = "lpmatrix")
set.seed(42)
beta_sim <- rmvn(n = 1000, beta, Vb) # simulate parameters
eta_p <- Xp %*% t(beta_sim)        # form linear predictor values
mu_p <- inv_link(m_rich)(eta_p)    # apply inverse link function

mean(mu_p[1, ]) # mean of posterior for the first observation in the data
quantile(mu_p[1, ], probs = c(0.025, 0.975))


## ----posterior-sim-model-hist, fig.height = 5---------------------------------
ggplot(tibble(richness = mu_p[587, ]), aes(x = richness)) +
    geom_histogram() + labs(title = "Posterior richness for obs #587")


## ----richness-fitted-samples, fig.height = 4.5--------------------------------
rich_post <- fitted_samples(m_rich, n = 1000, newdata = shrimp, seed = 42)
ggplot(filter(rich_post, row == 587), aes(x = fitted)) +
    geom_histogram() + labs(title = "Posterior richness for obs #587", x = "Richness")


## ----total-biomass-posterior-1------------------------------------------------
sp_new <- with(shrimp, expand.grid(x = seq_min_max(x, n = 100), y = seq_min_max(y, n = 100),
                                   year = 2007))
Xp <- predict(m_spt, newdata = sp_new, type = "lpmatrix")

## work out now which points are too far now
too_far <- exclude.too.far(sp_new$x, sp_new$y, shrimp$x, shrimp$y, dist = 0.1)

beta <- coef(m_spt)                  # vector of model parameters
Vb <- vcov(m_spt)                    # default is the bayesian covariance matrix
set.seed(42)
beta_sim <- rmvn(n = 1000, beta, Vb) # simulate parameters
eta_p <- Xp %*% t(beta_sim)          # form linear predictor values
mu_p <- inv_link(m_spt)(eta_p)       # apply inverse link function


## ----total-biomass-posterior-2, dependson = -1--------------------------------
mu_copy <- mu_p              # copy mu_p
mu_copy[too_far, ] <- NA     # set cells too far from data to be NA
total_biomass <- colSums(mu_copy, na.rm = TRUE)  # total biomass over the region

mean(total_biomass)
quantile(total_biomass, probs = c(0.025, 0.975))


## ----total-biomass-histogram, echo = FALSE------------------------------------
ggplot(tibble(biomass = total_biomass), aes(x = biomass)) +
    geom_histogram()


## ----biomass-fitted-samples-example-------------------------------------------
bio_post <- fitted_samples(m_spt, n = 1000,
                           newdata = sp_new[!too_far, ],
                           seed = 42) %>%
    group_by(draw) %>%
    summarise(total = sum(fitted),
              .groups = "drop_last")

with(bio_post, mean(total))
with(bio_post, quantile(total, probs = c(0.025, 0.975)))


## ----biomass-fitted-samples-plot, fig.width = 5, fig.height = 5---------------
ggplot(bio_post, aes(x = total)) +
    geom_histogram() +
    labs(x = "Total biomass")


## -----------------------------------------------------------------------------
spp_url <- "https://bit.ly/spp-gradient"
gradient <- read_csv(spp_url, col_types = "dd")
gradient


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

