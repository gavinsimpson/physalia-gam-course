## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = "svg", echo = TRUE, message = FALSE,
                      warning = FALSE,
                      fig.height = 6, fig.width = 1.777777 * 6)

library("gridGraphics")
library('here')
library('mgcv')
library('gratia')
library('ggplot2')
library('forcats')
library('purrr')
library('mvnfast')
library("tibble")
library('patchwork')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('dplyr')
library('sf')

## plot defaults
theme_set(theme_bw(base_size = 16, base_family = 'Fira Sans'))



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


## ----echo = FALSE, out.width = "90%"------------------------------------------
knitr::include_graphics("resources/miller-bayesian-gam-interpretation-fig.svg")


## ----setup-confint-example, fig = TRUE, fig.width = 11, fig.height = 5.5, results = "hide", echo = FALSE----
library(mgcv)
set.seed(0)
## fake some data...
f1 <- function(x) {exp(2 * x)}
f2 <- function(x) { 
  0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10 
}
f3 <- function(x) {x*0}

n<-200
sig2 <- 12
x0 <- rep(1:4,50)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
e <- rnorm(n, 0, sqrt(sig2))
y <- 2*x0 + f1(x1) + f2(x2) + f3(x3) + e
x0 <- factor(x0)

## fit and plot...
b <- gam(y ~ x0 + s(x1) + s(x2) + s(x3))

op <- par(mar = c(4,4,1,1) + 0.1)
layout(matrix(1:9, ncol = 3, byrow = TRUE))
curve(f1)
curve(f2)
curve(f3)
plot(b, shade=TRUE)
plot(b, shade = TRUE, seWithMean = TRUE) ## better coverage intervals
layout(1)
par(op)


## ----draw-coverage-bands-closup, echo = FALSE---------------------------------
f1_fun <- function(x) {
  y <- f1(x)
  y - mean(y)
}
p1 <- draw(b, select = c("s(x1)"), overall_uncertainty = FALSE) +
  labs(title = "s(x1): without uncertainty in constant") +
  geom_function(fun = f1_fun, colour = "red")
p2 <- draw(b, select = c("s(x1)")) +
  labs(title = "s(x1): with uncertainty in constant") +
  geom_function(fun = f1_fun, colour = "red")

p1 + p2 + plot_layout(ncol = 2)


## ----draw-coverage-bands-closup-x3, echo = FALSE------------------------------
f3_fun <- function(x) {
  y <- f3(x)
  y - mean(y)
}
p1 <- draw(b, select = c("s(x3)"), overall_uncertainty = FALSE) +
  labs(title = "s(x3): without uncertainty in constant") +
  geom_function(fun = f3_fun, colour = "red")
p2 <- draw(b, select = c("s(x3)")) +
  labs(title = "s(x3): with uncertainty in constant") +
  geom_function(fun = f3_fun, colour = "red")

p1 + p2 + plot_layout(ncol = 2)


## ----aic-models-setup, echo = FALSE-------------------------------------------
n <- 200
dat <- data_sim("eg1", n = n, scale = .15, dist = "poisson", seed = 3)
set.seed(22)
dat <- dat %>% mutate(x4 = runif(n, 0, 1), x5 = runif(n, 0, 1),
                      f4 = rep(0, n), f5 = rep(0, n))   ## spurious


## ----shrinkage-example-summary, results = "hide"------------------------------
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) +
           s(x4) + s(x5),
         data = dat,
         method = "REML",
         select = TRUE)

summary(m) # ==>


## ----shrinkage-example-summary, echo = FALSE----------------------------------
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) +
           s(x4) + s(x5),
         data = dat,
         method = "REML",
         select = TRUE)

summary(m) # ==>


## ----aic-models, dependson="aic-models-setup"---------------------------------
b0 <- gam(y ~ s(x0) + s(x1) + s(x2),
          data = dat, family = poisson, method = "REML")
b1 <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) + s(x4) + s(x5),
          data = dat, family = poisson, method = "REML", select = TRUE)
b2 <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) + s(x4) + s(x5),
          data = dat, family = poisson, method = "REML")


## ----aic-example, echo = TRUE, dependson = -1---------------------------------
AIC(b0, b1, b2)


## ----aic-chisq, echo = TRUE---------------------------------------------------
pchisq(2, 1, lower.tail = FALSE)


## ----ranefs-------------------------------------------------------------------
m_nlme <- lme(travel ~ 1, data = Rail, ~ 1 | Rail, method = "REML")

m_gam  <- gam(travel ~ s(Rail, bs = "re"), data = Rail, method = "REML")


## -----------------------------------------------------------------------------
head(Rail)


## -----------------------------------------------------------------------------
m_nlme <- lme(travel ~ 1, data = Rail, ~ 1 | Rail, method = "REML")

m_gam  <- gam(travel ~ s(Rail, bs = "re", k = 2), data = Rail, method = "REML")

unlist(c(fixef(m_nlme), ranef(m_nlme)))
coef(m_gam)


## ----variance-comp-nlme-------------------------------------------------------
m_nlme


## ----variance-comp-gam--------------------------------------------------------
variance_comp(m_gam)


## ----re-basis, fig.show = "hide"----------------------------------------------
pm <- penalty(m_gam, smooth = "s(Rail)")
draw(pm)


## ----re-basis, echo = FALSE, dev = "png", dpi = 360, fig.align = "center", fig.width = 6.5, fig.height = 6----
pm <- penalty(m_gam, smooth = "s(Rail)")
draw(pm)


## ----setup-rat-hormone-example, echo = FALSE----------------------------------
rats_url <- "https://bit.ly/rat-hormone"
rats <- read_table(rats_url, col_types = "dddddddddddd-")
# ignore the warning - it"s due to trailing white space at the ends of each
#   row in the file

rats <- rats %>%
    mutate(treatment = fct_recode(factor(group, levels = c(1, 2, 3)),
                                  Low = "1",
                                  High = "2",
                                  Control = "3"),
           treatment = fct_relevel(treatment, c("Control", "Low", "High")),
           subject = factor(subject))

plt_labs <- labs(y = "Head height (distance in pixels)",
                 x = "Age in days",
                 colour = "Treatment")

rat_plt <- ggplot(rats, aes(x = time, y = response,
                            group = subject, colour = treatment)) +
    geom_line() +
    facet_wrap(~ treatment, ncol = 3) +
    plt_labs


## ----plot-rat-data, echo = FALSE----------------------------------------------
rat_plt


## ----obs-per-rat, echo = FALSE------------------------------------------------
rats %>%
    na.omit() %>%
    count(subject) %>%
    count(n, name = "n_rats")


## ---- fig.align = "center", out.width = "95%", echo = FALSE-------------------
knitr::include_graphics("resources/lawton-et-al-hgam-locust-paper-fig.svg")


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
library('sf')
coast <- read_sf(here("data", "nl_coast.shp"))
ggplot(shrimp) +
  geom_point(aes(x = long, y = lat, size = shrimp), alpha = 0.5) +
  geom_sf(data = coast) +
  facet_wrap(~year, ncol = 5)


## ----fit-shrimp-space-time----------------------------------------------------
m_spt <- gam(shrimp ~ te(x, y, year, d = c(2,1),
                         bs = c('tp', 'cr'), k = c(20, 5)),
             data = shrimp,
             family = tw(),
             method = "REML")


## ----predict-newdata----------------------------------------------------------
new_year <- with(shrimp, tibble(year = seq(min(year), max(year), length.out = 100)))
pred <- predict(m_rich, newdata = new_year, se.fit = TRUE, type = 'link')
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
sp_new <- with(shrimp, expand.grid(x = evenly(x, n = 100), y = evenly(y, n = 100),
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
ti_new <- with(shrimp, expand.grid(x = mean(x), y = mean(y), year = evenly(year, n = 100)))

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
knitr::include_graphics("resources/miller-bayesian-gam-interpretation-fig.svg")


## ----richness-coefs-----------------------------------------------------------
sm_year <- get_smooth(m_rich, "s(year)") # extract the smooth object from model
idx <- gratia:::smooth_coef_indices(sm_year) # indices of the coefs for this smooth
idx

beta <- coef(m_rich)                     # vector of model parameters
beta[idx]                                # coefs for this smooth


## ----richness-vcov, results = "hide", dependson=-1----------------------------
Vb <- vcov(m_rich) # default is the bayesian covariance matrix
Vb


## ----richness-vcov-print, echo = FALSE, dependson=-1--------------------------
op <- options(width = 170)
Vb
options(op)


## ----richness-xp-matrix, dependson=-1-----------------------------------------
new_year <- with(shrimp, tibble(year = evenly(year, n = 100)))
Xp <- predict(m_rich, newdata = new_year, type = 'lpmatrix')
dim(Xp)


## ----richness-reduce-xp, dependson=-1-----------------------------------------
Xp <- Xp[, idx, drop = FALSE]
dim(Xp)


## ----richness-simulate-params, dependson=-1-----------------------------------
set.seed(42)
beta_sim <- rmvn(n = 20, beta[idx], Vb[idx, idx, drop = FALSE])
dim(beta_sim)


## ----richness-posterior-draws, fig.height = 5, fig.show = 'hide', dependson=-1----
sm_draws <- Xp %*% t(beta_sim)
dim(sm_draws)
matplot(sm_draws, type = 'l')


## ----richness-posterior-draws, fig.height = 5, fig.width = 5, echo = FALSE, results = 'hide'----
sm_draws <- Xp %*% t(beta_sim)
dim(sm_draws)
matplot(sm_draws, type = 'l')


## ----plot-posterior-smooths, fig.height = 5, dependson=-1---------------------
sm_post <- smooth_samples(m_rich, 's(year)', n = 20, seed = 42)
draw(sm_post)


## ----posterior-sim-model, dependson=-1----------------------------------------
beta <- coef(m_rich)   # vector of model parameters
Vb <- vcov(m_rich)     # default is the bayesian covariance matrix
Xp <- predict(m_rich, type = 'lpmatrix')
set.seed(42)
beta_sim <- rmvn(n = 1000, beta, Vb) # simulate parameters
eta_p <- Xp %*% t(beta_sim)        # form linear predictor values
mu_p <- inv_link(m_rich)(eta_p)    # apply inverse link function

mean(mu_p[1, ]) # mean of posterior for the first observation in the data
quantile(mu_p[1, ], probs = c(0.025, 0.975))


## ----posterior-sim-model-hist, fig.height = 5, dependson=-1-------------------
ggplot(tibble(richness = mu_p[587, ]), aes(x = richness)) +
    geom_histogram() + labs(title = "Posterior richness for obs #587")


## ----richness-fitted-samples, fig.height = 4.5, dependson=-1------------------
rich_post <- fitted_samples(m_rich, n = 1000, newdata = shrimp, seed = 42)
ggplot(filter(rich_post, row == 587), aes(x = fitted)) +
    geom_histogram() + labs(title = "Posterior richness for obs #587", x = "Richness")


## ----total-biomass-posterior-1, dependson=-1----------------------------------
sp_new <- with(shrimp, expand.grid(x = evenly(x, n = 100), y = evenly(y, n = 100),
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


## ----total-biomass-histogram, echo = FALSE, dependson=-1----------------------
ggplot(tibble(biomass = total_biomass), aes(x = biomass)) +
    geom_histogram()


## ----biomass-fitted-samples-example, dependson=-1-----------------------------
bio_post <- fitted_samples(m_spt, n = 1000,
                           newdata = sp_new[!too_far, ],
                           seed = 42) %>%
    group_by(draw) %>%
    summarise(total = sum(fitted),
              .groups = "drop_last")

with(bio_post, mean(total))
with(bio_post, quantile(total, probs = c(0.025, 0.975)))


## ----biomass-fitted-samples-plot, fig.width = 5, fig.height = 5, dependson=-1----
ggplot(bio_post, aes(x = total)) +
    geom_histogram() +
    labs(x = "Total biomass")


## -----------------------------------------------------------------------------
spp_url <- "https://bit.ly/spp-gradient"
gradient <- read_csv(spp_url, col_types = "dd")
gradient

