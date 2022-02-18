# Analyse the Global temperature anaomaly record with BRMS
pkgs <- c("mgcv", "brms", "ggplot2", "readr", "dplyr", "tidyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load temperature record
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn',
                    col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>%
    setNames(nm = c('Year', 'Temperature'))

# Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_line() + 
    geom_point() +
    labs(x = 'Year', y = expression(Temeprature ~ degree * C))
gtemp_plt

# fit the model, default everything; adjust cores to reflect the number of
#   phsyical cores on your machine
m <- brm(bf(Temperature ~ s(Year)), # <== formula wrapped in bf()
            data = gtemp,
            family = gaussian(),
            cores = 4, # <== adjust me
            seed = 17,
            control = list(adapt_delta = 0.95))

# Model summary
summary(m)

# The equivalent of draw.gam / plot.gam is
c_sm <- conditional_smooths(m)

# posterior predictive checks; can we reproduce the data with samples from
#   the posterior?
pp_check(m) # density plot

pp_check(m, type = "ecdf_overlay") # epirical cummulative density function

## What about priors?

# use `get_prior()` with the formula to get info on the defualt priors that will
#   be used by brms
get_prior(bf(Temperature ~ s(Year)),
          data = gtemp, family = gaussian())

## This includes flat priors on the linear basis `sYear_1` and something for
##  the Intercept.
##
## We can put a prior on the intercept and the linear term by setting a prior
##   on the `b` class of parameters (BRMS speak for common or garden fixed
##   effects *excluding the Intercept* [by default]). Here we set all fixed
##   effect terms to have student t priors, 3 df, scale of 2.5
pr <- prior(student_t(3, 0, 2.5), class = "b", coef = "sYear_1")

## To get BRMS to apply the above penalty to the intercept as well, we have to
##   change the formula to be `0 + Intercept` where 0 says drop the usual
##   intercept (which BRMS treats in a special way). Otherwise we get the
##   default prior on the intercept, which may be more efficient, but is flat
##   IIRC.
##
## Summary when using `prior = pr`
##
## * Temperature ~ s(Year)
##
##     - get t prior on sYear_1 but flat on intercept
##
## * Temperature ~ 0 + Intercept + s(Year)
##
##     - get t prior on sYear_1 and intercept

# Here we only sample from the prior
m_prior1 <- brm(bf(Temperature ~ s(Year)),
                prior = pr, # <== specify our prior
                data = gtemp,
                family = gaussian(),
                sample_prior = "only", # <== only prior samples
                cores = 4,
                seed = 17,
                control = list(adapt_delta = 0.95))

## If you want to look at draws from the prior for the smooth (think equivalent
##   of draw.gam but sampling coefficients for the smooths from the prior only,
##   no data fitting), then we can grab the prior draws for the smooth only with
##   (confusingly!!) `posterior_smooths()`:
p_year <- posterior_smooths(m_prior1, smooth = "s(Year)")

## then we need to massage into something we can plot as `p_year` is just a
##   matrix of draws, 1 row per draw, 1 column per data/observations
dim(p_year1)

## so massage away
pr_draws1 <- data.frame(t(p_year1)) %>%  # convert to a df transposed(!)
  setNames(seq_len(nrow(p_year1))) %>%   # need names for the draws 1:4000
  as_tibble() %>%                        # convert to a tibble
  add_column(Year = gtemp[["Year"]],     # add a column for the data: Year
             .before = 1L) %>%
  pivot_longer(-Year, names_to = "draw") # pivot from wide -> long
## look at what we did
pr_draws1

## Plot a sample of the draws
n_take <- 50                                       # how many draws to plot?
n <- nrow(gtemp)                                   # how many data?
set.seed(3)                                        # repeatible
pr_draws1 %>%
  filter(draw %in% sample(seq_len(n), n_take)) %>% # sample draws
  ggplot(aes(x = Year, y = value, group = draw)) + # plot
    geom_line(alpha = 0.5)                         # 1 line per draw

## can use the posterior predictive check plot pp_check() with draws from the
##   prior only; these then a *prior* predictive checks
pp_check(m_prior1)
## or
pp_check(m_prior1, type = "ecdf_overlay")

## Now lets change the priors for the variance of the random effect / wiggliness
##   of the smooth as well as the other terms
pr <- c(prior(student_t(10, 0, 2.5), class = "b"),
        prior(student_t(10, 0, 2), class = "b", coef = "sYear_1"),
        prior(student_t(10, 0, 2), class = "sds"))

# prior sampling, default everything
m_prior2 <- brm(bf(Temperature ~ 0 + Intercept + s(Year)), # <== alt version
                prior = pr, # <== set out new prior
                data = gtemp,
                family = gaussian(),
                sample_prior = "only", # <== sample from prior only
                cores = 4, # <== change me
                seed = 17,
                control = list(adapt_delta = 0.95))

## posterior (prior!!) smooths for this new prior
p_year2 <- posterior_smooths(m_prior2, smooth = "s(Year)")

pr_draws2 <- data.frame(t(p_year2)) %>%  # convert to a df transposed(!)
  setNames(seq_len(nrow(p_year2))) %>%   # need names for the draws 1:4000
  as_tibble() %>%                        # convert to a tibble
  add_column(Year = gtemp[["Year"]],     # add a column for the data: Year
             .before = 1L) %>%
  pivot_longer(-Year, names_to = "draw") # pivot from wide -> long
## check what we did
pr_draws2

# plot a sample of the draws
n_take <- 50                                       # how many draws to plot?
n <- nrow(gtemp)                                   # how many data?
set.seed(3)                                        # repeatible
pr_draws2 %>%
  filter(draw %in% sample(seq_len(n), n_take)) %>% # sample draws
  ggplot(aes(x = Year, y = value, group = draw)) + # plot
    geom_line(alpha = 0.5)                         # 1 line per draw

# conditional smooths is partial effect of smooth (but only the prior)
conditional_smooths(m_prior2)

## prior predictive checks
pp_check(m_prior2)
pp_check(m_prior2, type = "ecdf_overlay")

## plot observed vs simulated data from the prior
prior_pred <- predict(m_prior2) %>%
  as_tibble() %>%
  add_column(Temperature = gtemp$Temperature, .before = 1L)

# check what we did
prior_pred

# plot
prior_pred %>%
  ggplot(aes(x = Temperature, y = Estimate)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    coord_equal()

## Is this any good?

##Nope!!

## Note this is the mean (or the median) of the prior distribution of Y
## We could modify this to show the 95% quantiles of the draws, which would show
## the variation in the prior data samples
# plot
prior_pred %>%
  ggplot(aes(x = Temperature, y = Estimate)) +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5))

## Still not good

##   If you want draws from the prior distribution of Y, use
pr_draws_y <- predict(m_prior2, summary + FALSE)
## This is a matrix which you can wrangle as we did for pr_draws2 above if you
## want to extract and plot an actual draw from the prior against the observed
## Y.

## Task: alter the priors above to give prior draws that are more like the data
## 
## * copy the code for pr below
pr2 <-                                              # <== your prior code here

## modify the priors so that we get data that looks more like the data
## you could change the df and/or the scale of the student_t() priors, e.g.
##
##   student_t(3, 0, 2.5)
##
## will have the slightly larger scale (spread, SD) but fatter tails than the
## prior we used before. You can try normal/Gaussian priors too:
##
##   normal(0, 5)
##
## for exmaple would be a mean 0 Guassian prior with a SD of 5

## Change one or all of the priors - it seems reasonable to change them all the
## same way as the betas for the fully penalized smooth would be MVN(0, S^-)
##
## Then sample from your new prior
m_prior3 <- brm(bf(Temperature ~ 0 + Intercept + s(Year)), # <== alt version
                prior = pr2, # <== set our new prior
                data = gtemp,
                family = gaussian(),
                sample_prior = "only", # <== sample from prior only
                cores = 4, # <== change me
                seed = 17,
                control = list(adapt_delta = 0.95))

## get the prior draws for the smooth and plot

## use predict() to get the prior predictive sample (i.e. means/median of prior)

## plot the prior predictive values for Y against observed Temperature

## Are you getting better prior draws? Do they look more like the data?

##----------------------------------------------------------------------------##

## Alternative Task:
##
## If you'd rather just try out some more model fits, take one of the models
## we've worked on (best to use a simpler model from earlier in the week) and
## convert the mgcv model to a brms one and fit it.

## As an example
data(mcycle, package = "MASS")
m <- gam(accel ~ s(times), data = mcycle, method = "REML")

## convert m to be a BRMS model
m_brms <- brm(
    bf(_____ ~ _____),                          # <== change me
    data = -----,                               # <== change me
    family = gaussian(),
    cores = 4, # <== change me
    seed = 17,
    control = list(adapt_delta = 0.8)           # <== increase if divergences
  )

## as an extra challenge, figure out how to fit a distributional model with brms
## which would be the equivalent of:
m_distr <- gam(list(
                    accel ~ s(times), # <== eta for mean/mu
                    ~ s(times)        # eta for the "variance"
               ),
               data = mcycle,
               family = gaulss(),
               method = "REML")

## Read ?bf for help or take a look at:
## https://paul-buerkner.github.io/brms/articles/brms_distreg.html
## for examples