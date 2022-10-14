# Ruddy turnstone example from Zuur et al

# packages
pkgs <- c("mgcv", "lme4", "ggplot2", "readr", "dplyr", "forcats", "tidyr",
          "gratia", "purrr", "patchwork")

# load packages
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load data
ruddy_url <- "https://bit.ly/gam-ruddy"
ruddy <- read_table(ruddy_url)

ruddy <- read_table(here("data", "ruddy.txt"), col_types = "ddddcdddddd") %>%
    mutate(TideState = factor(TideState, levels = c("Low", "High")),
        FlockID = factor(FlockID))

# Variables are

# * `BirdID` subject identifier
# * `FlockID` group identifier
# * `HeadUps` vigilence counts per 30 second period of focal bird (response)
# * `FlockSize` number of birds in each flock
# * `TideState` categorical, `High` or `Low` tide
# * `NumPecks` number of attempts by focal bird to capture a single prey item
# * `TimeOfDay` proportion of daylight period elapsed at time of observation
# * `TimeHighTide` time to the nearest high tide
# * `DaysStartWinter` number of days since November 1,
# * `Temperature` temperature in degrees celcius

plts <- ruddy %>%
    select(Temperature, DaysStartWinter, TimeHighTide, TimeOfDay, FlockSize,
        NumPecks) %>%
    imap(~ ggplot(ruddy, aes(x = ..1, y = HeadUps)) +
        geom_point() +
        geom_smooth(method = "loess", formula = y ~ x) +
        labs(x = ..2))

wrap_plots(plts, ncol = 3)

# model
m1 <- gam(HeadUps ~ TideState +
    FlockSize +
    s(NumPecks) +
    s(TimeOfDay) +
    s(TimeHighTide) +
    s(Temperature) +
    s(FlockID, bs = "re"),
data = ruddy,
family = poisson(link = "log"),
method = "REML")

# summary
summary(m1)

# diagnostics
appraise(m1, method = "simulate")

rg_m1 <- rootogram(m1)
draw(rg_m1)
# some underestimations of the zeros

# plot the fitted functions
draw(m1, parametric = TRUE)

# negbin
m2 <- update(m1, . ~ ., family = nb())

# summary
summary(m2)

# diagnostics
appraise(m2, method = "simulate")

rg_m2 <- rootogram(m2)
draw(rg_m2)
# some underestimations of the zeros still

# model
m3 <- update(m1, . ~ ., family = tw())
appraise(m3)

AIC(m1, m2, m3)

# anyway, model selection
summary(m1)

# we have insignificant terms; what should we do? I would argue that from the
# point of view of inference, we should just leave the model as is. You (we)
# thought all these covariates had effects on the response, so we should control
# for them even if they turn out to be not statistically significant at some
# arbitray threshold the community has decided to use. Importantly, note that
# excluding the insignificant terms would be tantamount to saying that the
# effect of these variables is === 0! Clearly the estimates *data* disagree with
# that!

# *if* you are interested in prediction, you could prune these out and refit,
# continuing to minimise some estimator of out-of-sample prediction error. You
# could compare models using Cross Validation or applied to a test set for
# example. Or you could use AIC.

# *if* you really aren't sure what should be in a model, then another
# alternative, but again, this comes from more of a prediction point of view
# rather than an inference or estimation point of view, you could use
# `select = TRUE`.

# of these options for predictions, I would probably use the shrinkage
# penalties (`select = TRUE`). That said, for prediction, and if you have a
# large enough data set, I would do repeated split smapling of the data into
# test and train sets, fitting with `select = TRUE` and compare which terms
# are in the models for each set of training data.

# Alternatively, I would only consider fitting separate models if there was
# strong a priori evidence or discussion suggestions two different hypotheses
# as to controls on the response and I wanted to test those hypotheses. In
# which case I would the models corresponding to the known hyoptheses and
# compare their properties in terms of ability to replicate the data, predict
# new observations, etc.

# fit with extra shrinkage penalties
m1a <- update(m1, . ~ ., select = TRUE)

summary(m1a)
appraise(m1a, method = "simulate")
draw(m1a)

# compare the smooths of temperature
c_temp <- compare_smooths(m1, m1a, smooths = "s(Temperature)")
draw(c_temp)

# As a scientist, would you treat these two estimated effects as different even
# though one of them was significan and one not at the threshold of 0.05?