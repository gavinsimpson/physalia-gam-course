# Analyse the motor cycle experiment data with a QGAM
pkgs <- c("mgcv", "qgam", "ggplot2", "readr", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

## Load the mcycle data
data(mcycle, package = "MASS")

## Fit the standard GAM but to the 0.8 probability quantile of Y
## using an adaptive smoother
m_q <- qgam(accel ~ s(times, k = 20, bs = "ad"),
            data = mcycle,
            qu = 0.8)

## check the learning rate
check(m_q$calibr, nbin = 2)
## if the optimisation worked, this should have a minimum

## we can use draw() and appraise() as for a single quantile the fits are GAM
## objects with some extra things
draw(m_q)

## model diagnostics
appraise(m_q)

## predict to visualise the estimated quantile
## new data to predict at
new_df <- with(mcycle,
               tibble(times = seq(round(min(times)), round(max(times)),
                                  length.out = 200)))

## for a single QGAM the fits work just fine with fitted_values too
fv <- fitted_values(m_q, data = new_df)

## plot
fv %>%
  ggplot(aes(x = times)) +
    geom_point(data = mcycle, aes(x = times, y = accel)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(aes(y = fitted))

## Repeat the exercise but for a lower quantile, say the 0.2 probability
##   quantile
