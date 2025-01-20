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

# Can also
set.seed(6436)
cal <- tuneLearn(accel ~ s(times, k = 20, bs = "ad"),
  data = mcycle,
  qu = 0.8,
  lsig = seq(1, 3, length.out = 20), #<- sequence of values for learning rate
  control = list("progress" = "none")
)

layout(matrix(1:2, ncol = 2))
check(cal)
layout(1)

# summary as usual
summary(m_q)

## we can use draw() and appraise() as for a single quantile the fits are GAM
## objects with some extra things
draw(m_q)

## model diagnostics
appraise(m_q, use_worm = TRUE)

## predict to visualise the estimated quantile
## new data to predict at
new_df <- data_slice(mcycle, times = evenly(times, n = 200))

## for a single QGAM the fits work just fine with fitted_values too
fv <- fitted_values(m_q, data = new_df) # I broke this - FIXME!

## plot
fv |>
  ggplot() +
    geom_point(data = mcycle, aes(x = times, y = accel)) +
    geom_ribbon(aes(x = times, ymin = .lower_ci, ymax = .upper_ci),
    alpha = 0.2) +
    geom_line(aes(y = .fitted, x = times))

## Repeat the exercise but for a lower quantile, say the 0.2 probability
##   quantile

m_mq <- mqgam(accel ~ s(times, k = 20, bs = "ad"),
              data = mcycle,
              qu = c(0.2, 0.8))

invisible( qdo(m_mq, 0.2, plot, pages = 1) )

p_80 <- qdo(m_mq, 0.8, draw) +
  labs(subtitle = "0.8 quantile")
p_20 <- qdo(m_mq, 0.2, draw) +
  labs(subtitle = "0.2 quantile")

p_20 + p_80 + plot_layout(ncol = 2, nrow = 1)

sm_80 <- qdo(m_mq, 0.8, smooth_estimates)
sm_20 <- qdo(m_mq, 0.2, smooth_estimates)

sm_q <- sm_20 |>
  bind_rows(sm_80) |>
  mutate(.q = rep(c(0.2, 0.8), each = 100))

sm_q |>
  ggplot(aes(x = times, y = .estimate, colour = .q, group = .q)) +
  geom_line()
