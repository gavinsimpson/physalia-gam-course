## Fit distributional GAMs with mgcv

# packages
pkgs <- c("mgcv", "ggplot2", "readr", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# Analyse the motor cycle experiment data with a QGAM
# Load the mcycle data
data(mcycle, package = "MASS")

# Fit a GAM with linear predictors for the mean and the variance of accel
m_dist1 <- gam(list(accel ~ s(times, k = 20, bs = "ad"),
                          ~ s(times, k = 10)),
               data = mcycle,
               method = "REML", # <== IIRC REML is only option for these LSS
               family = gaulss())

## we can use draw() and appraise() as usual
## don't want the uncertainty in the intercept here
draw(m_dist1, overall_uncertainty = FALSE)

## model diagnostics
appraise(m_dist1)

## model summary
summary(m_dist1)

## looks OK, what about a standard Gaussian? This is a gaulss() with a constant
## in the variance
# Fit a GAM with linear predictors for the mean and the variance of accel
m_dist2 <- gam(list(accel ~ s(times, k = 20, bs = "ad"),
                          ~ 1),
               data = mcycle,
               method = "REML", # <== IIRC REML is only option for these LSS
               family = gaulss())

## Do we need a smooth in the variance or is linear OK?
m_dist3 <- gam(list(accel ~ s(times, k = 20, bs = "ad"),
                          ~ times + s(times, m = c(2,0))),
               data = mcycle,
               method = "REML", # <== IIRC REML is only option for these LSS
               family = gaulss())

## model summary
summary(m_dist3)

## AIC as another model guide
AIC(m_dist1, m_dist2, m_dist3)

## predict to visualise the estimated quantile
## new data to predict at
new_df <- with(mcycle,
               tibble(times = seq(round(min(times)), round(max(times)),
                                  length.out = 200)))

## sorry, doesn't work yet
fitted_values(m_dist1)

## so we go back to our recipe
pred <- predict(m_dist1, newdata = new_df, se.fit = TRUE, type = "link")
head(pred$fit) ## what ???

## we can handle this, get the inverse link functions
mu_ilink <- inv_link(m_dist1, parameter = "mu")
scale_ilink <- inv_link(m_dist1, parameter = "scale")

## now follow the recipe, slightly modified
pred <- pred %>%
  as.data.frame() %>%
  as_tibble() %>%
  setNames(c("fit_mu", "fit_scale", "se_mu", "se_scale")) %>%
  bind_cols(new_df)

pred

## now finish off by doing the mapping from the link scale to the resposne scale
crit <- gratia:::coverage_normal(0.95)
pred <- pred %>%
  mutate(fitted_mu = mu_ilink(fit_mu),
         fitted_sd = scale_ilink(fit_scale),
         lower_mu  = mu_ilink(fit_mu - (crit * se_mu)),
         upper_mu  = mu_ilink(fit_mu + (crit * se_mu)),
         lower_sd  = scale_ilink(fit_scale - (crit * se_scale)),
         upper_sd  = scale_ilink(fit_scale + (crit * se_scale))) %>%
  mutate(across(matches("_sd$"), .fns = ~ 1 / .x))

## plot
mu_plt <- pred %>%
  ggplot(aes(x = times)) +
    geom_point(data = mcycle, aes(x = times, y = accel)) +
    geom_ribbon(aes(ymin = lower_mu, ymax = upper_mu), alpha = 0.2) +
    geom_line(aes(y = fitted_mu)) +
    labs(y = "Acceleration", x = "Milliseconds after impact")

## have something to plot for "data" for the std dev plot
## take the absolute value of the response residual
res_data <- mcycle %>%
  mutate(abs_residual = abs(accel - fitted(m_dist1)[,1]))

sd_plt <- pred %>%
  ggplot(aes(x = times)) +
    geom_point(data = res_data, aes(x = times, y = abs_residual)) +
    geom_ribbon(aes(ymin = lower_sd, ymax = upper_sd), alpha = 0.2) +
    geom_line(aes(y = fitted_sd)) +
    labs(y = "Std. Deviation", x = "Milliseconds after impact")

mu_plt + sd_plt
