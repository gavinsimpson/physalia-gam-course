# South Pole CO2 example from slides

# packages
library("readr")
library("dplyr")
library("mgcv")
library("gratia")

south_url <- "https://bit.ly/gam-south-pole-co2"
south <- read_csv(south_url, col_types = "ddd")
south

# co2 is the CO2 concentration in PPM
# c.month is a monthly counter 1, 2, M, where M is the last n_year * 12 month
#   alternatively just think of this as a timestep variable
# month is a numeric indicator of month of year

# plot the data
ggplot(south, aes(x = c.month, y = co2)) +
  geom_line()

# Fit a naive GAM using the time step variable
# we need really large k to capture the within and between year variability
# in a single smooth
m_co2 <- gam(co2 ~ s(c.month, k = 300, bs = "cr"),
  data = south, method = "REML")

# model summary
summary(m_co2)

# let's predict from this naive model
# this requires a new data frame of values we want to predict at
# lets also forecast the next few years (36 months)
new_df <- with(south, tibble(c.month = 1:(nrow(south) + 36)))

# look at our data
new_df

# then we use the fitted_values function from gratia to generate predictions
fv <- fitted_values(m_co2, data = new_df, scale = "response")

# look at fv
fv

# plot our predictions
fv |>
  ggplot(aes(x = c.month, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line(data = south, aes(c.month, co2), col = "red") +
  geom_line(alpha = 0.4)

# linear extrapolation - be very careful extrapolating with splines...

# a better model would decompose into seasonal variation and long term trend
# ideally we'd do this with `year` and `month` variables but here we just have
# a time variable (`c.month` and a `month`` variable)
# Note how we set the knots
m2_co2 <- gam(co2 ~ s(month, bs = "cc") +
                s(c.month, bs = "cr", k = 300),
  data = south, method = "REML",
  knots = list(month = c(0.5, 12.5)))

# model summary
summary(m2_co2)

# plot these two smooths
draw(m2_co2, residuals = TRUE, rug = FALSE)

# compare their complexities
model_edf(m_co2, m2_co2)

# what about our forecasts
nr <- nrow(south)
new_df <- with(south,
               tibble(c.month = 1:(nr + 36),
                      month = rep(seq_len(12), length = nr + 36)))
fv2 <- fitted_values(m2_co2, data = new_df, scale = "response")

# and plot
fv2 |>
ggplot(aes(x = c.month, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line(data = south, aes(c.month, co2), col = 'red') +
  geom_line(alpha = 0.4)
