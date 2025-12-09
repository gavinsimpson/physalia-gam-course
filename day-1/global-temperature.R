# The HadCRUTv4 Global Temperature record example from the slides

# packages
library("readr")
library("ggplot2")
library("dplyr")
library("mgcv")
library("gratia")
library("marginaleffects")

# load the data from the web
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = "nnnnnnnnnnnn", col_names = FALSE) |>
  select(num_range("X", 1:2)) |>
  setNames(nm = c("Year", "Temperature"))

# plot
plt_labs <- labs(x = "Year", y = expression(Temperature ~ degree * C))
gtemp |>
  ggplot(
    aes(x = Year, y = Temperature)
  ) +
  geom_line() +
  geom_point() +
  plt_labs

# fit the Gaussian GAM
m_gtemp <- gam(
  Temperature ~ s(Year),
  data = gtemp,
  method = "REML",
  family = gaussian()
)

# model summary
summary(m_gtemp)
# or overview
m_gtemp |> overview()

# plot the estimate smooth
draw(m_gtemp, residuals = TRUE, rug = FALSE)

# plot the fitted smooth on the original data
#
# generate new data to predict at
newd <- m_gtemp |>
  data_slice(
    Year = evenly(Year, n = 200)
  )

# use fitted_values to get predictions on the response scale
fv_gtemp <- m_gtemp |>
  fitted_values(
    data = newd, scale = "response"
  )

# plot
fv_gtemp |>
ggplot(aes(x = Year, y = .fitted)) +
  geom_point(data = gtemp, aes(x = Year, y = Temperature)) +
  geom_ribbon(
    aes(
      ymin = .lower_ci, ymax = .upper_ci, x = Year
    ),
    alpha = 0.4,
    inherit.aes = FALSE,
    fill = "#fdb338"
  ) +
  geom_line(linewidth = 1, colour = "#025196") +
  plt_labs

# or with conditional_values
m_gtemp |>
  conditional_values(
    condition = "Year"
  ) |>
  draw() +
  geom_point(
    data = gtemp,
    aes(x = Year, y = Temperature)
  ) +
  plt_labs

# or with plot_predictions
m_gtemp |>
  plot_predictions(
    condition = "Year",
    points = 1
  ) +
  plt_labs
