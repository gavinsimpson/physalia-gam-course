# packages
pkgs <- c("mgcv", "lme4", "ggplot2", "readr", "dplyr", "tidyr",
          "gratia", "purrr", "patchwork", "hexbin")

# load packages
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

ocean <- read_csv(here("data", "chlorophyll.csv"), col_types = "dddddd") %>%
    rename(longitude = lon, latitude = lat, julian_day = jul.day,
    ocean_depth = bath, chl_sw = chl.sw)
# Variables are

# * `lon` longitude
# * `lat` latitude
# * `julian_day` day of year of observation
# * `ocean_depth` depth of the ocean in m
# * `chl` measured chlorophyll from bottle samples
# * `chl_sw` chlorophyll as estimated using the SeaWifs satellite

plts <- ocean %>%
    select(latitude, longitude, julian_day, ocean_depth) %>%
    imap(~ ggplot(ocean, aes(x = ..1, y = chl)) +
        geom_hex() +
        scale_fill_viridis_c(option = "plasma") +
        labs(x = ..2, y = "Chlorophyll"))

wrap_plots(plts, ncol = 3)

ocean %>%
    ggplot(aes(x = chl)) +
    geom_density()

# I don't believe that, the data are likely censored, but not much we can do
# in mgcv for that
ocean %>%
    summarise(min_chl = min(chl))

ocean %>%
    filter(chl > 0) %>%
    summarise(min_chl = min(chl))

# fit a simple model for space only
m <- gam(chl ~ s(latitude, longitude, bs = "sos"),
    family = tw(),
    data = chl, method = "REML")

# model summary
summary(m)

# plot the smooth
draw(m, n = 50)

# predict for new data
new_df <- data_slice(m, longitude = evenly(longitude, n = 50),
    latitude = evenly(latitude, n = 50))
ind <- too_far(new_df$longitude, new_df$latitude,
    ocean$longitude, ocean$latitude, dist = 0.1)

fv <- fitted_values(m, data = new_df) %>%
    mutate(fitted = if_else(ind, NA_real_, fitted))

# plot on the response scale
ggplot(chl, aes(x = longitude, y = latitude)) +
    geom_tile(data = fv, aes(fill = fitted)) +
    geom_point() +
    scale_fill_viridis_c(option = "plasma") +
    coord_map("orthographic", orientation = c(20, 0, 0)) +
    labs(fill = "Chlorophyll") +
    theme(legend.position = "top")
