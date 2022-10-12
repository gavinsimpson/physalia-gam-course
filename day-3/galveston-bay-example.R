# Galveston Bay example

# packages
pkgs <- c("ggplot2", "here", "readr", "mgcv", "gratia", "dplyr", "patchwork")
vapply(pkgs, library, logical(1L), logical.return = TRUE,
       character.only = TRUE)

# read in the data
galveston <- read_csv(here("data", "galveston.csv")) %>%
    mutate(datetime = as.POSIXct(paste(DATE, TIME),
                                 format = "%m/%d/%y %H:%M", tz = "CDT"),
           STATION_ID = factor(STATION_ID),
           DoY = as.numeric(format(datetime, format = "%j")),
           ToD = as.numeric(format(datetime, format = "%H")) +
               (as.numeric(format(datetime, format = "%M")) / 60))
galveston

# set boundary knots for the DoY smooth
knots <- list(DoY = c(0.5, 366.5))

# fit using BAM
m <- bam(MEASUREMENT ~
             s(ToD, k = 10) +
             s(DoY, k = 12, bs = "cc") +
             s(YEAR, k = 30) +
             s(LONGITUDE, LATITUDE, k = 100, bs = "ds", m = c(1, 0.5)) +
             ti(DoY, YEAR, bs = c("cc", "tp"), k = c(12, 15)) +
             ti(LONGITUDE, LATITUDE, ToD, d = c(2,1), bs = c("ds","tp"),
                m = list(c(1, 0.5), NA), k = c(20, 10)) +
             ti(LONGITUDE, LATITUDE, DoY, d = c(2,1), bs = c("ds","cc"),
                m = list(c(1, 0.5), NA), k = c(25, 12)) +
             ti(LONGITUDE, LATITUDE, YEAR, d = c(2,1), bs = c("ds","tp"),
                m = list(c(1, 0.5), NA), k = c(25, 15)),
         data = galveston,
         method = "fREML",
         knots = knots,
         nthreads = c(4, 4),
         discrete = TRUE)

## fit the simpler model
m_sub <- bam(MEASUREMENT ~
             s(ToD, k = 10) +
             s(DoY, k = 12, bs = "cc") +
             s(YEAR, k = 30) +
             s(LONGITUDE, LATITUDE, k = 100, bs = "ds", m = c(1, 0.5)) +
             ti(DoY, YEAR, bs = c("cc", "tp"), k = c(12, 15)),
         data = galveston,
         method = "fREML",
         knots = knots,
         nthreads = c(4, 4),
         discrete = TRUE)

# comapre fits
AIC(m, m_sub)

# Look at the complex model
summary(m)

# plot.gam output
plot(m, pages = 1, scheme = 2, shade = TRUE)

# draw() output
draw(m, scales = "free")

# spatial predictions over time
pdata <- data_slice(m, ToD = 12, DoY = 180,
                    YEAR = evenly(YEAR, by = 1),
                    LONGITUDE = evenly(LONGITUDE, n = 50),
                    LATITUDE  = evenly(LATITUDE, n = 50))
fv <- fitted_values(m, data = pdata)
# set fitted values to NA for grid points that are too far from the data
ind <- too_far(pdata$LONGITUDE, pdata$LATITUDE,
               galveston$LONGITUDE, galveston$LATITUDE, dist = 0.1)
fv <- fv %>%
  mutate(fitted = if_else(ind, NA_real_, fitted))

# plot the estimated spatial field over time
ggplot(fv, aes(x = LONGITUDE, y = LATITUDE)) +
  geom_raster(aes(fill = fitted)) + facet_wrap(~ YEAR, ncol = 12) +
  scale_fill_viridis(name = expression(degree*C), option = "plasma",
                     na.value = "transparent") +
  coord_quickmap() +
  scale_x_continuous(guide = guide_axis(n.dodge = 2,
                                        check.overlap = TRUE)) +
  theme(legend.position = "top")

# predict temperature time series at a particular location for 4 different
# periods of the year
ds <- data_slice(m, ToD = 12, DoY = c(1, 90, 180, 270),
                 YEAR = evenly(YEAR, n = 250),
                 LONGITUDE = -94.8751, LATITUDE  = 29.50866)
fv2 <- fitted_values(m, data = ds, scale = "response")

ggplot(fv2, aes(x = YEAR, y = fitted, group = factor(DoY))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
  geom_line() + facet_wrap(~ DoY, scales = "free_y") +
  labs(x = NULL, y = expression(Temperature ~ (degree * C)))
