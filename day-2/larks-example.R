# Portugese Larks example from slides

# packages
library("readr")
library("dplyr")
library("mgcv")
library("gratia")
library("patchwork")

# read in the data
larks_url <- "https://bit.ly/gam-larks"
larks <-  read_csv(larks_url, col_types = "ccdddd")

# Or as per the slides:
# larks <-  read_csv(here("data", "larks.csv"),
#                   col_types = "ccdddd")

# convert some variables to factors (two bird species), and scale the
# easting and northing variables to be in KM
# numerically we prefer smaller valued variables
larks <- larks %>%
  mutate(crestlark = factor(crestlark),
         linnet = factor(linnet),
         e = x / 1000,
         n = y / 1000)
larks

# plot the data
larks %>%
    ggplot(aes(x = e, y = n, colour = crestlark)) +
    geom_point(size = 0.5) + # add a point layer
    coord_fixed() +
    scale_colour_discrete(na.value = "#bbbbbb33") +
    labs(x = NULL, y = NULL)

# fit our model
# isotropic thin plate regression spline smooth of spatial coordinates
crest <- gam(crestlark ~ s(e, n, k = 100),
             data = larks,
             family = binomial,
             method = "REML")

# model summary
summary(crest)

# visualise the estimate smooth
# a grid of 75x75 is a bit more efficient than the default of 100x100
draw(crest, rug = FALSE, n = 75)

# model checking is a pain for 0/1 data, so aggregate to binomial counts
larks2 <- larks %>%
  mutate(crestlark = as.numeric(as.character(crestlark)),  # to numeric
         linnet = as.numeric(as.character(linnet)),
         tet_n = rep(1, nrow(larks)), # counter for how many grid cells we sum
         N = rep(1, nrow(larks)),     # number of obs, 1 per row currently
         N = if_else(is.na(crestlark), NA_real_, N)) %>% # set N to NA if no obs
  group_by(QUADRICULA) %>%                     # group by the larger grid square
  summarise(across(c(N, crestlark, linnet, tet_n, e, n),
                   ~ sum(., na.rm = TRUE))) %>%  # sum all needed variables
  mutate(e = e / tet_n, n = n / tet_n) # rescale to get avg E,N coords

# fit binomial GAM
crest2 <- gam(cbind(crestlark, N - crestlark) ~ s(e, n, k = 100),
    data = larks2,
    family = binomial,
    method = "REML")

# model summary
summary(crest2)

# compare the fits
draw(crest, n = 75, rug = FALSE) +
    draw(crest2, n = 75, rug = FALSE) +
    plot_layout(ncol = 2)
# look pretty similar

# let's look at some model diagnostics
appraise(crest2, method = "simulate")

# hmmm, some overdispersion - refit as a quasi binomial
crest3 <- gam(cbind(crestlark, N - crestlark) ~ s(e, n, k = 100),
    data = larks2,
    family = quasibinomial,
    method = "REML")

# look at the estimate dispersion parameter - this should be 1 for true
# binomial counts
crest3$scale # should be == 1

# model summary
summary(crest3)

# let's look at some model diagnostics
appraise(crest3, method = "simulate")

# compare the fits
draw(crest, n = 75, rug = FALSE) +
    draw(crest3, n = 75, rug = FALSE) +
    plot_layout(ncol = 2)

# generate fitted values
ds <- data_slice(crest3, e = evenly(e, n = 75), n = evenly(n, n = 75))
fv <- fitted_values(crest3, data = ds, scale = "response")

fv %>%
ggplot(aes(x = e, y = n, fill = fitted)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal()

# Repeat the modelling using the linnet data