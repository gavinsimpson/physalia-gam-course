# packages
pkgs <- c("here", "sf", "mgcv", "lme4", "ggplot2", "readr", "dplyr", "tidyr",
  "gratia", "patchwork")

# load packages
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
  quietly = TRUE)

# Load the shrimp data
shrimp <- read_csv("https://bit.ly/nl-shrimp",
  col_types = "dcdddddddddddd") %>%
  mutate(stratum = factor(stratum)) # to allow a random effect of stratum

# Plot the richness data
ggplot(shrimp) +
  geom_violin(aes(x = richness, y = factor(year))) +
  labs(x = "Number of species", y = "Year")

# A simple model for shrimp species richness
# completely ignoring the repeated measures on trawl locations / space
m_rich <- gam(richness ~ s(year),
  family = poisson,
  method = "REML",
  data = shrimp)

# plot the biomass data spatially - can't load this from the internet
# as a shapefile is actually many files
coast <- read_sf(here("data", "nl_coast.shp"))
ggplot(shrimp) +
  geom_point(aes(x = long, y = lat, size = shrimp), alpha = 0.5) +
  geom_sf(data = coast) +
  facet_wrap(~year, ncol = 5)

# A spatio temporal model for shrimp biomass
m_spt <- gam(shrimp ~ te(x, y, year, d = c(2, 1),
    bs = c("tp", "cr"), k = c(20, 5)),
  data = shrimp,
  family = tw(),
  method = "REML")

# data to predict at for thw two models
new_year <- with(shrimp, tibble(year = seq(min(year), max(year),
  length.out = 100)))

sp_new <- with(shrimp,
  expand.grid(x = evenly(x, n = 100), y = evenly(y, n = 100),
    year = unique(year))
)

# or with data_slice
ds_year <- data_slice(m_rich, year = evenly(year, n = 100))
ds_spt <- data_slice(m_spt,
  x = evenly(x, n = 100),
  y = evenly(y, n = 100),
  year = evenly(year, by = 1)
)
spt_too_far <- too_far(ds_spt$x, ds_spt$y, shrimp$x, shrimp$y, dist = 0.05)

# fitted values for shrimp richness
fv_rich <- fitted_values(m_rich, data = ds_year)

# plot fitted values of richness
fv_rich |>
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line(aes(y = .fitted)) +
  labs(y = "Species richness", x = NULL) +
  scale_x_continuous(breaks = 2005:2014)

# fitted values for biomass
fv_spt <- fitted_values(m_spt, data = ds_spt) |>
  mutate(.fitted = case_when(spt_too_far ~ NA_real_, TRUE ~ .fitted))

# plot the fitted values for biomass
fv_spt |>
  ggplot(aes(x = x, y = y, fill = .fitted)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma") +
  facet_wrap(~year, ncol = 5) +
  coord_equal() +
  labs(fill = "Biomass")

# decomposed version of the spatio temporal model for biomass
m_ti <- gam(shrimp ~
  ti(x, y, year, d = c(2, 1), bs = c("tp", "cr"), k = c(20, 5)) +
  s(x, y, bs = "tp", k = 20) +
  s(year, bs = "cr", k = 5),
data = shrimp, family = tw, method = "REML")

# We will now predict for the average f(year) trend

# 1. identify the names of the smooths
smooths(m_ti)

# 2. data to predict at
ds_ti <- data_slice(m_ti,
  x = mean(x), y = mean(y), year = evenly(year, n = 100))

# 3. predict
fv_spt2 <- fitted_values(m_ti, data = ds_ti,
  scale = "response", exclude = c("ti(x,y,year)", "s(x,y)"))

# 4. plot
fv_spt2 |>
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3) +
  geom_line(aes(y = .fitted)) +
  labs(y = "Biomass", x = NULL) +
  scale_x_continuous(breaks = 2005:2014)

## Posterior simulation

# posterior smoothers for richness model
sm_post <- smooth_samples(m_rich, "s(year)", n = 20, seed = 42)

# plot the draws
sm_post |>
  draw(alpha = 0.5)

# instead, use the Metropolis Hasting sampler and form a Bayesian interval
# takes a ~10 seconds to sample
sm_mh <- smooth_samples(m_rich, "s(year)", data = ds_year, n = 10000, seed = 42,
  method = "mh", burnin = 1500, thin = 2)

# plot 20 of the draws
sm_mh |>
  draw(n_samples = 20, seed = 2, alpha = 0.5)

# create the quantile interval from the posterior

# 1. a function to compute quantiles
quantile_fun <- function(x, probs = c(0.025, 0.975), ...) {
  tibble::tibble(
    .value = quantile(x, probs = probs, ...),
    .q = probs * 100
  )
}

# 2. add the row number to the data we are sampling at
ds_year <- ds_year |>
  mutate(.row = row_number())

# 3. compute the quantiles of the posterior for all values of year
#    and rearrange the data into long format ready to plot
q_int <- sm_mh |>
  group_by(.row) |>
  reframe(quantile_fun(.value)) |>
  pivot_wider(
    id_cols = .row, names_from = .q, values_from = .value,
    names_prefix = ".q"
  ) |>
  left_join(ds_year, by = join_by(.row == .row))

# 4. plot
sm_mh |>
  draw(n_samples = 20, seed = 2, alpha = 0.5) +
  geom_line(data = q_int, aes(x = year, y = .q2.5), inherit.aes = FALSE,
    colour = "yellow", linewidth = 1) +
  geom_line(data = q_int, aes(x = year, y = .q97.5), inherit.aes = FALSE,
    colour = "yellow", linewidth = 1)

## Posterior simulation for Richness

# To explore the uncertainty in the fitted values due to model uncertainty we
# can use `fitted_samples()`
rich_post <- fitted_samples(m_rich, n = 1000, data = shrimp, seed = 42)

# plot posterior distribution for a selected sample
rich_post |>
  filter(.row == 587) |>
  ggplot(aes(x = .fitted)) +
  geom_histogram() +
  labs(title = "Posterior richness for obs #587", x = "Richness")

# if you want to simulate new data from the model you can use
# `predicted_samples()` - note this doesn't include any uncertainty in the model
rich_ppred <- predicted_samples(m_rich, n = 1000, data = shrimp, seed = 42)

# plot posterior predictions for a selected sample
rich_ppred |>
  filter(.row == 587) |>
  ggplot(aes(x = .response)) +
  geom_histogram() +
  labs(title = "Posterior richness for obs #587", x = "Richness")

# if you want to include both the uncertainty in the model estimates and the
# sampling variation, use `posterior_samples()`
rich_ppost <- posterior_samples(m_rich, n = 1000, data = shrimp, seed = 42)

# plot posterior sample for a selected sample
rich_ppost |>
  filter(.row == 587) |>
  ggplot(aes(x = .response)) +
  geom_histogram() +
  labs(title = "Posterior richness for obs #587", x = "Richness")

# Note that, for speed, I have only taken 1000 draws here. Typically we might
# use more draws when doing this for real, especially if creating a quantile
# based interval

# Posterior prediction for uncertainty in 2007 biomass
ds_2007 <- data_slice(m_spt, x = evenly(x, n = 100), y = evenly(y, n = 100),
  year = 2007)
ds_2007_far <- too_far(ds_2007$x, ds_2007$y, shrimp$x, shrimp$y, dist = 0.05)

bio_post <- fitted_samples(m_spt, n = 1000,
  data = ds_2007[!ds_2007_far, ],
  seed = 42) |>
  group_by(.draw) |>
  summarise(total = sum(.fitted),
    .groups = "drop_last")

# note I used a smaller `dist` in `too_far()` than in the slides
with(bio_post, mean(total))
with(bio_post, quantile(total, probs = c(0.025, 0.975)))

# If you want samples from the full posterior, we use posterior_samples

bio_post_full <- posterior_samples(m_spt, n = 1000,
  data = ds_2007[!ds_2007_far, ],
  seed = 42) |>
  group_by(.draw) |>
  summarise(total = sum(.response), # <-- changed
    .groups = "drop_last")

# note I used a smaller `dist` in `too_far()` than in the slides
with(bio_post_full, mean(total))
with(bio_post_full, quantile(total, probs = c(0.025, 0.975)))

# because we are using the Gaussian approximation, the mean in `bio_post` and
# the one in `bio_post_full` are very close (increasing `n` will make them
# closer still), but we have a larger 95% interval because we include the
# sampling variation too

# If you want to use the Metropolis Hasting sampler, just set the `method`
# this takes about 10 seconds for 1000 draws
bio_post_full_mh <- posterior_samples(
  m_spt, n = 1000,
  data = ds_2007[!ds_2007_far, ],
  seed = 42,
  method = "mh"                     # <-- changed
) |>
  group_by(.draw) |>
  summarise(total = sum(.response), # <-- changed
    .groups = "drop_last")

# note I used a smaller `dist` in `too_far()` than in the slides
with(bio_post_full_mh, mean(total))
with(bio_post_full_mh, quantile(total, probs = c(0.025, 0.5, 0.975)))

# plot the posteriors for comparison
bio_posteriors <- bio_post |>
  bind_rows(bio_post_full, bio_post_full_mh) |>
  mutate(type = rep(c("Expectations", "Full posterior", "Full posterior MH"),
  each = 1000))

bio_posteriors |>
  ggplot(aes(x = total)) +
  geom_histogram() +
  facet_wrap(~ type) +
  labs(x = "Total biomass")
