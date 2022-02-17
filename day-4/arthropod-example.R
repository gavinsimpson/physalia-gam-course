# Reanalyse arthropod data from Seibold et al (2019) using some
# ideas from Daskalova, Phillimore & Meyers-Smith (2021)
#
# Refs
# Daskalova, G.N., Phillimore, A.B., Myers‐Smith, I.H., 2021. Accounting for
#   year effects and sampling error in temporal analyses of invertebrate
#   population and biodiversity change: a comment on Seibold et al. 2019.
#   *Insect Conserv. Divers.* 14, 149–154. https://doi.org/10.1111/icad.12468
# Seibold, S., Gossner, M.M., Simons, N.K., Blüthgen, N., Müller, J., Ambarlı,
#   D., Ammer, C., Bauhus, J., Fischer, M., Habel, J.C., Linsenmair, K.E.,
#   Nauss, T., Penone, C., Prati, D., Schall, P., Schulze, E.-D., Vogt, J.,
#   Wöllauer, S., Weisser, W.W., 2019. Arthropod decline in grasslands and
#   forests is associated with landscape-level drivers. Nature 574, 671–674.
#   https://doi.org/10.1038/s41586-019-1684-3

# packages
pkgs <- c("here", "readr", "mgcv", "gratia", "dplyr", "ggplot2")
vapply(pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##
## Option 1:
##
## If you have the jonitor package or can install it
##   install.packages("janitor")
## then follow this

# data - from Seibold et al Nature 2019
# https://doi.org/10.25829/bexis.25786-1.3.11
seibold <- read_csv2(here("data/arthropod.csv"),
                     col_types = "cccccdndddddcddnddcddcddcddnnnnnnnnnnnn")

# basic data wrangling so we can plot
seibold <- seibold %>%
  janitor::clean_names() %>%
  rename(year = collection_year,
         region = exploratory,
         habitat = habitat_type) %>%
  mutate(across(c(habitat, sampling_regime, region, plot_id_year, plot_id),
                as.factor))
## skip to basic plot ==>

## ALTERNATIVE if you don't have janitor - sorry!
## use
arthropod_url <- "https://bit.ly/arthropod-csv"
seibold <- read_csv2(arthropod_url,
                     col_types = "cccccdndddddcddnddcddcddcddnnnnnnnnnnnn")
# basic data wrangling so we can plot
seibold <- seibold %>%
  rename(year = collection_year,
         region = exploratory,
         habitat = habitat_type) %>%
  mutate(across(c(habitat, sampling_regime, region, plot_id_year, plot_id),
                as.factor))
## end alternative

# basic plot <== skip to here if you do have janitor and used option 1
seibold %>%
  ggplot(aes(x = year, y = abundance_identified, group = plot_id)) +
  geom_line(alpha = 0.75) +
  facet_grid(habitat ~ region) +
  scale_y_log10() + theme_bw() +
  labs(x = NULL, y = "Abundance")

# more wrangling
seibold <- seibold %>%
  mutate(year_f = factor(year), # year as a factor for ranef
         year_c = year - 2012) # centre year

# filter out the forest sites; many have fewer than 9 observations
seibold <- seibold %>%
  filter(habitat == "grassland")

# abundance --------------------------------------------------------------------

# generate a simple plot with an appropriate GAM for illustration of fitting to
# a single site -- this is just for the figure
seibold %>%
  filter(plot_id == "AEG1") %>%
  ggplot(aes(x = year, y = abundance_identified)) +
  geom_line(alpha = 0.5) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 9),
              method.args = list(family = nb())) +
  labs(y = "Abundance", x = NULL) +
  theme_bw()

# select data for fitting a GAM to a single site AEG1
site <- seibold %>%
  filter(plot_id == "AEG1")

# Negative binomial with theta (constant) estimated
m_site <- gam(abundance_identified ~ s(year), data = site,
              method = "REML", family = nb())

summary(m_site) # model summary

draw(m_site) # partial effect of smooth

appraise(m_site, method = "simulate") # model diagnostics

k.check(m_site) # check basis size

### Models now for all grassland sites -----------------------------------------

# use multiple threads for fitting --- set this to number of physical cores
ctrl <- gam.control(nthreads = 4)

# fit to all data
m1 <- gam(abundance_identified ~ s(year) + # smooth of year (trend)
            s(plot_id, bs = "re"), # random effect of plot
          data = seibold,
          method = "ML",
          family = nb(),
          control = ctrl)

# perhaps region-specific trends?
m2 <- gam(abundance_identified ~ region + # regional means
            s(year, by = region) + # region-specific trends
            s(plot_id, bs = "re"), # random intercept for plot
          data = seibold,
          method = "ML",
          family = nb(),
          control = ctrl)

# individual random trends?
# about 8 minutes on my laptop
system.time(
m3 <- gam(abundance_identified ~
            s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
          data = seibold,
          method = "ML",
          family = nb(),
          control = ctrl)
)

# region-specific trends, plus individual trends
# about 19 minutes on my laptop
system.time(
m4 <- gam(abundance_identified ~ region + # regional means fixef
            s(year, by = region) + # region-specific smooths
            s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
          data = seibold,
          method = "ML",
          family = nb(),
          control = ctrl)
)

# region-specific trends, plus individual trends, plus year-to-year effects
# about 24 minutes on my laptop
system.time(
m5 <- gam(abundance_identified ~ region + # regional means fixef
            s(year_f, bs = "re") + # year-to-year effects
            s(year, by = region) + # region-specific smooths
            s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
          data = seibold,
          method = "ML",
          family = nb(),
          control = ctrl)
)

# a linear mixed effect version
# about 2 seconds on my laptop
system.time(
m6 <- gam(abundance_identified ~ region + region:year + # regional means fixef
            s(year_f, bs = "re") + # year-to-year effects
            s(year, plot_id, bs = "re"), # plot-specific linear trends
          data = seibold,
          method = "ML",
          family = nb(),
          control = ctrl)
)

## AIC
AIC(m1, m2, m3, m4, m5, m6)

## refit models that perform the best, with reml
m4_ml <- m4
m5_ml <- m5
# region-specific trends, plus individual trends
# about 10 minutes on my laptop
system.time(
m4 <- gam(abundance_identified ~ region + # regional means fixef
            s(year, by = region) + # region-specific smooths
            s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
          data = seibold,
          method = "REML",
          family = nb(),
          control = ctrl)
)
# region-specific trends, plus individual trends, plus year-to-year effects
# about 25 minutes on my laptop
system.time(
m5 <- gam(abundance_identified ~ region + # regional means fixef
            s(year_f, bs = "re") + # year-to-year effects
            s(year, by = region) + # region-specific smooths
            s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
          data = seibold,
          method = "REML",
          family = nb(),
          control = ctrl)
)

# plot smooths
draw(m5)

# model diagnostics
appraise(m5, method = "simulate")

# basis size
k.check(m5)

# rootogram
rootogram(m5, max_count = 300) %>%
  draw()

## other covariates
seibold <- seibold %>%
  mutate(year_s = scale(year),
         landuse_intensity_s = scale(landuse_intensity)[,1],
         mean_winter_temperature_s = scale(mean_winter_temperature)[,1],
         precipitation_sum_growing_preriod_s =
           scale(precipitation_sum_growing_preriod)[,1],
         grassland_cover_1000_s = scale(grassland_cover_1000)[,1],
         arable_cover_1000_s = scale(arable_cover_1000)[,1])

## about 1.5 minutes
system.time(
m_cov <- gam(abundance_identified ~
               s(year_s) + # overall trend
               s(landuse_intensity_s) +
               s(mean_winter_temperature_s) +
               s(precipitation_sum_growing_preriod_s) +
               s(grassland_cover_1000_s) +
               s(arable_cover_1000_s) +
               ti(mean_winter_temperature_s,
                  precipitation_sum_growing_preriod_s) +
               ti(year_s, landuse_intensity_s) +
               ti(year_s, grassland_cover_1000_s) +
               ti(year_s, arable_cover_1000_s) +
               s(year_f, bs = "re") + # year-to-year effects
               s(plot_id, bs = "re"), # site specific mean abundance
             family = nb(),
             method = "REML",
             control = ctrl,
             data = seibold,
             select = TRUE)
  )

# plot the smooths
sms <- smooths(m_cov)

# plot the univariate smooths
draw(m_cov, select = sms[1:6])

# plot the tensor product interation smooths
draw(m_cov, select = sms[7:10], rug = FALSE)

# plot the ranefs
draw(m_cov, select = sms[11:12])

# plot the overall trend effect
draw(m_cov, select = "s(year_s)")

# DON'T RUN IN WEBINAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# put in site specific trends
# ~17 minutes 
system.time(
m_cov2 <- gam(abundance_identified ~
                s(year_s) +
                s(landuse_intensity_s) +
                s(mean_winter_temperature_s) +
                s(precipitation_sum_growing_preriod_s) +
                s(grassland_cover_1000_s) +
                s(arable_cover_1000_s) +
                ti(mean_winter_temperature_s,
                   precipitation_sum_growing_preriod_s) +
                ti(year_s, landuse_intensity_s) +
                ti(year_s, grassland_cover_1000_s) +
                ti(year_s, arable_cover_1000_s) +
                s(year_f, bs = "re") +
                s(year_s, plot_id, bs = "fs", k = 5), # <-- here
              family = nb(),
              method = "REML",
              control = ctrl,
              data = seibold)
  )
