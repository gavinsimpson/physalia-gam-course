# packages
pkgs <- c("mgcv", "lme4", "ggplot2", "readr", "dplyr", "tidyr",
          "gratia")

# load packages
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
    quietly = TRUE)

#
shrimp <- read_csv(here("data", "trawl_nl.csv"),
    col_types = "dcdddddddddddd") %>%
    mutate(stratum = factor(stratum))

ggplot(shrimp) +
  geom_violin(aes(x = richness, y = factor(year))) +
    labs(x = "Number of species", y = "Year")

#
m_rich <- gam(richness ~ s(year) + s(stratum, bs = "re"),
              family = poisson,
              method = "REML",
              data = shrimp)

#
draw(m_rich)

# How to predict for the last year
fv <- fitted_values(m_rich, data = data_slice(m_rich, year = 2015:2016))

# Better control on uncertainty
knots <- list(year = c(2004, 2005, 2014, 2016))
bs_1 <- gam(richness ~ s(year, bs = "bs", k = 10, m = c(3, 2)) +
        s(stratum, bs = "re"),
    family = poisson,
    method = "REML",
    data = shrimp, knots = knots)

bs <- basis(s(year, bs = "bs", k = 10, m = c(3, 2)), knots = knots,
    data = data_slice(m_rich, year = seq(2004, 2016, length = 100)))
draw(bs)

summary(bs_1)

draw(bs_1)

fitted_values(m_rich, data = data_slice(m_rich, year = 2015:2016),
    exclude = "s(stratum)")

fitted_values(bs_1, data = data_slice(m_rich, year = 2015:2016),
    exclude = "s(stratum)")

ds <- data_slice(m_rich, year = seq(2005, 2016, length = 100))
fv_bs1 <- fitted_values(bs_1, data = ds, exclude = "s(stratum)")

fv_bs1 %>%
    ggplot(aes(x = year, y = fitted)) +
    geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.2) +
    geom_line()

# change the penalty
bs_2 <- gam(richness ~ s(year, bs = "bs", k = 10, m = c(3, 1)) +
    s(stratum, bs = "re"),
family = poisson,
method = "REML",
data = shrimp, knots = knots)

fv_bs2 <- fitted_values(bs_2, data = ds, exclude = "s(stratum)")

fv_bs1 %>%
    bind_rows(fv_bs2) %>%
    mutate(penalty = factor(rep(c(2, 1), each = 100))) %>%
    ggplot(aes(x = year, y = fitted, group = penalty)) +
    geom_ribbon(aes(ymax = upper, ymin = lower, fill = penalty), alpha = 0.2) +
    geom_line(aes(colour = penalty))

# change the penalty
bs_3 <- gam(richness ~ s(year, bs = "bs", k = 10, m = c(3, 0)) +
    s(stratum, bs = "re"),
family = poisson,
method = "REML",
data = shrimp, knots = knots)

fv_bs3 <- fitted_values(bs_3, data = ds, exclude = "s(stratum)")

fv_bs1 %>%
    bind_rows(fv_bs2, fv_bs3) %>%
    mutate(penalty = factor(rep(c(2, 1, 0), each = 100))) %>%
    ggplot(aes(x = year, y = fitted, group = penalty)) +
    geom_ribbon(aes(ymax = upper, ymin = lower, fill = penalty), alpha = 0.2) +
    geom_line(aes(colour = penalty))

# mix the penalties
bs_4 <- gam(richness ~ s(year, bs = "bs", k = 10, m = c(3, 2, 1)) +
    s(stratum, bs = "re"),
family = poisson,
method = "REML",
data = shrimp, knots = knots)

fv_bs4 <- fitted_values(bs_4, data = ds, exclude = "s(stratum)")

fv_bs1 %>%
    bind_rows(fv_bs2, fv_bs3, fv_bs4) %>%
    mutate(penalty = factor(rep(c("2", "1", "0", "2,1"), each = 100))) %>%
    ggplot(aes(x = year, y = fitted, group = penalty)) +
    geom_ribbon(aes(ymax = upper, ymin = lower, fill = penalty), alpha = 0.2) +
    geom_line(aes(colour = penalty))
