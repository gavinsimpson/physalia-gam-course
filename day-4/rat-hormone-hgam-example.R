# Analyse the Rat hormone data from Fahrmeir et al (2013) Regression: Models,
#  Methods and Applications. Springer
pkgs <- c("mgcv", "lme4", "ggplot2", "readr", "dplyr", "forcats", "tidyr",
          "gratia")

# load data
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

rats_url <- "https://bit.ly/rat-hormone"
rats <- read_table(rats_url, col_types = "dddddddddddd-")
# ignore the warning - it"s due to trailing white space at the ends of each
#   row in the file

rats <- rats %>%
    mutate(treatment = fct_recode(factor(group, levels = c(3,2,1)),
                                  Control = "3",
                                  Low = "1",
                                  High = "2"),
           subject = factor(subject))

rats %>%
    na.omit() %>%
    count(subject) %>%
    count(n, name = "n_rats")

plt_labs <- labs(y = "Head height (distance in pixels)",
                 x = "Age in days",
                 colour = "Treatment")

ggplot(rats, aes(x = time, y = response,
                 group = subject, colour = treatment)) +
    geom_line() +
    facet_wrap(~ treatment, ncol = 3) +
    plt_labs

K <- 7

## This is Model G
m1_hgam <- gam(response ~ s(time, k = K) +
                  s(subject, bs = "re"),
               data = rats, method = "REML")

## this is Model S at the treatment level
m2_hgam <- gam(response ~ s(time, treatment, bs = "fs", k = K) +
                  s(subject, bs = "re"),
               data = rats, method = "REML")

## this is Model I at the treatment level
m3_hgam <- gam(response ~ treatment +
                 s(time, by = treatment, k = K) +
                 s(subject, bs = "re"),
               data = rats, method = "REML")

## this is Model GI at the treatment level
m4_hgam <- gam(response ~ treatment + s(time, k = K) +
                  s(time, by = treatment, k = K, m = 1) +
                  s(subject, bs = "re"),
               data = rats, method = "REML")

## this is Model GS at the treatment level
m5_hgam <- gam(response ~ s(time, k = K) +
                  s(time, treatment, bs = "fs", k = K) +
                  s(subject, bs = "re"),
               data = rats, method = "REML")

## This is Model GS as the rat/subject level (GSS if you will)
m6_hgam <- gam(response ~ s(time, k = K) +
                  s(time, treatment, bs = "fs", k = K) +
                  s(time, subject, bs = "fs", k = 4), # not enough data for more
               data = rats, method = "REML")

## This is Model GS at the subject/rat level (SGS if you will)
## no top level common smooth effect
m7_hgam <- gam(response ~ s(time, treatment, bs = "fs", k = K) +
                  s(time, subject, bs = "fs", k = 4), # not enough data for more
               data = rats, method = "REML")

## This is Model S smooths only at the lowest subject level
m8_hgam <- gam(response ~ s(time, subject, bs = "fs", k = 5), # not enough data for more
               data = rats, method = "REML")

AIC(m1_hgam, m2_hgam, m3_hgam, m4_hgam, m5_hgam, m6_hgam, m7_hgam, m8_hgam)
